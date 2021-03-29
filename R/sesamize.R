#' "fix" an RGChannelSet (for which IDATs may be unavailable) with Sesame
#' The input is an RGSet and the output is a sesamized minfi::GenomicRatioSet
#' 
#' \code{HDF5Array} package required.
#' 
#' @param rgSet an RGChannelSet, perhaps with colData of various flavors
#' @param naFrac maximum NA fraction for a probe before it gets dropped (1)
#' @param BPPARAM get parallel with MulticoreParam(n)
#' @param HDF5 is the rgSet HDF5-backed? if so, avoid eating RAM (perhaps)
#' @param HDF5SEdestination character(1) path to where the
#' HDF5-backed GenomicRatioSet will be stored
#' @param replace logical(1) passed to saveHDF5SummarizedExperiment
#' @note We employ BPREDO for a second chance if bplapply hits an error.
#' @return a sesamized GenomicRatioSet
#' @import BiocParallel
#' @importFrom SummarizedExperiment start
#' @importFrom SummarizedExperiment end
#' @importFrom SummarizedExperiment rowRanges
#' @importFrom SummarizedExperiment assays
#' @importFrom SummarizedExperiment assays<-
#' @importFrom SummarizedExperiment colData
#' @importFrom SummarizedExperiment colData<-
#' @importFrom S4Vectors metadata
#' @importFrom S4Vectors metadata<-
#'
#' @export 
sesamize <- function(
    rgSet, naFrac=1, BPPARAM=SerialParam(), HDF5=NULL,
    HDF5SEdestination=paste0(tempdir(check=TRUE), "/sesamize_HDF5_scratch"),
    replace=FALSE) {

    stopifnot(is(rgSet, "RGChannelSet"))

    pkgTest('minfi')
    pkgTest('SummarizedExperiment')
    samples <- colnames(rgSet)
    names(samples) <- samples

    if (is.null(HDF5)) {
        ## are we working on an HDF5-backed RGChannelSet?
        HDF5 <- (class(assays(rgSet)[[1]])[1] == "DelayedMatrix")
    }
    t1 =  bptry(bplapply(samples, function(sample) {
            message("Sesamizing ", sample, "...")
            sset <- RGChannelSet1ToSigSet(rgSet[,sample])
            sset <- dyeBiasCorrTypeINorm(noob(sset))
            SigSetToRatioSet(sset)}, BPPARAM=BPPARAM))
    lk = vapply(t1, inherits, logical(1), "bperror")  # second try?
    if (any(lk)) {
      t1 =  bptry(bplapply(samples, function(sample) {
            message("Sesamizing ", sample, "...")
            sset <- RGChannelSet1ToSigSet(rgSet[,sample])
            sset <- dyeBiasCorrTypeINorm(noob(sset))
            SigSetToRatioSet(sset)}, BPREDO=t1, BPPARAM=BPPARAM))
      }

    ratioSet <- do.call(
        SummarizedExperiment::cbind, t1
    )
    colnames(ratioSet) = colnames(rgSet)
    if (HDF5) {
        pkgTest('HDF5Array')
        #td <- paste(tempdir(check=TRUE), "sesamize_HDF5_scratch", sep="/")
        ratioSet <- HDF5Array::saveHDF5SummarizedExperiment(
            ratioSet, dir=HDF5SEdestination, replace=replace) #td, replace=TRUE)
    }

    ## mapping occurs first, SNPs get separated here
    ratioSet <- minfi::mapToGenome(ratioSet)
    
    ## keep only probes surviving naFrac
    kept <- seq_len(nrow(ratioSet))
    if (naFrac < 1) { 
        kept <- which((rowSums(is.na(minfi::getBeta(ratioSet))) / 
                           ncol(ratioSet)) <= naFrac)
        if (length(kept) < 1) 
            stop("No probes survived with naFrac <= ",naFrac,".")
    } 
    
    ## put back colData(), @processMethod, $SNPs
    mfst <- packageVersion(paste(minfi::annotation(ratioSet), collapse="anno."))
    ratioSet@preprocessMethod <- c(
        rg.norm="SeSAMe (type I)",
        p.value="SeSAMe (pOOBAH)",
        sesame=as.character(packageVersion("sesame")),
        minfi=as.character(packageVersion("minfi")),
        manifest=as.character(mfst))

    ## SNP not adjusted in minfi, so keep them that way
    metadata(ratioSet)$SNPs <- minfi::getSnpBeta(rgSet)
    assays(ratioSet)[["M"]] <- NULL 
    colData(ratioSet) <- colData(rgSet)
    
    return(ratioSet[kept, ])
}

platformMinfiToSm <- function(platform) {
    plf <- sub("HMEPIC", "EPIC", 
        sub("IlluminaHumanMethylation", "HM", 
            sub("k$", "", platform)))
    stopifnot(plf %in% c('EPIC','HM450','HM27'))
    plf
}

platformSmToMinfi <- function(platform) {
    plf <- sub(
        "HM", "IlluminaHumanMethylation",
        sub("EPIC", "HMEPIC",
            sub("(450|27)$", "\\1k", platform)))
    stopifnot(plf %in% c(
        'IlluminaHumanMethylationEPIC',
        'IlluminaHumanMethylation450k',
        'IlluminaHumanMethylation27k'
    ))
    plf
}

## reverse of chipAddressToSignal
SigSetToRGChannel <- function(sset, manifest = NULL, controls = NULL) {

    if (is.null(manifest)) {
        dfAddress <- sesameDataGet(paste0(sset@platform,'.address'))
        manifest <- dfAddress$ordering
        controls <- dfAddress$controls
    }

    SSRed <- NULL
    SSGrn <- NULL
    
    IIdf <- manifest[
        manifest$COLOR_CHANNEL=='Both', c('Probe_ID','U')]
    SSRed <- c(SSRed, setNames(sset@II[match(
        IIdf$Probe_ID, rownames(sset@II)), 'U'],
        as.character(IIdf$U)))
    SSGrn <- c(SSGrn, setNames(sset@II[match(
        IIdf$Probe_ID, rownames(sset@II)), 'M'],
        as.character(IIdf$U)))
    
    IRdf <- manifest[
        manifest$COLOR_CHANNEL=='Red', c('Probe_ID','M','U')]
    SSRed <- c(SSRed, setNames(sset@IR[match(
        IRdf$Probe_ID, rownames(sset@IR)), 'M'],
        as.character(IRdf$M)))
    SSRed <- c(SSRed, setNames(sset@IR[match(
        IRdf$Probe_ID, rownames(sset@IR)), 'U'],
        as.character(IRdf$U)))
    ## OOB signals
    SSGrn <- c(SSGrn, setNames(sset@oobG[match(
        IRdf$Probe_ID, rownames(sset@oobG)), 'M'],
        as.character(IRdf$M)))
    SSGrn <- c(SSGrn, setNames(sset@oobG[match(
        IRdf$Probe_ID, rownames(sset@oobG)), 'U'],
        as.character(IRdf$U)))
    
    IGdf <- manifest[
        manifest$COLOR_CHANNEL=='Grn', c('Probe_ID','M','U')]
    SSGrn <- c(SSGrn, setNames(sset@IG[match(
        IGdf$Probe_ID, rownames(sset@IG)), 'M'], 
        as.character(IGdf$M)))
    SSGrn <- c(SSGrn, setNames(sset@IG[match(
        IGdf$Probe_ID, rownames(sset@IG)), 'U'], 
        as.character(IGdf$U)))
    ## OOB signals
    SSRed <- c(SSRed, setNames(sset@oobR[match(
        IGdf$Probe_ID, rownames(sset@oobR)), 'M'],
        as.character(IGdf$M)))
    SSRed <- c(SSRed, setNames(sset@oobR[match(
        IGdf$Probe_ID, rownames(sset@oobR)), 'U'],
        as.character(IGdf$U)))
    
    ## controls
    if (!is.null(controls)) {
        control.names <- make.names(controls$Name, unique = TRUE)
        SSGrn <- c(SSGrn, setNames(sset@ctl[match(
            control.names, rownames(sset@ctl)),'G'], 
            as.character(controls$Address)))
        SSRed <- c(SSRed, setNames(sset@ctl[match(
            control.names, rownames(sset@ctl)),'R'], 
            as.character(controls$Address)))
    } ## else TODO controls obtained from manifest
    
    list(grn=SSGrn, red=SSRed)
}

## annotation, if not given is guessed
guessMinfiAnnotation <- function(platform, annotation = NA) {
    if (is.na(annotation)) {
        if (platform %in% c("HM450", "HM27")) {
            'ilmn12.hg19'
        } else { # EPIC
            'ilm10b4.hg19'
        }
    } else {
        annotation
    }
}

#' Convert sesame::SigSet to minfi::RGChannelSet
#' 
#' @param ssets a list of sesame::SigSet
#' @param BPPARAM get parallel with MulticoreParam(n)
#' @param annotation the minfi annotation string, guessed if not given
#' @return a minfi::RGChannelSet
#' @import BiocParallel
#' @examples
#'
#' sset <- sesameDataGet('EPIC.1.LNCaP')$sset
#' rgSet <- SigSetsToRGChannelSet(sset)
#'
#' @export 
SigSetsToRGChannelSet <- function(ssets, BPPARAM=SerialParam(), annotation=NA) {
    if (is(ssets, 'SigSet')) {
        ssets <- list(sample=ssets)
    }

    platform <- ssets[[1]]@platform
    annotation <- guessMinfiAnnotation(annotation)
    
    ss_all <- bplapply(ssets, SigSetToRGChannel, BPPARAM=BPPARAM)
    rgset <- minfi::RGChannelSet(
        Green=do.call(cbind, lapply(ss_all, function(ss) ss$grn)), 
        Red=do.call(cbind, lapply(ss_all, function(ss) ss$red)), 
        annotation=c(
            array=unname(platformSmToMinfi(platform)),
            annotation=annotation))
}

## helper: convert RGChannelSet of one sample
RGChannelSet1ToSigSet <- function(rgSet1, manifest = NULL, controls = NULL) {

    pkgTest('minfi')
    stopifnot(ncol(rgSet1) == 1)
    
    # chipaddress/rownames are automatically the same
    dm <- cbind(
        G=as.matrix(minfi::getGreen(rgSet1)),
        R=as.matrix(minfi::getRed(rgSet1)))
    
    colnames(dm) <- c('G','R') # just in case..
    if (is.null(manifest)) {
        attr(dm, 'platform') <- platformMinfiToSm(
            minfi::annotation(rgSet1)['array'])
        df_address <- sesameDataGet(paste0(
            attr(dm, 'platform'), '.address'))
        manifest <- df_address$ordering
        controls <- df_address$controls
    }

    pOOBAH(chipAddressToSignal(dm, manifest, controls))
}

#' Convert RGChannelSet (minfi) to a list of SigSet (SeSAMe)
#'
#' Notice the colData() and rowData() is lost. Most cases, rowData is empty
#' anyway.
#'
#' @param rgSet a minfi::RGChannelSet
#' @param BPPARAM get parallel with MulticoreParam(n)
#' @param manifest manifest file
#' @return a list of sesame::SigSet
#' @import BiocParallel
#' @examples
#'
#' if (require(FlowSorted.Blood.450k)) {
#'     rgSet <- FlowSorted.Blood.450k[,1:2]
#'     ssets <- RGChannelSetToSigSets(rgSet)
#' }
#' @export
RGChannelSetToSigSets <- function(
    rgSet, manifest=NULL, BPPARAM=SerialParam()) {

    pkgTest('minfi')
    samples <- colnames(rgSet)
    setNames(bplapply(
        samples, function(sample) {
        RGChannelSet1ToSigSet(rgSet[,sample], manifest=manifest)
    }, BPPARAM=BPPARAM), samples)
}

#' Convert one sesame::SigSet to minfi::RatioSet
#'
#' @param sset a sesame::SigSet
#' @param annotation minfi annotation string
#' @return a minfi::RatioSet
#' @examples
#'
#' sset <- sesameDataGet('EPIC.1.LNCaP')$sset
#' ratioSet <- SigSetToRatioSet(sset)
#' 
#' @export
SigSetToRatioSet <- function(sset, annotation = NA) {
    Beta <- as.matrix(getBetas(sset))
    CN <- as.matrix(log2(totalIntensities(sset))[rownames(Beta)])
    annotation <- guessMinfiAnnotation(sset@platform, annotation)
    platform <- platformSmToMinfi(sset@platform)
    minfi::RatioSet(Beta = Beta, CN = CN, annotation = c(
        array = unname(platform), annotation = annotation))
}

