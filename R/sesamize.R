#' "fix" an RGset (for which IDATs may be unavailable) with Sesame
#' The input is an RGSet and the output is a sesamized RGSet
#' 
#' @param x an RGChannelSet, perhaps with colData of various flavors
#' @param naFrac maximum NA fraction for a probe before it gets dropped (1)
#' @param parallel attempt to run in parallel? (This is a bad idea on laptops)
#' @return a sesamized GenomicRatioSet from the input RGChannelSet 
#' @import parallel
#' @importFrom S4Vectors metadata
#' @importFrom S4Vectors metadata<-
#' 
#' @examples
#' # Takes about two minutes to process 48 samples on my 48-core desktop
#' if (require(FlowSorted.CordBloodNorway.450k)) {
#'     sesamized <- sesamize(
#'         FlowSorted.CordBloodNorway.450k[,1:2])
#' } 
#' @export 
sesamize <- function(x, naFrac=1, parallel=FALSE) { 
    stopifnot(is(x, "RGChannelSet"))
    pkgTest('minfi')
    pkgTest('SummarizedExperiment')
    
    if (ncol(x) > 1) {
        nameses <- colnames(x)
        names(nameses) <- nameses
        if (parallel == TRUE) { 
            res <- do.call(
                SummarizedExperiment::cbind,
                mclapply(nameses, function(y) sesamize(x[,y])))
        } else { 
            res <- do.call(
                SummarizedExperiment::cbind, 
                lapply(nameses, function(y) sesamize(x[,y])))
        }
        res <- minfi::mapToGenome(res)
        kept <- which((rowSums(is.na(
            minfi::getBeta(res))) / ncol(res)) <= naFrac)
        if (length(kept) < 1) 
            stop("No probes survived with naFrac <= ",naFrac,".")
        mfst <- packageVersion(paste(minfi::annotation(x), collapse="anno."))
        res@preprocessMethod <- c(
            rg.norm="SeSAMe (type I)",
            p.value="SeSAMe (pOOBAH)",
            sesame=as.character(packageVersion("sesame")),
            minfi=as.character(packageVersion("minfi")),
            manifest=as.character(mfst))
        
        metadata(res)$SNPs <- minfi::getSnpBeta(x)
        SummarizedExperiment::assays(res)[["M"]] <- NULL 
        SummarizedExperiment::colData(res) <- SummarizedExperiment::colData(x)
        return(res[kept, ])
    } else {
        message("Sesamizing ", colnames(x), "...")
        dm <- cbind(G=as.matrix(minfi::getGreen(x)),
                    R=as.matrix(minfi::getRed(x)))
        colnames(dm) <- c("G", "R") # just in case...
        attr(dm, "platform") <- sub("HMEPIC", "EPIC", 
            sub("IlluminaHumanMethylation", "HM", 
                sub("k$", "", minfi::annotation(x)["array"])))
        sset <- chipAddressToSignal(dm) # see above for kludge
        Beta <- as.matrix(getBetas(dyeBiasCorrTypeINorm(noob(sset))))
        CN <- as.matrix(log2(totalIntensities(sset))[rownames(Beta)])
        minfi::RatioSet(Beta=Beta, CN=CN, annotation=minfi::annotation(x))
    }
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


.SigSetToRGChannelSet <- function(sset) {
    
    dfAddress <- sesameDataGet(paste0(sset@platform,'.address'))
    SSRed <- NULL
    SSGrn <- NULL
    
    IIdf <- dfAddress$ordering[
        dfAddress$ordering$COLOR_CHANNEL=='Both', c('Probe_ID','U')]
    SSRed <- c(SSRed, setNames(sset@II[match(
        IIdf$Probe_ID, rownames(sset@II)), 'U'],
        as.character(IIdf$U)))
    SSGrn <- c(SSGrn, setNames(sset@II[match(
        IIdf$Probe_ID, rownames(sset@II)), 'M'],
        as.character(IIdf$U)))
    
    IRdf <- dfAddress$ordering[
        dfAddress$ordering$COLOR_CHANNEL=='Red', c('Probe_ID','M','U')]
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
    
    IGdf <- dfAddress$ordering[
        dfAddress$ordering$COLOR_CHANNEL=='Grn', c('Probe_ID','M','U')]
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
    
    # controls
    control.names <- make.names(dfAddress$controls$Name, unique = TRUE)
    SSGrn <- c(SSGrn, setNames(sset@ctl[match(
        control.names, rownames(sset@ctl)),'G'], 
        as.character(dfAddress$controls$Address)))
    SSRed <- c(SSRed, setNames(sset@ctl[match(
        control.names, rownames(sset@ctl)),'R'], 
        as.character(dfAddress$controls$Address)))
    
    list(grn=SSGrn, red=SSRed)
}

#' Convert SeSAMe::SigSet to minfi::RGChannelSet
#' 
#' @param ssets a list of SigSet
#' @return a minfi::RGChannelSet
#' @examples
#'
#' sset <- sesameDataGet('EPIC.1.LNCaP')$sset
#' rgSet <- SigSetToRGChannelSet(sset)
#'
#' @export 
SigSetToRGChannelSet <- function(ssets) {
    if (is(ssets, 'SigSet')) {
        ssets <- list(sample=ssets)
    }

    platform <- ssets[[1]]@platform
    
    ss_all <- lapply(ssets, .SigSetToRGChannelSet)
    rgset <- minfi::RGChannelSet(
        Green=do.call(cbind, lapply(ss_all, function(ss) ss$grn)), 
        Red=do.call(cbind, lapply(ss_all, function(ss) ss$red)), 
        annotation=c(
            array=platformSmToMinfi(platform),
            annotation=NA)) # can be ilm10b4.hg19 or ilmn12.hg19
}

#' Convert RGChannelSet (minfi) to SigSet (SeSAMe)
#'
#' Notice the colData() and rowData() is lost. Most cases, rowData is empty
#' anyway.
#' 
#' @param rgSet a minfi::RGChannelSet
#' @return a list of sesame::SigSet
#' @examples
#'
#' if (require(FlowSorted.Blood.450k)) {
#'     rgSet <- FlowSorted.Blood.450k[,1:2]
#'     ssets <- RGChannelSetToSigSet(rgSet)
#' }
#' @export
RGChannelSetToSigSet <- function(rgSet) {

    pkgTest('minfi')
    rg_grn <- minfi::getGreen(rgSet)
    rg_red <- minfi::getRed(rgSet)
    samples <- colnames(rg_grn) # rg.red should be the same
    prb_addr_grn <- rownames(rg_grn)
    prb_addr_red <- rownames(rg_red)
    
    platform <- platformMinfiToSm(minfi::annotation(rgSet)[['array']])
    
    # Process mappings
    addr <- sesameDataGet(paste0(platform,'.address'))
    IIdf <- addr$ordering[
        addr$ordering$COLOR_CHANNEL=='Both', c('Probe_ID','U')]
    IItoSSU <- match(IIdf$U, as.numeric(prb_addr_red))
    IItoSSM <- match(IIdf$U, as.numeric(prb_addr_grn))
    
    IGdf <- addr$ordering[
        addr$ordering$COLOR_CHANNEL=='Grn', c('Probe_ID','M','U')]
    IGtoSSU <- match(IGdf$U, as.numeric(prb_addr_grn))
    IGtoSSM <- match(IGdf$M, as.numeric(prb_addr_grn))
    oobRtoSSU <- match(IGdf$U, as.numeric(prb_addr_red))
    oobRtoSSM <- match(IGdf$M, as.numeric(prb_addr_red))
    
    IRdf <- addr$ordering[
        addr$ordering$COLOR_CHANNEL=='Red', c('Probe_ID','M','U')]
    IRtoSSU <- match(IRdf$U, as.numeric(prb_addr_red))
    IRtoSSM <- match(IRdf$M, as.numeric(prb_addr_red))
    oobGtoSSU <- match(IRdf$U, as.numeric(prb_addr_grn))
    oobGtoSSM <- match(IRdf$M, as.numeric(prb_addr_grn))
    
    controlToSSG <- match(addr$controls$Address, as.numeric(prb_addr_red))
    controlToSSR <- match(addr$controls$Address, as.numeric(prb_addr_grn))
    
    # actually conversion
    setNames(lapply(samples, function(sample) {
        .II <- cbind(M=rg_grn[IItoSSM, sample], U=rg_red[IItoSSU, sample])
        rownames(.II) <- IIdf$Probe_ID
        
        .IG <- cbind(M=rg_grn[IGtoSSM, sample], U=rg_grn[IGtoSSU, sample])
        rownames(.IG) <- IGdf$Probe_ID
        
        .IR <- cbind(M=rg_red[IRtoSSM, sample], U=rg_red[IRtoSSU, sample])
        rownames(.IR) <- IRdf$Probe_ID
        
        .oobG <- cbind(
            M=rg_grn[oobGtoSSM, sample], 
            U=rg_grn[oobGtoSSU, sample])
        rownames(.oobG) <- IRdf$Probe_ID
        
        .oobR <- cbind(
            M=rg_red[oobRtoSSM, sample], 
            U=rg_red[oobRtoSSU, sample])
        rownames(.oobR) <- IGdf$Probe_ID
        
        .ctl <- data.frame(
            G=rg_grn[controlToSSG, sample], 
            R=rg_red[controlToSSR, sample], 
            col=addr$controls$Color_Channel,
            type=addr$controls$Type, stringsAsFactors = FALSE)
        rownames(.ctl) <- make.names(addr$controls$Name, unique = TRUE)
        
        sset <- new('SigSet', IG=.IG, IR=.IR, II=.II,
            oobG=.oobG, oobR=.oobR, ctl=.ctl,
            platform=platform)
        detectionPoobEcdf(sset)
    }), samples)
}

