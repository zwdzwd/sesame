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
#' \donttest{
#' # Takes about two minutes to process 48 samples on my 48-core desktop
#' if (require(FlowSorted.CordBloodNorway.450k)) {
#'     sesamized <- sesamize(
#'         FlowSorted.CordBloodNorway.450k[,1:2])
#' } 
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
        sub(
            "EPIC", "HMEPIC",
            sub("(450|27)$", "\\1k", platform)))
    stopifnot(plf %in% c(
        'IlluminaHumanMethylationEPIC',
        'IlluminaHumanMethylation450k',
        'IlluminaHumanMethylation27k'
    ))
    plf
}


SigSetToRGChannelSet <- function(sset) {
    
}

#' Conversion of RGChannelSet (Minfi) to SigSet (SeSAMe)
#' 
#' @param rgSet - a minfi::RGChannelSet
#' @return a sesame::SigSet
#' @export
RGChannelSetToSigSet <- function(rgSet) {
    
    rg_grn <- getGreen(rgSet)
    rg_red <- getRed(rgSet)
    samples <- colnames(rg_grn) # rg.red should be the same
    probe_address_grn <- rownames(rg_grn)
    probe_address_red <- rownames(rg_red)
    
    platform <- platformMinfiToSm(annotation(rgSet)[['array']])
    
    # Process mappings
    dfAddress <- sesameDataGet(paste0(platform,'.address'))
    IIdf <- dfAddress$ordering[
        dfAddress$ordering$COLOR_CHANNEL=='Both', c('Probe_ID','U')]
    IItoSSU <- match(IIdf$U, as.numeric(probe_address_red))
    IItoSSM <- match(IIdf$U, as.numeric(probe_address_grn))
    
    IGdf <- dfAddress$ordering[
        dfAddress$ordering$COLOR_CHANNEL=='Grn', c('Probe_ID','M','U')]
    IGtoSSU <- match(IGdf$U, as.numeric(probe_address_grn))
    IGtoSSM <- match(IGdf$M, as.numeric(probe_address_grn))
    oobRtoSSU <- match(IGdf$U, as.numeric(probe_address_red))
    oobRtoSSM <- match(IGdf$M, as.numeric(probe_address_red))
    
    IRdf <- dfAddress$ordering[
        dfAddress$ordering$COLOR_CHANNEL=='Red', c('Probe_ID','M','U')]
    IRtoSSU <- match(IRdf$U, as.numeric(probe_address_red))
    IRtoSSM <- match(IRdf$M, as.numeric(probe_address_red))
    oobGtoSSU <- match(IRdf$U, as.numeric(probe_address_grn))
    oobGtoSSM <- match(IRdf$M, as.numeric(probe_address_grn))
    
    controlToSSG <- match(dfAddress$controls$Address, as.numeric(probe_address_red))
    controlToSSR <- match(dfAddress$controls$Address, as.numeric(probe_address_grn))
    
    # actually conversion
    lapply(samples, function(sample) {
        .II <- cbind(M=rg_grn[IItoSSM, sample], U=rg_red[IItoSSU, sample])
        rownames(.II) <- IIdf$Probe_ID
        
        .IG <- cbind(M=rg_grn[IGtoSSM, sample], U=rg_grn[IGtoSSU, sample])
        rownames(.IG) <- IGdf$Probe_ID
        
        .IR <- cbind(M=rg_red[IRtoSSM, sample], U=rg_red[IRtoSSU, sample])
        rownames(.IR) <- IRdf$Probe_ID
        
        .oobG <- cbind(
            M=rg_grn[oobGtoSSM, sample], 
            U=rg_grn[oobGtoSSU, sample])
        
        .oobR <- cbind(
            M=rg_red[oobRtoSSM, sample], 
            U=rg_red[oobRtoSSU, sample])
        
        .ctl <- data.frame(
            G=rg_grn[controlToSSG, sample], 
            R=rg_red[controlToSSR, sample], 
            col=dfAddress$controls$Color_Channel, 
            type=dfAddress$controls$Type)
        rownames(.ctl) <- make.names(dfAddress$controls$Name, unique = TRUE)
        
        # TODO: pval
        new('SigSet', IG=.IG, IR=.IR, II=.II,
            oobG=.oobG, oobR=.oobR, ctl=.ctl,
            platform=platform)
    })
}

