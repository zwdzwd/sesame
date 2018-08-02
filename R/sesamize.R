#' "fix" an RGset (for which IDATs may be unavailable) with Sesame
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
#' 
#' # Takes about two minutes to process 48 samples on my 48-core desktop
#' if (require(FlowSorted.CordBloodNorway.450k)) {
#'     sesamized <- sesamize(
#'         FlowSorted.CordBloodNorway.450k[,1:2])
#' } 
#' 
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
        dm <- cbind(G=minfi::getGreen(x), R=minfi::getRed(x))
        attr(dm, "platform") <- sub(
            "IlluminaHumanMethylation", "HM", 
            sub("k$", "", minfi::annotation(x)["array"]))
        sset <- chipAddressToSignal(dm) # see above for kludge
        Beta <- as.matrix(getBetas(dyeBiasCorrTypeINorm(noob(sset))))
        CN <- as.matrix(log2(totalIntensities(sset))[rownames(Beta)])
        minfi::RatioSet(
            Beta=Beta, CN=CN, annotation=minfi::annotation(x))
    }
}


