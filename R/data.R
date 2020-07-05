

#' Retrieve manifest file from the supporting website
#' at http://zwdzwd.github.io/InfiniumAnnotation
#'
#' @param platform Infinium platform
#' @param refversion human reference version, irrelevant for mouse array
#' @param version manifest version, default to the latest/current.
#' @param probeType cg, ch or rs, default to all probes
#' @param designType I (Infinium-I) or II (Infinium-II), default to both
#' @return manifest file of requested probes
#' @examples
#'
#' mft <- getManifest('HM27', 'hg38')
#' 
#' @export
getManifest <- function(
    platform=c('EPIC','HM450','HM27'),
    refversion=c('hg19','hg38'),
    version="current",
    probeType=c('all','cg','ch','rs'),
    designType=c('all','I','II')) {

    platform <- match.arg(platform)
    refversion <- match.arg(refversion)
    probeType <- match.arg(probeType)
    designType <- match.arg(designType)
        
    download_path <-
        sprintf(
            'http://zwdzwd.io/InfiniumAnnotation/%s/%s/%s.%s.manifest.rds',
            version, platform, platform, refversion)

    cat("Retrieving manifest from ",download_path, "... ")
    mft <- readRDS(url(download_path))
    cat("Done.\n")
    if (probeType != 'all')
        mft <- mft[mft$probeType == probeType]

    if (designType != 'all')
        mft <- mft[mft$designType == designType]
    
    mft
}

#' Retrieve variant annotation file for explicit rs probes
#' from the supporting website
#' at http://zwdzwd.github.io/InfiniumAnnotation
#'
#' @param platform Infinium platform
#' @param refversion human reference version, irrelevant for mouse array
#' @param version manifest version, default to the latest/current.
#' @return variant annotation file of explicit rs probes
#' @examples
#'
#' annoS <- getVariantAnno_SNP('EPIC', 'hg38')
#' 
#' @export
getVariantAnno_SNP <- function(
    platform = c('EPIC'),
    refversion = c('hg19','hg38'),
    version = '20200704') {

    platform <- match.arg(platform)
    refversion <- match.arg(refversion)

    download_path <-
        sprintf(
            paste0(
                'http://zwdzwd.io/InfiniumAnnotation/',
                '%s/%s/%s.%s.snp_overlap_b151.rds'),
            version, platform, platform, refversion)

    cat("Retrieving SNP annotation from ",download_path, "... ")
    anno <- readRDS(url(download_path))
    cat("Done.\n")
    
    anno
}

#' Retrieve variant annotation file for Infinium-I probes
#' from the supporting website
#' at http://zwdzwd.github.io/InfiniumAnnotation
#'
#' @param platform Infinium platform
#' @param refversion human reference version, irrelevant for mouse array
#' @param version manifest version, default to the latest/current.
#' @return variant annotation file of infinium I probes
#' @examples
#'
#' annoI <- getVariantAnno_InfiniumI('EPIC', 'hg38')
#' 
#' @export
getVariantAnno_InfiniumI <- function(
    platform = c('EPIC'),
    refversion = c('hg19','hg38'),
    version = '20200704') {

    platform <- match.arg(platform)
    refversion <- match.arg(refversion)

    download_path <-
        sprintf(
            paste0(
                'http://zwdzwd.io/InfiniumAnnotation/',
                '%s/%s/%s.%s.typeI_overlap_b151.rds'),
            version, platform, platform, refversion)

    cat("Retrieving SNP annotation from ",download_path, "... ")
    anno <- readRDS(url(download_path))
    cat("Done.\n")
    
    anno
}
