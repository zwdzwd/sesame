
#' De-identify IDATs by removing SNP probes
#'
#' Mask SNP probe intensity mean by zero.
#'
#' @param path input IDAT file
#' @param out_path output IDAT file
#' @param snps SNP definition, if not given, default to SNP probes
#' @param mft sesame-compatible manifest if non-standard
#' @param secret a numeric integer, max 2^31-1
#' when not NULL randomize the signal intensities
#' when NULL (default), this sets all SNP intensities to zero.
#' @return NULL
#' @examples
#'
#' temp_out <- tempfile("test")
#' deIdentify(system.file(
#'     "extdata", "4207113116_A_Grn.idat", package = "sesameData"),
#'      temp_out)
#' unlink(temp_out)
#' @export
deIdentify <- function(path, out_path=NULL, snps=NULL, mft=NULL, secret=NULL) {

    res <- suppressWarnings(illuminaio::readIDAT(path))
    platform <- inferPlatform(res)

    if(is.null(out_path)) {
        pfx <- sub('.idat(.gz)?$','', path)
        if(grepl('_Grn$', pfx)) {
            out_path <- paste0(sub('_Grn$','',pfx), '_noid_Grn.idat')
        } else if (grepl('_Red$', pfx)) {
            out_path <- paste0(sub('_Red$','',pfx), '_noid_Red.idat')
        }
    }

    if (is.null(mft)) {
        mft <- sesameDataGet(paste0(platform, '.address'))$ordering
    }
    if (is.null(snps))
        snps <- grep("^rs", mft$Probe_ID, value=TRUE)
    mft <- mft[mft$Probe_ID %in% snps,]
    
    snpsTango <- na.omit(c(mft$M, mft$U))
    qt <- res$Quants
    snpsIdx <- match(snpsTango, rownames(qt))
    dt <- qt[,'Mean']
    if (is.null(secret)) {
        dt[snpsIdx] <- 0
    } else {
        set.seed(secret)
        dt[snpsIdx] <- sample(dt[snpsIdx])
    }
    
    if(grepl("\\.gz$", path)) {
        con <- gzfile(path, "rb")
    } else {
        con <- file(path, "rb")
    }
    
    con2 <- file(out_path, "wb")

    ## before Mean section
    writeBin(readBin(con, "raw", n = res$fields["Mean", 'byteOffset']), con2)

    ## write new Mean section
    writeBin(as.integer(dt), con2, size=2, endian='little')

    ## after Mean section
    a <- readBin(con, "raw", n = res$nSNPsRead*2) # skip by reading..., seek might not work for gzfile
    while (length(a <- readBin(con, 'raw', n=1))>0) writeBin(a, con2)

    close(con)
    close(con2)
}

#' Re-identify IDATs by restoring scrambled SNP intensities
#'
#' This requires a secret number that was used to de-idenitfy the IDAT
#'
#' @param path input IDAT file
#' @param secret the numeric number used to de-identify the IDAT
#' @param out_path output IDAT file
#' @param snps SNP definition, if not given, default to SNP probes
#' @param mft sesame-compatible manifest if non-standard
#' @return NULL
#' @examples
#'
#' temp_out <- tempfile("test")
#' reIdentify(system.file(
#'     "extdata", "4207113116_A_Grn.idat", package = "sesameData"),
#'      temp_out, 123)
#' unlink(temp_out)
#' @export
reIdentify <- function(path, secret, out_path=NULL, snps=NULL, mft=NULL) {

    res <- suppressWarnings(illuminaio::readIDAT(path))
    platform <- inferPlatform(res)

    if(is.null(out_path)) {
        pfx <- sub('.idat(.gz)?$','', path)
        if(grepl('_Grn$', pfx)) {
            out_path <- paste0(sub('_Grn$','',pfx), '_reid_Grn.idat')
        } else if (grepl('_Red$', pfx)) {
            out_path <- paste0(sub('_Red$','',pfx), '_reid_Red.idat')
        }
    }

    if (is.null(mft)) {
        mft <- sesameDataGet(paste0(platform, '.address'))$ordering
    }
    if (is.null(snps))
        snps <- grep("^rs", mft$Probe_ID, value=TRUE)
    mft <- mft[mft$Probe_ID %in% snps,]
    
    snpsTango <- na.omit(c(mft$M, mft$U))
    qt <- res$Quants
    snpsIdx <- match(snpsTango, rownames(qt))
    dt <- qt[,'Mean']
    idx <- seq_along(snpsIdx)
    set.seed(secret)
    dt[snpsIdx] <- dt[snpsIdx[match(idx, sample(idx))]]
    
    if(grepl("\\.gz$", path)) {
        con <- gzfile(path, "rb")
    } else {
        con <- file(path, "rb")
    }
    
    con2 <- file(out_path, "wb")

    ## before Mean section
    writeBin(readBin(con, "raw", n = res$fields["Mean", 'byteOffset']), con2)

    ## write new Mean section
    writeBin(as.integer(dt), con2, size=2, endian='little')

    ## after Mean section
    a <- readBin(con, "raw", n = res$nSNPsRead*2) # skip by reading..., seek might not work for gzfile
    while (length(a <- readBin(con, 'raw', n=1))>0) writeBin(a, con2)

    close(con)
    close(con2)
}
