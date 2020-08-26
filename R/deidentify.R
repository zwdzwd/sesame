
#' De-identify IDATs by removing SNP probes
#'
#' @param path input IDAT file
#' @param out_path output IDAT file
#' @param snps SNP definition, if not given, default to SNP probes
#' @param mft sesame-compatible manifest if non-standard
#' @return NULL
#' @examples
#'
#' temp_out <- tempfile("test")
#' deidentify(system.file(
#'     "extdata", "4207113116_A_Grn.idat", package = "sesameData"),
#'      temp_out)
#' unlink(temp_out)
#' @export
deidentify <- function(path, out_path=NULL, snps=NULL, mft=NULL) {

    res <- illuminaio::readIDAT(path)
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
    dt[snpsIdx] <- 0
    
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
