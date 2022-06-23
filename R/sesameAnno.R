

#' download Infinium manifest from the associated Github repository
#'
#' https://github.com/zhou-lab/InfiniumAnnotationV1
#' 
#' Please check the repo itself for what is available.
#' See also
#' http://zwdzwd.github.io/InfiniumAnnotation
#'
#' Unless return_path = TRUE, This function calls import function
#' depending on the resource name suffix. If the url ends with
#' .rds, it will use readRDS. If the url ends with .tsv.gz it will
#' use read_tsv. For all other cases, the function will return the
#' cached file name.
#'
#' This function replaces sesameAnno_getManifestDF.
#'
#' @param title the title of the resource
#' @param return_path return cached file path
#' @param version release version, default is the latest
#' @return tibble
#' @import BiocFileCache
#' @importFrom readr read_tsv
#' @importFrom readr cols
#' @importFrom readr col_integer
#' @importFrom readr col_character
#' @importFrom GenomeInfoDb Seqinfo
#' @examples
#' 
#' ## avoid testing since it depends on external host
#' if (FALSE) {
#' mapping <- sesameAnno_get("Mammal40/hg38.tsv.gz")
#' annoI <- sesameAnno_get("Anno/EPIC/EPIC.hg19.typeI_overlap_b151.rds")
#' mft <- sesameAnno_get("Anno/MM285/MM285.mm10.manifest.tsv.gz")
#' }
#' @export
sesameAnno_get <- function(title,
    return_path = FALSE,
    version = 1) {
    
    rpath <- paste(
        "https://github.com/zhou-lab",
        sprintf("InfiniumAnnotationV%d", version),
        "raw/main", title, sep="/")
    
    path <- bfcrpath(BiocFileCache(), rpath)
    if (return_path) { return(path); }
    if (endsWith(path, ".tsv.gz")) {
        read_tsv(path, col_types=cols(
            CpG_beg=col_integer(),
            CpG_end=col_integer(),
            address_A=col_integer(),
            address_B=col_integer(),
            .default=col_character()))
    } else if (endsWith(path, ".rds")) {
        readRDS(path)
    } else {
        return(path)
    }
}

#' Download additional annotation files
#'
#' see also
#' http://zwdzwd.github.io/InfiniumAnnotation
#'
#' This function acts similarly as sesameAnno_get except that it directly
#' download files without invoking BiocFileCache. This is needed in some
#' situation because BiocFileCache may change the file name and downstream
#' program may depend on the correct file names. It also lets you download
#' files in a cleaner way without routing through BiocFileCache
#'
#' @param title title of the annotation file
#' @param dest_dir download to this directory
#' @param version version number
#' @return annotation file
#' @examples
#'
#' ## avoid testing as this function uses external host
#' if (FALSE) {
#' sesameAnno_download("Test/3999492009_R01C01_Grn.idat", tempdir())
#' }
#' 
#' @export
sesameAnno_download <- function(
    title, dest_dir, version = 1) {

    stopifnot(!is.null(dest_dir))
    dest_file <- sprintf("%s/%s", dest_dir, title)
    dir.create(sprintf("%s/%s", dest_dir, dirname(title)),
        recursive = TRUE, showWarnings = FALSE)

    rpath <- paste(
        "https://github.com/zhou-lab",
        sprintf("InfiniumAnnotationV%d", version),
        "raw/main",
        title,
        sep="/")
    
    stopifnot(valid_url(rpath))
    download.file(rpath, dest_file, mode="wb")
    stopifnot(file.exists(dest_file) && file.info(dest_file)$size > 0)
    
    invisible(list(
        url = rpath,
        dest_dir = dest_dir,
        dest_file = dest_file,
        file_name = title))
}

guess_chrmorder <- function(chrms) {
    chrms1 <- chrms[!(chrms %in% c("chrX","chrY","chrM"))]
    paste0("chr",c(as.character(seq_len(max(as.integer(str_replace(
        sort(unique(chrms1)), "chr", "")), na.rm=TRUE))), c("X","Y","M")))
}

buildManifestGRanges <- function(
    platform, genome = NULL, version = 1,
    decoy = FALSE, columns = NULL) {

    platform <- sesameData_check_platform(platform)
    genome <- sesameData_check_genome(genome, platform)

    df <- sesameAnno_get(
        sprintf("%s/%s.tsv.gz", platform, genome), version=version)
    chrms <- df$CpG_chrm
    chrms <- chrms[!is.na(chrms)]
    if (genome %in% c("mm10","mm39","hg19","hg38")) {
        if (decoy) {
            chrms <- c(
                guess_chrmorder(chrms[!grepl("_", chrms)]),
                sort(unique(chrms[grepl("_", chrms)])))
        } else {
            chrms <- guess_chrmorder(chrms[!grepl("_", chrms)])
        }
    } else {
        chrms <- sort(unique(chrms))
    }
    df <- df[!is.na(df$CpG_chrm) & !is.na(df$CpG_beg) & !is.na(df$CpG_end),]
    df <- df[df$CpG_chrm %in% chrms,]
    gr <- GRanges(df$CpG_chrm,
        IRanges::IRanges(df$CpG_beg+1, df$CpG_end),
        strand = ifelse(df$mapFlag_A=="0", "+", "-"),
        seqinfo = Seqinfo(chrms))
    if (length(columns) > 0) {
        SummarizedExperiment::mcols(gr) <- df[,columns] }
    names(gr) <- df$Probe_ID
    sort(gr, ignore.strand = TRUE)
}

#' retrieve additional annotation files
#' @param title title of the annotation file
#' @param version version number
#' @param dest_dir if not NULL, download to this directory
#' @return annotation file
#' @examples
#' cat("Deprecated!")
#' 
#' @export
sesameData_getAnno <- function(
    title, version = 1, dest_dir = NULL) {
    .Deprecated("sesameAnno_get")
}

valid_url <- function(url_in,t=2){
    con <- url(url_in)
    check <- suppressWarnings(try(open.connection(
        con,open="rt",timeout=t),silent=TRUE)[1])
    suppressWarnings(try(close.connection(con),silent=TRUE))
    ifelse(is.null(check),TRUE,FALSE)
}
