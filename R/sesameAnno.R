
anno_base_default_version <- 1

#' download Infinium manifest from Github repositories
#'
#' @param platform Mammal40, MM285, EPIC, and HM450
#' @param genome hg38, mm10 etc.
#' @param version release version, default is the latest
#' @return tibble
#' @import BiocFileCache
#' @importFrom readr read_tsv
#' @importFrom readr cols
#' @importFrom readr col_integer
#' @importFrom readr col_character
#' @examples
#' 
#' ## avoid testing since it depends on external host
#' if (FALSE) {
#' mft <- sesameAnno_getManifestDF("Mammal40")
#' }
#' @export
sesameAnno_getManifestDF <- function(platform, genome=NULL,
    version = anno_base_default_version) {
    
    platform <- sesameData_check_platform(platform)
    genome <- sesameData_check_genome(genome, platform)

    rpath <- paste(
        "https://github.com/zhou-lab",
        sprintf("InfiniumAnnotationV%d", version),
        "raw/main",
        sprintf("%s/%s.tsv.gz", platform, genome),
        sep="/")
    
    read_tsv(bfcrpath(BiocFileCache(), rpath),
        col_types=cols(
            CpG_beg=col_integer(),
            CpG_end=col_integer(),
            address_A=col_integer(),
            address_B=col_integer(),
            .default=col_character()))
}

buildManifestGRanges <- function(
    platform, genome = NULL, version = anno_base_default_version,
    decoy = FALSE, columns = NULL) {

    platform <- sesameData_check_platform(platform)
    genome <- sesameData_check_genome(genome, platform)

    df <- sesameAnno_getManifestDF(platform, genome=genome, version=version)
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
        IRanges(df$CpG_beg+1, df$CpG_end),
        strand = ifelse(df$mapFlag_A=="0", "+", "-"),
        seqinfo = Seqinfo(chrms))
    if (length(columns) > 0) {
        mcols(gr) <- df[,columns] }
    names(gr) <- df$Probe_ID
    sort(gr, ignore.strand = TRUE)
}


#' Download additional annotation files
#'
#' From the Infinium annotation website associated github repo
#' e.g., 
#' https://github.com/zhou-lab/InfiniumAnnotationV1
#'
#' The default version number should always work. One need to
#' refer to the actual repo to see which one of the other versions
#' also work.
#' 
#' See also
#' http://zwdzwd.github.io/InfiniumAnnotation
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
    title, dest_dir, version = anno_base_default_version) {

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

#' Retrieve additional annotation Rds data
#'
#' From the Infinium annotation website associated github repo
#' e.g., 
#' https://github.com/zhou-lab/InfiniumAnnotationV1
#'
#' The default version number should always work. One need to
#' refer to the actual repo to see which one of the other versions
#' also work.
#' 
#' See also
#' http://zwdzwd.github.io/InfiniumAnnotation
#'
#' @param title title of the annotation file
#' @param version version number
#' @param dest_dir if not NULL, download to this directory
#' @return annotation file
#' @examples
#'
#' ## avoided testing as this function uses external host
#' if (FALSE) {
#' annoI <- anno_get("Anno/EPIC/EPIC.hg19.typeI_overlap_b151.rds")
#' }
#' 
#' @export
sesameAnno_get <- function(
    title, version = anno_base_default_version, dest_dir = NULL) {

    rpath <- paste(
        "https://github.com/zhou-lab",
        sprintf("InfiniumAnnotationV%d", version),
        "raw/main",
        title,
        sep="/")

    if (!valid_url(rpath)) {
        message(sprintf("File not available %s.", rpath))
        return(NULL)
    }
    stopifnot(endsWith(title, ".rds"))
    message("Retrieving annotation from ",rpath,
        "... ", appendLF = FALSE)
    anno <- readRDS(url(rpath))
    message("Done.")
    anno
}

#' retrieve additional annotation files
#' @param title title of the annotation file
#' @param version version number
#' @param dest_dir if not NULL, download to this directory
#' @return annotation file
#' @examples
#' sesameData_getAnno(NULL)
#' 
#' @export
sesameData_getAnno <- function(
    title, version = anno_base_default_version, dest_dir = NULL) {
    .Deprecated("anno_get")
}

valid_url <- function(url_in,t=2){
    con <- url(url_in)
    check <- suppressWarnings(try(open.connection(
        con,open="rt",timeout=t),silent=TRUE)[1])
    suppressWarnings(try(close.connection(con),silent=TRUE))
    ifelse(is.null(check),TRUE,FALSE)
}
