

#' download Infinium manifest from the associated Github repository
#'
#' Since most of the annotation is not essential to sesame functioning,
#' sesameData package no longer host the full manifest. This is
#' the command to use to retrieve the full manifest and other annotation
#' from the following Github host:
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
        df <- read_tsv(path)
        for (col in c("CpG_beg","CpG_end","address_A","address_B")) {
            if (col %in% colnames(df)) {
                df[[col]] <- as.integer(df[[col]])
            }
        }
        return(df)
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

#' Build manifest GRanges from tsv
#'
#' @param tsv a file path, a platform (e.g., EPIC), or
#' a tibble/data.frame object
#' @param genome a genome string, e.g., hg38, mm10
#' @param decoy consider decoy sequence in chromosome order
#' @param columns the columns to include in the GRanges
#' @importFrom SummarizedExperiment metadata<-
#' @return GRanges
#'
#' ## genome <- sesameAnno_buildManifestGRanges("HM450", "hg38")
#' ## or
#' ## download tsv directly from
#' ## http://zwdzwd.github.io/InfiniumAnnotation or
#' ## genome <- sesameAnno_buildManifestGRanges("downloaded_file")
#' @export
sesameAnno_buildManifestGRanges <- function(
    tsv, genome, decoy = FALSE, columns = NULL) {

    if (is.character(tsv)) { # file path
        if (file.exists(tsv)) {
            tsv <- sesameAnno_readManifestTSV(tsv)
        } else {
            platform <- tsv
            if (is.null(genome)) {
                genome <- sesameData_check_genome(NULL, platform)
            }
            tsv <- sesameAnno_get(sprintf("%s/%s.tsv.gz", platform, genome))
        }
    }

    for (col in c("CpG_beg","CpG_end","address_A","address_B")) {
        if (col %in% colnames(tsv)) {
            tsv[[col]] <- as.integer(tsv[[col]]) }}
    
    chrms <- tsv$CpG_chrm
    chrms <- chrms[!is.na(chrms)]
    if (genome %in% c("mm10","mm39","hg19","hg38")) {
        if (decoy) { chrms <- c(
            guess_chrmorder(chrms[!grepl("_", chrms)]),
            sort(unique(chrms[grepl("_", chrms)])))
        } else {
            chrms <- guess_chrmorder(chrms[!grepl("_", chrms)])
        }
    } else {
        chrms <- sort(unique(chrms))
    }
    chrms <- c(chrms, "*")
    idx <- is.na(tsv$CpG_chrm) | !(tsv$CpG_chrm %in% chrms)
    tsv$CpG_chrm[idx] <- "*"
    tsv$CpG_beg[idx] <- -1
    tsv$CpG_end[idx] <- 0
    
    gr <- GRanges(tsv$CpG_chrm, IRanges::IRanges(tsv$CpG_beg+1, tsv$CpG_end),
        strand = ifelse(
            is.na(tsv$mapFlag_A), "*", ifelse(tsv$mapFlag_A=="0", "+", "-")),
        seqinfo = Seqinfo(chrms))
    if (length(columns) > 0) {
        SummarizedExperiment::mcols(gr) <- tsv[,columns] }
    names(gr) <- tsv$Probe_ID
    metadata(gr)[["genome"]] <- genome
    message(sprintf("%d probes in GRanges.", length(gr)))
    message(sprintf("%d probes belong to chr*.", sum(seqnames(gr)=="*")))
    message(sprintf("%d probes on decoy chr.", sum(grepl("_", seqnames(gr)))))
    sort(gr, ignore.strand = TRUE)
}

create_mask <- function(df) {
    unmapped <- (
        is.na(df$mapAS_A) | df$mapAS_A < 35 |
        (!is.na(df$mapAS_B) & df$mapAS_B < 35))
    masks <- data.frame(
        Probe_ID = df$Probe_ID,
        nonunique = (
            (!unmapped) &
            (df$mapQ_A == 0 | (!is.na(df$mapQ_B) & df$mapQ_B == 0))),
        missing_target = (
            (!unmapped) &
            (is.na(df$target) | (df$target != "CG")) &
            grepl("^cg", df$Probe_ID)))
    masks$control <- grepl("^ctl", df$Probe_ID)
    masks$design_issue <- grepl("^uk", df$Probe_ID)
    masks$unmapped <- (unmapped & masks$control != 1 & masks$design_issue != 1)
    masks$low_mapq <- (
        (!is.na(df$mapQ_A)) &
        (df$mapQ_A < 30 | (!is.na(df$mapQ_B) & df$mapQ_B < 30)))
    masks$ref_issue <- (unmapped | masks$missing_target)
    masks[c("Probe_ID","unmapped","missing_target",
        "ref_issue","nonunique","low_mapq","control","design_issue")]
}

sesameAnno_readManifestTSV <- function(tsv_fn) {
    read_tsv(
        tsv_fn, col_types=cols(
            CpG_chrm = col_character(),
            CpG_beg = col_integer(),
            CpG_end = col_integer(),
            address_A = col_integer(), address_B = col_integer(),
            target = col_character(), nextBase = col_character(),
            channel = col_character(),
            Probe_ID = col_character(), mapFlag_A = col_integer(),
            mapChrm_A = col_character(),
            mapPos_A = col_integer(), mapQ_A = col_integer(),
            mapCigar_A = col_character(),
            AlleleA_ProbeSeq = col_character(),
            mapNM_A = col_character(), mapAS_A = col_integer(),
            mapYD_A = col_character(),
            mapFlag_B = col_integer(),
            mapChrm_B = col_character(), mapPos_B = col_integer(),
            mapQ_B = col_integer(), mapCigar_B = col_character(),
            AlleleB_ProbeSeq = col_character(),
            mapNM_B = col_character(), mapAS_B = col_integer(),
            mapYD_B = col_character(), type = col_character()))
}

#' Build sesame ordering address file from tsv
#'
#' @param tsv a platform name, a file path or a tibble/data.frame manifest file
#' @return a list of ordering and controls
#'
#' ## addr <- sesameAnno_buildAddressFile("HM450")
#' ## or
#' ## download tsv directly from
#' ## http://zwdzwd.github.io/InfiniumAnnotation or
#' ## addr <- sesameAnno_buildAddressFile("downloaded_file")
#' @export
sesameAnno_buildAddressFile <- function(tsv) {

    if (is.character(tsv)) { # string input
        if (startsWith(tsv, "https:")) { # url
            tsv <- sesameAnno_readManifestTSV(gzcon(url(tsv)))
        } else if (file.exists(tsv) && !dir.exists(tsv)) { # file path
            tsv <- sesameAnno_readManifestTSV(tsv)
        } else { # use some guess
            ## note there is a discrepancy in the
            ## control probe IDs, to update
            platform <- tsv
            genome <- sesameData_check_genome(NULL, platform)

            if (platform %in% c("EPICv2")) {
                tsv <- sesameAnno_get(sprintf("Anno/%s/%s.%s.manifest.tsv.gz", platform, platform, genome))
            } else {
                tsv <- sesameAnno_get(sprintf("%s/%s.tsv.gz", platform, genome))
            }
        }
    }
    
    ordering <- data.frame(
        Probe_ID = tsv$Probe_ID,
        M=tsv$address_B, U=tsv$address_A,
        col=factor(tsv$channel, levels=c("G","R")), mask=FALSE)
    ordering$mask <- create_mask(tsv)$ref_issue

    message(sprintf("%d probes masked", sum(ordering$mask)))
    message(sprintf("%d probes/rows in ordering", nrow(ordering)))
    message(sprintf("%d probes masked", sum(ordering$mask)))
    message(sprintf("%d red probes", sum(na.omit(ordering$col=="R"))))
    message(sprintf("%d grn probes", sum(na.omit(ordering$col=="G"))))
    ordering
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
