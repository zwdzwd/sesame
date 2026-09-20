
guess_chrmorder <- function(chrms) {
    chrms1 <- chrms[!(chrms %in% c("chrX","chrY","chrM"))]
    paste0("chr",c(as.character(seq_len(max(as.integer(str_replace(
        sort(unique(chrms1)), "chr", "")), na.rm=TRUE))), c("X","Y","M")))
}

#' Build manifest GRanges from tsv
#'
#' manifest tsv files can be downloaded from
#' http://zwdzwd.github.io/InfiniumAnnotation
#'
#' @param tsv a file path, a platform (e.g., EPIC), or
#' a tibble/data.frame object
#' @param genome a genome string, e.g., hg38, mm10
#' @param decoy consider decoy sequence in chromosome order
#' @param columns the columns to include in the GRanges
#' @importFrom SummarizedExperiment metadata<-
#' @importFrom BiocFileCache BiocFileCache
#' @importFrom BiocFileCache bfcrpath
#' @importFrom readr read_tsv
#' @importFrom readr cols
#' @importFrom readr col_integer
#' @importFrom readr col_character
#' @importFrom Seqinfo Seqinfo
#' @return GRanges
#' @examples
#' \dontrun{
#' tsv = sesameAnno_download("HM450.hg38.manifest.tsv.gz")
#' gr <- sesameAnno_buildManifestGRanges(tsv)
#' ## direct access
#' gr <- sesameAnno_buildManifestGRanges("HM450.hg38.manifest")
#' }
#' @export
sesameAnno_buildManifestGRanges <- function(
    tsv, genome = NULL, decoy = FALSE, columns = NULL) {

    if (is.character(tsv)) {
        tsv <- sesameAnno_readManifestTSV(tsv)
    }
    
    chrms <- tsv$CpG_chrm
    chrms <- chrms[!is.na(chrms)]
    if (!is.null(genome) && genome %in% c("mm10","mm39","hg19","hg38")) {
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

create_default_mask <- function(df) {
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

valid_url <- function(url_in,t=2){
    con <- url(url_in)
    check <- suppressWarnings(try(open.connection(
        con,open="rt",timeout=t),silent=TRUE)[1])
    suppressWarnings(try(close.connection(con),silent=TRUE))
    ifelse(is.null(check),TRUE,FALSE)
}

#' Read manifest file to a tsv format
#' 
#' @param tsv_fn tsv file path
#' @return a manifest as a tibble
#' @examples
#' \dontrun{
#' tsv = sesameAnno_download("HM450.hg38.manifest.tsv.gz")
#' mft <- sesameAnno_readManifestTSV(tsv)
#' ## direct access
#' mft <- sesameAnno_readManifestTSV("HM450.hg38.manifest")
#' }
#' @export
sesameAnno_readManifestTSV <- function(tsv_fn) {

    if (is.character(tsv_fn) && !file.exists(tsv_fn)) {
        ## if file doesn't exist, try remote
        tsv_fn <- expand_url(tsv_fn)
        if (!valid_url(tsv_fn)) {
            stop(sprintf("File %s cannot be found.", tsv_fn))
        }
        return(sesameAnno_readManifestTSV(gzcon(url(tsv_fn))))
    }

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
#' @examples
#' \dontrun{
#' tsv = sesameAnno_download("HM450.hg38.manifest.tsv.gz")
#' addr <- sesameAnno_buildAddressFile(tsv)
#' }
#' @export
sesameAnno_buildAddressFile <- function(tsv) {
    if (is.character(tsv)) {
        tsv <- sesameAnno_readManifestTSV(tsv)
    }
    ordering <- data.frame(
        Probe_ID = tsv$Probe_ID,
        M=tsv$address_B, U=tsv$address_A,
        col=factor(tsv$channel, levels=c("G","R")), mask=FALSE)
    ordering$mask <- create_default_mask(tsv)$ref_issue

    message(sprintf("%d probes masked", sum(ordering$mask)))
    message(sprintf("%d probes/rows in ordering", nrow(ordering)))
    message(sprintf("%d probes masked", sum(ordering$mask)))
    message(sprintf("%d red probes", sum(na.omit(ordering$col=="R"))))
    message(sprintf("%d grn probes", sum(na.omit(ordering$col=="G"))))
    ordering
}


#' Build heatmap annotations
#'
#' Builds probe-level CpG-context and Ensembl Regulatory Build annotations
#' for regional methylation heatmaps.
#'
#' CpG context is classified as CpG Island, CpG Shore, CpG Shelf, or
#' Open Sea using the UCSC cpgIslandExt table. Ensembl regulatory feature
#' types are kept intact and assigned to annotation tracks such that
#' different feature types share a track only when their intervals do not
#' overlap.
#'
#' Annotation calculations are restricted to primary chromosomes and,
#' where possible, to genomic features relevant to probes represented on
#' the requested array.
#'
#' @param platform array platform, e.g. MM285 or EPIC
#' @param genome genome build, e.g. hg38 or mm39
#' @param annotations whether to build heatmap annotations
#' @param cpg whether to build CpG annotations
#' @param regulatory whether to build Ensembl regulatory annotations
#' @param shore.width CpG shore width in base pairs
#' @param shelf.width CpG shelf width in base pairs
#' @param biomart.mirror Ensembl BioMart mirror
#' @return a data.frame containing probe-level heatmap annotations
#' @examples
#' \dontrun{
#' anno <- sesameAnno_buildHeatmapAnnotations(
#'     platform="MM285", genome="mm39")
#' }
#' @export
sesameAnno_buildHeatmapAnnotations <- function(
    platform, genome=NULL, annotations=TRUE, cpg=TRUE,
    regulatory=TRUE, shore.width=2000, shelf.width=2000,
    biomart.mirror="www") {
  
  ## Return immediately when annotations are not requested.
  if (!annotations) return(NULL)
  
  ## Require the genome explicitly so all annotation sources use the
  ## same genomic assembly.
  if (is.null(genome))
    stop("genome must be specified.")
  
  message(
    "Building heatmap annotations for ",
    platform, " (", genome, ")...")
  
  ## Define primary chromosomes. Alternate loci, patches, random contigs,
  ## and unplaced sequences are excluded from annotation calculations.
  main.chr <- switch(
    genome,
    mm10=paste0("chr", c(1:19, "X", "Y", "M")),
    mm39=paste0("chr", c(1:19, "X", "Y", "M")),
    hg19=paste0("chr", c(1:22, "X", "Y", "M")),
    hg38=paste0("chr", c(1:22, "X", "Y", "M")),
    stop("Unsupported genome: ", genome))
  
  ## Load the requested array manifest.
  message("Loading ", platform, " manifest...")
  
  tsv <- sesameAnno_readManifestTSV(
    sprintf("%s.%s.manifest", platform, genome))
  
  req <- c("Probe_ID", "CpG_chrm", "CpG_beg", "CpG_end")
  if (length(x <- setdiff(req, names(tsv))))
    stop("Missing columns: ", paste(x, collapse=", "))
  
  ## ---------------------------------------------------------------------
  ## Construct genomic ranges for array probes
  ## ---------------------------------------------------------------------
  
  message("Preparing array probe coordinates...")
  
  ## Only probes with complete coordinates on primary chromosomes
  ## participate in genomic overlap calculations.
  mapped <- complete.cases(
    tsv[,c("CpG_chrm", "CpG_beg", "CpG_end")]) &
    tsv$CpG_chrm %in% main.chr
  
  ## Manifest CpG coordinates use a 0-based beginning whereas GRanges
  ## uses 1-based coordinates, so CpG_beg is incremented by one.
  probes <- GenomicRanges::GRanges(
    tsv$CpG_chrm[mapped],
    IRanges::IRanges(
      tsv$CpG_beg[mapped]+1,
      tsv$CpG_end[mapped]))
  names(probes) <- tsv$Probe_ID[mapped]
  
  ## Preserve the original manifest rows corresponding to mapped probes.
  ## This permits overlap results to be written back into a table
  ## containing every manifest probe.
  probe.rows <- which(mapped)
  
  anno <- data.frame(
    Probe_ID=tsv$Probe_ID,
    stringsAsFactors=FALSE)
  
  message(
    length(probes), " of ", nrow(tsv),
    " probes retained with primary-chromosome coordinates.")
  
  ## Convert genomic ranges into human-readable chr:start-end labels.
  range_label <- function(gr) {
    paste0(
      GenomicRanges::seqnames(gr), ":",
      IRanges::start(gr), "-",
      IRanges::end(gr))
  }
  
  ## =====================================================================
  ## CpG Island / Shore / Shelf / Open Sea annotations
  ## =====================================================================
  
  if (cpg) {
    
    message("Building CpG-context annotations...")
    
    ## UCSC provides a precomputed cpgIslandExt table for each assembly.
    url <- sprintf(
      paste0(
        "https://hgdownload.soe.ucsc.edu/goldenPath/",
        "%s/database/cpgIslandExt.txt.gz"),
      genome)
    
    message("Downloading UCSC CpG island annotations...")
    
    ## The UCSC file is needed only while constructing the annotation,
    ## so it is downloaded to a temporary file.
    fn <- tempfile(fileext=".txt.gz")
    download.file(url, fn, mode="wb", quiet=TRUE)
    
    message("Reading UCSC CpG island annotations...")
    
    x <- read.table(
      gzfile(fn),
      sep="\t",
      quote="",
      comment.char="",
      stringsAsFactors=FALSE)
    
    ## Remove alternate/random/unplaced chromosomes before constructing
    ## genomic ranges.
    x <- x[x[[2]] %in% main.chr,,drop=FALSE]
    
    message(
      nrow(x),
      " CpG islands retained on primary chromosomes.")
    
    ## UCSC chromStart is 0-based. Adding one converts it to the
    ## corresponding 1-based GRanges interval.
    island <- GenomicRanges::GRanges(
      x[[2]],
      IRanges::IRanges(
        x[[3]]+1,
        x[[4]]))
    
    ## -----------------------------------------------------------------
    ## Construct CpG shores
    ## -----------------------------------------------------------------
    
    message(
      "Constructing CpG shores (",
      shore.width, " bp)...")
    
    ## Shores immediately flank CpG islands on both sides.
    left.shore <- IRanges::flank(
      island, shore.width,
      start=TRUE,
      ignore.strand=TRUE)
    
    right.shore <- IRanges::flank(
      island, shore.width,
      start=FALSE,
      ignore.strand=TRUE)
    
    ## Merge overlapping shore intervals generated by nearby islands.
    shore <- IRanges::reduce(
      c(left.shore, right.shore),
      ignore.strand=TRUE)
    
    ## Remove any island sequence from the shore intervals, enforcing:
    ##
    ##     CpG Island > CpG Shore
    ##
    shore <- IRanges::setdiff(
      shore,
      IRanges::reduce(
        island,
        ignore.strand=TRUE))
    
    ## -----------------------------------------------------------------
    ## Construct CpG shelves
    ## -----------------------------------------------------------------
    
    message(
      "Constructing CpG shelves (",
      shelf.width, " bp beyond shores)...")
    
    ## Shelves extend outward from the outer boundaries of the shores.
    left.shelf <- IRanges::flank(
      left.shore, shelf.width,
      start=TRUE,
      ignore.strand=TRUE)
    
    right.shelf <- IRanges::flank(
      right.shore, shelf.width,
      start=FALSE,
      ignore.strand=TRUE)
    
    shelf <- IRanges::reduce(
      c(left.shelf, right.shelf),
      ignore.strand=TRUE)
    
    ## Remove any sequence already classified as an island or shore,
    ## producing mutually exclusive CpG-context categories:
    ##
    ##     Island > Shore > Shelf > Open Sea
    ##
    shelf <- IRanges::setdiff(
      shelf,
      IRanges::reduce(
        c(island, shore),
        ignore.strand=TRUE))
    
    ## -----------------------------------------------------------------
    ## Assign CpG context to probes
    ## -----------------------------------------------------------------
    
    message("Mapping CpG-context annotations to array probes...")
    
    ## Every probe initially receives Open Sea. Shelf, Shore, and Island
    ## assignments subsequently overwrite this value in priority order.
    anno$CpG_Location_Type <- "Open Sea"
    anno$CpG_Location_Type_Genomic_Coordinates <- ""
    
    ## Assignment order is intentionally lowest-to-highest priority.
    cpg.features <- list(
      "CpG Shelf"=shelf,
      "CpG Shore"=shore,
      "CpG Island"=island)
    
    for (type in names(cpg.features)) {
      message("  Mapping ", type, "...")
      
      z <- cpg.features[[type]]
      
      ## Determine which probes overlap genomic intervals belonging
      ## to the current CpG-context category.
      h <- IRanges::findOverlaps(
        probes, z,
        ignore.strand=TRUE)
      
      if (!length(h)) {
        message("    No overlapping probes.")
        next
      }
      
      ## Translate mapped-probe indices back to complete manifest rows.
      q <- probe.rows[S4Vectors::queryHits(h)]
      
      anno$CpG_Location_Type[q] <- type
      anno$CpG_Location_Type_Genomic_Coordinates[q] <-
        range_label(z[S4Vectors::subjectHits(h)])
      
      message(
        "    ", length(unique(q)),
        " probe(s) overlap ", type, ".")
    }
    
    message("CpG-context annotations complete.")
  }
  
  ## =====================================================================
  ## Ensembl Regulatory Build annotations
  ## =====================================================================
  
  if (regulatory) {
    
    message("Building Ensembl Regulatory Build annotations...")
    
    ## Select the Ensembl regulatory dataset corresponding to the
    ## requested genome.
    species <- switch(
      genome,
      hg38="hsapiens",
      mm39="mmusculus",
      stop(
        "Current Ensembl Regulatory Build unavailable for ",
        genome))
    
    dataset <- paste0(
      species,
      "_regulatory_feature")
    
    message(
      "Connecting to Ensembl BioMart dataset ",
      dataset, "...")
    
    ## Connect to the Ensembl Functional Genomics BioMart. biomaRt may
    ## automatically try another mirror when the requested site is
    ## temporarily unavailable.
    mart <- biomaRt::useEnsembl(
      biomart="ENSEMBL_MART_FUNCGEN",
      dataset=dataset,
      mirror=biomart.mirror)
    
    ## -----------------------------------------------------------------
    ## Verify genome assembly
    ## -----------------------------------------------------------------
    
    message("Verifying Ensembl genome assembly...")
    
    ## Regulatory annotations and array coordinates must refer to the
    ## same assembly before genomic overlaps are calculated.
    ds <- biomaRt::listDatasets(mart)
    j <- match(dataset, ds$dataset)
    
    expected <- switch(
      genome,
      hg38="GRCh38",
      mm39="GRCm39")
    
    if (is.na(j) || !startsWith(ds$version[j], expected))
      stop("Ensembl Regulatory Build assembly mismatch.")
    
    message(
      "Using Ensembl assembly: ",
      ds$version[j])
    
    ## -----------------------------------------------------------------
    ## Download regulatory features
    ## -----------------------------------------------------------------
    
    message("Downloading Ensembl regulatory features...")
    
    ## Retrieve genomic coordinates, current Ensembl feature labels,
    ## and regulatory stable IDs. Feature labels are not hard-coded,
    ## allowing future Ensembl feature categories to be incorporated.
    reg <- biomaRt::getBM(
      attributes=c(
        "chromosome_name",
        "chromosome_start",
        "chromosome_end",
        "feature_type_name",
        "regulatory_stable_id"),
      mart=mart)
    
    n.reg.total <- nrow(reg)
    
    message(
      n.reg.total,
      " Ensembl regulatory features downloaded.")
    
    ## Normalize Ensembl chromosome names to the chr-prefixed naming
    ## convention used by the array manifest.
    chr <- ifelse(
      grepl("^chr", reg$chromosome_name),
      reg$chromosome_name,
      paste0("chr", reg$chromosome_name))
    
    ## -----------------------------------------------------------------
    ## Restrict to primary chromosomes
    ## -----------------------------------------------------------------
    
    message("Filtering regulatory features to primary chromosomes...")
    
    keep <- chr %in% main.chr
    reg <- reg[keep,,drop=FALSE]
    chr <- chr[keep]
    
    reg.gr <- GenomicRanges::GRanges(
      chr,
      IRanges::IRanges(
        reg$chromosome_start,
        reg$chromosome_end))
    
    ## Store Ensembl feature category and stable ID directly on each
    ## regulatory genomic interval.
    S4Vectors::mcols(reg.gr)$type <- reg$feature_type_name
    S4Vectors::mcols(reg.gr)$id <- reg$regulatory_stable_id
    
    n.reg.primary <- length(reg.gr)
    
    message(
      n.reg.primary,
      " regulatory features retained on primary chromosomes.")
    
    ## -----------------------------------------------------------------
    ## Restrict to features relevant to this array
    ## -----------------------------------------------------------------
    
    message(
      "Filtering regulatory features to those overlapping ",
      platform, " probes...")
    
    ## Features overlapping no array probe can never appear in the
    ## probe-level annotation table. Removing them before category-level
    ## overlap calculations substantially reduces unnecessary work.
    reg.gr <- IRanges::subsetByOverlaps(
      reg.gr, probes,
      ignore.strand=TRUE)
    
    n.reg.probe <- length(reg.gr)
    
    message(
      n.reg.probe,
      " regulatory features overlap ",
      platform, " probes.")
    
    ## -----------------------------------------------------------------
    ## Remove exact duplicate regulatory features
    ## -----------------------------------------------------------------
    
    message("Removing duplicate regulatory features...")
    
    ## Chromosome + start + end + feature type defines a duplicate here.
    ## Different feature types occupying identical coordinates are kept
    ## because they represent biologically distinct annotations.
    key <- paste(
      GenomicRanges::seqnames(reg.gr),
      IRanges::start(reg.gr),
      IRanges::end(reg.gr),
      S4Vectors::mcols(reg.gr)$type,
      sep=":")
    
    reg.gr <- reg.gr[!duplicated(key)]
    
    message(
      length(reg.gr),
      " unique probe-relevant regulatory features retained.")
    
    ## -----------------------------------------------------------------
    ## Identify regulatory feature categories
    ## -----------------------------------------------------------------
    
    ## Regulatory tracks are constructed using complete feature TYPES,
    ## not individual intervals. Every interval belonging to a category
    ## such as Enhancer or Promoter therefore remains in one track.
    feature.types <- sort(
      unique(
        S4Vectors::mcols(reg.gr)$type))
    
    features <- lapply(
      feature.types,
      function(type)
        reg.gr[
          S4Vectors::mcols(reg.gr)$type == type])
    
    names(features) <- feature.types
    
    message(
      length(feature.types),
      " regulatory feature type(s) identified: ",
      paste(feature.types, collapse=", "))
    
    ## -----------------------------------------------------------------
    ## Determine overlaps between complete feature categories
    ## -----------------------------------------------------------------
    
    message(
      "Calculating overlaps between regulatory feature types...")
    
    ## type.overlap[i,j] is TRUE when at least one interval belonging to
    ## feature type i overlaps at least one interval belonging to type j.
    ##
    ## Within-category overlaps do not matter here because an entire
    ## category must remain together in a single annotation track.
    n.types <- length(feature.types)
    
    type.overlap <- matrix(
      FALSE,
      nrow=n.types,
      ncol=n.types,
      dimnames=list(
        feature.types,
        feature.types))
    
    if (n.types > 1) {
      for (i in seq_len(n.types-1)) {
        for (j in (i+1):n.types) {
          
          message(
            "  Checking ",
            feature.types[i], " vs ",
            feature.types[j], "...")
          
          overlap <- any(
            IRanges::overlapsAny(
              features[[i]],
              features[[j]],
              ignore.strand=TRUE))
          
          type.overlap[i,j] <- overlap
          type.overlap[j,i] <- overlap
        }
      }
    }
    
    ## -----------------------------------------------------------------
    ## Assign complete feature types to heatmap tracks
    ## -----------------------------------------------------------------
    
    message(
      "Assigning regulatory feature types to annotation tracks...")
    
    ## Each feature type is placed into the first existing track that
    ## contains no category with which it overlaps. If no compatible
    ## track exists, a new Regulatory_Feature_n track is created.
    type.tracks <- list()
    
    for (type in feature.types) {
      
      available <- which(vapply(
        type.tracks,
        function(existing.types)
          !any(
            type.overlap[
              type,
              existing.types]),
        logical(1)))
      
      if (!length(available)) {
        type.tracks[[length(type.tracks)+1]] <- type
      } else {
        i <- available[1]
        type.tracks[[i]] <- c(
          type.tracks[[i]],
          type)
      }
    }
    
    ## Verify that every regulatory feature type occurs exactly once
    ## across the dynamically generated tracks.
    stopifnot(
      identical(
        sort(unlist(type.tracks)),
        sort(feature.types)))
    
    ## Verify that categories sharing the same annotation track do not
    ## overlap one another anywhere in the probe-relevant intervals.
    for (track.types in type.tracks) {
      if (length(track.types) > 1) {
        stopifnot(
          !any(
            type.overlap[
              track.types,
              track.types,
              drop=FALSE]))
      }
    }
    
    message(
      length(feature.types),
      " feature type(s) assigned to ",
      length(type.tracks),
      " annotation track(s).")
    
    for (i in seq_along(type.tracks)) {
      message(
        "  Regulatory_Feature_", i, ": ",
        paste(type.tracks[[i]], collapse=", "))
    }
    
    ## -----------------------------------------------------------------
    ## Map regulatory tracks onto array probes
    ## -----------------------------------------------------------------
    
    message("Mapping regulatory annotation tracks to array probes...")
    
    for (i in seq_along(type.tracks)) {
      
      type.col <- paste0(
        "Regulatory_Feature_", i)
      
      coord.col <- paste0(
        type.col,
        "_Genomic_Coordinates")
      
      anno[[type.col]] <- ""
      anno[[coord.col]] <- ""
      
      ## Collect all regulatory intervals belonging to the complete
      ## feature categories assigned to this heatmap track.
      track.types <- type.tracks[[i]]
      
      message(
        "  Mapping Regulatory_Feature_", i,
        " (", paste(track.types, collapse=", "), ")...")
      
      track <- reg.gr[
        S4Vectors::mcols(reg.gr)$type %in%
          track.types]
      
      ## Determine which probes overlap any regulatory feature
      ## represented in this annotation track.
      h <- IRanges::findOverlaps(
        probes, track,
        ignore.strand=TRUE)
      
      if (!length(h)) {
        message("    No overlapping probes.")
        next
      }
      
      qh <- S4Vectors::queryHits(h)
      sh <- S4Vectors::subjectHits(h)
      
      ## Group all regulatory-feature hits belonging to each probe.
      ## A probe can overlap more than one regulatory interval.
      by.probe <- split(sh, qh)
      
      ## Translate mapped-probe indices back into complete manifest
      ## annotation-table row numbers.
      out.rows <- probe.rows[
        as.integer(
          names(by.probe))]
      
      ## Store the unique Ensembl regulatory category or categories
      ## overlapping each probe.
      anno[[type.col]][out.rows] <- vapply(
        by.probe,
        function(k)
          paste(
            unique(
              S4Vectors::mcols(track)$type[k]),
            collapse=";"),
        character(1))
      
      ## Store the corresponding genomic regulatory intervals.
      track.coords <- range_label(track)
      
      anno[[coord.col]][out.rows] <- vapply(
        by.probe,
        function(k)
          paste(
            unique(
              track.coords[k]),
            collapse=";"),
        character(1))
      
      message(
        "    ", length(by.probe),
        " probe(s) annotated.")
    }
    
    message("Ensembl Regulatory Build annotations complete.")
    
    ## -----------------------------------------------------------------
    ## Preserve regulatory-build metadata
    ## -----------------------------------------------------------------
    
    ## Store information describing how the regulatory annotation was
    ## constructed without duplicating it across every probe row.
    attr(anno, "regulatory_tracks") <-
      length(type.tracks)
    
    attr(anno, "regulatory_track_types") <-
      type.tracks
    
    attr(anno, "regulatory_feature_types") <-
      feature.types
    
    attr(anno, "regulatory_type_overlap") <-
      type.overlap
    
    attr(anno, "ensembl_assembly") <-
      ds$version[j]
    
    attr(anno, "regulatory_features_downloaded") <-
      n.reg.total
    
    attr(anno, "regulatory_features_primary") <-
      n.reg.primary
    
    attr(anno, "regulatory_features_probe_relevant") <-
      n.reg.probe
  }
  
  ## Record the genome assembly used for the complete annotation table.
  attr(anno, "genome") <- genome
  
  message(
    "Heatmap annotations complete for ",
    platform, " (", genome, ").")
  
  anno
}

#' Annotate a data.frame using manifest
#' 
#' Annotation source: https://zwdzwd.github.io/InfiniumAnnotation
#' e.g., EPICv2.hg38.manifest
#' 
#' @param df input data frame with Probe_ID as a column
#' @param probe_id the Probe_ID column name, default to "Probe_ID" or
#' rownames
#' @param platform which array platform, guess from probe ID if not given
#' @param genome the genome build, use default if not given
#' @return a new data.frame with manifest attached
#' @examples
#' \dontrun{
#' df <- data.frame(Probe_ID = c("cg00101675_BC21", "cg00116289_BC21"))
#' sesameAnno_attachManifest(df)
#' }
#' @export
sesameAnno_attachManifest <- function(
    df, probe_id="Probe_ID", platform=NULL, genome=NULL) {

    if (is.numeric(df)) {
        if (is.matrix(df)) {
            df <- cbind(Probe_ID=rownames(df), as.data.frame(df))
        } else {
            df <- data.frame(Probe_ID = names(df), beta = df)
        }
    }
    stopifnot(is(df, "data.frame"))
    stopifnot(probe_id %in% colnames(df))

    if (is.null(platform)) {
        platform <- inferPlatformFromProbeIDs(df[[probe_id]]) }

    genome <- sesameData_check_genome(genome, platform)

    mft <- sesameAnno_readManifestTSV(
        sprintf("%s.%s.manifest", platform, genome))
    if (platform %in% c("HM27","HM450")) {
        mft_probeid <- "probeID"
    } else {
        mft_probeid <- "Probe_ID"
    }
    cbind(df, as.data.frame(mft)[match(df[[probe_id]], mft[[mft_probeid]]),])
}

expand_url <- function(url,
    base = "https://github.com/zhou-lab/InfiniumAnnotationData/raw/main/") {
    ## input can be :
    ## EPIC.hg38.manifest.tsv.gz
    ## /Test/3999492009_R01C01_Grn.idat
    ## https://github.com/zhou-lab/InfiniumAnnotationData/raw/main/Anno/EPIC/EPIC.hg19.manifest.tsv.gz
    if (!any(endsWith(url, c("rds","tsv.gz","mask.cm")))) {
        url <- sprintf("%s.tsv.gz", url)
    }
    if (!grepl("http", url)) {
        if (!grepl("/", url)) {
            platform <- strsplit(url, "\\.")[[1]][1]
            ## Lean, versioned release files (ordering / coord / mask.cm / snp)
            ## live in zhou-lab/InfiniumAnnotation/<platform>/; the larger
            ## annotation tables (manifest, gene, ...) in
            ## InfiniumAnnotationData/Anno/<platform>/.
            if (grepl("\\.(ordering|coord|snp)\\.|\\.mask\\.cm$", url)) {
                base <- "https://github.com/zhou-lab/InfiniumAnnotation/raw/main/"
                url <- sprintf("%s/%s", platform, url)
            } else {
                url <- sprintf("Anno/%s/%s", platform, url)
            }
        }
        url <- sprintf("%s/%s", base, url)
    }
    url
}

#' Download SeSAMe annotation files
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
#' @param url url or title of the annotation file
#' @param destfile download to this file, a temp file if unspecified
#' @return the path to downloaded file
#' @importFrom utils download.file
#' @examples
#'
#' \dontrun{
#' ## avoid testing as this function uses external host
#' sesameAnno_download("Test/3999492009_R01C01_Grn.idat")
#' sesameAnno_download("EPIC.hg38.manifest.tsv.gz")
#' sesameAnno_download("EPIC.hg38.snp.tsv.gz")
#' }
#' 
#' @export
sesameAnno_download <- function(
    url, destfile = tempfile(basename(url))) {
    url <- expand_url(url)
    download.file(url, destfile=destfile)
    destfile
}
