#' Visualize Region
#'
#' The function takes a genomic coordinate (chromosome, start and end) and a
#' beta value matrix (probes on the row and samples on the column). It plots
#' the beta values as a heatmap for all probes falling into the genomic region.
#' If `draw=TRUE` the function returns the plotted grid graphics object.
#' Otherwise, the selected beta value matrix is returned.
#' `cluster.samples=TRUE/FALSE` controls whether hierarchical clustering is
#' applied to the subset beta value matrix.
#'
#' @param chrm chromosome
#' @param beg begin of the region
#' @param end end of the region
#' @param betas beta value matrix (row: probes, column: samples)
#' @param platform EPIC, HM450, or MM285
#' @param genome hg38, mm10, ..., will infer if not given.
#' For additional mapping, download the GRanges object from
#' http://zwdzwd.github.io/InfiniumAnnotation
#' and provide the following argument
#' ..., genome = sesameAnno_buildManifestGRanges("downloaded_file"),...
#' to this function.
#' @param draw draw figure or return betas
#' @param cluster.samples whether to cluster samples
#' @param na.rm remove probes with all NA.
#' @param nprobes.max maximum number of probes to plot
#' @param txn.types default to protein_coding, use NULL for all
#' @param txn.font.size transcript name font size
#' @param ... additional options, see assemble_plots
#' @return graphics or a matrix containing the captured beta values
#' @import grid
#' @importMethodsFrom IRanges subsetByOverlaps
#' @importFrom GenomicRanges GRanges
#' @examples
#' betas <- sesameDataGet('HM450.76.TCGA.matched')$betas
#' visualizeRegion('chr20', 44648623, 44652152, betas, 'HM450')
#' @export
visualizeRegion <- function(chrm, beg, end, betas, platform = NULL,
    genome = NULL, draw = TRUE, cluster.samples = FALSE,
    na.rm = FALSE, nprobes.max = 1000, txn.types = "protein_coding",
    txn.font.size = 6, ...) {

    if (is.null(dim(betas))) { betas <- as.matrix(betas) }
    platform <- sesameData_check_platform(platform, rownames(betas))
    genome <- sesameData_check_genome(genome, platform)

    reg <- GRanges(chrm, IRanges::IRanges(beg, end))
    
    genomeInfo <- sesameData_getGenomeInfo(genome)
    txns <- subsetByOverlaps(genomeInfo$txns, reg)

    probes <- sesameData_getManifestGRanges(platform, genome = genome)
    probes <- subsetByOverlaps(probes, reg)
    probes <- probes[names(probes) %in% rownames(betas)]
    if (na.rm) { probes <- probes[apply(
        betas[names(probes), ], 1, function(x) !all(is.na(x)))] }
    
    if (length(probes) == 0) {
        stop("No probe overlap region ", sprintf(
            '%s:%d-%d', chrm, beg, end)) }
    if (length(probes) > nprobes.max) {
        stop(sprintf('Too many probes (%d). Shrink region?', length(probes))) }

    ## plot transcripts
    plt.txns <- plotTranscripts(txns, reg, beg, end,
        txn.types = txn.types, txn.font.size = txn.font.size)
    plt.mapLines <- plotMapLines(probes, beg, end)
    plt.cytoband <- plotCytoBand(chrm, beg, end, genomeInfo)
    
    ## clustering
    betas <- betas[names(probes),,drop=FALSE]
    if (cluster.samples) {
        betas <- column.cluster(betas[names(probes),,drop=FALSE])$mat }

    if (draw) { assemble_plots(betas, txns, probes,
        plt.txns, plt.mapLines, plt.cytoband, ...)
    } else { return(betas); }
}

#' Visualize Gene
#'
#' Visualize the beta value in heatmaps for a given gene. The function takes
#' a gene name which is taken from the UCSC refGene. It searches all the
#' transcripts for the given gene and optionally extend the span by certain
#' number of base pairs. The function also takes a beta value matrix with
#' sample names on the columns and probe names on the rows. The function can
#' also work on different genome builds (default to hg38, can be hg19).
#'
#' @param gene_name gene name
#' @param betas beta value matrix (row: probes, column: samples)
#' @param platform HM450, EPIC, or MM285 (default)
#' @param upstream distance to extend upstream
#' @param dwstream distance to extend downstream
#' @param genome hg19, hg38, or mm10 (default)
#' @param ... additional options, see visualizeRegion, assemble_plots
#' @import grid
#' @return None
#' @examples
#' betas <- sesameDataGet('HM450.76.TCGA.matched')$betas
#' visualizeGene('ADA', betas, 'HM450')
#' @export
visualizeGene <- function(gene_name, betas,
    platform = NULL, genome = NULL,
    upstream = 2000, dwstream = 2000, ...) {

    if (is.null(dim(betas))) { betas <- as.matrix(betas); }
    platform <- sesameData_check_platform(platform, rownames(betas))
    genome <- sesameData_check_genome(genome, platform)

    txns <- sesameData_getGenomeInfo(genome)$txns
    target.txns <- txns[GenomicRanges::mcols(txns)$gene_name == gene_name]
    stopifnot(length(target.txns) > 0)
    target.strand <- as.character(GenomicRanges::strand(target.txns[[1]][1]))
    if (target.strand == '+') {
        pad.start <- upstream
        pad.end <- dwstream
    } else {
        pad.start <- dwstream
        pad.end <- upstream
    }

    merged.exons <- GenomicRanges::reduce(unlist(target.txns))
    visualizeRegion(
        as.character(GenomicRanges::seqnames(merged.exons[1])),
        min(GenomicRanges::start(merged.exons)) - pad.start,
        max(GenomicRanges::end(merged.exons)) + pad.end,
        betas, platform = platform, genome = genome, ...)
}

#' Visualize Region that Contains the Specified Probes
#'
#' Visualize the beta value in heatmaps for the genomic region containing
#' specified probes. The function works only if specified probes can be
#' spanned by a single genomic region. The region can cover more probes
#' than specified. Hence the plotting heatmap may encompass more probes.
#' The function takes as input a string vector of probe IDs (cg/ch/rs-numbers).
#' if draw is FALSE, the function returns the subset beta value matrix
#' otherwise it returns the grid graphics object.
#' 
#' @param probeNames probe names
#' @param betas beta value matrix (row: probes, column: samples)
#' @param platform HM450, EPIC or MM285 (default)
#' @param genome hg19, hg38 or mm10 (default)
#' @param upstream distance to extend upstream
#' @param dwstream distance to extend downstream
#' @param ... additional options, see visualizeRegion and assemble_plots
#' @return None
#' @examples
#' betas <- sesameDataGet('HM450.76.TCGA.matched')$betas
#' visualizeProbes(c('cg22316575', 'cg16084772', 'cg20622019'), betas, 'HM450')
#' @export
visualizeProbes <- function(
    probeNames, betas,
    platform = NULL, genome = NULL,
    upstream = 1000, dwstream = 1000, ...) {

    if (is.null(dim(betas))) { betas <- as.matrix(betas); }
    platform <- sesameData_check_platform(platform, rownames(betas))
    genome <- sesameData_check_genome(genome, platform)

    probes <- sesameData_getManifestGRanges(platform, genome)
    probeNames <- probeNames[probeNames %in% names(probes)]

    if (length(probeNames)==0)
        stop('Probes specified are not well mapped.')

    target.probes <- probes[probeNames]

    regBeg <- min(GenomicRanges::start(target.probes)) - upstream
    regEnd <- max(GenomicRanges::end(target.probes)) + dwstream
    
    visualizeRegion(
        as.character(GenomicRanges::seqnames(
            target.probes[1])), regBeg, regEnd,
        betas, platform = platform, genome = genome, ...)
}
