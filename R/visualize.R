#' Visualize Gene
#'
#' Visualize the beta value in heatmaps for a given gene. The function takes
#' a gene name which is taken from the UCSC refGene. It searches all the
#' transcripts for the given gene and optionally extend the span by certain
#' number of base pairs. The function also takes a beta value matrix with
#' sample names on the columns and probe names on the rows. The function can
#' also work on different genome builds (default to hg38, can be hg19).
#'
#' @param geneName gene name
#' @param betas beta value matrix (row: probes, column: samples)
#' @param platform HM450, EPIC, or MM285 (default)
#' @param upstream distance to extend upstream
#' @param dwstream distance to extend downstream
#' @param genome hg19, hg38, or mm10 (default)
#' @param ... additional options, see visualizeRegion
#' @import grid
#' @return None
#' @examples
#' betas <- sesameDataGet('HM450.76.TCGA.matched')$betas
#' visualizeGene('ADA', betas, 'HM450')
#' @export
visualizeGene <- function(
    geneName, betas, platform = c('EPIC','HM450','MM285'),
    upstream = 2000, dwstream = 2000,
    genome = c('hg38', 'hg19', 'mm10'), ...) {

    if (is.null(dim(betas))) { betas <- as.matrix(betas); }
    platform <- sesameData_check_platform(platform, rownames(betas))
    genome <- sesameData_check_genome(genome, platform)
    
    pkgTest('GenomicRanges')
    target.txns <- sesameData_getTranscriptsByGene(geneName, genome)
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
#' @param ... additional options, see visualizeRegion
#' @return None
#' @import wheatmap
#' @examples
#' betas <- sesameDataGet('HM450.76.TCGA.matched')$betas
#' visualizeProbes(c('cg22316575', 'cg16084772', 'cg20622019'), betas, 'HM450')
#' @export
visualizeProbes <- function(
    probeNames, betas,
    platform = c('EPIC', 'HM450', 'MM285'),
    genome = c('hg38','hg19','mm10'),
    upstream = 1000, dwstream = 1000, ...) {

    if (is.null(dim(betas))) { betas <- as.matrix(betas); }
    platform <- sesameData_check_platform(platform, rownames(betas))
    genome <- sesameData_check_genome(genome, platform)
    
    platform <- match.arg(platform)
    genome <- match.arg(genome)
    
    pkgTest('GenomicRanges')
    probes <- sesameDataGet(paste0(
        platform, '.probeInfo'))[[paste0('mapped.probes.',genome)]]
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
#' @param plt.beg begin of the region
#' @param plt.end end of the region
#' @param betas beta value matrix (row: probes, column: samples)
#' @param platform EPIC, HM450, or MM285
#' @param genome hg38, hg19, or mm10
#' @param draw draw figure or return betas
#' @param cluster.samples whether to cluster samples
#' @param nprobes.max maximum number of probes to plot
#' @param na.rm remove probes with all NA.
#' @param ... additional options, see plot_assemble
#' @return graphics or a matrix containing the captured beta values
#' @import grid
#' @importMethodsFrom IRanges subsetByOverlaps
#' @importMethodsFrom IRanges IRanges
#' @importFrom GenomicRanges GRanges
#' @examples
#' betas <- sesameDataGet('HM450.76.TCGA.matched')$betas
#' visualizeRegion('chr20', 44648623, 44652152, betas, 'HM450')
#' @export
visualizeRegion <- function(
    chrm, plt.beg, plt.end, betas,
    platform = NULL, genome = NULL,
    draw = TRUE, cluster.samples = FALSE,
    nprobes.max = 1000, na.rm = FALSE, ...) {

    if (is.null(dim(betas))) { betas <- as.matrix(betas) }
    platform <- sesameData_check_platform(platform, rownames(betas))
    genome <- sesameData_check_genome(genome, platform)

    reg <- GRanges(chrm, IRanges(plt.beg, plt.end))
    
    pkgTest('GenomicRanges')
    txns <- sesameDataGet(paste0('genomeInfo.',genome))$txns
    txns <- subsetByOverlaps(txns, reg)

    probes <- sesameData_getManifestGRanges(platform, genome)
    probes <- subsetByOverlaps(probes, reg)
    probes <- probes[names(probes) %in% rownames(betas)]
    if (na.rm) { probes <- probes[apply(
        betas[names(probes), ], 1, function(x) !all(is.na(x)))] }
    
    if (length(probes) == 0) {
        stop("No probe overlap region ", sprintf(
            '%s:%d-%d', chrm, plt.beg, plt.end)) }
    if (length(probes) > nprobes.max) {
        stop(sprintf('Too many probes (%d). Shrink region?', length(probes))) }

    ## plot transcripts
    plt.txns <- plotTranscripts(txns, reg, plt.beg, plt.end)
    plt.mapLines <- plotMapLines(probes, plt.beg, plt.end)
    plt.cytoband <- plotCytoBand(chrm, plt.beg, plt.end, genome)
    
    ## clustering
    pkgTest('wheatmap')
    betas <- betas[names(probes),,drop=FALSE]
    if (cluster.samples) {
        betas <- column.cluster(betas[names(probes),,drop=FALSE])$mat }

    if (draw) {
        plot_assemble(
            betas, txns, probes, plt.txns, plt.mapLines, plt.cytoband, ...)
    } else {
        betas
    }
}

## plot chromosome of genomic ranges and cytobands
#' @importFrom grDevices gray.colors
plotCytoBand <- function(
    chrom, plt.beg, plt.end, genome) {

    cytoBand <- sesameDataGet(paste0('genomeInfo.',genome))$cytoBand

    ## set cytoband color
    cytoBand2col <- setNames(
        gray.colors(7, start=0.9,end=0),
        c('stalk', 'gneg', 'gpos25', 'gpos50', 'gpos75', 'gpos100'))
    cytoBand2col['acen'] <- 'red'
    cytoBand2col['gvar'] <- cytoBand2col['gpos75']

    ## chromosome range
    cytoBand.target <- cytoBand[cytoBand$chrom == chrom,]
    chromEnd <- max(cytoBand.target$chromEnd)
    chromBeg <- min(cytoBand.target$chromStart)
    chromWid <- chromEnd - chromBeg
    bandColor <- cytoBand2col[as.character(cytoBand.target$gieStain)]

    pltx0 <- (c(plt.beg, plt.end)-chromBeg)/chromWid
    gList(
        grid.text(
            sprintf("%s:%d-%d", chrom, plt.beg, plt.end), 0, 0.9,
            just = c('left','bottom'), draw = FALSE),
        grid.rect(
            vapply(
                cytoBand.target$chromStart,
                function(x) (x-chromBeg)/chromWid, 1),
            0.35,
            (cytoBand.target$chromEnd - cytoBand.target$chromStart)/chromWid,
            0.5, gp = gpar(fill = bandColor, col = bandColor),
            just = c('left','bottom'), draw = FALSE),
        grid.segments(
            x0 = pltx0, y0 = 0.3,
            x1 = c(0,1), y1 = 0.1, draw = FALSE, gp = gpar(lty='dotted')),
        grid.segments(
            x0 = pltx0, y0 = 0.3,
            x1 = pltx0, y1 = 0.85, draw = FALSE))
}
