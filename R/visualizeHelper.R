exonToCDS <- function(exons, cdsStart, cdsEnd) {
    if (is.na(cdsStart) || is.na(cdsEnd) || cdsEnd <= cdsStart) {
        return(NULL); } 
    cds <- exons[(
        (GenomicRanges::start(exons) < cdsEnd) &
            (GenomicRanges::end(exons) > cdsStart))]
    GenomicRanges::start(cds) <- pmax(
        GenomicRanges::start(cds), cdsStart)
    GenomicRanges::end(cds) <- pmin(
        GenomicRanges::end(cds), cdsEnd)
    cds
}

plotTranscript1 <- function(txn, reg, i, beg, end,
    isoformHeight, padHeight, txn.font.size) {

    txn_name <- names(txn)[1]; exons <- txn[[1]]
    meta <- as.data.frame(GenomicRanges::mcols(txn))
    plt.width <- end - beg
    txn.beg <- max(beg, min(GenomicRanges::start(exons))-2000)
    txn.end <- min(end, max(GenomicRanges::end(exons))+2000)
    exons <- subsetByOverlaps(exons, reg)
    txn.strand <- as.character(GenomicRanges::strand(exons[1]))
    lined <- (c(txn.beg, txn.end)-beg) / plt.width # direction is in arrow ends
    
    y.bot <- (i-1) * isoformHeight + padHeight
    y.bot.exon <- y.bot + padHeight
    y.hei <- isoformHeight - 2 * padHeight

    ## transcript name
    g <- gList(grid.text(sprintf('%s (%s)', meta$gene_name, txn_name),
        x=mean(lined), y=y.bot + y.hei + padHeight * 0.5,
        just=c('center','bottom'),
        gp = gpar(fontsize = txn.font.size), draw=FALSE))

    ## plot transcript line
    g <- gList(g, gList(grid.lines(
        x=lined, y=y.bot+y.hei/2, arrow=arrow(length=unit(0.06, "inches"),
            ends=ifelse(txn.strand == "+", "last", "first")), draw=FALSE)))

    g <- gList(g, gList(grid.lines(x=c(0,1), y=y.bot+y.hei/2,
        gp=gpar(lty='dotted'), draw=FALSE)))

    ## plot exons
    g <- gList(g, gList(
        grid.rect((GenomicRanges::start(exons)-beg)/plt.width,
            y.bot + y.hei/2 - y.hei/3, GenomicRanges::width(exons)/plt.width,
            y.hei/3*2, gp=gpar(fill='grey10', lwd=0),
            just=c('left','bottom'), draw=FALSE)))

    ## plot cds
    cds <- exonToCDS(exons, as.integer(meta$cdsStart), as.integer(meta$cdsEnd))
    if (length(cds) > 0) {
        g <- gList(g, gList(
            grid.rect((GenomicRanges::start(cds)-beg)/plt.width,
                y.bot + y.hei/2 - y.hei/6, GenomicRanges::width(cds)/plt.width,
                y.hei/6*2, gp=gpar(fill='red', lwd=0),
                just=c('left','bottom'), draw=FALSE)))
    }
    g
}

## helper function to plot transcript
plotTranscripts <- function(
    txns, reg, beg, end,
    txn.types = c("protein_coding"), txn.font.size = 6) {

    if (!is.null(txn.types)) {
        txns <- txns[
            GenomicRanges::mcols(txns)$transcript_type %in% txn.types] }
    
    if (length(txns) == 0) {
        return(gList(
            grid.rect(0,0.1,1,0.8, just = c('left','bottom'), draw=FALSE),
            grid.text('No transcript found', x=0.5, y=0.5, draw=FALSE)))
    }
    
    isoformHeight <- 1/length(txns)
    padHeight <- isoformHeight*0.2

    do.call(gList, lapply(seq_along(txns), function(i) {
        plotTranscript1(txns[i], reg, i, beg, end,
            isoformHeight, padHeight, txn.font.size)
    }))
}

plotMapLines <- function(probes, beg, end) {
    nprobes <- length(probes)
    x00 <- ((GenomicRanges::start(probes) - beg) / (end - beg))
    y0 <- rep(0.5, length.out=length(probes))
    x1 <- ((seq_len(nprobes) - 0.5)/nprobes)
    y1 <- rep(0, length.out=nprobes)
    x0 <- c(x00, x00)
    x1 <- c(x1, x00)
    y0 <- c(y0, rep(0.5, length.out=length(probes)))
    y1 <- c(y1, rep(1, length.out=length(probes)))
    grid.segments(x0, y0, x1, y1, draw=FALSE)
}

plotCytoBand <- function(
    chrom, beg, end, genomeInfo) {

    cytoBand <- genomeInfo$cytoBand

    ## set cytoband color
    requireNamespace("pals")
    cytoBand2col <- setNames(
        pals::ocean.gray(10)[seq(9,3)],
        c('stalk', 'gneg', 'gpos25', 'gpos50', 'gpos75', 'gpos100'))
    cytoBand2col['acen'] <- 'red'
    cytoBand2col['gvar'] <- cytoBand2col['gpos75']

    ## chromosome range
    cytoBand.target <- cytoBand[cytoBand$chrom == chrom,]
    chromEnd <- max(cytoBand.target$chromEnd)
    chromBeg <- min(cytoBand.target$chromStart)
    chromWid <- chromEnd - chromBeg
    bandColor <- cytoBand2col[as.character(cytoBand.target$gieStain)]

    pltx0 <- (c(beg, end)-chromBeg)/chromWid
    gList(
        grid.text( # coordinate name
            sprintf("%s:%d-%d", chrom, beg, end), 0, 0.9,
            just = c('left','bottom'), draw = FALSE),
        ## cytoband box
        grid.rect(0, 0.35, 1, 0.35, just = c("left", "bottom"),
            gp = gpar(col = "black", lwd=2, lty="solid"), draw = FALSE),
        grid.rect( # cytoband
            vapply(cytoBand.target$chromStart,
                function(x) (x-chromBeg)/chromWid, 1),
            0.35,
            (cytoBand.target$chromEnd - cytoBand.target$chromStart)/chromWid,
            0.35, gp = gpar(fill = bandColor, col = bandColor),
            just = c('left','bottom'), draw = FALSE),
        grid.segments( # sentinel bar
            x0 = pltx0, y0 = 0.1, x1 = pltx0, y1 = 0.9,
            gp = gpar(col = "red"), draw = FALSE))
}

#' assemble plots
#'
#' @param betas beta value
#' @param txns transcripts GRanges
#' @param probes probe GRanges
#' @param plt.txns transcripts plot objects
#' @param plt.mapLines map line plot objects
#' @param plt.cytoband cytoband plot objects
#' @param heat.height heatmap height (auto inferred based on rows)
#' @param mapLine.height height of the map lines
#' @param show.probeNames whether to show probe names
#' @param show.samples.n number of samples to show (default: all)
#' @param show.sampleNames whether to show sample names
#' @param sample.name.fontsize sample name font size
#' @param dmin data min
#' @param dmax data max
#' @return a grid object
assemble_plots <- function(
    betas, txns, probes, plt.txns, plt.mapLines, plt.cytoband,
    heat.height = NULL, mapLine.height = 0.2,
    show.probeNames = TRUE, show.samples.n = NULL,
    show.sampleNames = TRUE, sample.name.fontsize = 10,
    dmin = 0, dmax = 1) {
    
    if (is.null(show.samples.n)) { show.samples.n <- ncol(betas); }
    if (is.null(heat.height) && length(txns) > 0) {
        heat.height <- 10 / length(txns); }
    w <- WGrob(plt.txns, name = 'txn')
    w <- w + WGrob(plt.mapLines, Beneath(pad=0, height=mapLine.height))
    w <- w + WHeatmap(
        t(betas), Beneath(height = heat.height),
        name = 'betas',
        cmp = CMPar(dmin=dmin, dmax=dmax),
        xticklabels = show.probeNames,
        xticklabel.rotat = 45,
        yticklabels = show.sampleNames,
        yticklabel.fontsize = sample.name.fontsize,
        yticklabels.n = show.samples.n,
        xticklabels.n = length(probes))
    w <- w + WGrob(plt.cytoband, TopOf('txn', height=0.15))
    w
}
