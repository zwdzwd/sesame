## helper function to plot transcript
plotTranscripts <- function(
    target.txns, target.region, plt.beg, plt.end) {

    if (length(target.txns) == 0) {
        return(gList(
            grid.rect(0,0.1,1,0.8, just = c('left','bottom'), draw=FALSE),
            grid.text('No transcript found', x=0.5, y=0.5, draw=FALSE)))
    }

    plt.width <- plt.end - plt.beg
    isoformHeight <- 1/length(target.txns)
    padHeight <- isoformHeight*0.2

    do.call(gList, lapply(seq_along(target.txns), function(i) {
        txn <- target.txns[[i]]
        txn.name <- names(target.txns)[i]

        txn.beg <- max(plt.beg, min(GenomicRanges::start(txn))-2000)
        txn.end <- min(plt.end, max(GenomicRanges::end(txn))+2000)

        txn <- subsetByOverlaps(txn, target.region)

        txn.strand <- as.character(GenomicRanges::strand(txn[1]))
        if (txn.strand == '+') {
            line.direc <- c(txn.beg-plt.beg, txn.end-plt.beg) / plt.width
        } else {
            line.direc <- c(txn.end-plt.beg, txn.beg-plt.beg) / plt.width
        }
        
        y.bot <- (i-1)*isoformHeight+padHeight
        y.bot.exon <- y.bot+padHeight
        y.hei <- isoformHeight-2*padHeight
        y.hei.exon <- isoformHeight-4*padHeight

        g <- gList(
            ## plot transcript name
            grid.text(sprintf('%s (%s)', txn.name,
                GenomicRanges::mcols(target.txns)$gene_name[i]),
                x=mean(line.direc), y=y.bot+y.hei+padHeight*0.5,
                just=c('center','bottom'), draw=FALSE),
            
            ## plot transcript line
            grid.lines(
                x=line.direc, y=y.bot+y.hei/2, arrow=arrow(), draw=FALSE),

            grid.lines(
                x=c(0,1), y=y.bot+y.hei/2,
                gp=gpar(lty='dotted'), draw=FALSE),

            ## plot exons
            grid.rect(
            (GenomicRanges::start(txn)-plt.beg)/plt.width, y.bot.exon, 
            GenomicRanges::width(txn)/plt.width, y.hei.exon,
            gp=gpar(fill='red',col='red'),
            just=c('left','bottom'), draw=FALSE))

        ## plot cds
        cdsEnd <- as.integer(GenomicRanges::mcols(txn)$cdsEnd[1])
        cdsStart <- as.integer(GenomicRanges::mcols(txn)$cdsStart[1])
        txnCds <- txn[(GenomicRanges::start(txn) < cdsEnd) &
                          (GenomicRanges::end(txn) > cdsStart)]
        
        GenomicRanges::start(txnCds) <- pmax(
            GenomicRanges::start(txnCds), cdsStart)
        
        GenomicRanges::end(txnCds) <- pmin(
            GenomicRanges::end(txnCds), cdsEnd)

        if (cdsEnd > cdsStart && length(txnCds) > 0) {
            g <- gList(g, gList(grid.rect(
            (GenomicRanges::start(txnCds)-plt.beg)/plt.width, y.bot,
            GenomicRanges::width(txnCds)/plt.width,
            y.hei, gp=gpar(fill='grey'),
            just=c('left','bottom'), draw=FALSE)))
        }
        g
    }))
}

plotMapLines <- function(probes, plt.beg, plt.end) {
    nprobes <- length(probes)
    grid.segments(
    ((GenomicRanges::start(probes) - plt.beg) / (plt.end - plt.beg)), 1,
    ((seq_len(nprobes) - 0.5)/nprobes), 0, draw=FALSE)
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
#' @param show.probeNames whether to show probe names
#' @param show.samples.n number of samples to show (default: all)
#' @param show.sampleNames whether to show sample names
#' @param sample.name.fontsize sample name font size
#' @param dmin data min
#' @param dmax data max
#' @return a grid object
plot_assemble <- function(
    betas, txns, probes, plt.txns, plt.mapLines, plt.cytoband,
    heat.height = NULL, 
    show.probeNames = TRUE, show.samples.n = NULL,
    show.sampleNames = TRUE, sample.name.fontsize = 10,
    dmin = 0, dmax = 1) {
    
    if (is.null(show.samples.n)) { show.samples.n <- ncol(betas); }
    if (is.null(heat.height)) { heat.height <- 10 / length(txns); }
    w <- WGrob(plt.txns, name = 'txn')
    w <- w + WGrob(plt.mapLines, Beneath(pad=0, height=0.15))
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
    w <- w + WGrob(plt.cytoband, TopOf('txn', height=0.25))
    w
}
