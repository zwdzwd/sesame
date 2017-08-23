
#' visualize gene
#'
#' @param geneName gene name
#' @param betas beta value matrix (row: probes, column: samples)
#' @param platform hm450 (default) or EPIC
#' @param upstream distance to extend upstream
#' @param dwstream distance to extend downstream
#' @param refversion hg19 or hg38
#' @param ... additional options, see visualizeRegion
#' @import grid
#' @export
visualizeGene <- function(geneName, betas, platform='EPIC',
                          upstream=2000, dwstream=2000,
                          refversion='hg38', ...) {

  if (is.null(dim(betas))) {
    betas <- as.matrix(betas);
  }
  
  pkgTest('GenomicRanges')

  gene2txn <- getBuiltInData(paste0('UCSC.refGene.gene2txn.', refversion))
  if (!(geneName %in% names(gene2txn))) {
    stop('Gene ', geneName, ' not found in this reference.');
  }
  txns <- getBuiltInData(paste0('UCSC.refGene.txns.', refversion))

  target.txns <- txns[gene2txn[[geneName]]]
  target.strand <- as.character(GenomicRanges::strand(target.txns[[1]][1]))
  if (target.strand == '+') {
    pad.start <- upstream
    pad.end <- dwstream
  } else {
    pad.start <- dwstream
    pad.end <- upstream
  }

  merged.exons <- GenomicRanges::reduce(GenomicRanges::unlist(target.txns))
  visualizeRegion(as.character(GenomicRanges::seqnames(merged.exons[1])),
                  min(GenomicRanges::start(merged.exons)) - pad.start,
                  max(GenomicRanges::end(merged.exons)) + pad.end,
                  betas, platform=platform, refversion=refversion, ...)
}

#' visualize probes
#'
#' @param probeNames probe names
#' @param betas beta value matrix (row: probes, column: samples)
#' @param platform hm450 (default) or EPIC
#' @param refversion hg19 or hg38
#' @param upstream distance to extend upstream
#' @param dwstream distance to extend downstream
#' @param ... additional options, see visualizeRegion
#' @export
visualizeProbes <- function(probeNames, betas, platform='EPIC', refversion='hg38', upstream=1000, dwstream=1000, ...) {
  pkgTest('GenomicRanges')

  probes <- getBuiltInData(paste0('mapped.probes.', refversion), platform=platform)
  probeNames <- probeNames[probeNames %in% names(probes)]
  if (length(probeNames)==0)
    stop('Probes specified are not well mapped.')
  target.probes <- probes[probeNames]
  regBeg <- min(GenomicRanges::start(target.probes)) - upstream
  regEnd <- max(GenomicRanges::end(target.probes)) + dwstream
  visualizeRegion(as.character(GenomicRanges::seqnames(target.probes[1])), regBeg, regEnd,
                  betas, platform=platform, refversion=refversion, ...)
}

#' visualize region
#'
#' @param chrm chromosome
#' @param plt.beg begin of the region
#' @param plt.end end of the region
#' @param betas beta value matrix (row: probes, column: samples)
#' @param platform hm450 (default) or EPIC
#' @param refversion hg19 or hg38
#' @param draw draw figure or return betas
#' @param heat.height heatmap height (auto inferred based on rows)
#' @param show.sampleNames whether to show sample names
#' @param show.probeNames whether to show probe names
#' @param show.samples.n number of samples to show (default: all)
#' @param cluster.samples whether to cluster samples
#' @param sample.name.fontsize sample name font size
#' @param dmin data min
#' @param dmax data max
#' @param na.rm remove probes with all NA.
#' @import grid
#' @export
visualizeRegion <- function(chrm, plt.beg, plt.end, betas, platform='EPIC', refversion='hg38', sample.name.fontsize=10, heat.height=NULL, draw=TRUE, show.sampleNames=TRUE, show.samples.n=NULL, show.probeNames=TRUE, cluster.samples=FALSE, na.rm=FALSE, dmin=0, dmax=1) {

  pkgTest('GenomicRanges')
  
  if (is.null(dim(betas))) {
    betas <- as.matrix(betas);
  }

  txns <- getBuiltInData(paste0('UCSC.refGene.txns.', refversion))
  txn2gene <- getBuiltInData(paste0('UCSC.refGene.txn2gene.', refversion))
  probes <- getBuiltInData(paste0('mapped.probes.', refversion), platform=platform)

  target.region <- GenomicRanges::GRanges(chrm, IRanges::IRanges(plt.beg, plt.end))
  target.txns <- GenomicRanges::subsetByOverlaps(txns, target.region)

  plt.width <- plt.end-plt.beg
  probes <- GenomicRanges::subsetByOverlaps(probes, target.region)
  probes <- probes[names(probes) %in% rownames(betas)]
  if (na.rm)
    probes <- probes[apply(betas[names(probes), ], 1, function(x) !all(is.na(x)))]
  nprobes <- length(probes)

  if (is.null(show.samples.n))
    show.samples.n <- dim(betas)[2]
  
  if (nprobes == 0)
    stop("No probe overlap region ", sprintf('%s:%d-%d', chrm, plt.beg, plt.end))

  if (nprobes > 1000) {
    stop('Too many probes. Consider smaller region?')
  }

  isoformHeight <- 1/length(target.txns)
  padHeight <- isoformHeight*0.2

  ## plot transcripts
  if (length(target.txns) > 0) {
    plt.txns <- do.call(gList, lapply(1:length(target.txns), function(i) {
      txn <- target.txns[[i]]
      txn.name <- names(target.txns)[i]

      txn.beg <- max(plt.beg, min(GenomicRanges::start(txn))-2000)
      txn.end <- min(plt.end, max(GenomicRanges::end(txn))+2000)

      txn <- GenomicRanges::subsetByOverlaps(txn, target.region)

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
        grid.text(sprintf('%s (%s)', txn.name, txn2gene[[txn.name]][1]),
                  x=mean(line.direc), y=y.bot+y.hei+padHeight*0.5,
                  just=c('center','bottom'), draw=FALSE),
        
        ## plot transcript line
        grid.lines(x=line.direc, y=y.bot+y.hei/2, arrow=arrow(), draw=FALSE),
        grid.lines(x=c(0,1), y=y.bot+y.hei/2, gp=gpar(lty='dotted'), draw=FALSE),

        ## plot exons
        grid.rect((GenomicRanges::start(txn)-plt.beg)/plt.width, y.bot.exon, 
                  GenomicRanges::width(txn)/plt.width, y.hei.exon, gp=gpar(fill='red',col='red'),
                  just=c('left','bottom'), draw=FALSE))

      ## plot cds
      cdsEnd <- as.integer(GenomicRanges::mcols(txn)$cdsEnd[1])
      cdsStart <- as.integer(GenomicRanges::mcols(txn)$cdsStart[1])
      txnCds <- txn[(GenomicRanges::start(txn) < cdsEnd) & (GenomicRanges::end(txn) > cdsStart)]
      GenomicRanges::start(txnCds) <- pmax(GenomicRanges::start(txnCds), cdsStart)
      GenomicRanges::end(txnCds) <- pmin(GenomicRanges::end(txnCds), cdsEnd)

      if (cdsEnd > cdsStart && length(txnCds) > 0) {
        g <- gList(g, gList(
          grid.rect((GenomicRanges::start(txnCds)-plt.beg)/plt.width, y.bot,
                    GenomicRanges::width(txnCds)/plt.width, y.hei, gp=gpar(fill='grey'),
                    just=c('left','bottom'), draw=FALSE)))
      }
      g
    }))
  } else {
    plt.txns <- gList(
      grid.rect(0,0.1,1,0.8, just = c('left','bottom'), draw=FALSE),
      grid.text('No transcript found', x=0.5, y=0.5, draw=FALSE))
  }

  plt.chromLine <- grid.lines(x=c(0, 1), y=c(1,1), draw=FALSE)
  plt.mapLines <- grid.segments((GenomicRanges::start(probes)-plt.beg) / plt.width,
                                1, ((1:nprobes-0.5)/nprobes), 0, draw=FALSE)
  if (draw) {
    pkgTest('wheatmap')
    if (cluster.samples) {
      betas <- betas[,wheatmap::column.cluster(betas[names(probes),,drop=FALSE])$column.clust$order]
    }
    w <- wheatmap::WGrob(plt.txns, name='txn') +
      wheatmap::WGrob(plt.mapLines, wheatmap::Beneath(pad=0, height=0.15)) +
        wheatmap::WHeatmap(t(betas[names(probes),,drop=FALSE]), wheatmap::Beneath(height=heat.height), name='betas', cmp=wheatmap::CMPar(dmin=dmin, dmax=dmax), xticklabels=show.probeNames, xticklabel.rotat=45, yticklabels=show.sampleNames, yticklabel.fontsize=sample.name.fontsize, yticklabels.n=show.samples.n, xticklabels.n=nprobes)
    w <- w + wheatmap::WGrob(
      plotCytoBand(chrm, plt.beg, plt.end, refversion=refversion),
      wheatmap::TopOf('txn', height=0.25))
    w
  } else {
    if (cluster.samples) {
      pkgTest('wheatmap')
      betas <- betas[,wheatmap::column.cluster(betas[names(probes),,drop=FALSE])$column.clust$order]
    } else {
      betas <- betas[names(probes),,drop=FALSE]
    }
    betas
  }
}

## plot chromosome of genomic ranges and cytobands
plotCytoBand <- function(chrom, plt.beg, plt.end, refversion='hg38') {

  cytoBand <- getBuiltInData(paste0('cytoBand.', refversion))
  
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
    grid.text(sprintf("%s:%d-%d", chrom, plt.beg, plt.end), 0, 0.9,
              just=c('left','bottom'), draw=FALSE),
    grid.rect(sapply(cytoBand.target$chromStart, function(x) (x-chromBeg)/chromWid), 0.35,
              (cytoBand.target$chromEnd - cytoBand.target$chromStart)/chromWid, 0.5,
              gp=gpar(fill=bandColor, col=bandColor),
              just=c('left','bottom'), draw=FALSE),
    grid.segments(x0=pltx0, y0=0.3, x1=c(0,1), y1=0.1, draw=FALSE, gp=gpar(lty='dotted')),
    grid.segments(x0=pltx0, y0=0.3, x1=pltx0, y1=0.85, draw=FALSE))
}
