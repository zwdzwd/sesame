#' perform copy number segmentation
#'
#' @param sset \code{SignalSet}
#' @param ssets.normal \code{SignalSet} for normalization
#' @param refversion hg19 or hg38
#' @return segmentation
#' @export
cn.segmentation <- function(sset, ssets.normal=NULL, refversion='hg19') {

  ## retrieve chromosome info and probe coordinates
  chrominfo <- getBuiltInData(paste0(refversion, '.chrominfo'))
  probe.coords <- getBuiltInData(paste0(sset$platform, '.mapped.probes.', refversion))

  ## fix bin coordinates
  bin.coords <- getBinCoordinates(chrominfo, probe.coords)

  if (is.null(ssets.normal)) {
    ssets.normal <- getBuiltInData('cn.normals', sset$platform);
  }

  ## segmentation
  probe.signals <- probeSignals(sset, ssets.normal)
  bin.signals <- binSignals(probe.signals, bin.coords, probe.coords)
  seg <- list()
  seg$seg.signals <- segmentBins(bin.signals, bin.coords)
  seg$chrominfo <- chrominfo
  seg$bin.coords <- bin.coords
  seg$bin.signals <- bin.signals

  seg
}

## Left-right merge bins
leftRightMerge1 <- function(chrom.windows, min.probes.per.bin=20) {
  while (dim(chrom.windows)[1] > 0 && min(chrom.windows[,'probes']) < min.probes.per.bin) {
    min.window <- which.min(chrom.windows$probes)
    merge.left <- FALSE
    merge.right <- FALSE
    if (min.window > 1 && 
          chrom.windows[min.window,'start']-1 == chrom.windows[min.window-1, 'end']) {
      merge.left <- TRUE
      }
    if (min.window < dim(chrom.windows)[1] &&
          chrom.windows[min.window, 'end']+1 == chrom.windows[min.window+1, 'start']) {
      merge.right <- TRUE
      }
    
    if (merge.left && merge.right) {
      if (chrom.windows[min.window-1,'probes'] < chrom.windows[min.window+1,'probes']) {
        merge.right <- FALSE
        }
      }
    
    if (merge.left) {
      chrom.windows[min.window-1, 'end'] <- chrom.windows[min.window, 'end']
      chrom.windows[min.window-1, 'probes'] <- chrom.windows[min.window-1, 'probes'] + chrom.windows[min.window, 'probes']
      chrom.windows <- chrom.windows[-min.window,]
    } else if (merge.right) {
      chrom.windows[min.window+1, 'start'] <- chrom.windows[min.window, 'start']
      chrom.windows[min.window+1, 'probes'] <- chrom.windows[min.window+1, 'probes'] + chrom.windows[min.window, 'probes']
      chrom.windows <- chrom.windows[-min.window,]
    } else {
      chrom.windows <- chrom.windows[-min.window,]
    }
  }
  chrom.windows
}

#' Get bin coordinates
#'
#' requires GenomicRanges, IRanges
#' 
#' @param chrominfo chromosome information object
#' @param probe.coords probe coordinates
#' @return bin.coords
#' @export
getBinCoordinates <- function(chrominfo, probe.coords) {
  pkgTest('GenomicRanges')
  pkgTest('IRanges')
  
  tiles <- sort(GenomicRanges::tileGenome(chrominfo$seqinfo, tilewidth=50000, cut.last.tile.in.chrom = T))
  tiles <- sort(c(setdiff(tiles[seq(1, length(tiles), 2)], chrominfo$gap), 
                  setdiff(tiles[seq(2, length(tiles), 2)], chrominfo$gap)))
  GenomicRanges::values(tiles)$probes <- GenomicRanges::countOverlaps(tiles, probe.coords)
  chrom.windows <- as.data.frame(split(sort(tiles), GenomicRanges::seqnames(tiles))[[2]])

  bin.coords <- do.call(rbind, lapply(split(tiles, GenomicRanges::seqnames(tiles)),
                                      function(chrom.tiles) leftRightMerge1(as.data.frame(sort(chrom.tiles)))))
  bin.coords <- sort(GenomicRanges::GRanges(seqnames=bin.coords$seqnames, IRanges::IRanges(start=bin.coords$start, end=bin.coords$end), seqinfo = GenomicRanges::seqinfo(tiles)))
  chr.cnts <- table(as.vector(GenomicRanges::seqnames(bin.coords)))
  chr.names <- as.vector(GenomicRanges::seqnames(bin.coords))
  names(bin.coords) <- paste(as.vector(GenomicRanges::seqnames(bin.coords)), formatC(unlist(lapply(GenomicRanges::seqnames(bin.coords)@lengths, seq_len)), width=nchar(max(chr.cnts)), format='d', flag='0'), sep='-')
  bin.coords
}

#' Compute probe signals
#'
#' @param sset.target \code{SignalSet} for target
#' @param ssets.normal \code{SignalSet}s from normal samples
#' @return probe.signals
#' @export
probeSignals <- function(sset.target, ssets.normal) {
  target.intens <- sset.target$totalIntensities()
  normal.intens <- sapply(ssets.normal, function(sset) {
    sset$totalIntensities()
  })

  fit <- lm(y~., data=data.frame(y=target.intens, X=normal.intens))
  probe.signals <- setNames(log2(target.intens / pmax(predict(fit), 1)), names(target.intens))
  probe.signals
}

#' Bin signals from probe signals
#'
#' require GenomicRanges
#' 
#' @param probe.signals probe signals
#' @param bin.coords bin coordinates
#' @param probe.coords probe coordinates
#' @export
binSignals <- function(probe.signals, bin.coords, probe.coords) {
  pkgTest('GenomicRanges')
  ov <- GenomicRanges::findOverlaps(probe.coords, bin.coords)
  if (.hasSlot(ov, 'queryHits')) {
      .bins <- names(bin.coords)[ov@subjectHits]
      .probe.signals <- probe.signals[names(probe.coords)[ov@queryHits]]
    } else {
      .bins <- names(bin.coords)[ov@to]
      .probe.signals <- probe.signals[names(probe.coords)[ov@from]]
    }

  bin.signals <- sapply(split(.probe.signals, .bins), median, na.rm=TRUE)
}

#' Segment bins
#'
#' require DNAcopy
#' @param bin.signals bin signals (input)
#' @param bin.coords bin coordinates
#' @return segment signal data frame
#' @export
segmentBins <- function(bin.signals, bin.coords) {

  pkgTest('DNAcopy')
  
  ## make input data frame
  cna <- DNAcopy::CNA(genomdat=bin.signals,
                      chrom=as.character(GenomicRanges::seqnames(bin.coords)),
                      maploc=as.integer((start(bin.coords) + end(bin.coords))/2),
                      data.type='logratio')

  seg <- DNAcopy::segment(x=cna, min.width = 5, 
                          nperm= 10000, alpha = 0.001, undo.splits = 'sdundo',
                          undo.SD= 2.2, verbose=0)
  summary <- DNAcopy::segments.summary(seg)
  pval <- DNAcopy::segments.p(seg)
  seg.signals <- cbind(summary, pval[,c('pval','lcl','ucl')])
  seg.signals$chrom <- as.character(seg.signals$chrom)
  seg.signals
}

#' Visualize segments
#'
#' require ggplot2, scales
#' @param seg segments
#' @param to.plot chromosome to plot (by default plot all chromosomes)
#' @export
visualize.segments <- function(seg, to.plot=NULL) {

  pkgTest('ggplot2')
  pkgTest('scales')
  
  chrominfo <- seg$chrominfo
  bin.coords <- seg$bin.coords
  bin.signals <- seg$bin.signals
  seg.signals <- seg$seg.signals
  total.length <- sum(as.numeric(chrominfo$seqinfo@seqlengths), na.rm=T)
  # skip chromosome too small (e.g, chrM)
  if (is.null(to.plot))
    to.plot <- (chrominfo$seqinfo@seqlengths > total.length*0.01)
  seqlen <- as.numeric(chrominfo$seqinfo@seqlengths[to.plot])
  seq.names <- chrominfo$seqinfo@seqnames[to.plot]
  total.length <- sum(seqlen, na.rm=T)
  seqcumlen <- cumsum(seqlen)
  seqcumlen <- seqcumlen[-length(seqcumlen)]
  seqstart <- setNames(c(0,seqcumlen), seq.names)

  bin.coords <- bin.coords[as.vector(GenomicRanges::seqnames(bin.coords)) %in% seq.names]
  GenomicRanges::values(bin.coords)$bin.mids <- (GenomicRanges::start(bin.coords) + GenomicRanges::end(bin.coords))/2
  GenomicRanges::values(bin.coords)$bin.x <- seqstart[as.character(GenomicRanges::seqnames(bin.coords))] + bin.coords$bin.mids

  ## plot bin
  p <- ggplot2::qplot(bin.coords$bin.x / total.length, bin.signals, color=bin.signals, alpha=I(0.7))
  ## plot segment
  p <- p + ggplot2::geom_segment(ggplot2::aes(x=(seqstart[seg.signals$chrom] + seg.signals$loc.start) / total.length, xend=(seqstart[seg.signals$chrom] + seg.signals$loc.end) / total.length, y=seg.signals$seg.mean, yend=seg.signals$seg.mean), size=1.5, color='blue')
  ## chromosome boundary
  p <- p + ggplot2::geom_vline(xintercept=seqstart[-1]/total.length, alpha=I(0.5))
  
  ## chromosome label
  p <- p + ggplot2::scale_x_continuous(labels=seq.names, breaks=(seqstart+seqlen/2)/total.length) + ggplot2::theme(axis.text.x = ggplot2::element_text(angle=90, hjust=0.5))
  
  ## styling
  p <- p + ggplot2::scale_colour_gradient2(limits=c(-0.3,0.3), low='red', mid='grey', high='green', oob=scales::squish) + ggplot2::xlab('') + ggplot2::ylab('') + ggplot2::theme(legend.position="none")
  p
}
