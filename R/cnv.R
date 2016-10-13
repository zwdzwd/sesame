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

  ## segmentation
  probe.signals <- probeSignals(sset.target, ssets.normal)
  bin.signals <- binSignals(probe.signals, bin.coords, probe.coords)
  segs <- segmentBins(bin.signals)
  segs$chrominfo <- chrominfo
  segs$bin.coords <- bin.coords

  segs
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
      if (is.na(chrom.windows[min.window-1,'probes']) || is.na(chrom.windows[min.window+1,'probes'])) {
        message(min.window)
        browser()
      }
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
#' @param chrominfo chromosome information object
#' @param probe.coords probe coordinates
#' @return bin.coords
#' @export
getBinCoordinates <- function(chrominfo, probe.coords) {
  tiles <- sort(tileGenome(chrominfo$seqinfo, tilewidth=50000, cut.last.tile.in.chrom = T))
  tiles <- sort(c(setdiff(tiles[seq(1, length(tiles), 2)], chrominfo$gap), 
                  setdiff(tiles[seq(2, length(tiles), 2)], chrominfo$gap)))
  values(tiles)$probes <- countOverlaps(tiles, probe.coords)
  chrom.windows <- as.data.frame(split(sort(tiles), seqnames(tiles))[[2]])

  bin.coords <- do.call(rbind, lapply(split(tiles, seqnames(tiles)), function(chrom.tiles) leftRightMerge1(as.data.frame(sort(chrom.tiles)))))
  bin.coords <- sort(GRanges(seqnames=bin.coords$seqnames, IRanges(start=bin.coords$start, end=bin.coords$end), seqinfo = seqinfo(tiles)))
  names(bin.coords) <- paste(seqnames(bin.coords), 
                             formatC(unlist(lapply(table(seqnames(bin.coords)), seq_len)), 
                                     width=nchar(max(table(seqnames(bin.coords)))),
                                     format='d', flag='0'), sep='-')
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
#' @param probe.signals
#' @param bin.coords bin coordinates
#' @param probe.coords probe coordinates
#' @export
binSignals <- function(probe.signals, bin.coords, probe.coords) {
  ov <- findOverlaps(probe.coords, bin.coords)
  .bins <- names(bin.coords)[ov@subjectHits]
  .probe.signals <- probe.signals[names(probe.coords)[ov@queryHits]]

  bin.signals <- sapply(split(.probe.signals, .bins), median, na.rm=TRUE)
}

#' Segment bins
#'
#' @param bin.signals bin signals (input)
#' @return segment data frame
#' @export
segmentBins <- function(bin.signals) {
  ## make input data frame
  cna <- DNAcopy::CNA(genomdat=bin.signals,
                      chrom=as.character(seqnames(bin.coords)),
                      maploc=(start(bin.coords) + end(bin.coords))/2,
                      data.type='logratio')

  seg <- DNAcopy::segment(x=cna, min.width = 5, 
                          nperm= 10000, alpha = 0.001, undo.splits = 'sdundo',
                          undo.SD= 2.2, verbose=0)
  summary <- segments.summary(seg)
  pval <- segments.p(seg)
  segs <- cbind(summary, pval[,c('pval','lcl','ucl')])
  segs$chrom <- as.character(segs$chrom)
  segs
}

#' Visualize segments
#'
#' @param segs
#' @param to.plot chromosome to plot (by default plot all chromosomes)
#' @export
visualize.segments <- function(segs, to.plot=NULL) {

  chrominfo <- segs$chrominfo
  bin.coords <- segs$bin.coords
  library(scales)
  total.length <- sum(as.numeric(chrominfo$seqinfo@seqlengths), na.rm=T)
  # skip chromosome too small (e.g, chrM)
  if (is.null(to.plot))
    to.plot <- (chrominfo$seqinfo@seqlengths > total.length*0.01)
  seqlen <- as.numeric(chrominfo$seqinfo@seqlengths[to.plot])
  seqnames <- chrominfo$seqinfo@seqnames[to.plot]
  total.length <- sum(seqlen, na.rm=T)
  seqcumlen <- cumsum(seqlen)
  seqcumlen <- seqcumlen[-length(seqcumlen)]
  seqstart <- setNames(c(0,seqcumlen), seqnames)
  
  bin.coords <- bin.coords[seqnames(bin.coords) %in% seqnames]
  values(bin.coords)$bin.mids <- (start(bin.coords) + end(bin.coords))/2
  values(bin.coords)$bin.x <- seqstart[as.character(seqnames(bin.coords))] + bin.coords$bin.mids

  ## plot bin
  p <- qplot(bin.coords$bin.x / total.length, bin.signals, color=bin.signals, alpha=I(0.7))
  ## plot segment
  p <- p + geom_segment(aes(x=(seqstart[segs$chrom] + segs$loc.start) / total.length,
                            xend=(seqstart[segs$chrom] + segs$loc.end) / total.length,
                            y=segs$seg.mean, yend=segs$seg.mean), size=1.5, color='blue')
  ## chromosome boundary
  p <- p + geom_vline(xintercept=seqstart[-1]/total.length, alpha=I(0.5))
  
  ## chromosome label
  p <- p + scale_x_continuous(labels=seqnames, breaks=(seqstart+seqlen/2)/total.length) + theme(axis.text.x = element_text(angle=90, hjust=0.5))
  
  ## styling
  p <- p + scale_colour_gradient2(limits=c(-0.3,0.3), low='red', mid='grey', high='green', oob=squish) + xlab('') + ylab('') + theme(legend.position="none")
  p
}
