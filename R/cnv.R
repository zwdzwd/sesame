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
#' @param chrominfo chromosome information object
#' @param probe.coords probe coordinates
#' @return bin.coords
#' @import GenomicRanges
#' @importFrom IRanges IRanges
#' @export
getBinCoordinates <- function(chrominfo, probe.coords) {
  tiles <- sort(tileGenome(chrominfo$seqinfo, tilewidth=50000, cut.last.tile.in.chrom = T))
  tiles <- sort(c(setdiff(tiles[seq(1, length(tiles), 2)], chrominfo$gap), 
                  setdiff(tiles[seq(2, length(tiles), 2)], chrominfo$gap)))
  values(tiles)$probes <- countOverlaps(tiles, probe.coords)
  chrom.windows <- as.data.frame(split(sort(tiles), seqnames(tiles))[[2]])

  bin.coords <- do.call(rbind, lapply(split(tiles, seqnames(tiles)),
                                      function(chrom.tiles) leftRightMerge1(as.data.frame(sort(chrom.tiles)))))
  bin.coords <- sort(GRanges(seqnames=bin.coords$seqnames, IRanges(start=bin.coords$start, end=bin.coords$end), seqinfo = seqinfo(tiles)))
  chr.cnts <- table(as.vector(seqnames(bin.coords)))
  chr.names <- as.vector(seqnames(bin.coords))
  names(bin.coords) <- paste(as.vector(seqnames(bin.coords)),
                             formatC(unlist(lapply(seqnames(bin.coords)@lengths, seq_len)),
                                     width=nchar(max(chr.cnts)),
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
#' @param probe.signals probe signals
#' @param bin.coords bin coordinates
#' @param probe.coords probe coordinates
#' @import GenomicRanges
#' @export
binSignals <- function(probe.signals, bin.coords, probe.coords) {
  ov <- findOverlaps(probe.coords, bin.coords)
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
#' @param bin.signals bin signals (input)
#' @param bin.coords bin coordinates
#' @return segment signal data frame
#' @import DNAcopy
#' @export
segmentBins <- function(bin.signals, bin.coords) {

  ## make input data frame
  cna <- DNAcopy::CNA(genomdat=bin.signals,
                      chrom=as.character(seqnames(bin.coords)),
                      maploc=as.integer((start(bin.coords) + end(bin.coords))/2),
                      data.type='logratio')

  seg <- DNAcopy::segment(x=cna, min.width = 5, 
                          nperm= 10000, alpha = 0.001, undo.splits = 'sdundo',
                          undo.SD= 2.2, verbose=0)
  summary <- segments.summary(seg)
  pval <- segments.p(seg)
  seg.signals <- cbind(summary, pval[,c('pval','lcl','ucl')])
  seg.signals$chrom <- as.character(seg.signals$chrom)
  seg.signals
}

#' Visualize segments
#'
#' @param seg segments
#' @param to.plot chromosome to plot (by default plot all chromosomes)
#' @importFrom scales squish
#' @import ggplot2 
#' @export
visualize.segments <- function(seg, to.plot=NULL) {

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

  bin.coords <- bin.coords[as.vector(seqnames(bin.coords)) %in% seq.names]
  values(bin.coords)$bin.mids <- (start(bin.coords) + end(bin.coords))/2
  values(bin.coords)$bin.x <- seqstart[as.character(seqnames(bin.coords))] + bin.coords$bin.mids

  ## plot bin
  p <- qplot(bin.coords$bin.x / total.length, bin.signals, color=bin.signals, alpha=I(0.7))
  ## plot segment
  p <- p + geom_segment(aes(x=(seqstart[seg.signals$chrom] + seg.signals$loc.start) / total.length,
                            xend=(seqstart[seg.signals$chrom] + seg.signals$loc.end) / total.length,
                            y=seg.signals$seg.mean, yend=seg.signals$seg.mean), size=1.5, color='blue')
  ## chromosome boundary
  p <- p + geom_vline(xintercept=seqstart[-1]/total.length, alpha=I(0.5))
  
  ## chromosome label
  p <- p + scale_x_continuous(labels=seq.names, breaks=(seqstart+seqlen/2)/total.length) + theme(axis.text.x = element_text(angle=90, hjust=0.5))
  
  ## styling
  p <- p + scale_colour_gradient2(limits=c(-0.3,0.3), low='red', mid='grey', high='green', oob=squish) + xlab('') + ylab('') + theme(legend.position="none")
  p
}
