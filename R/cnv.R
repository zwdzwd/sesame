#' Perform copy number segmentation
#'
#' Perform copy number segmentation using the signals in the signal set.
#' The function takes a \code{SigDF} for the target sample and a set of
#' normal \code{SigDF} for the normal samples. An optional arguments specifies
#' the version of genome build that the inference will operate on. The function
#' outputs an object of class \code{CNSegment} with signals for the segments (
#' seg.signals), the bin coordinates (
#' bin.coords) and bin signals (bin.signals).
#'
#' @param sdf \code{SigDF}
#' @param sdfs.normal a list of \code{SigDF}s for normalization, if not given,
#' use the stored normal data from sesameData. However, we do recommend using
#' a matched copy number normal dataset for normalization.
#' @param refversion hg19 or hg38
#' @return an object of \code{CNSegment}
#' @examples
#'
#' sesameDataCache("EPIC") # in case not done yet
#'
#' sdf <- sesameDataGet('EPIC.1.SigDF')
#' sdfs.normal <- sesameDataGet('EPIC.5.SigDF.normal')
#' seg <- cnSegmentation(sdf, sdfs.normal)
#' 
#' @export
cnSegmentation <- function(sdf, sdfs.normal=NULL, refversion=c('hg19','hg38')) {

    stopifnot(is(sdf, "SigDF"))
    pkgTest('GenomicRanges')
    refversion <- match.arg(refversion)

    if (is.null(sdfs.normal)) {
        if (sdfPlatform(sdf) == "EPIC") {
            sdfs.normal = sesameDataGet("EPIC.5.SigDF.normal")
        } else {
            stop(sprintf("For %s, please provide the sdfs.normal argument.",
                sdfPlatform(sdf)))
        }
    }
    
    ## retrieve chromosome info and probe coordinates
    options(sesameData_use_alternative = TRUE) # TODO: REMOVE THIS AFTER EH UPDATE
    seqLength <- sesameDataGet(paste0('genomeInfo.', refversion))$seqLength
    gapInfo <- sesameDataGet(paste0('genomeInfo.', refversion))$gapInfo
    options(sesameData_use_alternative = FALSE) # TODO: REMOVE THIS AFTER EH UPDATE
    probe.coords <- sesameDataGet(paste0(
        sdfPlatform(sdf), '.probeInfo'))[[paste0('mapped.probes.', refversion)]]
    
    ## extract intensities
    target.intens <- totalIntensities(sdf)
    normal.intens <- do.call(cbind, lapply(sdfs.normal, function(sdf) {
        totalIntensities(sdf) }))

    ## find overlapping probes
    pb <- intersect(rownames(normal.intens), names(target.intens))
    pb <- intersect(names(probe.coords), pb)
    target.intens <- target.intens[pb]
    normal.intens <- normal.intens[pb,]
    probe.coords <- probe.coords[pb]
    
    ## normalize probes intensities
    fit <- lm(y~., data=data.frame(y=target.intens, X=normal.intens))
    probe.signals <- setNames(log2(target.intens / pmax(predict(fit), 1)), pb)

    ## bin signals
    ## fix bin coordinates, TODO: this is too time-consuming
    bin.coords <- getBinCoordinates(seqLength, gapInfo, probe.coords)
    bin.signals <- binSignals(probe.signals, bin.coords, probe.coords)

    ## segmentation
    seg <- structure(list(
        seg.signals = segmentBins(bin.signals, bin.coords),
        bin.coords = bin.coords,
        bin.signals = bin.signals), class='CNSegment')
    
    seg
}

## Left-right merge bins
## THIS IS VERY SLOW, SHOULD OPTIMIZE
leftRightMerge1 <- function(chrom.windows, min.probes.per.bin=20) {
    
    while (
        dim(chrom.windows)[1] > 0 &&
            min(chrom.windows[,'probes']) < min.probes.per.bin) {
        
        min.window <- which.min(chrom.windows$probes)
        merge.left <- FALSE
        merge.right <- FALSE
        if (min.window > 1 && 
            chrom.windows[min.window,'start']-1 ==
            chrom.windows[min.window-1, 'end']) {
            merge.left <- TRUE
        }
        if (min.window < dim(chrom.windows)[1] &&
            chrom.windows[min.window, 'end']+1 ==
            chrom.windows[min.window+1, 'start']) {
            merge.right <- TRUE
        }
        
        if (merge.left && merge.right) {
            if (chrom.windows[min.window-1,'probes'] <
                chrom.windows[min.window+1,'probes']) {
                merge.right <- FALSE
            }
        }
        
        if (merge.left) {
            chrom.windows[min.window-1, 'end'] <-
                chrom.windows[min.window, 'end']

            chrom.windows[min.window-1, 'probes'] <-
                chrom.windows[min.window-1, 'probes'] +
                    chrom.windows[min.window, 'probes']
            
            chrom.windows <- chrom.windows[-min.window,]
        } else if (merge.right) {
            chrom.windows[min.window+1, 'start'] <-
                chrom.windows[min.window, 'start']
            chrom.windows[min.window+1, 'probes'] <-
                chrom.windows[min.window+1, 'probes'] +
                    chrom.windows[min.window, 'probes']
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
#' @param seqLength chromosome information object
#' @param gapInfo chromosome gap information
#' @param probe.coords probe coordinates
#' @return bin.coords
getBinCoordinates <- function(seqLength, gapInfo, probe.coords) {

    pkgTest('GenomicRanges')
    pkgTest('IRanges')

    tiles <- sort(GenomicRanges::tileGenome(
        seqLength, tilewidth=50000, cut.last.tile.in.chrom = TRUE))
    
    tiles <- sort(c(
        GenomicRanges::setdiff(tiles[seq(1, length(tiles), 2)], gapInfo), 
        GenomicRanges::setdiff(tiles[seq(2, length(tiles), 2)], gapInfo)))
    
    GenomicRanges::values(tiles)$probes <-
        GenomicRanges::countOverlaps(tiles, probe.coords)
    
    bin.coords <- do.call(rbind, lapply(
        split(tiles, as.vector(GenomicRanges::seqnames(tiles))),
        function(chrom.tiles)
            leftRightMerge1(GenomicRanges::as.data.frame(
                GenomicRanges::sort(chrom.tiles)))))
    
    bin.coords <- GenomicRanges::sort(
        GenomicRanges::GRanges(
            seqnames = bin.coords$seqnames,
            IRanges::IRanges(start = bin.coords$start, end = bin.coords$end),
            seqinfo = GenomicRanges::seqinfo(tiles)))
    
    chr.cnts <- table(as.vector(GenomicRanges::seqnames(bin.coords)))
    chr.names <- as.vector(GenomicRanges::seqnames(bin.coords))

    names(bin.coords) <- paste(
        as.vector(GenomicRanges::seqnames(bin.coords)),
        formatC(unlist(lapply(
            GenomicRanges::seqnames(bin.coords)@lengths, seq_len)),
                width=nchar(max(chr.cnts)), format='d', flag='0'), sep='-')

    bin.coords
}

#' Bin signals from probe signals
#'
#' require GenomicRanges
#' 
#' @param probe.signals probe signals
#' @param bin.coords bin coordinates
#' @param probe.coords probe coordinates
#' @import methods
#' @return bin signals
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

    bin.signals <- vapply(split(.probe.signals, .bins), median, 1, na.rm=TRUE)
    bin.signals
}

#' Segment bins using DNAcopy
#'
#' @import DNAcopy
#' @param bin.signals bin signals (input)
#' @param bin.coords bin coordinates
#' @return segment signal data frame
segmentBins <- function(bin.signals, bin.coords) {

    pkgTest('DNAcopy')
    bin.coords <- bin.coords[names(bin.signals)]
    
    ## make input data frame
    maplocs <- as.integer((
        GenomicRanges::start(bin.coords) + GenomicRanges::end(bin.coords))/2)
    
    cna <- DNAcopy::CNA(
        genomdat = bin.signals,
        chrom = as.character(GenomicRanges::seqnames(bin.coords)),
        maploc = maplocs,
        data.type = 'logratio')

    seg <- DNAcopy::segment(
        x=cna, min.width = 5, 
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
#' The function takes a \code{CNSegment} object obtained from cnSegmentation
#' and plot the bin signals and segments (as horizontal lines).
#'
#' require ggplot2, scales
#' @param seg a \code{CNSegment} object
#' @param to.plot chromosome to plot (by default plot all chromosomes)
#' @return plot graphics
#' @examples
#'
#' sesameDataCache("EPIC") # in case not done yet
#' sdf <- sesameDataGet('EPIC.1.SigDF')
#' sdfs.normal <- sesameDataGet('EPIC.5.SigDF.normal')
#' seg <- cnSegmentation(sdf, sdfs.normal)
#' 
#' visualizeSegments(seg)
#' 
#' @export
visualizeSegments <- function(seg, to.plot=NULL) {

    stopifnot(is(seg, "CNSegment"))
    pkgTest('ggplot2')
    pkgTest('scales')
    pkgTest('GenomicRanges')
    
    bin.coords <- seg$bin.coords
    bin.seqinfo <- GenomicRanges::seqinfo(bin.coords)
    bin.signals <- seg$bin.signals
    seg.signals <- seg$seg.signals
    total.length <- sum(as.numeric(bin.seqinfo@seqlengths), na.rm=TRUE)
    
    ## skip chromosome too small (e.g, chrM)
    if (is.null(to.plot))
        to.plot <- (bin.seqinfo@seqlengths > total.length*0.01)

    seqlen <- as.numeric(bin.seqinfo@seqlengths[to.plot])
    seq.names <- bin.seqinfo@seqnames[to.plot]
    total.length <- sum(seqlen, na.rm=TRUE)
    seqcumlen <- cumsum(seqlen)
    seqcumlen <- seqcumlen[-length(seqcumlen)]
    seqstart <- setNames(c(0,seqcumlen), seq.names)

    bin.coords <- bin.coords[as.vector(
        GenomicRanges::seqnames(bin.coords)) %in% seq.names]
    bin.signals <- bin.signals[names(bin.coords)]

    GenomicRanges::values(bin.coords)$bin.mids <-
        (GenomicRanges::start(bin.coords) + GenomicRanges::end(bin.coords))/2

    GenomicRanges::values(bin.coords)$bin.x <-
        seqstart[as.character(GenomicRanges::seqnames(
            bin.coords))] + bin.coords$bin.mids

    ## plot bin
    p <- ggplot2::qplot(
        bin.coords$bin.x / total.length,
        bin.signals, color=bin.signals, alpha=I(0.8))

    ## plot segment
    seg.beg <- (seqstart[seg.signals$chrom] +
                    seg.signals$loc.start) / total.length
    seg.end <- (seqstart[seg.signals$chrom] +
                    seg.signals$loc.end) / total.length
    p <- p +
        ggplot2::geom_segment(
            ggplot2::aes(
                x = seg.beg, xend = seg.end,
                y = seg.signals$seg.mean, yend=seg.signals$seg.mean),
            size=1.5, color='blue')

    ## chromosome boundary
    p <- p + ggplot2::geom_vline(
        xintercept=seqstart[-1]/total.length, alpha=I(0.5))
    
    ## chromosome label
    p <- p + ggplot2::scale_x_continuous(
        labels=seq.names, breaks=(seqstart+seqlen/2)/total.length) +
            ggplot2::theme(
                axis.text.x = ggplot2::element_text(angle=90, hjust=0.5))
    
    ## styling
    p <- p + ggplot2::scale_colour_gradient2(
        limits=c(-0.3,0.3), low='red', mid='grey', high='green',
        oob=scales::squish) +
            ggplot2::xlab('') + ggplot2::ylab('') +
                ggplot2::theme(legend.position="none")
    
    p
}
