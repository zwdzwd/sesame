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
#' assembly
#' @param genomeInfo the genomeInfo files. The default is retrieved from
#' sesameData. Alternative genomeInfo files can be found at
#' https://github.com/zhou-lab/GenomeInfo
#' @param probeCoords the probe coordinates in the corresponding genome
#' if NULL (default), then the default genome assembly is used.
#' Default genome is given by, e.g., sesameData_check_genome(NULL, "EPIC")
#' For additional mapping, download the GRanges object from
#' http://zwdzwd.github.io/InfiniumAnnotation
#' and provide the following argument
#' ..., probeCoords = sesameAnno_buildManifestGRanges("downloaded_file"),...
#' to this function.
#' @param tilewidth tile width for smoothing
#' @param return.probe.signals return probe-level instead of bin-level signal
#' @param verbose print more messages
#' @return an object of \code{CNSegment}
#' @examples
#'
#' sesameDataCache()
#' 
#' ## sdf <- sesameDataGet('EPIC.1.SigDF')
#' ## sdfs.normal <- sesameDataGet('EPIC.5.SigDF.normal')
#' ## seg <- cnSegmentation(sdf, sdfs.normal)
#' ## probeSignal <- cnSegmentation(sdf, return.probe.signals=TRUE)
#'
#' @export
cnSegmentation <- function(
    sdf, sdfs.normal=NULL, genomeInfo=NULL,
    probeCoords=NULL, tilewidth=50000, verbose = FALSE,
    return.probe.signals = FALSE) {
    
    stopifnot(is(sdf, "SigDF"))
    platform <- sdfPlatform(sdf, verbose = verbose)
    
    if (is.null(sdfs.normal)) {
        if (platform == "EPIC") {
            sdfs.normal <- sesameDataGet("EPIC.5.SigDF.normal")
        } else {
            stop(sprintf(
                "Please provide sdfs.normal=. No default for %s", platform))
        }
    }
    
    if (is.null(genomeInfo)) { # genome/chromosome info
        genome <- sesameData_check_genome(NULL, platform)
        genomeInfo <- sesameData_getGenomeInfo(genome)
    }
    if (is.null(probeCoords)) { # probe coordinates
        genome <- sesameData_check_genome(NULL, platform)
        probeCoords <- sesameData_getManifestGRanges(platform, genome = genome)
    }
    seqLength <- genomeInfo$seqLength
    gapInfo <- genomeInfo$gapInfo
    
    ## extract intensities
    target.intens <- totalIntensities(sdf)
    normal.intens <- do.call(cbind, lapply(sdfs.normal, function(sdf) {
        totalIntensities(sdf) }))

    ## find overlapping probes
    target.intens <- na.omit(target.intens)
    pb <- intersect(rownames(normal.intens), names(target.intens))
    pb <- intersect(names(probeCoords), pb)
    target.intens <- target.intens[pb]
    normal.intens <- normal.intens[pb,]
    probeCoords <- probeCoords[pb]
    
    ## normalize probes intensities
    fit <- lm(y~., data=data.frame(y=target.intens, X=normal.intens))
    probe.signals <- setNames(log2(target.intens / pmax(predict(fit), 1)), pb)
    if (return.probe.signals) { return(probe.signals); }

    ## bin signals
    ## fix bin coordinates, TODO: this is too time-consuming
    bin.coords <- getBinCoordinates(seqLength, gapInfo,
        tilewidth = tilewidth, probeCoords)
    bin.signals <- binSignals(probe.signals, bin.coords, probeCoords)

    ## segmentation
    structure(list(
        seg.signals = segmentBins(bin.signals, bin.coords),
        bin.coords = bin.coords,
        bin.signals = bin.signals), class='CNSegment')
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
#' @param probeCoords probe coordinates
#' @param tilewidth tile width for smoothing
#' @return bin.coords
getBinCoordinates <- function(
    seqLength, gapInfo, tilewidth=50000, probeCoords) {

    tiles <- sort(GenomicRanges::tileGenome(
        seqLength, tilewidth=tilewidth, cut.last.tile.in.chrom = TRUE))
    
    tiles <- sort(c(
        GenomicRanges::setdiff(tiles[seq(1, length(tiles), 2)], gapInfo), 
        GenomicRanges::setdiff(tiles[seq(2, length(tiles), 2)], gapInfo)))
    
    GenomicRanges::values(tiles)$probes <-
        GenomicRanges::countOverlaps(tiles, probeCoords)
    
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
#' @param probeCoords probe coordinates
#' @importFrom methods .hasSlot
#' @return bin signals
binSignals <- function(probe.signals, bin.coords, probeCoords) {
    ov <- GenomicRanges::findOverlaps(probeCoords, bin.coords)
    if (.hasSlot(ov, 'queryHits')) {
        .bins <- names(bin.coords)[ov@subjectHits]
        .probe.signals <- probe.signals[names(probeCoords)[ov@queryHits]]
    } else {
        .bins <- names(bin.coords)[ov@to]
        .probe.signals <- probe.signals[names(probeCoords)[ov@from]]
    }

    bin.signals <- vapply(split(.probe.signals, .bins), median, 1, na.rm=TRUE)
    bin.signals
}

#' Segment bins using DNAcopy
#'
#' @param bin.signals bin signals (input)
#' @param bin.coords bin coordinates
#' @return segment signal data frame
segmentBins <- function(bin.signals, bin.coords) {

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
        x = cna, min.width = 5,
        nperm = 10000, alpha = 0.001, undo.splits = 'sdundo',
        undo.SD = 2.2, verbose=0)
    
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
#' @importFrom GenomicRanges start
#' @importFrom GenomicRanges end
#' @importFrom GenomicRanges seqnames
#' @importFrom GenomicRanges seqinfo
#' @return plot graphics
#' @examples
#'
#' sesameDataCache()
#' ## sdf <- sesameDataGet('EPIC.1.SigDF')
#' ## sdfs.normal <- sesameDataGet('EPIC.5.SigDF.normal')
#' ## seg <- cnSegmentation(sdf, sdfs.normal)
#' ## visualizeSegments(seg)
#'
#' sesameDataGet_resetEnv()
#' 
#' @export
visualizeSegments <- function(seg, to.plot=NULL) {

    stopifnot(is(seg, "CNSegment"))
    bin.coords <- seg$bin.coords
    bin.seqinfo <- seqinfo(bin.coords)
    bin.signals <- seg$bin.signals
    sigs <- seg$seg.signals
    total.length <- sum(as.numeric(bin.seqinfo@seqlengths), na.rm=TRUE)
    
    ## skip chromosome too small (e.g, chrM)
    if (is.null(to.plot)) {
        to.plot <- (bin.seqinfo@seqlengths > total.length*0.01) }

    seqlen <- as.numeric(bin.seqinfo@seqlengths[to.plot])
    seq.names <- bin.seqinfo@seqnames[to.plot]
    totlen <- sum(seqlen, na.rm=TRUE)
    seqcumlen <- cumsum(seqlen)
    seqstart <- setNames(c(0,seqcumlen[-length(seqcumlen)]), seq.names)
    bin.coords <- bin.coords[as.vector(seqnames(bin.coords)) %in% seq.names]
    bin.signals <- bin.signals[names(bin.coords)]

    GenomicRanges::values(bin.coords)$bin.mids <-
        (start(bin.coords) + end(bin.coords)) / 2
    GenomicRanges::values(bin.coords)$bin.x <-
        seqstart[as.character(seqnames(bin.coords))] + bin.coords$bin.mids

    ## plot bin
    p <- ggplot2::qplot(bin.coords$bin.x / totlen,
        bin.signals, color=bin.signals, alpha=I(0.8))

    ## plot segment
    seg.beg <- (seqstart[sigs$chrom] + sigs$loc.start) / totlen
    seg.end <- (seqstart[sigs$chrom] + sigs$loc.end) / totlen
    p <- p + ggplot2::geom_segment(ggplot2::aes(x = seg.beg, xend = seg.end,
        y = sigs$seg.mean, yend=sigs$seg.mean), size=1.0, color='blue')

    ## chromosome boundary
    p <- p + ggplot2::geom_vline(xintercept=seqstart[-1]/totlen,
        linetype="dotted", alpha=I(0.5))

    ## chromosome label
    p <- p + ggplot2::scale_x_continuous(
        labels=seq.names, breaks=(seqstart+seqlen/2)/totlen) +
        ggplot2::theme(axis.text.x=ggplot2::element_text(angle=90, hjust=0.5))

    p <- p + ggplot2::scale_colour_gradient2(
        limits=c(-0.3,0.3), low='red', mid='grey', high='green',
        oob=scales::squish) + ggplot2::xlab('') + ggplot2::ylab('') +
        ggplot2::theme(legend.position="none")
    p
}
