cnv_normal_default <- function(platform) {
    if (platform == "EPICv2") {
        sdfs.normal <- sesameDataGet("EPICv2.8.SigDF")
        sdfs.normal[c("GM12878_206909630042_R08C01","GM12878_206909630040_R03C01")]
    } else if (platform == "EPIC") {
        sdfs.normal <- sesameDataGet("EPIC.5.SigDF.normal")
    } else {
        stop(sprintf(
            "Please provide sdfs.normal=. No default for %s", platform))
    }
}

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
#' \dontrun{
#' sdfs <- sesameDataGet('EPICv2.8.SigDF')
#' sdf <- sdfs[["K562_206909630040_R01C01"]]
#' seg <- cnSegmentation(sdf)
#' seg <- cnSegmentation(sdf, return.probe.signals=TRUE)
#' visualizeSegments(seg)
#' }
#' 
#' @export
cnSegmentation <- function(
    sdf, sdfs.normal=NULL, genomeInfo=NULL,
    probeCoords=NULL, tilewidth=50000, verbose = FALSE,
    return.probe.signals = FALSE) {
    
    stopifnot(is(sdf, "SigDF"))
    platform <- sdfPlatform(sdf, verbose = verbose)
    
    if (is.null(sdfs.normal)) {
        sdfs.normal <- cnv_normal_default(platform)
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
    if (return.probe.signals) {
        probeCoords$cnv <- probe.signals;
        return(probeCoords[seqnames(probeCoords) != "*"]); }

    ## bin signals
    ## fix bin coordinates, TODO: this is too time-consuming
    bin.coords <- getBinCoordinates(seqLength, gapInfo,
        tilewidth = tilewidth, probeCoords)
    bin.signals <- binSignals(probe.signals, bin.coords, probeCoords)

    ## segmentation
    structure(list(
        seg.signals = segmentBins(bin.signals, bin.coords),
        bin.coords = bin.coords,
        bin.signals = bin.signals,
        genomeInfo = genomeInfo), class='CNSegment')
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

cnv_plot_extra <- function(
    seg, genes.to.label, seq.names, seqstart, totlen, p) {

    ## cytoband
    cband <- seg$genomeInfo$cytoBand
    cband <- cband[cband$chrom %in% seq.names,]
    xmins <- (cband$chromStart + seqstart[as.character(cband$chrom)]) / totlen
    xmaxs <- (cband$chromEnd + seqstart[as.character(cband$chrom)]) / totlen
    requireNamespace("pals")
    cband2col <- setNames(
        pals::ocean.gray(10)[seq(9,3)],
        c('stalk', 'gneg', 'gpos25', 'gpos50', 'gpos75', 'gpos100'))
    cband2col['acen'] <- 'red'
    cband2col['gvar'] <- cband2col['gpos75']
    band_color <- cband2col[as.character(cband$gieStain)]
    band_loc <- min(-2, min(seg$bin.signals, na.rm=TRUE)) - 0.6
    p <- p + ggplot2::geom_rect(ggplot2::aes(
        xmin = xmins, xmax = xmaxs, ymin = band_loc-0.5, ymax = band_loc,
        fill=names(band_color))) +
        ggplot2::scale_fill_manual(values=band_color, guide="none")
    
    ## chromosome boundary
    p <- p + ggplot2::geom_vline(xintercept = c(0,seqstart[-1]/totlen, 1),
        linetype="dotted", alpha=I(0.5))
    p <- p + ggplot2::geom_segment(ggplot2::aes(
        x = c(0,seqstart[-1]/totlen, 1), xend = c(0,seqstart[-1]/totlen, 1),
        y = band_loc-0.7, yend = band_loc+0.2), linewidth=0.6, color="black")
    
    ## gene labels
    merged.exons <- lapply(genes.to.label, function(x) {
        target.txns <- seg$genomeInfo$txns[GenomicRanges::mcols(
            seg$genomeInfo$txns)$gene_name == x]
        GenomicRanges::reduce(unlist(target.txns))
    })
    chrm <- vapply(merged.exons, function(x) as.character(
        GenomicRanges::seqnames(x)[1]), character(1))
    pos <- vapply(merged.exons, function(x)
        min(GenomicRanges::start(x)), numeric(1))
    label_xpos <- (seqstart[chrm] + pos) / totlen
    label_ypos <- max(seg$bin.signals, na.rm=TRUE) * 1.1
    p <- p + ggplot2::geom_text(aes(
        label_xpos, label_ypos+0.6, label=genes.to.label), vjust=0)
    p <- p + annotate('point', x = label_xpos, y = label_ypos,
        shape = 25, size = 6, color = '#F2AE30', fill = '#F2AE30')
    p <- p + ggplot2::geom_segment(ggplot2::aes(
        x = label_xpos, xend = label_xpos,
        y = band_loc, yend = label_ypos), color="#F2AE30")
    p
}

#' Visualize segments
#'
#' The function takes a \code{CNSegment} object obtained from cnSegmentation
#' and plot the bin signals and segments (as horizontal lines).
#'
#' require ggplot2, scales
#' @param seg a \code{CNSegment} object
#' @param to.plot chromosome to plot (by default plot all chromosomes)
#' @param genes.to.label gene(s) to label
#' @importFrom GenomicRanges start
#' @importFrom GenomicRanges end
#' @importFrom GenomicRanges seqnames
#' @importFrom GenomicRanges seqinfo
#' @return plot graphics
#' @examples
#'
#' sesameDataCache()
#' \dontrun{
#' sdfs <- sesameDataGet('EPICv2.8.SigDF')
#' sdf <- sdfs[["K562_206909630040_R01C01"]]
#' seg <- cnSegmentation(sdf)
#' seg <- cnSegmentation(sdf, return.probe.signals=TRUE)
#' visualizeSegments(seg)
#' visualizeSegments(seg, to.plot=c("chr9","chr22"))
#' visualizeSegments(seg, genes.to.label=c("ABL1","BCR"))
#' }
#'
#' sesameDataGet_resetEnv()
#' 
#' @export
visualizeSegments <- function(seg, to.plot=NULL, genes.to.label = NULL) {

    stopifnot(is(seg, "CNSegment"))
    bin.coords <- seg$bin.coords
    bin.seqinfo <- seqinfo(bin.coords)
    bin.signals <- seg$bin.signals
    sigs <- seg$seg.signals
    total.length <- sum(as.numeric(bin.seqinfo@seqlengths), na.rm=TRUE)
    
    ## skip chromosome too small (e.g, chrM)
    if (is.null(to.plot)) {
        to.plot <- (bin.seqinfo@seqlengths > total.length*0.01)
    } else {
        to.plot <- seqnames(bin.seqinfo) %in% to.plot
    }
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
    p <- ggplot2::ggplot() + ggplot2::geom_point(
        ggplot2::aes(bin.coords$bin.x / totlen,
            bin.signals, color = bin.signals, alpha = I(0.8)))

    ## plot segment
    seg.beg <- (seqstart[sigs$chrom] + sigs$loc.start) / totlen
    seg.end <- (seqstart[sigs$chrom] + sigs$loc.end) / totlen
    p <- p + ggplot2::geom_segment(ggplot2::aes(x = seg.beg, xend = seg.end,
        y = sigs$seg.mean, yend=sigs$seg.mean), linewidth=1.0, color='blue')

    ## chromosome label
    p <- p + ggplot2::scale_x_continuous(
        labels=seq.names, breaks=(seqstart+seqlen/2)/totlen) +
        ggplot2::theme(axis.text.x=ggplot2::element_text(angle=90, hjust=0.5))

    ## color legend
    p <- p + ggplot2::scale_colour_gradient2(
        limits=c(-0.3,0.3), low='red', mid='grey', high='green',
        oob=scales::squish, guide=guide_legend(title="Log2 Signal Ratio")) +
        ggplot2::xlab('') + ggplot2::ylab('')

    p <- cnv_plot_extra(seg, genes.to.label, seq.names, seqstart, totlen, p)
    
    p + ggplot2::theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())
}
