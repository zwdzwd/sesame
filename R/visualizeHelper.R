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
#' @param platform methylation array platform
#' @param genome genome assembly
#' @param heat.height heatmap height (auto inferred based on rows)
#' @param mapLine.height height of the map lines
#' @param show.probeNames whether to show probe names
#' @param show.loci whether to show genomic loci instead of probe names
#' @param show.samples.n number of samples to show (default: all)
#' @param show.sampleNames whether to show sample names
#' @param sample.name.fontsize sample name font size
#' @param show.cgi whether to show CpG island-context annotations
#' @param show.regulatory whether to show Ensembl regulatory annotations
#' @param show.beta.legend whether to show the beta-value legend
#' @param show.annotation.legend whether to show annotation legends
#' @param dmin data min
#' @param dmax data max
#' @return a grid object
assemble_plots <- function(
    betas, txns, probes, plt.txns, plt.mapLines, plt.cytoband,
    platform, genome,
    heat.height = NULL, mapLine.height = 0.2,
    show.probeNames = TRUE, show.loci = FALSE, show.samples.n = NULL,
    show.sampleNames = TRUE, sample.name.fontsize = 10,
    show.cgi = TRUE, show.regulatory = TRUE,
    show.beta.legend = TRUE, show.annotation.legend = TRUE,
    dmin = 0, dmax = 1) {
    
    if (is.null(show.samples.n)) { show.samples.n <- ncol(betas); }


  ## -----------------------------------------------------------------------
  ## Heatmap annotations
  ## -----------------------------------------------------------------------
  ##
  ## CpG-context and Ensembl regulatory annotations are cached separately.
  ## This allows either annotation type to be requested independently and
  ## avoids requiring biomaRt when only CpG-context annotations are needed.
  
  anno <- NULL
  
  if (show.cgi || show.regulatory) {
    
    message("Heatmap annotations requested.")
    message("  Checking BiocFileCache...")
    
    bfc <- BiocFileCache::BiocFileCache(ask=FALSE)
    
    cgi.anno <- NULL
    regulatory.anno <- NULL
    
    ## -------------------------------------------------------------------
    ## CpG-context annotations
    ## -------------------------------------------------------------------
    
    if (show.cgi) {
      
      title <- paste(
        platform, genome,
        "cgi", "heatmapAnnotations",
        sep=".")
      
      message("  CpG resource: ", title)
      
      cached <- BiocFileCache::bfcquery(
        bfc,
        title,
        field="rname",
        exact=TRUE)
      
      if (nrow(cached)) {
        
        message("  Found CpG annotations in BiocFileCache.")
        
        cache.path <- BiocFileCache::bfcrpath(
          bfc,
          rids=cached$rid[1])
        
        e <- new.env(parent=emptyenv())
        load(cache.path, envir=e)
        cgi.anno <- e$heatmapAnnotations
        
      } else {
        
        message("  CpG annotations not found in BiocFileCache.")
        message("  Building CpG-context annotations...")
        
        cgi.anno <- sesameAnno_buildHeatmapAnnotations(
          platform=platform,
          genome=genome,
          cpg=TRUE,
          regulatory=FALSE)
        
        ## Capture the temporary UCSC source cache record, then remove
        ## that bookkeeping attribute from the final annotation object.
        cpg.source.rid <- attr(
          cgi.anno,
          "cpg_source_rid")
        
        attr(
          cgi.anno,
          "cpg_source_rid") <- NULL
        
        message("  Saving CpG annotations to BiocFileCache.")
        
        cache.path <- BiocFileCache::bfcnew(
          bfc,
          rname=title,
          ext=".rda",
          fname="exact")
        
        heatmapAnnotations <- cgi.anno
        
        save(
          heatmapAnnotations,
          file=cache.path)
        
        rm(heatmapAnnotations)
        
        ## The finished CpG RDA is now persistent, so the raw UCSC
        ## source file is no longer required in BiocFileCache.
        if (!is.null(cpg.source.rid) &&
            !is.na(cpg.source.rid)) {
          
          message(
            "  Removing temporary UCSC source from BiocFileCache.")
          
          BiocFileCache::bfcremove(
            bfc,
            rids=cpg.source.rid)
        }
      }
    }
    
    ## -------------------------------------------------------------------
    ## Ensembl regulatory annotations
    ## -------------------------------------------------------------------
    
    if (show.regulatory) {
      
      title <- paste(
        platform, genome,
        "regulatory", "heatmapAnnotations",
        sep=".")
      
      message("  Regulatory resource: ", title)
      
      cached <- BiocFileCache::bfcquery(
        bfc,
        title,
        field="rname",
        exact=TRUE)
      
      if (nrow(cached)) {
        
        message(
          "  Found regulatory annotations in BiocFileCache.")
        
        cache.path <- BiocFileCache::bfcrpath(
          bfc,
          rids=cached$rid[1])
        
        e <- new.env(parent=emptyenv())
        load(cache.path, envir=e)
        regulatory.anno <- e$heatmapAnnotations
        
      } else {
        
        message(
          "  Regulatory annotations not found in BiocFileCache.")
        message(
          "  Building Ensembl regulatory annotations...")
        
        regulatory.anno <- sesameAnno_buildHeatmapAnnotations(
          platform=platform,
          genome=genome,
          cpg=FALSE,
          regulatory=TRUE)
        
        message(
          "  Saving regulatory annotations to BiocFileCache.")
        
        cache.path <- BiocFileCache::bfcnew(
          bfc,
          rname=title,
          ext=".rda",
          fname="exact")
        
        heatmapAnnotations <- regulatory.anno
        
        save(
          heatmapAnnotations,
          file=cache.path)
        
        rm(heatmapAnnotations)
      }
    }
    
    ## -------------------------------------------------------------------
    ## Combine requested annotations in memory for plotting
    ## -------------------------------------------------------------------
    
    if (!is.null(cgi.anno))
      anno <- cgi.anno
    
    if (!is.null(regulatory.anno)) {
      
      if (is.null(anno)) {
        
        anno <- regulatory.anno
        
      } else {
        
        ## Align regulatory annotations to the CpG annotation table by
        ## Probe_ID. The two resources remain independent on disk and
        ## are combined only for the current plot.
        i <- match(
          anno$Probe_ID,
          regulatory.anno$Probe_ID)
        
        regulatory.cols <- setdiff(
          names(regulatory.anno),
          "Probe_ID")
        
        anno[regulatory.cols] <-
          regulatory.anno[
            i,
            regulatory.cols,
            drop=FALSE]
      }
    }
    
    ## Restrict the platform-wide annotation table to the probes displayed
    ## in this region while preserving their plotted order.
    anno <- anno[
      match(names(probes), anno$Probe_ID),
      ,
      drop=FALSE]
    
    message(
      "  Using annotations for ",
      nrow(anno),
      " plotted probe(s).")
    
  } else {
    
    message("Heatmap annotations not requested.")
  }
  

    if (is.null(heat.height) && length(txns) > 0) {
        heat.height <- 10 / length(txns); }
    w <- WGrob(plt.txns, name = 'txn')
    w <- w + WGrob(plt.mapLines, Beneath(pad=0, height=mapLine.height))
    
    ## Probe names take precedence when both probe names and loci are requested.
    if (show.loci && !show.probeNames) {
      rownames(betas) <- paste0(
        as.character(GenomicRanges::seqnames(probes)),
        ":",
        GenomicRanges::start(probes))
    }
    
    w <- w + WHeatmap(
        t(betas), Beneath(height = heat.height),
        name = 'betas',
        cmp = CMPar(dmin=dmin, dmax=dmax),
        xticklabels = (show.probeNames || show.loci) && is.null(anno),
        xticklabel.rotat = 45,
        yticklabels = show.sampleNames,
        yticklabel.fontsize = sample.name.fontsize,
        yticklabels.n = show.samples.n,
        xticklabels.n = length(probes))

    ## -----------------------------------------------------------------------
    ## Annotations
    ## -----------------------------------------------------------------------
    
    if (!is.null(anno)) {
      
      ## Select which annotation tracks to display.
      ## CpG context and Ensembl regulatory annotations can be enabled
      ## independently.
      cols <- character()
      
      if (show.cgi && "CpG_Location_Type" %in% names(anno)) {
        cols <- c(cols, "CpG_Location_Type")
      }
      
      if (show.regulatory) {
        cols <- c(
          cols,
          grep(
            "^Regulatory_Feature_[0-9]+$",
            names(anno),
            value=TRUE))
      }
      
      ## Fixed CpG-context colors.
      cpg.colors <- c(
        "CpG Island"="#009E73",
        "CpG Shore"="#D55E00",
        "CpG Shelf"="#E69F00",
        "Open Sea"="#0072B2")
      
      ## Fixed Ensembl regulatory-feature colors.
      regulatory.colors <- c(
        "Enhancer"="#8E197D",
        "Promoter"="#49429A",
        "CTCF Binding Site"="#E7298A",
        "Open chromatin"="#66A61E",
        "EMAR"="#E6AB02")
      
      ## Draw each selected annotation as a compact one-row heatmap
      ## directly beneath the beta-value heatmap.
      for (i in seq_along(cols)) {
        
        x <- anno[[cols[i]]]
        
        ## Empty strings indicate that a probe has no regulatory
        ## annotation in this track. Convert them to NA so they appear
        ## white rather than becoming a separate annotation category.
        x[x == ""] <- NA_character_
        
        mat <- matrix(x, nrow=1)
        colnames(mat) <- rownames(betas)
        rownames(mat) <- gsub("_", " ", cols[i])
        
        ## Use the appropriate fixed categorical color mapping.
        cmap <- if (cols[i] == "CpG_Location_Type")
          cpg.colors
        else
          regulatory.colors
        
        w <- w + WHeatmap(
          mat,
          Beneath(),
          name=paste0("anno", i),
          cmp=CMPar(
            label2color=cmap,
            na.color="white"),
          xticklabels=(show.probeNames || show.loci) &&
            i == length(cols),
          xticklabel.rotat=45,
          xticklabel.pad=0.9,
          xticklabels.n=length(probes),
          yticklabels=TRUE,
          yticklabels.n=1)
      }
    }
    
    ## -----------------------------------------------------------------------
    ## Beta-value legend
    ## -----------------------------------------------------------------------
    ##
    ## The beta-value legend is controlled independently from the annotation
    ## tracks and annotation legends.
    
    if (show.beta.legend) {
      
      w <- w +
        WLegendV(
          x="betas",
          RightOf(
            "betas",
            h.scale="betas",
            h.scale.proportional=TRUE,
            pad=.02),
          n.text=2,
          name="betalegend",
          decreasing=TRUE) +
        
        WLabel(
          x="Beta Values",
          TopOf("betalegend"),
          name="betaleglab",
          fontsize=)
    }
    
    ## -----------------------------------------------------------------------
    ## Annotation legends
    ## -----------------------------------------------------------------------
    ##
    ## Annotation legends are shown only for annotation types that are
    ## currently enabled. They can also be disabled altogether using
    ## `show.annotation.legend`.
    
    if (!is.null(anno) && show.annotation.legend &&
        (show.cgi || show.regulatory)) {
      
      ## Keep track of the most recently created legend so that the next
      ## annotation legend can be positioned immediately beneath it.
      last.legend <- NULL
      
      ## -------------------------------------------------------------------
      ## CpG-context legend
      ## -------------------------------------------------------------------
      
      if (show.cgi && "CpG_Location_Type" %in% cols) {
        
        ## Find which displayed annotation heatmap contains CpG context.
        ## This is usually anno1, but determining it from `cols` keeps
        ## the code correct if the displayed tracks change.
        cpg.anno <- paste0(
          "anno",
          match("CpG_Location_Type", cols))
        
        ## If the beta legend exists, place the CpG legend underneath it.
        ## Otherwise, begin the annotation-legend column directly to the
        ## right of the beta heatmap.
        if (show.beta.legend) {
          
          w <- w + WLabel(
            x="CpG Location Type",
            BottomRightOf(
              x="betalegend",
              just=c("center", "top"),
              v.pad=-.1),
            name="cpgleglab",
            fontsize=sample.name.fontsize)
          
        } else {
          
          w <- w + WLabel(
            x="CpG Location Type",
            RightOf(
              "betas",
              h.scale="betas",
              h.scale.proportional=TRUE,
              pad=.02),
            name="cpgleglab",
            fontsize=sample.name.fontsize)
        }
        
        w <- w + WLegendV(
          x=cpg.anno,
          name="cpgleg",
          Beneath(
            "cpgleglab",
            v.scale=cpg.anno,
            v.scale.proportional=TRUE),
          label.fontsize=sample.name.fontsize)
        
        last.legend <- "cpgleg"
      }
      
      ## -------------------------------------------------------------------
      ## Unified Ensembl Regulatory Build legend
      ## -------------------------------------------------------------------
      
      if (show.regulatory) {
        
        ## Identify all regulatory annotation columns in the annotation
        ## table, independent of which Regulatory_Feature_n track each
        ## category was assigned to.
        regulatory.cols <- grep(
          "^Regulatory_Feature_[0-9]+$",
          names(anno),
          value=TRUE)
        
        ## Identify the first displayed regulatory annotation heatmap so the
        ## regulatory legend can use the same vertical scale as an annotation row.
        regulatory.displayed.cols <- grep(
          "^Regulatory_Feature_[0-9]+$",
          cols,
          value=TRUE)
        
        regulatory.anno <- paste0(
          "anno",
          match(regulatory.displayed.cols[1], cols))
        
        ## Collect every regulatory category represented in this region.
        regulatory.types <- unique(unlist(
          anno[regulatory.cols],
          use.names=FALSE))
        
        ## Remove probes without regulatory annotations.
        regulatory.types <- sort(
          regulatory.types[
            !is.na(regulatory.types) &
              regulatory.types != ""])
        
        if (length(regulatory.types)) {
          
          ## Create a zero-height heatmap containing every regulatory
          ## category. It is used only as a common color-mapping source
          ## for WLegendV and does not add a visible annotation row.
          reg.legend.mat <- matrix(
            regulatory.types,
            nrow=1)
          
          w <- w + WHeatmap(
            reg.legend.mat,
            Beneath(height=0),
            name="reglegend.source",
            cmp=CMPar(
              label2color=regulatory.colors,
              na.color="white"),
            xticklabels=FALSE,
            yticklabels=FALSE)
          
          ## Position the regulatory legend beneath the CpG legend when
          ## the CpG legend exists. Otherwise, position it beneath the
          ## beta legend, or directly beside the heatmap when neither
          ## preceding legend exists.
          if (!is.null(last.legend)) {
            
            w <- w + WGrob(
              grid::gList(
                grid::textGrob(
                  "Ensembl Regulatory\nBuild Feature",
                  gp=grid::gpar(
                    fontsize=sample.name.fontsize,
                    lineheight=.8))),
              Beneath(last.legend, pad=.03),
              name="ensregleg")
            
          } else if (show.beta.legend) {
            
            w <- w + WGrob(
              grid::gList(
                grid::textGrob(
                  "Ensembl Regulatory\nBuild Feature",
                  gp=grid::gpar(
                    fontsize=sample.name.fontsize,
                    lineheight=.8))),
              Beneath(
                "betalegend",
                pad=.18,
                v.scale=regulatory.anno,
                v.scale.proportional=TRUE),
              name="ensregleg")
            
          } else {
            
            w <- w + WGrob(
              grid::gList(
                grid::textGrob(
                  "Ensembl Regulatory\nBuild ",
                  gp=grid::gpar(
                    fontsize=sample.name.fontsize,
                    lineheight=.8))),
              RightOf(
                "betas",
                h.scale="betas",
                h.scale.proportional=TRUE,
                pad=.02),
              name="ensregleg")
          }
          
          ## Draw one unified regulatory legend regardless of how many
          ## Regulatory_Feature_n annotation rows are displayed.
          w <- w + WLegendV(
            x="reglegend.source",
            Beneath(
              "ensregleg",
              v.scale=regulatory.anno,
              v.scale.proportional=TRUE),
            name="ensreg",
            label.fontsize=sample.name.fontsize)
          
          last.legend <- "ensreg"
        }
      }
    }
    
    w <- w + WGrob(plt.cytoband, TopOf('txn', height=0.15))
    w
}
