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
    grid.segments(
    ((GenomicRanges::start(probes) - beg) / (end - beg)), 1,
    ((seq_len(nprobes) - 0.5)/nprobes), 0, draw=FALSE)
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
            sprintf("%s:%d-%d", chrom, beg, end), 0, 0.6,
            just = c('left','bottom'), draw = FALSE),
        grid.rect(
            vapply(cytoBand.target$chromStart,
                function(x) (x-chromBeg)/chromWid, 1),
            0.35,
            (cytoBand.target$chromEnd - cytoBand.target$chromStart)/chromWid,
            0.2, gp = gpar(fill = bandColor, col = bandColor),
            just = c('left','bottom'), draw = FALSE),
        grid.rect(0, 0.35, 1, 0.2, just = c("left", "bottom"), draw = FALSE),
        grid.segments( # projection lines
            x0 = pltx0, y0 = 0.3, x1 = c(0,1), y1 = 0.1, draw = FALSE),
        grid.segments( # sentinel bar
            x0 = pltx0, y0 = 0.3, x1 = pltx0, y1 = 0.6,
            gp = gpar(col = "red"), draw = FALSE))
}

#this is a function to convert probes to loci and also de-duplicate loci (by combining probeID with the genomic loci)
#this is made for the assemble_plots function (the same code was more or less being used twice, so to reduce lines of code, it was made into one function)
rename_probes_to_loci_and_de_duplicate_if_needed function(probes.list.or.betas, probes) {
  #converting probes to genomic loci in probes.list.or.betas 'rownames' //START
  probes.list.or.betas <- GenomicRanges::as.data.frame(probes.list.or.betas)
  
  probes.list.or.betas["probes"]<- rownames(probes.list.or.betas)
  
  probes.df <- GenomicRanges::as.data.frame(probes)
  probes.df["genomic.loci"] <- paste0(probes.df$seqnames, ":", probes.df$start, "-", probes.df$end)
  probes.df["probes"]<- row.names(probes.df)
  probes.df
  
  probes.list.or.betas  <- merge(probes.df[c("genomic.loci", "probes")], probes.list.or.betas, by = "probes")
  
  #below 'if/else statement' will determine if any genomic loci coordinates are repeated, if so, de-duplicate by paste0'ing together probe_ID and genomic loci by an '_' in between and setting that to the row.name, instead of just the genomic.loci
  #for duplicates only, the paste0'd (combined probe_ID and genomic loci coordinates) will be their heatmap "x-axis", otherwise only genomic loci coordinate will be used
  if (any(duplicated(probes.list.or.betas$genomic.loci))) {
    
    probes.list.or.betas[which(duplicated(probes.list.or.betas$genomic.loci)), "redundant_loci"] <- c('TRUE')
    rownames(probes.list.or.betas)[which(probes.list.or.betas$redundant_loci == TRUE)] <- paste0(probes.list.or.betas[which(probes.list.or.betas$redundant_loci == TRUE),"probes"], "_", probes.list.or.betas[which(probes.list.or.betas$redundant_loci == TRUE),"genomic.loci"])
    rownames(probes.list.or.betas)[which(is.na(probes.list.or.betas$redundant_loci))] <- probes.list.or.betas[which(is.na(probes.list.or.betas$redundant_loci)), "genomic.loci"]
    probes.list.or.betas$redundant_loci <- NULL
    
  } else {
    
    row.names(probes.list.or.betas) <- probes.list.or.betas$genomic.loci
  }
  probes.list.or.betas$probes <- NULL
  probes.list.or.betas$genomic.loci <- NULL
  probes.list.or.betas <- as.matrix(probes.list.or.betas)
  #converting probes to genomic loci in probes.list.or.betas 'rownames' //END
  probes.list.or.betas
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
assemble_plots <- function(
    betas, txns, probes, plt.txns, plt.mapLines, plt.cytoband,
    heat.height = NULL, 
    show.probeNames = TRUE, 
    replace.probes.with.loci_and_show.regulatory.features = TRUE,
    show.samples.n = NULL,
    show.sampleNames = TRUE, sample.name.fontsize = 10,
    dmin = 0, dmax = 1) {
#replace.probes.with.loci_and_show.regulatory.features = TRUE is only really meant to be used with (mouse) mm39, since the regulatory features were obtained from Ensembl 108 (mouse) mm39
 if (replace.probes.with.loci_and_show.regulatory.features) {
      #START - Addition by Pratik - bring in regulatory features from mft
      mft <- sesameDataGet(sprintf('%s.%s.manifest', platform, genome))
      probe.list <- data.frame(mft[which(names(mft) %in% rownames(betas))])
      rownames(probe.list) <- probe.list$probeID
      probe.list <- probe.list[row.names(betas),]
      
      cpg.and.regulatory.feature.probe.list <- rename_probes_to_loci_and_de_duplicate_if_needed(probe.list, probes)
      #END - Addition by Pratik - Dynamic loci rearrangement with clustering on heatmaps
        
#START - PRATIK - add colors for features
#tol.rainbow color pallete is color-blind friendly
regulatory_features_4_colors <- pals::tol.rainbow(n=12)[1:4] 
regulatory_features_ctcf_colors <- pals::tol.rainbow(n=12)[5]

#from unique(mft$regulatory_feature)[1:4]
names(regulatory_features_4_colors) <- c("predicted enhancer", "predicted promoter", "open chromatin region", "TF binding site")
#from unique(mft$CTCF_binding_site)[2]
names(regulatory_features_ctcf_colors) <- "CTCF binding site"
      
cpg_island_color <- pals::tol.rainbow(n=12)[6]
#below from unique(mft$CpG_Island)[2]
names(cpg_island_color) <- "CpG Island"
#END - PRATIK - add colors for features
      
      betas <- rename_probes_to_loci_and_de_duplicate_if_needed(betas, probes)
        }
     }
    
  
  
  if (is.null(show.samples.n)) { show.samples.n <- ncol(betas); }
    if (is.null(heat.height) && length(txns) > 0) {
        heat.height <- 10 / length(txns); }
    w <- WGrob(plt.txns, name = 'txn')
    w <- w + WGrob(plt.mapLines, Beneath(pad=0, height=0.1))
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
    if (replace.probes.with.loci_and_show.regulatory.features) {
    w <- w + WLegendV(x = "betas", RightOf("betas"), n.text = 2, "betalegend", decreasing = TRUE)  +
             WColorBarH(probe.list$CpG_Island, Beneath('betas', v.scale.proportional = TRUE), cmp=CMPar(brewer.name= 'Set2'), label = "CpG Island", label.side = 'l', "CpG") +
             WColorBarH(probe.list$regulatory_feature, Beneath(v.scale.proportional = TRUE), cmp=CMPar(brewer.name= 'Paired', label2color = mycolors), label = "Genomic Feature Type", label.side = 'l', 'GFT') +
             WLabel(x= "Beta Values", TopOf("betalegend", pad=0.1), fontsize = 15) +
             WLegendV(x= "CpG", BottomRightOf("betas",just = c('center', 'center'), h.pad = .10, v.pad = -.3), height = rel(la.size), label.fontsize = sample.name.fontsize, yticklabel.pad=0.05) +
             WLegendV(x= "GFT",Beneath(), label.fontsize = sample.name.fontsize, yticklabel.pad=0.05, height = .3)
     w
  }
}
