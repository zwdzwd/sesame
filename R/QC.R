#' An S4 class to hold QC statistics
#'
#' @slot stat a list to store qc stats
#' @return sesameQC object
setClass("sesameQC", representation(stat="list", group="list"))

setMethod("as.data.frame", signature="sesameQC",
    definition = function(x, ...) as.data.frame(x@stat))

#' Convert data frame to sesameQC object
#'
#' The function convert a data frame back to a list of sesameQC objects
#' 
#' @param df a publicQC data frame
#' @return a list sesameQC objects
#' @importFrom methods new
dataFrame2sesameQC <- function(df) {
    ## df <- sesameDataGet("MM285.publicQC")
    ## qcs <- dataFrame2sesameQC(df[1:2,])
    groups <- c(
        .setGroup_detection(),
        .setGroup_numProbes(),
        .setGroup_intensity(),
        .setGroup_channel(),
        .setGroup_dyeBias(),
        .setGroup_betas())
    groups <- groups[vapply(groups, function(g) {
        if (all(names(g) %in% colnames(df))) {
            TRUE } else { FALSE }}, logical(1))]
    lapply(seq_len(nrow(df)), function(i) {
        new("sesameQC", group=groups, stat=df[i,])})
}

#' Convert a list of sesameQC to data frame
#'
#' @param qcs sesameQCs
#' @param cols QC columns, use NULL to report all
#' @return a data frame
#' @examples
#' sdf <- sesameDataGet("EPIC.1.SigDF")
#' qcs <- sesameQC_calcStats(sdf, "detection")
#' sesameQCtoDF(qcs)
#' @export
sesameQCtoDF <- function(qcs, cols=c("frac_dt_cg","RGdistort", "RGratio")) {
    if (is.list(qcs)) {
        df <- bind_rows(lapply(qcs, sesameQCtoDF, cols=cols))
    } else if (is(qcs, "sesameQC")) {
        df <- as.data.frame(qcs@stat)
    }
    if (!is.null(cols) && any(colnames(df) %in% cols)) {
        df <- df[, colnames(df) %in% cols, drop=FALSE]
    }
    if (!is.null(names(qcs))) {
        df <- cbind(IDAT=names(qcs), df)
    }
    df
}

## ' The display method for sesameQC
## '
## ' The function outputs the number of probes in each category and the first
## ' few signal measurements.
## '
## ' @param object object to be displayed
## ' @return None
## ' @rdname show-methods
## ' @aliases show,sesameQC-method
setMethod("show", "sesameQC", function(object)  {
    s <- object@stat
    g <- object@group
    
    cat("\n")
    for (gname in names(g)) {
        cat("=====================\n")
        cat("|", gname,"\n")
        cat("=====================\n")
        g1 <- g[[gname]]
        for (metric in names(g1)) {
            s_display <- s[[metric]]
            if (startsWith(metric, "frac_")) {
                s_display <- sprintf("%1.1f %%", s[[metric]] * 100)
            }
            if (is.integer(s_display)) {
                s_display <- sprintf("%d", s_display)
            }
            if (is.numeric(s_display)) {
                s_display <- sprintf("%1.2f", s_display)
            }
            cat(g1[metric], ":", s_display, sprintf("(%s)", metric))
            if (paste0("rank_", metric) %in% names(s)) {
                cat(sprintf(" - Rank %1.1f%% (N=%d)",
                    s[[paste0("rank_", metric)]]*100, s$rankN))
            }
            cat("\n")
        }
        cat("\n")
    }
})

#' Get stat numbers from an sesameQC object
#' 
#' @param qc a sesameQC object
#' @param stat_names which stat(s) to retrieve, default to all.
#' @param drop whether to drop to a string when stats_names has
#' only one element.
#' @return a list of named stats to be retrieved
#' @examples 
#' sdf <- sesameDataGet("EPIC.1.SigDF")
#' qc <- sesameQC_calcStats(sdf, "detection")
#' sesameQC_getStats(qc, "frac_dt")
#' @export
sesameQC_getStats <- function(qc, stat_names = NULL, drop = TRUE) {
    if(is.null(stat_names)) {
        stat_names <- names(qc@stat)
    }
    if (length(stat_names) == 1 && drop) {
        qc@stat[[stat_names]]
    } else {
        qc@stat[stat_names]
    }
}

#' This function compares the input sample with public data.
#' Only overlapping metrics will be compared.
#'
#' @param qc a sesameQC object
#' @param publicQC public QC statistics, filtered from e.g.: EPIC.publicQC,
#' MM285.publicQC and Mammal40.publicQC
#' @param platform EPIC, MM285 or Mammal40, used when publicQC is not given
#' @return a sesameQC
#' @importFrom methods new
#' @examples
#'
#' sesameDataCache() # if not done yet
#' sdf <- sesameDataGet('EPIC.1.SigDF')
#' sesameQC_rankStats(sesameQC_calcStats(sdf, "intensity"))
#' 
#' @export
sesameQC_rankStats <- function(qc, publicQC=NULL, platform="EPIC") {

    if (is.null(publicQC)) {
        publicQC <- sesameDataGet(sprintf("%s.publicQC", platform))
    }
    s <- qc@stat; g <- qc@group
    metrics <- intersect(names(qc@stat), colnames(publicQC))
    if (length(metrics) == 0) { return(qc); }
    ranks <- lapply(metrics, function(mt) {
        ecdf(publicQC[[mt]])(qc@stat[[mt]])
    })
    names(ranks) <- paste0("rank_", metrics)
    s <- c(s, ranks)
    s$rankN <- nrow(publicQC)
    new("sesameQC", stat=s, group=g)
}

#' Calculate QC statistics
#' 
#' It is a function to call one or multiple
#' sesameQC_calcStats functions
#'
#' currently supporting: detection, intensity, numProbes, channel,
#' dyeBias, betas
#'
#' @param sdf a SigDF object
#' @param funs a sesameQC_calcStats_* function or a list of them
#' default to all functions. One can also use a string such as
#' "detection" or c("detection", "intensity") to reduce typing
#' @return a sesameQC object
#' @importFrom methods new
#' @importFrom methods is
#' @examples
#' sesameDataCache() # if not done yet
#' sdf <- sesameDataGet('EPIC.1.SigDF')
#' sesameQC_calcStats(sdf)
#' sesameQC_calcStats(sdf, "detection")
#' sesameQC_calcStats(sdf, c("detection", "channel"))
#' ## retrieve stats as a list
#' sesameQC_getStats(sesameQC_calcStats(sdf, "detection"))
#' ## or as data frames
#' as.data.frame(sesameQC_calcStats(sdf, "detection"))
#' 
#' @export
sesameQC_calcStats <- function(sdf, funs = NULL) {
    if (is.null(funs)) {
        funs <- c(
            sesameQC_calcStats_detection,
            sesameQC_calcStats_intensity,
            sesameQC_calcStats_numProbes,
            sesameQC_calcStats_channel,
            sesameQC_calcStats_dyeBias,
            sesameQC_calcStats_betas)
    }
    if (!is(funs, "list")) { funs <- c(funs) }
    qc <- new("sesameQC")

    for (func in funs) {
        if (is.character(func)) {
            func <- get(paste0("sesameQC_calcStats_",func))
            stopifnot(is(func, "function"))
        }
        qc <- func(sdf, qc = qc)
    }
    qc
}

.setGroup_detection <- function() {
    list("Detection" = c(
        num_dtna    = "N. Probes w/ Missing Raw Intensity  ",
        frac_dtna   = "% Probes w/ Missing Raw Intensity   ",
        num_dt      = "N. Probes w/ Detection Success      ",
        frac_dt     = "% Detection Success                 ",
        num_dt_mk   = "N. Detection Succ. (after masking)  ",
        frac_dt_mk  = "% Detection Succ. (after masking)   ",
        num_dt_cg   = "N. Probes w/ Detection Success (cg) ",
        frac_dt_cg  = "% Detection Success (cg)            ",
        num_dt_ch   = "N. Probes w/ Detection Success (ch) ",
        frac_dt_ch  = "% Detection Success (ch)            ",
        num_dt_rs   = "N. Probes w/ Detection Success (rs) ",
        frac_dt_rs  = "% Detection Success (rs)            "))
}
sesameQC_calcStats_detection <- function(sdf, qc = NULL) {

    g1 <- .setGroup_detection()
    group_nm <- names(g1)[1]
    if (is.null(qc)) { s <- list(); g <- list()
    } else { s <- qc@stat; g <- qc@group }
    if (group_nm %in% names(g)) { return(qc); }
    g[[group_nm]] <- g1[[group_nm]]

    pvals0 <- pOOBAH(sdf, return.pval = TRUE)
    pvals <- na.omit(pvals0)
    s$num_dtna <- sum(is.na(pvals0))
    s$frac_dtna <- s$num_dtna / length(pvals0)
    s$num_dt <- sum(pvals <= 0.05)
    s$frac_dt <- s$num_dt / length(pvals)
    idx_mk <- !is.na(pvals0) & !sdf$mask
    s$num_dt_mk <- sum(pvals0[idx_mk] <= 0.05)
    s$frac_dt_mk <- s$num_dt_mk / sum(idx_mk)
    for (pt in c('cg','ch','rs')) {
        p1 <- pvals[grep(paste0('^', pt), names(pvals))]
        s[[paste0('num_dt_', pt)]] <- sum(p1 <= 0.05)
        s[[paste0('frac_dt_', pt)]] <- sum(p1 <= 0.05) / length(p1)
    }
    new("sesameQC", stat=s, group=g)
}

.setGroup_numProbes <- function() {
    list("Number of Probes" = c(
        num_probes    = "N. Probes         ",
        num_probes_II = "N. Inf.-II Probes ",
        num_probes_IR = "N. Inf.-I (Red)   ",
        num_probes_IG = "N. Inf.-I (Grn)   ",
        num_probes_cg = "N. Probes (CG)    ",
        num_probes_ch = "N. Probes (CH)    ",
        num_probes_rs = "N. Probes (RS)    "))
}
sesameQC_calcStats_numProbes <- function(sdf, qc = NULL) {

    g1 <- .setGroup_numProbes()
    group_nm <- names(g1)[1]
    if (is.null(qc)) { s <- list(); g <- list()
    } else { s <- qc@stat; g <- qc@group }
    if (group_nm %in% names(g)) { return(qc); }
    g[[group_nm]] <- g1[[group_nm]]
    
    s$num_probes <- nrow(sdf)
    s$num_probes_II <- nrow(InfII(sdf))
    s$num_probes_IR <- nrow(InfIR(sdf))
    s$num_probes_IG <- nrow(InfIG(sdf))
    s$num_probes_cg <- sum(startsWith(sdf$Probe_ID,"cg"))
    s$num_probes_ch <- sum(startsWith(sdf$Probe_ID,"ch"))
    s$num_probes_rs <- sum(startsWith(sdf$Probe_ID,"rs"))
    new("sesameQC", stat=s, group=g)
}

.setGroup_intensity <- function() {
    list("Signal Intensity" = c(
        mean_intensity    = "Mean sig. intensity         ",
        mean_intensity_MU = "Mean sig. intensity (M+U)   ",
        mean_ii           = "Mean sig. intensity (Inf.II)",
        mean_inb_grn      = "Mean sig. intens.(I.Grn IB) ",
        mean_inb_red      = "Mean sig. intens.(I.Red IB) ",
        mean_oob_grn      = "Mean sig. intens.(I.Grn OOB)",
        mean_oob_red      = "Mean sig. intens.(I.Red OOB)",
        na_intensity_M    = "N. NA in M (all probes)     ",
        na_intensity_U    = "N. NA in U (all probes)     ",
        na_intensity_ig   = "N. NA in raw intensity (IG) ",
        na_intensity_ir   = "N. NA in raw intensity (IR) ",
        na_intensity_ii   = "N. NA in raw intensity (II) "))
}
sesameQC_calcStats_intensity <- function(sdf, qc = NULL) {

    g1 <- .setGroup_intensity()
    group_nm <- names(g1)[1]
    if (is.null(qc)) { s <- list(); g <- list()
    } else { s <- qc@stat; g <- qc@group }
    if (group_nm %in% names(g)) { return(qc); }
    g[[group_nm]] <- g1[[group_nm]]

    dG <- InfIG(sdf); dR <- InfIR(sdf); d2 <- InfII(sdf)
    s$mean_intensity <- meanIntensity(sdf) # excluding type-I out-of-band
    s$mean_intensity_MU <- mean(totalIntensities(sdf), na.rm=TRUE) # M + U
    s$mean_ii <- mean(c(d2$UG,d2$UR), na.rm = TRUE)
    s$mean_inb_grn <- mean(c(dG$MG, dG$UG), na.rm = TRUE)
    s$mean_inb_red <- mean(c(dR$MR, dR$UR), na.rm = TRUE)
    s$mean_oob_grn <- mean(c(dR$MG, dR$UG), na.rm = TRUE)
    s$mean_oob_red <- mean(c(dG$MR, dG$UR), na.rm = TRUE)
    mu <- signalMU(sdf)
    s$na_intensity_M <- sum(is.na(mu$M))
    s$na_intensity_U <- sum(is.na(mu$U))
    s$na_intensity_ig <- sum(is.na(c(dG$MG, dG$MR, dG$UG, dG$UR)))
    s$na_intensity_ir <- sum(is.na(c(dR$MG, dR$MR, dR$UG, dR$UR)))
    s$na_intensity_ii <- sum(is.na(c(d2$UG, d2$UR)))
    
    new("sesameQC", stat=s, group=g)
}

.setGroup_channel <- function() {
    list("Color Channel" = c(
        InfI_switch_R2R = "N. Inf.I Probes Red -> Red ",
        InfI_switch_G2G = "N. Inf.I Probes Grn -> Grn ",
        InfI_switch_R2G = "N. Inf.I Probes Red -> Grn ",
        InfI_switch_G2R = "N. Inf.I Probes Grn -> Red "))
}
sesameQC_calcStats_channel <- function(sdf, qc = NULL) {

    g1 <- .setGroup_channel()
    group_nm <- names(g1)[1]
    if (is.null(qc)) { s <- list(); g <- list()
    } else { s <- qc@stat; g <- qc@group }
    if (group_nm %in% names(g)) { return(qc); }
    g[[group_nm]] <- g1[[group_nm]]
    
    res <- inferInfiniumIChannel(sdf, summary = TRUE)
    for (nm in names(res)) {
        s[[paste0('InfI_switch_', nm)]] <- unname(res[nm])
    }
    new("sesameQC", stat=s, group=g)
}

.setGroup_dyeBias <- function() {
    list("Dye Bias" = c(
        medR      = "Median Inf.I Intens. Red           ",
        medG      = "Median Inf.I Intens. Grn           ",
        topR      = "Median of Top 20 Inf.I Intens. Red ",
        topG      = "Median of Top 20 Inf.I Intens. Grn ",
        RGratio   = "Ratio of Red-to-Grn median Intens. ",
        RGdistort = "Ratio of Top vs. Global R/G Ratios "))
}
sesameQC_calcStats_dyeBias <- function(sdf, qc = NULL) {

    g1 <- .setGroup_dyeBias()
    group_nm <- names(g1)[1]
    if (is.null(qc)) { s <- list(); g <- list()
    } else { s <- qc@stat; g <- qc@group }
    if (group_nm %in% names(g)) { return(qc); }
    g[[group_nm]] <- g1[[group_nm]]
    
    t1 <- InfI(sdf)
    intens <- totalIntensities(sdf)
    s$medR <- median(sort(intens[t1[t1$col == "R", "Probe_ID"]]))
    s$medG <- median(sort(intens[t1[t1$col == "G", "Probe_ID"]]))
    s$topR <- median(tail(sort(intens[t1[t1$col == "R", "Probe_ID"]]), n=20))
    s$topG <- median(tail(sort(intens[t1[t1$col == "G", "Probe_ID"]]), n=20))
    s$RGratio <- s$medR / s$medG
    s$RGdistort <- (s$topR / s$topG) / (s$medR / s$medG)

    new("sesameQC", stat=s, group=g)
}

.setGroup_betas <- function() {
    list("Beta Value" = c(
        mean_beta      = "Mean Beta           ",
        median_beta    = "Median Beta         ",
        frac_unmeth    = "% Beta < 0.3        ",
        frac_meth      = "% Beta > 0.7        ",
        num_na         = "N. is.na(Beta)      ",
        frac_na        = "% is.na(Beta)       ",
        mean_beta_cg   = "Mean Beta (CG)      ",
        median_beta_cg = "Median Beta (CG)    ",
        frac_unmeth_cg = "% Beta < 0.3 (CG)   ",
        frac_meth_cg   = "% Beta > 0.7 (CG)   ",
        num_na_cg      = "N. is.na(Beta) (CG) ",
        frac_na_cg     = "% is.na(Beta) (CG)  ",
        mean_beta_ch   = "Mean Beta (CH)      ",
        median_beta_ch = "Median Beta (CH)    ",
        frac_unmeth_ch = "% Beta < 0.3 (CH)   ",
        frac_meth_ch   = "% Beta > 0.7 (CH)   ",
        num_na_ch      = "N. is.na(Beta) (CH) ",
        frac_na_ch     = "% is.na(Beta) (CH)  ",
        mean_beta_rs   = "Mean Beta (RS)      ",
        median_beta_rs = "Median Beta (RS)    ",
        frac_unmeth_rs = "% Beta < 0.3 (RS)   ",
        frac_meth_rs   = "% Beta > 0.7 (RS)   ",
        num_na_rs      = "N. is.na(Beta) (RS) ",
        frac_na_rs     = "% is.na(Beta) (RS)  "))
}
sesameQC_calcStats_betas <- function(sdf, qc = NULL) {

    g1 <- .setGroup_betas()
    group_nm <- names(g1)[1]
    if (is.null(qc)) { s <- list(); g <- list()
    } else { s <- qc@stat; g <- qc@group }
    if (group_nm %in% names(g)) { return(qc); }
    g[[group_nm]] <- g1[[group_nm]]
    
    betas <- getBetas(pOOBAH(noob(dyeBiasNL(sdf))))
    s$mean_beta <- mean(betas, na.rm = TRUE)
    s$median_beta <- median(betas, na.rm = TRUE)
    s$frac_unmeth <- sum(betas < 0.3, na.rm = TRUE)/sum(!is.na(betas))
    s$frac_meth <- sum(betas > 0.7, na.rm = TRUE)/sum(!is.na(betas))
    s$num_na <- sum(is.na(betas))
    s$frac_na <- sum(is.na(betas)) / length(betas)

    for (pt in c('cg','ch','rs')) {
        b1 <- betas[grep(paste0('^', pt), names(betas))]
        s[[paste0('mean_beta_', pt)]] <- mean(b1, na.rm = TRUE)
        s[[paste0('median_beta_', pt)]] <- median(b1, na.rm = TRUE)
        s[[paste0('frac_unmeth_', pt)]] <-
            sum(b1 < 0.3, na.rm = TRUE) / sum(!is.na(b1))
        s[[paste0('frac_meth_', pt)]] <-
            sum(b1 > 0.7, na.rm = TRUE) / sum(!is.na(b1))
        s[[paste0('num_na_', pt)]] <- sum(is.na(b1))
        s[[paste0('frac_na_', pt)]] <- sum(is.na(b1)) / length(b1)
    }
    new("sesameQC", stat=s, group=g)
}

#' Plot red-green QQ-Plot using Infinium-I Probes
#'
#' @param sdf a \code{SigDF}
#' @param main plot title
#' @param ... additional options to qqplot
#' @return create a qqplot
#' @import graphics
#' @examples
#' sesameDataCache() # if not done yet
#' sdf <- sesameDataGet('EPIC.1.SigDF')
#' sesameQC_plotRedGrnQQ(sdf)
#' @export
sesameQC_plotRedGrnQQ <- function(sdf, main="R-G QQ Plot", ...) {
    dG <- InfIG(noMasked(sdf)); dR <- InfIR(noMasked(sdf))
    m <- max(c(dR$MR,dR$UR,dG$MG,dG$UG), na.rm=TRUE)
    
    qqplot(
        c(dR$MR, dR$UR), c(dG$MG, dG$UG),
        xlab = 'Infinium-I Red Signal', ylab = 'Infinium-I Grn Signal',
        main = main, cex = 0.5,
        xlim = c(0,m), ylim = c(0,m), ...)
    graphics::abline(0,1,lty = 'dashed')
}

#' Plot betas distinguishing different Infinium chemistries
#'
#' @param sdf SigDF
#' @param prep prep codes to step through
#' @param legend_pos legend position (default: top)
#' @param mar margin of layout when showing steps of prep
#' @param main main title in plots
#' @param ... additional options to plot
#' @import graphics
#' @return create a density plot
#' @examples
#' sdf <- sesameDataGet("EPIC.1.SigDF")
#' sesameQC_plotBetaByDesign(sdf, prep="DB")
#' @export
sesameQC_plotBetaByDesign <- function(
    sdf, prep=NULL, legend_pos="top", mar=c(3,3,1,1), main="", ...) {

    if (!is.null(prep)) {
        par(mfrow=c(nchar(prep)+1,1), mar=mar)
        for (n in c(0,seq_len(nchar(prep)))) {
            sesameQC_plotBetaByDesign(
                prepSesame(sdf, substr(prep,1,n)),
                prep = NULL, legend_pos = legend_pos,
                main=sprintf("%s %s", main, substr(prep,1,n)), ...) }
        return(invisible(NULL)) }
    
    dA <- density(na.omit(getBetas(sdf)))
    dR <- density(na.omit(getBetas(InfIR(sdf))))
    dG <- density(na.omit(getBetas(InfIG(sdf))))
    d2 <- density(na.omit(getBetas(InfII(sdf))))
    
    plot(dA, main=main, ylim=c(0, max(dA$y, dR$y, dG$y, d2$y)), ...)
    lines(dR, col='red')
    lines(dG, col='darkgreen')
    lines(d2, col='blue')
    legend(legend_pos, legend=c(
        "All","Infinium-I Red","Infinium-I Grn","Infinium-II"),
        col=c("black","red","darkgreen","blue"), lty="solid")
}

#' Plot Total Signal Intensities vs Beta Values
#' This plot is helpful in revealing the extent of signal background
#' and dye bias.
#'
#' @param sdf a \code{SigDF}
#' @param mask whether to remove probes that are masked
#' @param intens.range plot range of signal intensity
#' @param use_max to use max(M,U) or M+U
#' @param pal color palette, whiteturbo, whiteblack, whitejet
#' @param ... additional arguments to smoothScatter
#' @return create a total signal intensity vs beta value plot
#' @import graphics
#' @examples
#' sesameDataCache() # if not done yet
#' sdf <- sesameDataGet('EPIC.1.SigDF')
#' sesameQC_plotIntensVsBetas(sdf)
#' @export
sesameQC_plotIntensVsBetas <- function(
    sdf, mask=TRUE, use_max=FALSE, intens.range=c(5,15),
    pal="whiteturbo", ...) {

    if (use_max) {
        df <- signalMU(sdf)
        intens <- setNames(pmax(df$M,df$U), df$Probe_ID)
    } else {
        intens <- totalIntensities(sdf, mask=mask)
    }

    requireNamespace("KernSmooth")
    smoothScatter(log2(intens), getBetas(sdf, mask=mask)[names(intens)],
        xlab='Total Intensity (Log2(M+U))',
        ylab=expression(paste(beta, " (DNA methylation Level)")),
        nrpoints=0, 
        colramp=palgen(pal), xlim=intens.range, ...)
    graphics::abline(h=0.5, lty='dashed')
    ## plot envelope lines
    x <- c(seq(1,100,by=1), seq(101,10000,by=100))
    dG <- InfIG(sdf); dR <- InfIR(sdf)
    bG <- median(c(dR$MG, dR$UG), na.rm=TRUE)
    bR <- median(c(dG$MR, dG$UR), na.rm=TRUE)
    if (use_max) {
        lines(log2(pmax(bG, x + bR)), (0 + bG) / (0 + bG + x + bR), col='blue')
        lines(log2(pmax(x + bG, bR)), (x + bG) / (x + bG + 0 + bR), col='blue')
        lines(log2(pmax(bR, x + bR)), (0 + bR) / (x + bR + 0 + bR), col='red')
        lines(log2(pmax(x + bR, bR)), (x + bR) / (x + bR + 0 + bR), col='red')
        lines(log2(pmax(bG, x + bG)), (0 + bG) / (x + bG + 0 + bG), col='green')
        lines(log2(pmax(x + bG, bG)), (x + bG) / (x + bG + 0 + bG), col='green')
    } else {
        lines(log2(x + bG + bR), (0 + bG) / (0 + bG + x + bR), col='blue')
        lines(log2(x + bG + bR), (x + bG) / (x + bG + 0 + bR), col='blue')
        lines(log2(x + bR + bR), (0 + bR) / (x + bR + 0 + bR), col='red')
        lines(log2(x + bR + bR), (x + bR) / (x + bR + 0 + bR), col='red')
        lines(log2(x + bG + bG), (0 + bG) / (x + bG + 0 + bG), col='green')
        lines(log2(x + bG + bG), (x + bG) / (x + bG + 0 + bG), col='green')
    }
}

#' Bar plots for sesameQC
#'
#' By default, it plots median_beta_cg, median_beta_ch, RGratio,
#' RGdistort, frac_dt
#'
#' @param qcs a list of SigDFs
#' @param keys optional, other key to plot, instead of the default
#' keys can be found in the parenthesis of the print output of each
#' sesameQC output.
#' @return a bar plot comparing different QC metrics
#' @examples
#' sesameDataCache() # if not done yet
#' sdfs <- sesameDataGet("EPIC.5.SigDF.normal")[1:2]
#' sesameQC_plotBar(lapply(sdfs, sesameQC_calcStats, "detection"))
#' @import ggplot2
#' @importFrom methods is
#' @export
sesameQC_plotBar <- function(qcs, keys = NULL) {
    if (is(qcs, "sesameQC")) { qcs <- list(qcs); }
    df <- do.call(rbind, lapply(qcs, function(x) as.data.frame(x@stat)))

    ## set display names
    g <- qcs[[1]]@group
    display_nms <- do.call(c, lapply(names(g), function(gn) { setNames(
        sprintf("%s | %s", gn, str_trim(g[[gn]])), names(g[[gn]]))}))

    if (is.null(keys)) {
        keys <- c("frac_dt", "mean_intensity",
            "median_beta_cg", "median_beta_ch",
            "RGratio", "RGdistort")
        df <- df[, keys[keys %in% colnames(df)], drop=FALSE]
    }
    if (ncol(df) == 0) { stop("There is no QC metrics to plot") }

    df$sample_name <- names(qcs) # infer sample name
    
    plt <- NULL
    for (x in colnames(df)) {
        if (x == "sample_name") { next; }
        p <- ggplot(df) +
            geom_bar(aes_string("sample_name", x), stat="identity") +
            ylab("") + xlab("") + ggtitle(display_nms[x]) +
            theme(axis.text.x = element_text(angle=-90, vjust=0.5, hjust=0))
        ## customization for the most important
        if (x == "frac_dt") {
            p <- p + scale_y_continuous(labels=scales::percent)
        } else if (x == "median_beta_cg" || x == "median_beta_ch") {
            p <- p + ylim(c(0,1))
        }
                
        if (is.null(plt)) {
            plt <- wheatmap::WGG(p)
        } else {
            plt <- plt + wheatmap::WGG(p, Beneath(pad=0))
        }
    }
    plt
}

#' Plot SNP heatmap
#'
#' @param sdfs beta value matrix, row: probes; column: samples
#' @param cluster show clustered heatmap
#' @param filter.nonvariant whether to filter nonvariant (range < 0.3)
#' @return a grid graphics object
#' @examples
#'
#' sdfs <- sesameDataGet("EPIC.5.SigDF.normal")[1:2]
#' plt <- sesameQC_plotHeatSNPs(sdfs, filter.nonvariant = FALSE)
#' @export
sesameQC_plotHeatSNPs <- function(
    sdfs, cluster = TRUE, filter.nonvariant = TRUE) {
    
    afs <- openSesame(sdfs, func = getAFs, mask = FALSE)
    if (cluster) {
        afs <- both.cluster(afs)$mat
    }
    if (filter.nonvariant) {
        rg <- apply(afs, 1, function(x) {
            max(x, na.rm=TRUE) - min(x, na.rm=TRUE)})
        afs <- afs[rg > 0.3,]
    }
    stopifnot(nrow(afs) > 0)
    wheatmap::WHeatmap(afs, xticklabels = TRUE,
        cmp=CMPar(stop.points=c("white", "yellow", "red"), dmin=0, dmax=1)) +
        WCustomize(mar.bottom=0.15)
}

## #' Retrieve stats of public data of similar nature
## #' e.g., tissue, FFPE vs non-FFPE, etc.
## #' @param platform HM450, EPIC, MM285 etc.
## #' @param tissue optional, blood, buccal, saliva, etc.
## #' @param samplePrep optional, fresh, FF, etc.
## #' @return a data frame of public data QC stats
## #'
## #' @examples
## #' publicQC <- sesameQC_publicQC()
## #' 
## #' @export
## sesameQC_publicQC <- function(platform, tissue=NULL, samplePrep=NULL) {
##     df <- sesame
##     df <- sesameDataGet('detection.stats')
##     if (!is.null(platform) && platform %in% df$Platform) {
##         df <- df[df$Platform == platform,]
##     }
##     if (!is.null(tissue) && tissue %in% df$Tissue) {
##         df <- df[df$Tissue==tissue,]
##     }
##     if (!is.null(samplePrep) && samplePrep %in% df$SamplePrep) {
##         df <- df[df$SamplePrep==samplePrep,]
##     }
##     stopifnot(nrow(df) >= 5) # stop if there are too few number of samples
##     df
## }
