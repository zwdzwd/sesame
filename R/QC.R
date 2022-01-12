#' An S4 class to hold QC statistics
#'
#' @slot stat a list to store qc stats
setClass("sesameQC", representation(stat="list", group="list"))

#' Convert sesameQC to data.frame
#'
#' @param x a sesameQC object
#' @param ... additional argument to as.data.frame
#' @return a data.frame
#' @rdname as.data.frame-methods
#' @aliases as.data.frame,sesameQC-method
#' @examples
#' as.data.frame(new("sesameQC"))
setMethod("as.data.frame", signature="sesameQC",
    definition = function(x, ...) as.data.frame(x@stat))

#' The display method for sesameQC
#'
#' The function outputs the number of probes in each category and the first
#' few signal measurements.
#'
#' @param object object to be displayed
#' @return None
#' @rdname show-methods
#' @aliases show,sesameQC-method
#' @examples
#' new("sesameQC")
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
                s_display <- sprintf("%1.1f", s_display)
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

#' This function compares the input sample with public data.
#' Only overlapping metrics will be compared.
#'
#' @param sdf a raw (unprocessed) \code{SigDF}
#' @param publicQC output of sesameQC_publicQC, optional
#' @return a sesameQC
#' @examples
#'
#' sesameDataCache("EPIC") # if not done yet
#' sdf <- sesameDataGet('EPIC.1.SigDF')
#' sesameQC_rankStats(sesameQC_calcStats_intens(sdf))
#' 
#' @export
sesameQC_rankStats <- function(qc, publicQC=NULL) {
    if (is.null(publicQC)) {
        publicQC <- sesameQC_publicQC(platform=sdfPlatform(sdf))
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

#' A convenience function to call one or multiple
#' sesameQC_calcStats functions
#'
#' @param sdf a SigDF object
#' @param funs a sesameQC_calcStats_* function or a list of them
#' default to sesameQC_calcStats_detection
#' @examples
#' sesameDataCache("EPIC") # if not done yet
#' sdf <- sesameDataGet('EPIC.1.SigDF')
#' sesameQC_calcStats(sdf)
#' @export
sesameQC_calcStats <- function(sdf, funs = NULL, use_all = FALSE) {
    if (is.null(funs)) {
        if (use_all) {
            funs <- c(
                sesameQC_calcStats_detection,
                sesameQC_calcStats_numProbes,
                sesameQC_calcStats_intens,
                sesameQC_calcStats_channel,
                sesameQC_calcStats_dyeBias,
                sesameQC_calcStats_betas)
        } else {
            funs <- c(sesameQC_calcStats_detection)
        }
    }
    if (is(funs, "function")) { funs <- c(funs) }
    qc <- new("sesameQC")
    for (func in funs) {
        qc <- func(sdf, qc = qc)
    }
    qc
}

#' Generate summary numbers that indicative of experiment quality
#' based on number of probes.
#' 
#' Please provide a raw SigDF(before any preprocessing). Usually
#' directly from readIDATpair.
#' 
#' @param sdf a \code{SigDF} object
#' @param qc existing sesameQC object to add to (optional)
#' @return a sesameQC
#' @examples
#' sesameDataCache("EPIC") # if not done yet
#' sdf <- sesameDataGet('EPIC.1.SigDF')
#' sesameQC_calcStats_numProbes(sdf)
#' @export
sesameQC_calcStats_numProbes <- function(sdf, qc = NULL) {

    group_nm <- "Number of Probes"
    if (is.null(qc)) { s <- list(); g <- list()
    } else { s <- qc@stat; g <- qc@group }
    if (group_nm %in% names(g)) { return(qc); }

    g[[group_nm]] <- c(
        num_probes    = "N. Probes         ",
        num_probes_II = "N. Inf.-II Probes ",
        num_probes_IR = "N. Inf.-I (Red)   ",
        num_probes_IG = "N. Inf.-I (Grn)   ",
        num_probes_cg = "N. Probes (CG)    ",
        num_probes_ch = "N. Probes (CH)    ",
        num_probes_rs = "N. Probes (RS)    ")
    
    s$num_probes = nrow(sdf)
    s$num_probes_II = nrow(InfII(sdf))
    s$num_probes_IR = nrow(InfIR(sdf))
    s$num_probes_IG = nrow(InfIG(sdf))
    new("sesameQC", stat=s, group=g)
}

#' Generate summary numbers that indicative of experiment quality
#' based on intensity.
#' 
#' Please provide a raw SigDF(before any preprocessing). Usually
#' directly from readIDATpair.
#' 
#' @param sdf a \code{SigDF} object
#' @param qc existing sesameQC object to add to (optional)
#' @return a sesameQC
#' @examples
#' sesameDataCache("EPIC") # if not done yet
#' sdf <- sesameDataGet('EPIC.1.SigDF')
#' sesameQC_calcStats_intens(sdf)
#' @export
sesameQC_calcStats_intens <- function(sdf, qc = NULL) {

    group_nm <- "Signal Intensity"
    if (is.null(qc)) { s <- list(); g <- list()
    } else { s <- qc@stat; g <- qc@group }
    if (group_nm %in% names(g)) { return(qc); }

    g[[group_nm]] <- c(
        mean_intensity    = "Mean sig. intensity         ",
        mean_intensity_MU = "Mean sig. intensity (M+U)   ",
        mean_ii           = "Mean sig. intensity (Inf.II)",
        mean_inb_grn      = "Mean sig. intens.(I.Grn IB) ",
        mean_inb_red      = "Mean sig. intens.(I.Red IB) ",
        mean_oob_grn      = "Mean sig. intens.(I.Grn OOB)",
        mean_oob_red      = "Mean sig. intens.(I.Red OOB)")

    dG <- InfIG(sdf); dR <- InfIR(sdf); d2 <- InfII(sdf)
    s$mean_intensity <- meanIntensity(sdf) # excluding type-I out-of-band
    s$mean_intensity_MU <- mean(totalIntensities(sdf), na.rm=TRUE) # M + U
    s$mean_ii <- mean(c(d2$UG,d2$UR), na.rm = TRUE)
    s$mean_inb_grn <- mean(c(dG$MG, dG$UG), na.rm = TRUE)
    s$mean_inb_red <- mean(c(dR$MR, dR$UR), na.rm = TRUE)
    s$mean_oob_grn <- mean(c(dR$MG, dR$UG), na.rm = TRUE)
    s$mean_oob_red <- mean(c(dG$MR, dG$UR), na.rm = TRUE)
    
    new("sesameQC", stat=s, group=g)
}

#' Generate summary numbers that indicative of experiment quality
#' based on channel.
#' 
#' Please provide a raw SigDF(before any preprocessing). Usually
#' directly from readIDATpair.
#' 
#' @param sdf a \code{SigDF} object
#' @param qc existing sesameQC object to add to (optional)
#' @return a sesameQC
#' @examples
#' sesameDataCache("EPIC") # if not done yet
#' sdf <- sesameDataGet('EPIC.1.SigDF')
#' sesameQC_calcStats_channel(sdf)
#' @export
sesameQC_calcStats_channel <- function(sdf, qc = NULL) {

    group_nm <- "Color Channel"
    if (is.null(qc)) { s <- list(); g <- list()
    } else { s <- qc@stat; g <- qc@group }
    if (group_nm %in% names(g)) { return(qc); }

    g[[group_nm]] <- c(
        InfI_switch_R2R = "N. Inf.I Probes Red > Red ",
        InfI_switch_G2G = "N. Inf.I Probes Grn > Grn ",
        InfI_switch_R2G = "N. Inf.I Probes Red > Grn ",
        InfI_switch_G2R = "N. Inf.I Probes Grn > Red ")

    res <- inferInfiniumIChannel(sdf, summary = TRUE)
    for (nm in names(res)) {
        s[[paste0('InfI_switch_', nm)]] <- unname(res[nm])
    }
    new("sesameQC", stat=s, group=g)
}

#' Quantify deviation of dye bias in the high signal range from the
#' global median
#'
#' Positive value indicates augmentation of high-end dye bias over
#' low-end. negative value represents high-end dye bias contradicts
#' that at low-end (a distorted dye bias). Negative distortion score
#' (< -1) suggests low experiment quality. 0 suggests a consistent
#' dye bias at high and low-end.
#'
#' @param sdf a \code{SigDF}
#' @param qc existing sesameQC object to add to (optional)
#' @return a sesameQC
#' @examples
#' sdf <- sesameDataGet('EPIC.1.SigDF')
#' sesameQC_calcStats_dyeBias(sdf)
#' @export
sesameQC_calcStats_dyeBias <- function(sdf, qc = NULL) {

    group_nm <- "Dye Bias"
    if (is.null(qc)) { s <- list(); g <- list()
    } else { s <- qc@stat; g <- qc@group }
    if (group_nm %in% names(g)) { return(qc); }

    g[[group_nm]] <- c(
        medR      = "Median Inf.I Intens. Red           ",
        medG      = "Median Inf.I Intens. Grn           ",
        topR      = "Median of Top 20 Inf.I Intens. Red ",
        topG      = "Median of Top 20 Inf.I Intens. Grn ",
        RGratio   = "Ratio of Red-to-Grn median intens. ",
        RGdistort = "Ratio of top vs. global R/G ratios ")
    
    t1 <- InfI(sdf)
    intens <- totalIntensities(sdf)
    s$medR <- median(sort(intens[t1[t1$col == "R", "Probe_ID"]]))
    s$medG <- median(sort(intens[t1[t1$col == "G", "Probe_ID"]]))
    s$topR <- median(tail(sort(intens[t1[t1$col == "R", "Probe_ID"]]), n=20))
    s$topG <- median(tail(sort(intens[t1[t1$col == "G", "Probe_ID"]]), n=20))
    s$RGratio <- s$medR / s$medG
    s$RGdistort <- log(s$topR / s$topG) / log(s$medR / s$medG)

    new("sesameQC", stat=s, group=g)
}

#' Generate summary numbers that indicative of experiment quality
#' based on detection.
#' 
#' Please provide a raw SigDF(before any preprocessing). Usually
#' directly from readIDATpair.
#' 
#' @param sdf a \code{SigDF} object
#' @param qc existing sesameQC object to add to (optional)
#' @return a sesameQC
#' @examples
#' sesameDataCache("EPIC") # if not done yet
#' sdf <- sesameDataGet('EPIC.1.SigDF')
#' sesameQC_calcStats_detection(sdf)
#' @export
sesameQC_calcStats_detection <- function(sdf, qc = NULL) {

    group_nm <- "Detection"
    if (is.null(qc)) { s <- list(); g <- list()
    } else { s <- qc@stat; g <- qc@group }
    if (group_nm %in% names(g)) { return(qc); }

    g[[group_nm]] <- c(
        num_dt      = "N. Probes w/ Detection Success      ",
        frac_dt     = "% Detection Success                 ",
        num_dt_cg   = "N. Probes w/ Detection Success (CG) ",
        frac_dt_cg  = "% Detection Success (CG)            ",
        num_dt_ch   = "N. Probes w/ Detection Success (CH) ",
        frac_dt_ch  = "% Detection Success (CH)            ",
        num_dt_rs   = "N. Probes w/ Detection Success (RS) ",
        frac_dt_rs  = "% Detection Success (RS)            ")

    pvals <- pOOBAH(sdf, return.pval = TRUE)
    s$num_dt <- sum(pvals <= 0.05)
    s$frac_dt <- s$num_dt / length(pvals)
    for (pt in c('cg','ch','rs')) {
        p1 <- pvals[grep(paste0('^', pt), names(pvals))]
        s[[paste0('num_dt_', pt)]] <- sum(p1 <= 0.05)
        s[[paste0('frac_dt_', pt)]] <- sum(p1 <= 0.05) / length(p1)
    }
    new("sesameQC", stat=s, group=g)
}

#' Generate summary numbers that indicative of experiment quality
#' based on betas.
#' 
#' Please provide a raw SigDF(before any preprocessing). Usually
#' directly from readIDATpair.
#' 
#' @param sdf a \code{SigDF} object
#' @param qc existing sesameQC object to add to (optional)
#' @return a sesameQC
#' @examples
#' sesameDataCache("EPIC") # if not done yet
#' sdf <- sesameDataGet('EPIC.1.SigDF')
#' sesameQC_calcStats_betas(sdf)
#' @export
sesameQC_calcStats_betas <- function(sdf, qc = NULL) {

    group_nm <- "Number of Probes"
    if (is.null(qc)) { s <- list(); g <- list()
    } else { s <- qc@stat; g <- qc@group }
    if (group_nm %in% names(g)) { return(qc); }

    g[[group_nm]] <- c(
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
        frac_na_rs     = "% is.na(Beta) (RS)  ")

    betas <- getBetas(pOOBAH(noob(dyeBiasNL(sdf))))
    s$mean_beta <- mean(betas, na.rm = TRUE)
    s$median_beta <- median(betas, na.rm = TRUE)
    s$frac_unmeth <- sum(betas < 0.3, na.rm = TRUE)/sum(!is.na(betas))*100
    s$frac_meth <- sum(betas > 0.7, na.rm = TRUE)/sum(!is.na(betas))*100
    s$num_na <- sum(is.na(betas))
    s$frac_na <- sum(is.na(betas)) / length(betas)

    for (pt in c('cg','ch','rs')) {
        b1 <- betas[grep(paste0('^', pt), names(betas))]
        s[[paste0('mean_beta_', pt)]] <- mean(b1, na.rm = TRUE)
        s[[paste0('median_beta_', pt)]] <- median(b1, na.rm = TRUE)
        s[[paste0('frac_unmeth_', pt)]] <-
            sum(b1 < 0.3, na.rm = TRUE) / sum(!is.na(b1)) * 100
        s[[paste0('frac_meth_', pt)]] <-
            sum(b1 > 0.7, na.rm = TRUE) / sum(!is.na(b1)) * 100
        s[[paste0('num_na_', pt)]] <- sum(is.na(b1))
        s[[paste0('frac_na_', pt)]] <- sum(is.na(b1)) / length(b1)
    }
    new("sesameQC", stat=s, group=g)
}

#' Plot red-green QQ-Plot using Infinium-I Probes
#'
#' @param sdf a \code{SigDF}
#' @return create a qqplot
#' @examples
#' sesameDataCache("EPIC") # if not done yet
#' sdf <- sesameDataGet('EPIC.1.SigDF')
#' sesameQC_plotRedGrnQQ(sdf)
#' @import graphics
#' @export
sesameQC_plotRedGrnQQ <- function(sdf) {
    dG <- InfIG(noMasked(sdf)); dR <- InfIR(noMasked(sdf))
    m <- max(c(dR$MR,dR$UR,dG$MG,dG$UG), na.rm=TRUE)
    
    qqplot(
        c(dR$MR, dR$UR), c(dG$MG, dG$UG),
        xlab = 'Infinium-I Red Signal', ylab = 'Infinium-I Grn Signal',
        main = 'Red-Green QQ-Plot', cex = 0.5,
        xlim = c(0,m), ylim = c(0,m))
    abline(0,1,lty = 'dashed')
}

#' Plot Total Signal Intensities vs Beta Values
#' This plot is helpful in revealing the extent of signal background
#' and dye bias.
#'
#' @param sdf a \code{SigDF}
#' @param mask whether to remove probes that are masked
#' @param intens.range plot range of signal intensity
#' @param use_max to use max(M,U) or M+U
#' @param ... additional arguments to smoothScatter
#' @return create a total signal intensity vs beta value plot
#' @examples
#' sesameDataCache("EPIC") # if not done yet
#' sdf <- sesameDataGet('EPIC.1.SigDF')
#' sesameQC_plotIntensVsBetas(sdf)
#' @import graphics
#' @importFrom grDevices colorRampPalette
#' @export
sesameQC_plotIntensVsBetas <- function(
    sdf, mask=TRUE, use_max=FALSE, intens.range=c(5,15), ...) {

    if (use_max) {
        df <- signalMU(sdf)
        intens <- setNames(pmax(df$M,df$U), df$Probe_ID)
    } else {
        intens <- totalIntensities(sdf, mask=mask)
    }
    smoothScatter(log2(intens), getBetas(sdf, mask=mask)[names(intens)],
        xlab='Total Intensity (Log2(M+U))',
        ylab=expression(paste(beta, " (DNA methylation Level)")),
        nrpoints=0, 
        colramp=colorRampPalette(c("white","white","lightblue",
            "blue","green","yellow","orange","red","darkred"),
            space = "Lab"), xlim=intens.range, ...)
    abline(h=0.5, lty='dashed')
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

## #' Bar plot of probe detection success rate
## #'
## #' @param sdfs a list of SigDFs
## #' @return a bar plot comparing probe success rates
## #' @examples
## #' sesameDataCache("EPIC") # if not done yet
## #' sdfs <- sesameDataGet("EPIC.5.SigDF.normal")
## #' sesameQC_plotBarDetection(sdfs)
## #' @import ggplot2
## #' @importFrom wheatmap WGG
## #' @export
## sesameQC_plotBarDetection <- function(sdfs) {
##     sample_name <- mean_intensity <- mean_intensity_total <- NULL
##     x <- sesameQC_calcStats(sdfs, sesameQC_calcStats_detection)

##     p1 <- ggplot(x) +
##         geom_bar(aes(sample_name, num_na_cg), stat='identity') +
##         xlab('Sample') + ylab('N. detection failure') +
##         theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1))
    
##     p2 <- ggplot(x) +
##         geom_bar(aes(sample_name, (1-frac_na_cg)*100), stat='identity') +
##         xlab('Sample') + ylab('Detection success (%)') +
##         theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1))
    
##     WGG(p1) + WGG(p2, RightOf())
## }

## #' Bar plot of signal intensities
## #'
## #' @param sdfs a list of SigDFs
## #' @return a bar plot comparing signal intensities
## #' @examples
## #' sesameDataCache("EPIC") # if not done yet
## #' sdfs <- sesameDataGet("EPIC.5.SigDF.normal")
## #' sesameQC_plotBarIntens(sdfs)
## #' @import ggplot2
## #' @importFrom wheatmap WGG
## #' @export
## sesameQC_plotBarIntens <- function(sdfs) {
##     sample_name <- mean_intensity <- mean_intensity_total <- NULL
##     x <- sesameQC_calcStats(sdfs, sesameQC_calcStats_intens)

##     p1 <- ggplot(x) +
##         geom_bar(aes(sample_name, mean_intensity), stat='identity') +
##         xlab('Sample') + ylab('Mean Intensity') +
##         ylim(0, max(x$mean_intensity)*1.2) +
##         theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1))

##     p2 <- ggplot(x) +
##         geom_bar(aes(sample_name, mean_intensity_total), stat='identity') +
##         xlab('Sample') + ylab('Mean M+U Intensity') +
##         ylim(0, max(x$mean_intensity_total)*1.2) +
##         theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1))

##     WGG(p1) + WGG(p2, RightOf())
## }

#' Bar plot of color channel switch
#'
#' @param sdfs a list of SigDFs
#' @return a bar plot comparing color channel switches
#' @examples
#' sesameDataCache("EPIC") # if not done yet
#' sdfs <- sesameDataGet("EPIC.5.SigDF.normal")
#' sesameQC_plotBar(lapply(sdfs, sesameQC_calcStats_detection))
#' @import ggplot2
#' @importFrom wheatmap WGG
#' @export
sesameQC_plotBar <- function(qcs) {
    if (is(qcs, "sesameQC")) { qcs <- list(qcs); }
    df <- do.call(cbind, lapply(qcs, as.data.frame))

    plt <- NULL
    for (x in colnames(df)) {
        if (x == "sample_name") { next; }
        p <- ggplot(df) + geom_bar(aes_string("sample_name", x), stat="identity") + ylab(x)
        if (is.null(plt)) {
            plt <- WGG(p)
        } else {
            plt <- plt + WGG(p, Beneath())
        }
    }
}

## sesameQC_plotHeatSNPs <- function(betas) {
##     vafs <- betas[grep('^rs', rownames(betas)),]
##     WHeatmap(vafs)
## }

#' Retrieve stats of public data of similar nature
#' e.g., tissue, FFPE vs non-FFPE, etc.
#' @param platform HM450, EPIC, MM285 etc.
#' @param tissue optional, blood, buccal, saliva, etc.
#' @param samplePrep optional, fresh, FF, etc.
#' @return a data frame of public data QC stats
#'
#' @examples
#' publicQC <- sesameQC_publicQC()
#' 
#' @export
sesameQC_publicQC <- function(platform = NULL, tissue=NULL, samplePrep=NULL) {
    df <- sesameDataGet('detection.stats')
    if (!is.null(platform) && platform %in% df$Platform) {
        df <- df[df$Platform == platform,]
    }
    if (!is.null(tissue) && tissue %in% df$Tissue) {
        df <- df[df$Tissue==tissue,]
    }
    if (!is.null(samplePrep) && samplePrep %in% df$SamplePrep) {
        df <- df[df$SamplePrep==samplePrep,]
    }
    stopifnot(nrow(df) >= 5) # stop if there are too few number of samples
    df
}
