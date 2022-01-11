#' An S4 class to hold QC statistics
#'
#' @slot stat a list to store qc stats
setClass("sesameQC", representation(stat="list"))

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
#' print(new("sesameQC"))
setMethod("show", "sesameQC", function(object)  {
    s <- object@stat
    cat("\n")
    if ("num_probes" %in% names(s)) {
        cat('=======================\n')
        cat('=   Number of Probes  =\n')
        cat('=======================\n')
        cat('No. probes (num_probes_all)       ', s$num_probes_all, '\n')
        cat(sprintf(
            'No. probes: (num_probes_II)        %d (%1.3f%%)\n',
            s$num_probes_II, s$num_probes_II / s$num_probes_all * 100))
        cat(sprintf(
            'No. probes: (num_probes_IR)        %d (%1.3f%%)\n',
            s$num_probes_IR, s$num_probes_IR / s$num_probes_all * 100))
        cat(sprintf(
            'No. probes: (num_probes_IG)        %d (%1.3f%%)\n',
            s$num_probes_IG, s$num_probes_IG / s$num_probes_all * 100))
        cat('\n')
    }

    if ("mean_intensity" %in% names(s)) {
        cat('=======================\n')
        cat('=      Intensities    =\n')
        cat('=======================\n')
        cat('M/U (mean_intensity):          ', s$mean_intensity, '\n')
        cat('M+U (mean_intensity_total):    ', s$mean_intensity_total,'\n')

        cat('\n-- Infinium II --\n')
        cat('Mean Intensity (mean_ii):      ', s$mean_ii, '\n')
        
        cat('\n-- Infinium I (Red) -- \n')
        cat('Mean Intensity (in-band):      ', s$mean_inb_red, '\n')
        cat('Mean Intensity (out-of-band):  ', s$mean_oob_red, '\n')
        
        cat('\n-- Infinium I (Grn) -- \n')
        cat('Mean Intensity (in-band):      ', s$mean_inb_grn, '\n')
        cat('Mean Intensity (out-of-band):  ', s$mean_oob_grn, '\n')
        cat('\n')
    }

    if ("InfI_switch_R2R" %in% names(s)) {
        cat('===================\n')
        cat('=  Color Channel  =\n')
        cat('===================\n')
        cat('\n-- Infinium I (Red) -- \n')
        cat('No. Probes Consistent Channel: ', s$InfI_switch_R2R, '\n')
        cat('No. Porbes Swapped Channel:    ', s$InfI_switch_R2G, '\n')
        
        cat('\n-- Infinium I (Grn) -- \n')
        cat('No. Probes Consistent Channel: ', s$InfI_switch_G2G, '\n')
        cat('No. Probes Swapped Channel:    ', s$InfI_switch_G2R, '\n')
        cat('\n')
    }

    if ("medR" %in% names(s)) {
        cat('=================\n')
        cat('=    Dye bias   =\n')
        cat('=================\n')
        cat('\n-- Infinium I (Red) -- \n')
        cat('All Probe Median Signal (medR): ', s$medR, '\n')
        cat('Top 20 Probe Median (topR):     ', s$topR, '\n')
        cat('\n-- Infinium I (Grn) -- \n')
        cat('All Probe Median Signal (medG): ', s$medG, '\n')
        cat('Top 20 Probe Median (topG):     ', s$topG, '\n')
        cat('\n-- Dye Bias --\n')
        cat('Median R/G Ratio (RGratio):     ', s$RGratio, '\n')
        cat('Top R/G vs global (RGdistort):  ', s$RGdistort, '\n')
        cat('\n')
    }

    if("num_na" %in% names(s)) {
        cat('===================\n')
        cat('=    Detection    =\n')
        cat('===================\n')
        cat('No. probes w/ NA (num_na):     ',
            sprintf('%d (%1.3f%%)\n', s$num_na, s$frac_na*100))
        cat('No. nondetection (num_nondt):  ',
            sprintf('%d (%1.3f%%)\n', s$num_nondt, s$frac_nondt*100))

        for (pt in c('cg','ch','rs')) {
            cat(sprintf('\n-- %s probes --\n', pt))
            cat('No. Probes:                    ',
                s[[paste0('num_probes_', pt)]], '\n')
            cat('No. Probes with NA:            ',
                sprintf('%d (%1.3f%%)\n',
                    s[[paste0('num_na_', pt)]],
                    s[[paste0('frac_na_', pt)]]*100))
        }
    }

    if("mean_beta" %in% names(s)) {
        cat('=======================\n')
        cat('=      Beta Values    =\n')
        cat('=======================\n')
        cat('Mean Betas:                    ', s$mean_beta, '\n')
        cat('Median Betas:                  ', s$median_beta, '\n')
        cat(sprintf('%% Unmethylated (Beta < 0.3): %1.3f%%\n', s$Frac_Unmeth))
        cat(sprintf('%% Methylated (Beta > 0.7):   %1.3f%%\n', s$Frac_Meth))

        for (pt in c('cg','ch','rs')) {
            cat(sprintf('\n-- %s probes --\n', pt))
            cat('Mean Betas:                    ',
                s[[paste0('mean_beta_', pt)]], '\n')
            cat('Median Betas:                  ',
                s[[paste0('median_beta_', pt)]], '\n')
            cat(sprintf(
                '%% Unmethylated (Beta < 0.3):    %1.3f%%\n',
                s[[paste0('frac_unmeth_', pt)]]))
            cat(sprintf(
                '%% Methylated (Beta > 0.7):      %1.3f%%\n',
                s[[paste0('frac_meth_', pt)]]))
        }
        cat('\n')
    }
})

#' Calculate QC stats on multiple SigDFs
#'
#' @param sdfs a list of SigDFs
#' @param funs a sesameQC_calcStats_* function or a list of them
#' default to sesameQC_calcStats_detection
#' @return a data frame
#' @examples
#' sdfs <- sesameDataGet("EPIC.5.SigDF.normal")
#' sesameQC_calcStats(sdfs)
#' sesameQC_calcStats(sdfs, sesameQC_calcStats_numProbes)
#' sesameQC_calcStats(sdfs,
#'     c(sesameQC_calcStats_numProbes, sesameQC_calcStats_channel))
#' @export
sesameQC_calcStats <- function(sdfs, funs = NULL) {

    if (is.null(funs)) {
        funs <- c(sesameQC_calcStats_detection)
    }

    ## return a data frame on a list of SigDF
    stopifnot(is(sdfs,"list") && is(sdfs[[1]],"SigDF"))
    if (is(funs, "function")) { funs <- c(funs) }
    
    qc <- do.call(cbind, lapply(funs, function(func) {
        do.call(rbind, lapply(sdfs, function(x) {
            as.data.frame(func(x)@stat)
        }))
    }))
    qc$sample_name <- names(sdfs)
    qc
}

## #' Generate summary numbers that indicative of experiment quality
## #' Please provide a raw SigDF(before any preprocessing). Usually
## #' directly from readIDATpair
## #' 
## #' @param sdf a \code{SigDF} object
## #' @return a sesameQC class object
## #' @examples
## #' sesameDataCache("EPIC") # if not done yet
## #' sdf <- sesameDataGet('EPIC.1.SigDF')
## #' sesameQC_calcStats(sdf)
## #' @export
## sesameQC_calcStats <- function(sdf) {
##     if (is(sdf,"list") && is(sdf[[1]],"SigDF")) { # if given a list of SigDF
##         qc <- do.call(rbind, lapply(
##             sdf, function(x) as.data.frame(sesameQC_calcStats(x))))
##         qc$sample_name <- names(sdf)
##         return(qc)
##     }
##     qc <- list()
##     qc <- c(qc, sesameQC_calcStats_numProbes(sdf))
##     qc <- c(qc, sesameQC_calcStats_intens(sdf))
##     qc <- c(qc, sesameQC_calcStats_channel(sdf))
##     qc <- c(qc, sesameQC_calcStats_detection(sdf))
##     qc <- c(qc, sesameQC_calcStats_betas(sdf))
##     invisible(qc)
## }

#' Generate summary numbers that indicative of experiment quality
#' based on number of probes.
#' 
#' Please provide a raw SigDF(before any preprocessing). Usually
#' directly from readIDATpair.
#' 
#' @param sdf a \code{SigDF} object
#' @return a sesameQC
#' @examples
#' sesameDataCache("EPIC") # if not done yet
#' sdf <- sesameDataGet('EPIC.1.SigDF')
#' sesameQC_calcStats_numProbes(sdf)
#' @export
sesameQC_calcStats_numProbes <- function(sdf) {

    s <- list(
        num_probes <- nrow(sdf),
        num_probes_all <- nrow(sdf),
        num_probes_II <- nrow(InfII(sdf)),
        num_probes_IR <- nrow(InfIR(sdf)),
        num_probes_IG <- nrow(InfIG(sdf)))
    new("sesameQC", stat=s)
}

#' Generate summary numbers that indicative of experiment quality
#' based on intensity.
#' 
#' Please provide a raw SigDF(before any preprocessing). Usually
#' directly from readIDATpair.
#' 
#' @param sdf a \code{SigDF} object
#' @return a sesameQC
#' @examples
#' sesameDataCache("EPIC") # if not done yet
#' sdf <- sesameDataGet('EPIC.1.SigDF')
#' sesameQC_calcStats_intens(sdf)
#' @export
sesameQC_calcStats_intens <- function(sdf) {

    s <- list()
    dG <- InfIG(sdf); dR <- InfIR(sdf); d2 <- InfII(sdf)
    s$mean_ii <- mean(c(d2$UG,d2$UR), na.rm = TRUE)
    s$mean_intensity <- meanIntensity(sdf) # excluding type-I out-of-band
    s$mean_intensity_total <- mean(totalIntensities(sdf), na.rm=TRUE) # M + U
    s$mean_inb_grn <- mean(c(dG$MG, dG$UG), na.rm = TRUE)
    s$mean_inb_red <- mean(c(dR$MR, dR$UR), na.rm = TRUE)
    s$mean_oob_red <- mean(c(dG$MR, dG$UR), na.rm = TRUE)
    s$mean_oob_grn <- mean(c(dR$MG, dR$UG), na.rm = TRUE)
    new("sesameQC", stat=s)
}

#' Generate summary numbers that indicative of experiment quality
#' based on channel.
#' 
#' Please provide a raw SigDF(before any preprocessing). Usually
#' directly from readIDATpair.
#' 
#' @param sdf a \code{SigDF} object
#' @return a sesameQC
#' @examples
#' sesameDataCache("EPIC") # if not done yet
#' sdf <- sesameDataGet('EPIC.1.SigDF')
#' sesameQC_calcStats_channel(sdf)
#' @export
sesameQC_calcStats_channel <- function(sdf) {
    
    s <- list()
    res <- inferInfiniumIChannel(sdf, summary = TRUE)
    for (nm in names(res)) {
        s[[paste0('InfI_switch_', nm)]] <- unname(res[nm])
    }
    new("sesameQC", stat=s)
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
#' @return a sesameQC
#' @examples
#' sdf <- sesameDataGet('EPIC.1.SigDF')
#' sesameQC_calcStats_dyeBias(sdf)
#' @export
sesameQC_calcStats_dyeBias <- function(sdf) {
    t1 <- InfI(sdf)
    intens <- totalIntensities(sdf)
    s <- list()
    s$medR <- median(sort(intens[t1[t1$col == "R", "Probe_ID"]]))
    s$medG <- median(sort(intens[t1[t1$col == "G", "Probe_ID"]]))
    s$topR <- median(tail(sort(intens[t1[t1$col == "R", "Probe_ID"]]), n=20))
    s$topG <- median(tail(sort(intens[t1[t1$col == "G", "Probe_ID"]]), n=20))
    s$RGratio <- s$medR / s$medG
    s$RGdistort <- log(s$topR / s$topG) / log(s$medR / s$medG)
    new("sesameQC", stat=s)
}

#' Generate summary numbers that indicative of experiment quality
#' based on detection.
#' 
#' Please provide a raw SigDF(before any preprocessing). Usually
#' directly from readIDATpair.
#' 
#' @param sdf a \code{SigDF} object
#' @return a sesameQC
#' @examples
#' sesameDataCache("EPIC") # if not done yet
#' sdf <- sesameDataGet('EPIC.1.SigDF')
#' sesameQC_calcStats_detection(sdf)
#' @export
sesameQC_calcStats_detection <- function(sdf) {
    
    pvals <- pOOBAH(sdf, return.pval = TRUE)
    s <- list()
    ## should include quality mask, for num_na
    s$num_na <- sum(pvals > 0.05)
    s$frac_na <- s$num_na / length(pvals)
    s$num_nondt <- sum(pvals > 0.05)
    s$frac_nondt <- s$num_nondt / length(pvals)
    for (pt in c('cg','ch','rs')) {
        p1 <- pvals[grep(paste0('^', pt), names(pvals))]
        s[[paste0('num_probes_', pt)]] <- length(p1)
        s[[paste0('num_na_', pt)]] <- sum(p1 > 0.05)
        s[[paste0('frac_na_', pt)]] <- sum(p1 > 0.05) / length(p1)
    }
    new("sesameQC", stat=s)
}

#' Generate summary numbers that indicative of experiment quality
#' based on betas.
#' 
#' Please provide a raw SigDF(before any preprocessing). Usually
#' directly from readIDATpair.
#' 
#' @param sdf a \code{SigDF} object
#' @return a sesameQC
#' @examples
#' sesameDataCache("EPIC") # if not done yet
#' sdf <- sesameDataGet('EPIC.1.SigDF')
#' sesameQC_calcStats_betas(sdf)
#' @export
sesameQC_calcStats_betas <- function(sdf) {
    
    betas <- getBetas(pOOBAH(noob(dyeBiasNL(sdf))))
    s <- list()
    s$mean_beta <- mean(betas, na.rm = TRUE)
    s$median_beta <- median(betas, na.rm = TRUE)
    s$frac_unmeth <- sum(betas < 0.3, na.rm = TRUE)/sum(!is.na(betas))*100
    s$frac_meth <- sum(betas > 0.7, na.rm = TRUE)/sum(!is.na(betas))*100

    for (pt in c('cg','ch','rs')) {
        b1 <- betas[grep(paste0('^', pt), names(betas))]
        s[[paste0('mean_beta_', pt)]] <- mean(b1, na.rm = TRUE)
        s[[paste0('median_beta_', pt)]] <- median(b1, na.rm = TRUE)
        s[[paste0('frac_unmeth_', pt)]] <-
            sum(b1 < 0.3, na.rm = TRUE) / sum(!is.na(b1)) * 100
        s[[paste0('frac_meth_', pt)]] <-
            sum(b1 > 0.7, na.rm = TRUE) / sum(!is.na(b1)) * 100
    }
    new("sesameQC", stat=s)
}

## #' Print sesameQC object
## #'
## #' @param x a sesameQC object
## #' @param ... extra parameter for print
## #' @return print sesameQC result on screen
## #' @examples
## #' sesameDataCache("EPIC") # if not done yet
## #' sdf <- sesameDataGet('EPIC.1.SigDF')
## #' sesameQC_printStats(sesameQC_calcStats(sdf))
## #' @export
## sesameQC_printStats <- function(x, ...) {
## }


## #' Coerce a sesameQC into a dataframe
## #'
## #' @param x          a sesameQC object
## #' @param row.names  see as.data.frame
## #' @param optional   see as.data.frame
## #' @param ...        see as.data.frame
## #' @return           a data.frame
## #' @examples
## #' sesameDataCache("EPIC") # if not done yet
## #' sdf <- sesameDataGet('EPIC.1.SigDF')
## #' qc <- sesameQC(sdf)
## #' df <- as.data.frame(qc)
## #' @method as.data.frame sesameQC
## #' @export
## as.data.frame.sesameQC <- function(
##     x, row.names = NULL, optional = FALSE, ...) {

##     class(x) <- NULL;
##     as.data.frame(x)
## }

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

#' This function compares the input sample with public data
#' in terms of detection success rate.
#'
#' @param sdf a raw (unprocessed) \code{SigDF}
#' @param publicQC output of sesameQC_publicQC, optional
#' @return a fraction to represent the rank of the test stat
#' @examples
#'
#' sesameDataCache("EPIC") # if not done yet
#' sdf <- sesameDataGet('EPIC.1.SigDF')
#' ranks <- sesameQC_rankDetectionSuccess(sdf)
#' 
#' @export
sesameQC_rankDetectionSuccess <- function(sdf, publicQC=NULL) {
    if (is.null(publicQC)) {
        publicQC <- sesameQC_publicQC(platform=sdfPlatform(sdf))
    }
    
    pvals <- pOOBAH(sdf, return.pval = TRUE)
    x <- 1-ecdf(publicQC$frac_nondt)(sum(pvals > 0.05) / length(pvals))
    message(
        sprintf("Probe detection rate beats %1.2f%% of %d public %s data.",
            x * 100, nrow(publicQC), sdfPlatform(sdf)))
    invisible(x)
}

#' This function compares the input sample with public data
#' in terms of mean signal intensity.
#'
#' @param sdf a raw (unprocessed) \code{SigDF}
#' @param publicQC output of sesameQC_publicQC, optional
#' @return a fraction to represent the rank of the test stat
#' @examples
#'
#' sesameDataCache("EPIC") # if not done yet
#' sdf <- sesameDataGet('EPIC.1.SigDF')
#' sesameQC_rankMeanIntensity(sdf)
#' 
#' @export
sesameQC_rankMeanIntensity <- function(sdf, publicQC=NULL) {
    if (is.null(publicQC)) {
        publicQC <- sesameQC_publicQC(platform=sdfPlatform(sdf))
    }
    
    x <- meanIntensity(sdf)
    x <- ecdf(publicQC$mean_intensity)(x)
    message(sprintf(
        "Mean signal intensity beats %1.2f%% of %d public %s data.",
        x * 100, nrow(publicQC), sdfPlatform(sdf)))
    invisible(x)
}
