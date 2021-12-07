
#' This function looks at public data of similar nature
#' e.g., tissue, FFPE vs non-FFPE, etc to evaluate the quality
#' of the target data quality
#'
#' @param sdf a raw (unprocessed) \code{SigDF}
#' @param tissue A string (blood,buccal and saliva)
#' @param samplePrep FFPE, fresh, frozen
#' @param raw to return the raw comparison table
#' @return three numbers:
#' 1. The number of public samples compared.
#' 2. The fraction of public samples with more nondetection, and
#' 3. The fraction of public samples with lower mean intensity
#' 4. The higher the fraction, the better the sample.
#' @examples
#'
#' sesameDataCache("EPIC") # if not done yet
#' sdf <- sesameDataGet('EPIC.1.SigDF')
#' ranks <- qualityRank(sdf)
#' 
#' @export
qualityRank <- function(
    sdf,
    tissue=NULL,
    samplePrep=NULL,
    raw = FALSE) {

    df <- sesameDataGet('detection.stats')
    df <- df[df$Platform == sdfPlatform(sdf),]
    if (!is.null(tissue))
        df <- df[df$Tissue==tissue,]
    if (!is.null(samplePrep))
        df <- df[df$SamplePrep==samplePrep,]
    stopifnot(nrow(df) >= 5) # stop if there are too few number of samples
    if (raw) return (df);

    mean_intensity <- meanIntensity(sdf)
    pvals <- pOOBAH(sdf, return.pval = TRUE)
    frac_nondt <- sum(pvals > 0.05) / length(pvals)
    list(
        n_sample_compared = nrow(df),
        rank_probe_success_rate = 1-ecdf(df$frac_nondt)(frac_nondt),
        rank_mean_intensity = ecdf(df$mean_intensity)(mean_intensity))
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
#' sesameDataCache("EPIC")
#' sdf <- # if not done yet
#' sdf <- sesameDataGet('EPIC.1.SigDF')
#' sesamePlotIntensVsBetas(sdf)
#' @import graphics
#' @import KernSmooth
#' @importFrom grDevices colorRampPalette
#' @export
sesamePlotIntensVsBetas <- function(
    sdf, mask=TRUE, use_max=FALSE, intens.range=c(5,15), ...) {

    if (use_max) {
        df = signalMU(sdf)
        intens = setNames(pmax(df$M,df$U), df$Probe_ID)
    } else {
        intens = totalIntensities(sdf, mask=mask)
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
    dG = InfIG(sdf); dR = InfIR(sdf)
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

#' Generate summary numbers that indicative of experiment quality
#' Please provide a raw SigDF(before any preprocessing). Usually
#' directly from readIDATpair
#' 
#' @param sdf a \code{SigDF} object
#' @return a sesameQC class object
#' @examples
#' sesameDataCache("EPIC") # if not done yet
#' sdf <- sesameDataGet('EPIC.1.SigDF')
#' sesameQC(sdf)
#' @export
sesameQC <- function(sdf) {

    qc <- structure(data.frame(), class='sesameQC')
    ## number of type II probes
    qc$num_probes_II <- nrow(InfII(sdf))
    ## number of type I (red channel) probes
    qc$num_probes_IR <- nrow(InfIR(sdf))
    ## number of type I (grn channel) probes
    qc$num_probes_IG <- nrow(InfIG(sdf))
    ## number of all probes
    qc$num_probes_all <- qc$num_probes_II +
        qc$num_probes_IR + qc$num_probes_IG
    dG = InfIG(sdf); dR = InfIR(sdf); d2 = InfII(sdf)
    qc$mean_ii <- mean(c(d2$UG,d2$UR), na.rm = TRUE)
    
    qc$mean_intensity <- meanIntensity(sdf) # excluding type-I out-of-band
    qc$mean_intensity_total <- mean(totalIntensities(sdf), na.rm=TRUE) # M + U
    qc$mean_inb_grn <- mean(c(dG$MG, dG$UG), na.rm = TRUE)
    qc$mean_inb_red <- mean(c(dR$MR, dR$UR), na.rm = TRUE)
    qc$mean_oob_red <- mean(c(dG$MR, dG$UR), na.rm = TRUE)
    qc$mean_oob_grn <- mean(c(dR$MG, dR$UG), na.rm = TRUE)

    res <- inferInfiniumIChannel(sdf, summary = TRUE)
    for (nm in names(res)) {
        qc[[paste0('InfI_switch_', nm)]] <- unname(res[nm])
    }

    pvals = pOOBAH(sdf, return.pval=TRUE)
    betas <- getBetas(addMask(
        dyeBiasNL(noob(sdf)), pvals > 0.05))

    qc$num_probes <- length(betas)
    qc$num_na <- sum(is.na(betas))
    num_ok <- qc$num_probes - qc$num_na
    qc$frac_na <- qc$num_na / qc$num_probes
    qc$num_nondt <- sum(pOOBAH(sdf, return.pval=TRUE) > 0.05)
    qc$frac_nondt <- qc$num_nondt / qc$num_probes
    qc$mean_beta <- mean(betas, na.rm = TRUE)
    qc$median_beta <- median(betas, na.rm = TRUE)
    qc$frac_unmeth <- sum(betas < 0.3, na.rm = TRUE) / num_ok * 100
    qc$frac_meth <- sum(betas > 0.7, na.rm = TRUE) / num_ok * 100

    for (pt in c('cg','ch','rs')) {
        betas1 <- betas[grep(paste0('^', pt), names(betas))]
        num_probes <- length(betas1)
        num_na <- sum(is.na(betas1))
        num_ok <- num_probes - num_na
        frac_na <- num_na / num_probes
        qc[[paste0('num_probes_', pt)]] <- num_probes
        qc[[paste0('num_na_', pt)]] <- num_na
        qc[[paste0('frac_na_', pt)]] <- frac_na
        qc[[paste0('mean_beta_', pt)]] <- mean(betas1, na.rm = TRUE)
        qc[[paste0('median_beta_', pt)]] <- median(betas1, na.rm = TRUE)
        qc[[paste0('frac_unmeth_', pt)]] <-
            sum(betas1 < 0.3, na.rm = TRUE) / num_ok * 100
        qc[[paste0('frac_meth_', pt)]] <-
            sum(betas1 > 0.7, na.rm = TRUE) / num_ok * 100
    }
    
    qc
}

#' Print sesameQC object
#'
#' @param x a sesameQC object
#' @param ... extra parameter for print
#' @return print sesameQC result on screen
#' @examples
#' sesameDataCache("EPIC") # if not done yet
#' sdf <- sesameDataGet('EPIC.1.SigDF')
#' sesameQC(sdf)
#' @export
print.sesameQC <- function(x, ...) {

    cat('\n')
    cat('=======================\n')
    cat('=      Intensities    =\n')
    cat('=======================\n')
    cat('No. probes (num_probes_all)       ', x$num_probes_all, '\n')
    cat('mean (M/U) (mean_intensity):      ', x$mean_intensity, '\n')
    cat('mean (M+U) (mean_intensity_total):', x$mean_intensity_total, '\n')

    cat('\n-- Infinium II --\n')
    cat(sprintf(
        'No. probes: (num_probes_II)        %d (%1.3f%%)\n',
        x$num_probes_II, x$num_probes_II / x$num_probes_all * 100))
    cat('Mean Intensity (mean_ii):         ', x$mean_ii, '\n')
    
    cat('\n-- Infinium I (Red) -- \n')
    cat(sprintf(
        'No. probes: (num_probes_IR)        %d (%1.3f%%)\n',
        x$num_probes_IR, x$num_probes_IR / x$num_probes_all * 100))
    cat('No. Probes Consistent Channel:    ', x$InfI_switch_R2R, '\n')
    cat('No. Porbes Swapped Channel:       ', x$InfI_switch_R2G, '\n')
    cat('No. Probes Low Intensity:         ', x$InfI_switch_FailedR, '\n')
    cat('Mean Intensity (in-band):         ', x$mean_inb_red, '\n')
    cat('Mean Intensity (out-of-band):     ', x$mean_oob_red, '\n')
    
    cat('\n-- Infinium I (Grn) -- \n')
    cat(sprintf('No. probes:                     %d (%1.3f%%)\n',
        x$num_probes_IG, x$num_probes_IG / x$num_probes_all * 100))
    cat('No. Probes Consistent Channel:    ', x$InfI_switch_G2G, '\n')
    cat('No. Probes Swapped Channel:       ', x$InfI_switch_G2R, '\n')
    cat('No. Probes Low Intensity:         ', x$InfI_switch_FailedG, '\n')
    cat('Mean Intensity (in-band):         ', x$mean_inb_grn, '\n')
    cat('Mean Intensity (out-of-band):     ', x$mean_oob_grn, '\n')

    cat('\n')
    cat('=======================\n')
    cat('=    Non-detection    =\n')
    cat('=======================\n')
    cat('No. probes:                       ', x$num_na, '\n')
    cat('No. probes w/ NA (num_na):        ',
        sprintf('%d (%1.3f%%)\n', x$num_na, x$frac_na*100))
    cat('No. nondetection (num_nondt):     ',
        sprintf('%d (%1.3f%%)\n', x$num_nondt, x$frac_nondt*100))

    
    cat('\n')
    cat('=======================\n')
    cat('=      Beta Values    =\n')
    cat('=======================\n')
    cat('Mean Betas:                       ', x$mean_beta, '\n')
    cat('Median Betas:                     ', x$median_beta, '\n')
    cat(sprintf('%% Unmethylated (Beta < 0.3):    %1.3f%%\n', x$Frac_Unmeth))
    cat(sprintf('%% Methylated (Beta > 0.7):      %1.3f%%\n', x$Frac_Meth))

    for (pt in c('cg','ch','rs')) {
        cat(sprintf('\n-- %s probes --\n', pt))
        cat('No. Probes:                       ',
            x[[paste0('num_probes_', pt)]], '\n')
        cat('No. Probes with NA:               ',
            sprintf('%d (%1.3f%%)\n',
                x[[paste0('num_na_', pt)]], x[[paste0('frac_na_', pt)]]*100))
        cat('Mean Betas:                       ',
            x[[paste0('mean_beta_', pt)]], '\n')
        cat('Median Betas:                     ',
            x[[paste0('median_beta_', pt)]], '\n')
        cat(sprintf(
            '%% Unmethylated (Beta < 0.3):       %1.3f%%\n',
            x[[paste0('frac_unmeth_', pt)]]))
        cat(sprintf(
            '%% Methylated (Beta > 0.7):         %1.3f%%\n',
            x[[paste0('frac_meth_', pt)]]))
    }

    cat('\n')
}


#' Coerce a sesameQC into a dataframe
#'
#' @param x          a sesameQC object
#' @param row.names  see as.data.frame
#' @param optional   see as.data.frame
#' @param ...        see as.data.frame
#' @return           a data.frame
#' @examples
#' sesameDataCache("EPIC") # if not done yet
#' sdf <- sesameDataGet('EPIC.1.SigDF')
#' qc <- sesameQC(sdf)
#' df <- as.data.frame(qc)
#' @method as.data.frame sesameQC
#' @export
as.data.frame.sesameQC <- function(
    x, row.names = NULL, optional = FALSE, ...) {

    class(x) <- NULL;
    as.data.frame(x)
}

#' This function looks at public data of similar nature
#' e.g., tissue, FFPE vs non-FFPE, etc to evaluate the quality
#' of the target data quality
#'
#' @param sdf a raw (unprocessed) \code{SigDF}
#' @param tissue A string (blood,buccal and saliva)
#' @param samplePrep FFPE, fresh, frozen
#' @param raw to return the raw comparison table
#' @return three numbers:
#' 1. The number of public samples compared.
#' 2. The fraction of public samples with more nondetection, and
#' 3. The fraction of public samples with lower mean intensity
#' 4. The higher the fraction, the better the sample.
#' @examples
#'
#' sesameDataCache("EPIC") # if not done yet
#' sdf <- sesameDataGet('EPIC.1.SigDF')
#' ranks <- qualityRank(sdf)
#' 
#' @export
qualityRank <- function(
    sdf,
    tissue=NULL,
    samplePrep=NULL,
    raw = FALSE) {

    df <- sesameDataGet('detection.stats')
    df <- df[df$Platform == sdfPlatform(sdf),]
    if (!is.null(tissue))
        df <- df[df$Tissue==tissue,]
    if (!is.null(samplePrep))
        df <- df[df$SamplePrep==samplePrep,]
    stopifnot(nrow(df) >= 5) # stop if there are too few number of samples
    if (raw) return (df);

    mean_intensity <- meanIntensity(sdf)
    pvals <- pOOBAH(sdf, return.pval = TRUE)
    frac_nondt <- sum(pvals > 0.05) / length(pvals)
    list(
        n_sample_compared = nrow(df),
        rank_probe_success_rate = 1-ecdf(df$frac_nondt)(frac_nondt),
        rank_mean_intensity = ecdf(df$mean_intensity)(mean_intensity))
}

