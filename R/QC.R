
#' This function looks at public data of similar nature
#' e.g., tissue, FFPE vs non-FFPE, etc to evaluate the quality
#' of the target data quality
#'
#' @param sset a raw (unprocessed) \code{SigSet}
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
#' sset <- sesameDataGet('EPIC.1.LNCaP')$sset
#' ranks <- qualityRank(sset)
#' 
#' @export
qualityRank <- function(
    sset,
    tissue=NULL,
    samplePrep=NULL,
    raw = FALSE) {

    df <- sesameDataGet('detection.stats')
    df <- df[df$Platform == sset@platform,]
    if (!is.null(tissue))
        df <- df[df$Tissue==tissue,]
    if (!is.null(samplePrep))
        df <- df[df$SamplePrep==samplePrep,]
    stopifnot(nrow(df) >= 5) # stop if there are too few number of samples
    if (raw) return (df);

    mean_intensity <- meanIntensity(sset)
    pvals <- pOOBAH(sset, return.pval = TRUE)
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
#' @param sset a \code{SigSet}
#' @param mask whether to remove probes that are masked
#' @param intens.range plot range of signal intensity
#' @return create a total signal intensity vs beta value plot
#' @examples
#' sesameDataCache("EPIC")
#' sset <- # if not done yet
#' sset <- sesameDataGet('EPIC.1.LNCaP')$sset
#' sesamePlotIntensVsBetas(sset)
#' @import graphics
#' @importFrom grDevices colorRampPalette
#' @export
sesamePlotIntensVsBetas <- function(sset, mask=TRUE, intens.range=c(5,15)) {
    
    intens <- totalIntensities(sset)
    smoothScatter(log2(intens), getBetas(sset, mask=mask)[names(intens)],
        xlab='Total Intensity (Log2(M+U))',
        ylab=expression(paste(beta, " (DNA methylation Level)")),
        nrpoints=0, 
        colramp=colorRampPalette(c("white","white","lightblue",
            "blue","green","yellow","orange","red","darkred"),
            space = "Lab"), xlim=intens.range)
    abline(h=0.5, lty='dashed')
    ## plot envelope lines
    x <- c(seq(1,100,by=1), seq(101,10000,by=100))
    bG <- median(oobG(sset))
    bR <- median(oobR(sset))
    lines(log2(x + bG + bR), (0 + bG) / (0 + bG + x + bR), col='blue')
    lines(log2(x + bG + bR), (x + bG) / (x + bG + 0 + bR), col='blue')
    lines(log2(x + bR + bR), (0 + bR) / (x + bR + 0 + bR), col='red')
    lines(log2(x + bR + bR), (x + bR) / (x + bR + 0 + bR), col='red')
    lines(log2(x + bG + bG), (0 + bG) / (x + bG + 0 + bG), col='green')
    lines(log2(x + bG + bG), (x + bG) / (x + bG + 0 + bG), col='green')
}

#' Plot red-green QQ-Plot using Infinium-I Probes
#'
#' @param sset a \code{SigSet}
#' @return create a qqplot
#' @examples
#' sesameDataCache("EPIC")
#' sset <- # if not done yet
#' sset <- sesameDataGet('EPIC.1.LNCaP')$sset
#' sesamePlotRedGrnQQ(sset)
#' @import graphics
#' @export
sesamePlotRedGrnQQ <- function(sset) {
    m = max(max(slot(sset, 'IR'), na.rm=TRUE),
        max(slot(sset, "IG"), na.rm=TRUE))
    qqplot(
        slot(sset, 'IR'), slot(sset, 'IG'),
        xlab = 'Infinium-I Red Signal', ylab = 'Infinium-I Grn Signal',
        main = 'Red-Green QQ-Plot', cex = 0.5,
        xlim = c(0,m), ylim = c(0,m))
    abline(0,1,lty = 'dashed')
}

#' Generate summary numbers that indicative of experiment quality
#' Please provide a raw sigset (before any preprocessing). Usually
#' directly from readIDATpair
#' 
#' @param sset a \code{SigSet} object
#' @return a sesameQC class object
#' @examples
#' sesameDataCache("EPIC") # if not done yet
#' sset <- sesameDataGet('EPIC.1.LNCaP')$sset
#' sesameQC(sset)
#' @export
sesameQC <- function(sset) {

    qc <- structure(data.frame(), class='sesameQC')
    ## number of type II probes
    qc$num_probes_II <- nrow(II(sset))
    ## number of type I (red channel) probes
    qc$num_probes_IR <- nrow(IR(sset))
    ## number of type I (grn channel) probes
    qc$num_probes_IG <- nrow(IG(sset))
    ## number of all probes
    qc$num_probes_all <- qc$num_probes_II +
        qc$num_probes_IR + qc$num_probes_IG
    qc$mean_ii <- mean(II(sset), na.rm = TRUE)
    
    qc$mean_intensity <- meanIntensity(sset) # excluding type-I out-of-band
    qc$mean_intensity_total <- mean(totalIntensities(sset), na.rm=TRUE) # M + U
    qc$mean_inb_grn <- mean(IG(sset))
    qc$mean_inb_red <- mean(IR(sset))
    qc$mean_oob_grn <- mean(oobG(sset))
    qc$mean_oob_red <- mean(oobR(sset))

    res <- inferTypeIChannel(sset, summary = TRUE)
    for (nm in names(res)) {
        qc[[paste0('InfI_switch_', nm)]] <- unname(res[nm])
    }

    pvals = pOOBAH(sset, return.pval=TRUE)
    betas <- getBetas(addMask(
        dyeBiasNL(noob(sset)), pvals > 0.05))

    qc$num_probes <- length(betas)
    qc$num_na <- sum(is.na(betas))
    num_ok <- qc$num_probes - qc$num_na
    qc$frac_na <- qc$num_na / qc$num_probes
    qc$num_nondt <- sum(pvals > 0.05)
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

## ## inferences
## if (sset@platform %in% c("EPIC","HM450")) {
##     qc$sex <- inferSex(sset)
##     qc$ethnicity <- inferEthnicity(sset)
##     qc$GCT <- bisConversionControl(sset)
##     ## assuming no preprocessing
##     qc$age <- predictAgeHorvath353(getBetas(sset, mask=FALSE))
## } else {
##     qc$sex <- "na"
##     qc$ethnicity <- "na"
##     qc$GCT <- "na"
##     qc$age <- "na"
## }

#' Print sesameQC object
#'
#' @param x a sesameQC object
#' @param ... extra parameter for print
#' @return print sesameQC result on screen
#' @examples
#' sesameDataCache("EPIC") # if not done yet
#' sset <- sesameDataGet('EPIC.1.LNCaP')$sset
#' sesameQC(sset)
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
    cat('=======================\n')
    cat('=      Inferences     =\n')
    cat('=======================\n')
    cat('Sex:                           ', x$sex, '\n')
    cat('Ethnicity:                     ', x$ethnicity, '\n')
    cat('Age:                           ', x$age, '\n')
    cat('Bisulfite Conversion (GCT):    ', x$GCT, '\n')
    
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
#' sset <- sesameDataGet('EPIC.1.LNCaP')$sset
#' qc <- sesameQC(sset)
#' df <- as.data.frame(qc)
#' @method as.data.frame sesameQC
#' @export
as.data.frame.sesameQC <- function(
    x, row.names = NULL, optional = FALSE, ...) {

    class(x) <- NULL;
    as.data.frame(x)
}
