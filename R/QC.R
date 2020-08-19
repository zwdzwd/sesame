
#' Generate summary numbers that indicative of experiment quality
#' Please provide a raw sigset (before any preprocessing). Usually
#' directly from readIDATpair
#' 
#' @param sset a \code{SigSet} object
#' @return a sesameQC class object
#' @examples
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
    qc$mean_intensity_total <- mean(totalIntensities(sset)) # M + U
    qc$mean_inb_grn <- mean(IG(sset))
    qc$mean_inb_red <- mean(IR(sset))
    qc$mean_oob_grn <- mean(oobG(sset))
    qc$mean_oob_red <- mean(oobR(sset))

    ## inferences
    if (sset@platform %in% c("EPIC","HM450")) {
        qc$sex <- inferSex(sset)
        qc$ethnicity <- inferEthnicity(sset)
        qc$GCT <- bisConversionControl(sset)
    } else {
        qc$sex <- "na"
        qc$ethnicity <- "na"
        qc$GCT <- "na"
    }
     

    res <- inferTypeIChannel(sset, summary = TRUE)
    for (nm in names(res)) {
        qc[[paste0('InfI_switch_', nm)]] <- unname(res[nm])
    }

    betas <- getBetas(qualityMask(detectionMask(
        dyeBiasCorrTypeINorm(noob(sset)))))

    qc$num_probes <- length(betas)
    qc$num_na <- sum(is.na(betas))
    num_ok <- qc$num_probes - qc$num_na
    qc$frac_na <- qc$num_na / qc$num_probes
    qc$num_nondt <- sum(pval(sset) > 0.05)
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
    ## assuming no preprocessing
    qc$age <- predictAgeHorvath353(getBetas(sset, mask=FALSE))
    
    qc
}

#' Print sesameQC object
#'
#' @param x a sesameQC object
#' @param ... extra parameter for print
#' @return print sesameQC result on screen
#' @examples
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
#' sset <- sesameDataGet('EPIC.1.LNCaP')$sset
#' qc <- sesameQC(sset)
#' df <- as.data.frame(qc)
#' @exportS3Method as.data.frame sesameQC
as.data.frame.sesameQC <- function(
    x, row.names = NULL, optional = FALSE, ...) {

    class(x) <- NULL;
    as.data.frame(x)
}

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
#' # to be added
#' 
#' @export
qualityRank <- function(
    sset,
    tissue=NULL,
    samplePrep=NULL,
    raw = FALSE) {

    df <- sesameDataGet('detection.stats')
    df <- subset(df, Platform == sset@platform)
    if (!is.null(tissue))
        df <- subset(df, Tissue==tissue)
    if (!is.null(samplePrep))
        df <- subset(df, SamplePrep==samplePrep)
    stopifnot(nrow(df) < 5) # stop if there are too few number of samples
    if (raw) return (df);

    mean_intensity <- meanIntensity(sset)
    na_nondt <- sum(sset@extra$pvals[['pOOBAH']] > 0.05)
    c(
        n_compared = nrow(df),
        rank_nondetection = 1-ecdf(df$na_nondt)(na_nondt),
        rank_meanintensity = ecdf(df$mean_intensity)(mean_intensity))
}

