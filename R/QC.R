
#' Generate summary numbers that indicative of experiment quality
#'
#' @param sset a \code{SigSet} object
#' @param betas processed beta values
#' @return a sesameQC class object
#' @examples
#' sset <- sesameDataGet('EPIC.1.LNCaP')$sset
#' sesameQC(sset)
#' @export
sesameQC <- function(sset, betas = NULL) {

    qc <- structure(data.frame(), class='sesameQC')
    qc$num_probes_II <- nrow(II(sset))
    qc$num_probes_IR <- nrow(IR(sset))
    qc$num_probes_IG <- nrow(IG(sset))
    qc$num_probes_all <- qc$num_probes_II +
        qc$num_probes_IR + qc$num_probes_IG
    qc$mean_ii <- mean(II(sset), na.rm = TRUE)
    
    qc$mean_intensity <- meanIntensity(sset) # excluding type-I out-of-band
    qc$mean_intensity_total <- mean(totalIntensities(sset)) # M + U
    qc$mean_inb_grn <- mean(IG(sset))
    qc$mean_inb_red <- mean(IR(sset))
    qc$mean_oob_grn <- mean(oobG(sset))
    qc$mean_oob_red <- mean(oobR(sset))
    qc$sex <- inferSex(sset)
    qc$ethnicity <- inferEthnicity(sset)
    qc$GCT <- bisConversionControl(sset)
    res <- inferTypeIChannel(sset, summary = TRUE)
    for (nm in names(res)) {
        qc[[paste0('InfI_switch_', nm)]] <- unname(res[nm])
    }

    if (is.null(betas)) {
        betas <- getBetas(sset)
    }

    qc$num_probes <- length(betas)
    qc$num_na <- sum(is.na(betas))
    num_ok <- qc$num_probes - qc$num_na
    qc$frac_na <- qc$num_na / qc$num_probes * 100
    qc$mean_beta <- mean(betas, na.rm = TRUE)
    qc$median_beta <- median(betas, na.rm = TRUE)
    qc$frac_unmeth <- sum(betas < 0.3, na.rm = TRUE) / num_ok * 100
    qc$frac_meth <- sum(betas > 0.7, na.rm = TRUE) / num_ok * 100

    for (pt in c('cg','ch','rs')) {
        betas1 <- betas[grep(paste0('^', pt), names(betas))]
        num_probes <- length(betas1)
        num_na <- sum(is.na(betas1))
        num_ok <- num_probes - num_na
        frac_na <- num_na / num_probes * 100
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
    qc$age <- predictAgeHorvath353(betas)
    
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
    cat('No. probes                     ', x$num_probes_all, '\n')
    cat('mean (M/U) (in-band InfI):     ', x$mean_intensity, '\n')
    cat('mean (M+U) (in-band InfI):     ', x$mean_intensity_total, '\n')

    cat('\n-- Infinium II --\n')
    cat(sprintf('No. probes:                    %d (%1.3f%%)\n',
        x$num_probes_II, x$num_probes_II / x$num_probes_all * 100))
    cat('Mean Intensity:                ', x$mean_ii, '\n')
    
    cat('\n-- Infinium I (Red) -- \n')
    cat(sprintf('No. probes:                     %d (%1.3f%%)\n',
        x$num_probes_IR, x$num_probes_IR / x$num_probes_all * 100))
    cat('No. Probes Consistent Channel: ', x$InfI_switch_R2R, '\n')
    cat('No. Porbes Swapped Channel:    ', x$InfI_switch_R2G, '\n')
    cat('No. Probes Low Intensity:      ', x$InfI_switch_FailedR, '\n')
    cat('Mean Intensity (in-band):      ', x$mean_inb_red, '\n')
    cat('Mean Intensity (out-of-band):  ', x$mean_oob_red, '\n')
    
    cat('\n-- Infinium I (Grn) -- \n')
    cat(sprintf('No. probes:                     %d (%1.3f%%)\n',
        x$num_probes_IG, x$num_probes_IG / x$num_probes_all * 100))
    cat('No. Probes Consistent Channel: ', x$InfI_switch_G2G, '\n')
    cat('No. Probes Swapped Channel:    ', x$InfI_switch_G2R, '\n')
    cat('No. Probes Low Intensity:      ', x$InfI_switch_FailedG, '\n')
    cat('Mean Intensity (in-band):      ', x$mean_inb_grn, '\n')
    cat('Mean Intensity (out-of-band):  ', x$mean_oob_grn, '\n')
    
    cat('\n')
    cat('=======================\n')
    cat('=      Beta Values    =\n')
    cat('=======================\n')
    cat('No. probes:                    ', x$num_probes, '\n')
    cat('No. probes with NA:            ',
        sprintf('%d (%1.3f%%)\n', x$num_na, x$frac_na))
    cat('Mean Betas:                    ', x$mean_beta, '\n')
    cat('Median Betas:                  ', x$median_beta, '\n')
    cat(sprintf('%% Unmethylated (Beta < 0.3):    %1.3f%%\n', x$Frac_Unmeth))
    cat(sprintf('%% Methylated (Beta > 0.7):      %1.3f%%\n', x$Frac_Meth))

    for (pt in c('cg','ch','rs')) {
        cat(sprintf('\n-- %s probes --\n', pt))
        cat('No. Probes:                    ',
            x[[paste0('num_probes_', pt)]], '\n')
        cat('No. Probes with NA:            ',
            sprintf('%d (%1.3f%%)\n',
                x[[paste0('num_na_', pt)]], x[[paste0('frac_na_', pt)]]))
        cat('Mean Betas:                    ',
            x[[paste0('mean_beta_', pt)]], '\n')
        cat('Median Betas:                  ',
            x[[paste0('median_beta_', pt)]], '\n')
        cat(sprintf('%% Unmethylated (Beta < 0.3):    %1.3f%%\n',
            x[[paste0('frac_unmeth_', pt)]]))
        cat(sprintf('%% Methylated (Beta > 0.7):      %1.3f%%\n',
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
#' @export
as.data.frame.sesameQC <- function(
    x, row.names = NULL, optional = FALSE, ...) {

    class(x) <- NULL;
    as.data.frame(x)
}
