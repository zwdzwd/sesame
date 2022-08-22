#' get normalization control signal
#' 
#' get normalization control signal from SigDF. 
#' The function optionally takes mean for each channel.
#' 
#' @param sdf a SigDF
#' @param average whether to average
#' @param verbose print more messages
#' @return a data frame of normalization control signals
normControls <- function(sdf, average = FALSE, verbose = FALSE) {
    df <- controls(sdf)
    df <- df[grep('norm(_|\\.)', tolower(df$Type)),]

    ## stop if no control probes
    if (nrow(df) == 0)
        stop("No normalization control probes found!")

    if (sdfPlatform(sdf, verbose = verbose) == 'HM27') {
        df$channel <- ifelse(grepl('norm\\.green', tolower(df$Type)), 'G', 'R')
    } else {
        df$channel <- ifelse(grepl('norm_(c|g)', tolower(df$Type)), 'G', 'R')
    }
    
    if (average) {
        ## this is the old fashioned way
        c(G=mean(df[df$channel=='G','UG'], na.rm=TRUE),
            R=mean(df[df$channel=='R','UR'], na.rm=TRUE))
    } else {
        df
    }
}

#' Correct dye bias in by linear scaling.
#'
#' The function takes a \code{SigDF} as input and scale both the Grn and Red
#' signal to a reference (ref) level. If the reference level is not given, it
#' is set to the mean intensity of all the in-band signals. The function
#' returns a \code{SigDF} with dye bias corrected.
#' 
#' @param sdf a \code{SigDF}
#' @param ref reference signal level
#' @return a normalized \code{SigDF}
#' @examples
#' sesameDataCache() # if not done yet
#' sdf <- sesameDataGet('EPIC.1.SigDF')
#' sdf.db <- dyeBiasCorr(sdf)
#' @export
dyeBiasCorr <- function(sdf, ref=NULL) {

    stopifnot(is(sdf, "SigDF"))
    if (is.null(ref)) {
        ref <- meanIntensity(sdf)
    }

    normctl <- normControls(sdf, average=TRUE)
    fR <- ref/normctl['R']
    fG <- ref/normctl['G']

    sdf$MG <- sdf$MG * fG
    sdf$UG <- sdf$UG * fG
    sdf$MR <- sdf$MR * fR
    sdf$UR <- sdf$UR * fR

    sdf
}

#' Correct dye bias using most balanced sample as the reference
#'
#' The function chose the reference signal level from a list of \code{SigDF}.
#' The chosen sample has the smallest difference in Grn and Red signal
#' intensity as measured using the normalization control probes. In practice,
#' it doesn't matter which sample is chosen as long as the reference level
#' does not deviate much. The function returns a list of \code{SigDF}s with
#' dye bias corrected.
#' 
#' @param sdfs a list of normalized \code{SigDF}s
#' @return a list of normalized \code{SigDF}s
#' @examples
#' sesameDataCache() # if not done yet
#' sdfs <- sesameDataGet('HM450.10.SigDF')[1:2]
#' sdfs.db <- dyeBiasCorrMostBalanced(sdfs)
#' @export
dyeBiasCorrMostBalanced <- function(sdfs) {

    normctls <- vapply(sdfs, normControls, numeric(2), average=TRUE)
    most.balanced <- which.min(abs(normctls['G',] / normctls['R',] - 1))
    ref <- mean(normctls[,most.balanced], na.rm=TRUE)
    lapply(sdfs, function(sdf) dyeBiasCorr(sdf, ref))
}

maskIG <- function(sdf) { # mask IG if the grn channel fails completely
    sdf$mask[sdf$col=="G"] <- TRUE
    sdf
}

#' Dye bias correction by matching green and red to mid point
#'
#' This function compares the Type-I Red probes and Type-I Grn probes and
#' generates and mapping to correct signal of the two channels to the middle.
#' The function takes one single \code{SigDF} and returns a \code{SigDF}
#' with dye bias corrected.
#' 
#' @param sdf a \code{SigDF}
#' @param mask include masked probes in Infinium-I probes. No big difference is
#' noted in practice. More probes are generally better.
#' @param verbose print more messages
#' @return a \code{SigDF} after dye bias correction.
#' @importFrom preprocessCore normalize.quantiles.use.target
#' @importFrom stats approx
#' @examples
#' sesameDataCache() # if not done yet
#' sdf <- sesameDataGet('EPIC.1.SigDF')
#' sdf.db <- dyeBiasNL(sdf)
#' @export
dyeBiasNL <- function(sdf, mask = TRUE, verbose = FALSE) {

    stopifnot(is(sdf, "SigDF"))
    rgdistort <- sesameQC_calcStats(sdf, "dyeBias")@stat$RGdistort
    if (is.na(rgdistort) || rgdistort >10) {
        return(maskIG(sdf)); }
    
    ## we use all Inf-I probes so we capture the entire support range
    if (mask) { dG <- InfIG(sdf); dR <- InfIR(sdf)
    } else { dG <- InfIG(noMasked(sdf)); dR <- InfIR(noMasked(sdf)) }
    IG0 <- c(dG$MG, dG$UG); IR0 <- c(dR$MR, dR$UR)
    
    maxIG <- max(IG0, na.rm = TRUE); minIG <- min(IG0, na.rm = TRUE)
    maxIR <- max(IR0, na.rm = TRUE); minIR <- min(IR0, na.rm = TRUE)

    if (maxIG <= 0 || maxIR <= 0) { return(sdf); }
    IR1 <- sort(as.numeric(IR0))
    IR2 <- sort(as.vector(normalize.quantiles.use.target(
        matrix(IR1), as.vector(IG0))))
    IRmid <- (IR1 + IR2) / 2.0
    maxIRmid <- max(IRmid); minIRmid <- min(IRmid)
    fitfunRed <- function(data) {
        insupp    <- data <= maxIR & data >= minIR & (!is.na(data))
        oversupp  <- data > maxIR & (!is.na(data))
        undersupp <- data < minIR & (!is.na(data))
        data[insupp] <- approx(x=IR1, y=IRmid, xout=data[insupp], ties=mean)$y
        data[oversupp]  <- data[oversupp] - maxIR + maxIRmid
        data[undersupp] <- minIRmid/ minIR * data[undersupp]
        data
    }

    IG1 <- sort(as.numeric(IG0))
    IG2 <- sort(as.vector(normalize.quantiles.use.target(
        matrix(IG1), as.vector(IR0))))
    IGmid <- (IG1 + IG2) / 2.0
    maxIGmid <- max(IGmid); minIGmid <- min(IGmid)
    fitfunGrn <- function(data) {
        insupp    <- data <= maxIG & data >= minIG & (!is.na(data))
        oversupp  <- data > maxIG & (!is.na(data))
        undersupp <- data < minIG & (!is.na(data))
        data[insupp] <- approx(x=IG1, y=IGmid, xout=data[insupp], ties=mean)$y
        data[oversupp]  <- data[oversupp] - maxIG + maxIGmid
        data[undersupp] <- minIGmid/ minIG * data[undersupp]
        data
    }

    sdf$MR <- fitfunRed(sdf$MR); sdf$UR <- fitfunRed(sdf$UR)
    sdf$MG <- fitfunGrn(sdf$MG); sdf$UG <- fitfunGrn(sdf$UG)
    sdf
}

#' @rdname dyeBiasNL
#' @export
#' @examples
#' sdf <- sesameDataGet('EPIC.1.SigDF')
#' sdf <- dyeBiasCorrTypeINorm(sdf)
dyeBiasCorrTypeINorm <- dyeBiasNL

#' Correct dye bias in by linear scaling.
#'
#' The function takes a \code{SigDF} as input and scale both the Grn and Red
#' signal to a reference (ref) level. If the reference level is not given, it
#' is set to the mean intensity of all the in-band signals. The function
#' returns a \code{SigDF} with dye bias corrected.
#' 
#' @param sdf a \code{SigDF}
#' @param ref reference signal level
#' @return a normalized \code{SigDF}
#' @examples
#' sesameDataCache() # if not done yet
#' sdf <- sesameDataGet('EPIC.1.SigDF')
#' sdf.db <- dyeBiasL(sdf)
#' @export
dyeBiasL <- function(sdf, ref=NULL) {

    stopifnot(is(sdf, "SigDF"))
    if (is.null(ref)) {
        ref <- meanIntensity(sdf)
    }

    muR <- signalMU(InfIR(sdf))
    fR <- ref / median(muR$M, muR$U, na.rm = TRUE)
    muG <- signalMU(InfIG(sdf))
    fG <- ref / median(muG$M, muG$U, na.rm = TRUE)
    
    sdf$MG <- sdf$MG * fG
    sdf$UG <- sdf$UG * fG
    sdf$MR <- sdf$MR * fR
    sdf$UR <- sdf$UR * fR

    sdf
}

## the following three functions are retired since 1.9.1
## dyeBiasCorrTypeINormMpU, dyeBiasCorrTypeINormG2R, dyeBiasCorrTypeINormR2G

