#' get normalization control signal
#' 
#' get normalization control signal from SigDF. 
#' The function optionally takes mean for each channel.
#' 
#' @param sdf a SigDF
#' @param average whether to average
#' @return a data frame of normalization control signals
#' @examples 
#' sdf <- readIDATpair(file.path(system.file(
#'     'extdata','',package='sesameData'), '4207113116_B'))
#' 
#' df.ctl <- getNormCtls(sdf)
#' 
#' @export
getNormCtls <- function(sdf, average = FALSE) {
    df <- controls(sdf)
    df <- df[grep('norm(_|\\.)', tolower(rownames(df))),]

    ## stop if no control probes
    if (nrow(df) == 0)
        stop("No normalization control probes found!")

    if (sdfPlatform(sdf) == 'HM27') {
        df$channel <- ifelse(grepl(
            'norm\\.green', tolower(rownames(df))), 'G', 'R')
    } else {
        df$channel <- ifelse(grepl(
            'norm_(c|g)', tolower(rownames(df))), 'G', 'R')
    }
    
    if (average) {
        c(
            G=mean(df[df$channel=='G','G'], na.rm=TRUE), 
            R=mean(df[df$channel=='R','R'], na.rm=TRUE))
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
#' sesameDataCache("EPIC") # if not done yet
#' sdf <- sesameDataGet('EPIC.1.SigDF')
#' sdf.db <- dyeBiasCorr(sdf)
#' @export
dyeBiasCorr <- function(sdf, ref=NULL) {

    stopifnot(is(sdf, "SigDF"))
    if (is.null(ref)) {
        ref <- meanIntensity(sdf)
    }
    
    normctl <- getNormCtls(sdf, average=TRUE)
    fR <- ref/normctl['R']
    fG <- ref/normctl['G']

    sdf$MG = sdf$MG * fG
    sdf$UG = sdf$UG * fG
    sdf$MR = sdf$MR * fR
    sdf$UR = sdf$UR * fR

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
#' sesameDataCache("HM450") # if not done yet
#' sdfs <- sesameDataGet('HM450.10.SigDF')
#' sdfs.db <- dyeBiasCorrMostBalanced(sdfs)
#' @export
dyeBiasCorrMostBalanced <- function(sdfs) {

    normctls <- vapply(sdfs, getNormCtls, numeric(2), average=TRUE)
    most.balanced <- which.min(abs(normctls['G',] / normctls['R',] - 1))
    ref <- mean(normctls[,most.balanced], na.rm=TRUE)
    lapply(sdfs, function(sdf) dyeBiasCorr(sdf, ref))
}

#' Dye bias correction by matching green and red to mid point
#'
#' This function compares the Type-I Red probes and Type-I Grn probes and
#' generates and mapping to correct signal of the two channels to the middle.
#' The function takes one single \code{SigDF} and returns a \code{SigDF}
#' with dye bias corrected.
#' 
#' @param sdf a \code{SigDF}
#' @return a \code{SigDF} after dye bias correction.
#' @importFrom preprocessCore normalize.quantiles.use.target
#' @importFrom stats approx
#' @examples
#' sesameDataCache("EPIC") # if not done yet
#' sdf <- sesameDataGet('EPIC.1.SigDF')
#' sdf.db <- dyeBiasCorrTypeINorm(sdf)
#' @export
dyeBiasCorrTypeINorm <- function(sdf) {

    stopifnot(is(sdf, "SigDF"))

    ## we use all Inf-I probes so we capture the entire support range
    dG = InfIG(noMasked(sdf)); dR = InfIR(noMasked(sdf))
    IG0 <- c(dG$MG, dG$UG)
    IR0 <- c(dR$MR, dR$UR)
    
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

    sdf$MR = fitfunRed(sdf$MR)
    sdf$UR = fitfunRed(sdf$UR)
    sdf$MG = fitfunGrn(sdf$MG)
    sdf$UG = fitfunGrn(sdf$UG)
    sdf
}

#' @rdname dyeBiasCorrTypeINorm
#' @export
#' @examples
#' sdf <- sesameDataGet('EPIC.1.SigDF')
#' sdf <- dyeBiasNL(sdf)
dyeBiasNL <- dyeBiasCorrTypeINorm

## the following three functions are retired since 1.9.1
## dyeBiasCorrTypeINormMpU, dyeBiasCorrTypeINormG2R, dyeBiasCorrTypeINormR2G


#' Plot red-green QQ-Plot using Infinium-I Probes
#'
#' @param sdf a \code{SigDF}
#' @return create a qqplot
#' @examples
#' sesameDataCache("EPIC")  # if not done yet
#' sdf <- sesameDataGet('EPIC.1.SigDF')
#' sesamePlotRedGrnQQ(sdf)
#' @import graphics
#' @export
sesamePlotRedGrnQQ <- function(sdf) {
    dG = InfIG(noMasked(sdf)); dR = InfIR(noMasked(sdf))
    m = max(c(dR$MR,dR$UR,dG$MG,dG$UG), na.rm=TRUE)
    
    qqplot(
        c(dR$MR, dR$UR), c(dG$MG, dG$UG),
        xlab = 'Infinium-I Red Signal', ylab = 'Infinium-I Grn Signal',
        main = 'Red-Green QQ-Plot', cex = 0.5,
        xlim = c(0,m), ylim = c(0,m))
    abline(0,1,lty = 'dashed')
}

#' Quantify how much dye bias in high signal range deviates from the
#' global median
#'
#' Positive value indicates augmentation of high-end dye bias over
#' low-end. negative value represents high-end dye bias contradicts
#' that at low-end (a distorted dye bias). Negative distortion score
#' (< -1) suggests low experiment quality. 0 suggests a consistent
#' dye bias at high and low-end.
#'
#' @param sdf a \code{SigDF}
#' @return a numeric score
#' @examples
#' sdf <- sesameDataGet('EPIC.1.SigDF')
#' dyeBiasDistortion(sdf)
#' @export
dyeBiasDistortion = function(sdf) {
    t1 = InfI(sdf)
    intens = totalIntensities(sdf)
    intens[t1[t1$col == "G", "Probe_ID"]]
    medR = median(sort(intens[t1[t1$col == "R", "Probe_ID"]]))
    medG = median(sort(intens[t1[t1$col == "G", "Probe_ID"]]))
    topR = median(tail(sort(intens[t1[t1$col == "R", "Probe_ID"]]), n=20))
    topG = median(tail(sort(intens[t1[t1$col == "G", "Probe_ID"]]), n=20))
    log(topR / topG) / log(medR / medG)
}
