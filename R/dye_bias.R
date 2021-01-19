#' get normalization control signal
#' 
#' get normalization control signal from SigSet. 
#' The function optionally takes mean for each channel.
#' 
#' @param   sset a SigSet
#' @param   average whether to average
#' @return a data frame of normalization control signals
#' @examples 
#' sset <- readIDATpair(file.path(system.file(
#'     'extdata','',package='sesameData'), '4207113116_B'))
#' 
#' df.ctl <- getNormCtls(sset)
#' 
#' @export
getNormCtls <- function(sset, average = FALSE) {
    df <- ctl(sset)
    df <- df[grep('norm(_|\\.)', tolower(rownames(df))),]

    ## stop if no control probes
    if (nrow(df) == 0)
        stop("No normalization control probes found!")

    if (sset@platform == 'HM27') {
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
#' The function takes a \code{SigSet} as input and scale both the Grn and Red
#' signal to a reference (ref) level. If the reference level is not given, it
#' is set to the mean intensity of all the in-band signals. The function
#' returns a \code{SigSet} with dye bias corrected.
#' 
#' @param sset a \code{SigSet}
#' @param ref reference signal level
#' @return a normalized \code{SigSet}
#' @examples
#' sset <- sesameDataGet('EPIC.1.LNCaP')$sset
#' sset.db <- dyeBiasCorr(sset)
#' @export
dyeBiasCorr <- function(sset, ref=NULL) {

    stopifnot(is(sset, "SigSet"))
    if (is.null(ref)) {
        ref <- meanIntensity(sset)
    }
    
    normctl <- getNormCtls(sset, average=TRUE)
    fR <- ref/normctl['R']
    fG <- ref/normctl['G']

    IG(sset) <- matrix(
        c(fG*IG(sset)[,'M'], fG*IG(sset)[,'U']),
        nrow=nrow(IG(sset)), ncol=ncol(IG(sset)), dimnames=dimnames(IG(sset)))
    
    IR(sset) <- matrix(
        c(fR*IR(sset)[,'M'], fR*IR(sset)[,'U']),
        nrow=nrow(IR(sset)), ncol=ncol(IR(sset)), dimnames=dimnames(IR(sset)))
    
    II(sset) <- matrix(
        c(fG*II(sset)[,'M'], fR*II(sset)[,'U']),
        nrow=nrow(II(sset)), ncol=ncol(II(sset)), dimnames=dimnames(II(sset)))
    
    ctl(sset)$G <- fG*ctl(sset)$G
    ctl(sset)$R <- fR*ctl(sset)$R
    oobG(sset) <- fG*oobG(sset)
    oobR(sset) <- fR*oobR(sset)
    sset
}

#' Correct dye bias using most balanced sample as the reference
#'
#' The function chose the reference signal level from a list of \code{SigSet}.
#' The chosen sample has the smallest difference in Grn and Red signal
#' intensity as measured using the normalization control probes. In practice,
#' it doesn't matter which sample is chosen as long as the reference level
#' does not deviate much. The function returns a list of \code{SigSet}s with
#' dye bias corrected.
#' 
#' @param ssets a list of normalized \code{SigSet}s
#' @return a list of normalized \code{SigSet}s
#' @examples
#' ssets <- sesameDataGet('HM450.10.TCGA.BLCA.normal')
#' ssets.db <- dyeBiasCorrMostBalanced(ssets)
#' @export
dyeBiasCorrMostBalanced <- function(ssets) {

    normctls <- vapply(ssets, getNormCtls, numeric(2), average=TRUE)
    most.balanced <- which.min(abs(normctls['G',] / normctls['R',] - 1))
    ref <- mean(normctls[,most.balanced], na.rm=TRUE)
    lapply(ssets, function(sset) dyeBiasCorr(sset, ref))
}

#' Dye bias correction by matching green and red to mid point
#'
#' This function compares the Type-I Red probes and Type-I Grn probes and
#' generates and mapping to correct signal of the two channels to the middle.
#' The function takes one single \code{SigSet} and returns a \code{SigSet}
#' with dye bias corrected.
#' 
#' @param sset a \code{SigSet}
#' @return a \code{SigSet} after dye bias correction.
#' @importFrom preprocessCore normalize.quantiles.use.target
#' @importFrom stats approx
#' @examples
#' sset <- sesameDataGet('EPIC.1.LNCaP')$sset
#' sset.db <- dyeBiasCorrTypeINorm(sset)
#' @export
dyeBiasCorrTypeINorm <- function(sset) {

    stopifnot(is(sset, "SigSet"))

    IG0 <- IGpass(sset); IR0 <- IRpass(sset)
    
    maxIG <- max(IG0, na.rm = TRUE); minIG <- min(IG0, na.rm = TRUE)
    maxIR <- max(IR0, na.rm = TRUE); minIR <- min(IR0, na.rm = TRUE)

    if (maxIG <= 0 || maxIR <= 0) { return(sset); }

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

    ## fit type II
    II(sset)[,'U'] <- fitfunRed(II(sset)[,'U'])
    II(sset)[,'M'] <- fitfunGrn(II(sset)[,'M'])

    ## IR
    IR(sset) <- fitfunRed(IR(sset))
    ## IG
    IG(sset) <- fitfunGrn(IG(sset))

    ## fit control
    if (nrow(ctl(sset)) > 0) {
        ctl(sset)[,'R'] <- fitfunRed(ctl(sset)[,'R'])
        ctl(sset)[,'G'] <- fitfunGrn(ctl(sset)[,'G'])
    }

    ## fit oob
    oobR(sset) <- fitfunRed(oobR(sset))
    oobG(sset) <- fitfunGrn(oobG(sset))

    sset
}

## the following three functions are retired since 1.9.1
## dyeBiasCorrTypeINormMpU, dyeBiasCorrTypeINormG2R, dyeBiasCorrTypeINormR2G

