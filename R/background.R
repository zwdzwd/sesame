

#' SCRUB background correction
#'
#' This function takes a \code{SigDF} and returns a modified \code{SigDF}
#' with background subtracted. scrub subtracts residual background using
#' background median
#'
#' This function is meant to be used after noob.
#'
#' @param sdf a \code{SigDF}
#' @return a new \code{SigDF} with noob background correction
#' @examples
#' sdf <- sesameDataGet('EPIC.1.SigDF')
#' sdf.nb <- noob(sdf)
#' sdf.nb.scrub <- scrub(sdf.nb)
#' @export
scrub <- function(sdf) {
    bG <- median(oobG(noMasked(sdf)), na.rm=TRUE)
    bR <- median(oobR(noMasked(sdf)), na.rm=TRUE)
    sdf$MG <- pmax(sdf$MG - bG, 1)
    sdf$MR <- pmax(sdf$MR - bR, 1)
    sdf$UG <- pmax(sdf$UG - bG, 1)
    sdf$UR <- pmax(sdf$UR - bR, 1)
    sdf
}

noobSub <- function(sig, bg) {
    e <- MASS::huber(bg)
    mu <- e$mu
    sigma <- e$s
    alpha <- pmax(MASS::huber(sig)$mu-mu, 10)
    normExpSignal(mu, sigma, alpha, sig)
}

#' SCRUB background correction
#'
#' This function takes a \code{SigDF} and returns a modified \code{SigDF}
#' with background subtracted. scrubSoft subtracts residual background using a
#' noob-like procedure.
#'
#' This function is meant to be used after noob.
#'
#' @param sdf a \code{SigDF}
#' @return a new \code{SigDF} with noob background correction
#' @examples
#' sdf <- sesameDataGet('EPIC.1.SigDF')
#' sdf.nb <- noob(sdf)
#' sdf.nb.scrubSoft <- scrubSoft(sdf.nb)
#' @export
scrubSoft <- function(sdf) {
    bgR <- oobR(noMasked(sdf))
    bgG <- oobG(noMasked(sdf))

    sdf$MG <- noobSub(sdf$MG, bgG)
    sdf$MR <- noobSub(sdf$MR, bgR)
    sdf$UG <- noobSub(sdf$UG, bgG)
    sdf$UR <- noobSub(sdf$UR, bgR)
    sdf
}

#' Noob background subtraction
#'
#' The function takes a \code{SigDF} and returns a modified \code{SigDF}
#' with background subtracted. Background was modelled in a normal distribution
#' and true signal in an exponential distribution. The Norm-Exp deconvolution
#' is parameterized using Out-Of-Band (oob) probes. For species-specific
#' processing, one should call inferSpecies on SigDF first. Multi-mapping
#' probes are excluded.
#'
#' When combine.neg = TRUE, background will be parameterized by both
#' negative control and out-of-band probes.
#' 
#' @param sdf a \code{SigDF}
#' @param combine.neg whether to combine negative control probe.
#' @param offset offset
#' @return a new \code{SigDF} with noob background correction
#' @import stats
#' @examples
#' sdf <- sesameDataGet('EPIC.1.SigDF')
#' sdf.nb <- noob(sdf)
#' @export
noob <- function(sdf, combine.neg = TRUE, offset=15) {

    stopifnot(is(sdf, "SigDF"))
    nmk <- sdf[!(sdf$Probe_ID %in% nonuniqMask(sdfPlatform(sdf))),]
    bgG <- oobG(nmk)
    bgR <- oobR(nmk)
    if (combine.neg) {
        neg <- negControls(sdf)
        bgG <- c(bgG, neg$G)
        bgR <- c(bgR, neg$R)
    }
    
    ## if not enough out-of-band signal
    if (sum(bgG > 0, na.rm=TRUE) < 100 || sum(bgR > 0, na.rm=TRUE) < 100) {
        return(sdf)
    }
    bgR[bgR == 0] <- 1; bgG[bgG == 0] <- 1
    ## cap at 10xIQR, this is to proof against multi-mapping probes
    bgR <- bgR[bgR < median(bgR, na.rm=TRUE) + 10*IQR(bgR, na.rm=TRUE)]
    bgG <- bgG[bgG < median(bgG, na.rm=TRUE) + 10*IQR(bgG, na.rm=TRUE)]

    ## foreground
    ibG <- c(InfIG(nmk)$MG, InfIG(nmk)$UG, InfII(nmk)$UG)
    ibR <- c(InfIR(nmk)$MR, InfIR(nmk)$UR, InfII(nmk)$UR)
    ibG[ibG == 0] <- 1 # set signal to 1 if 0
    ibR[ibR == 0] <- 1 # set signal to 1 if 0

    ## grn channel
    fitG <- backgroundCorrectionNoobFit(ibG, bgG)
    sdf$MG <- normExpSignal(fitG$mu, fitG$sigma, fitG$alpha, sdf$MG) + 15
    sdf$UG <- normExpSignal(fitG$mu, fitG$sigma, fitG$alpha, sdf$UG) + 15

    ## red channel
    fitR <- backgroundCorrectionNoobFit(ibR, bgR)
    sdf$MR <- normExpSignal(fitR$mu, fitR$sigma, fitR$alpha, sdf$MR) + 15
    sdf$UR <- normExpSignal(fitR$mu, fitR$sigma, fitR$alpha, sdf$UR) + 15
    
    sdf
}

## Noob background correction for one channel
## ib array of in-band signal
## oob array of out-of-band-signal
## ctl control probe signals
## offset padding for normalized signal
## return normalized in-band signal
backgroundCorrectionNoobFit <- function(ib, bg) {
    e <- MASS::huber(bg)
    mu <- e$mu
    sigma <- e$s
    alpha <- pmax(MASS::huber(ib)$mu-mu, 10)
    list(mu = mu, sigma = sigma, alpha = alpha)
}


## the following is adapted from Limma
## normal-exponential deconvolution (conditional expectation of
## xs|xf; WEHI code)
normExpSignal <- function(mu, sigma, alpha, x)  {

    sigma2 <- sigma * sigma

    if (alpha <= 0)
        stop("alpha must be positive")
    if (sigma <= 0)
        stop("sigma must be positive")
    
    mu.sf <- x - mu - sigma2/alpha
    signal <- mu.sf + sigma2 * exp(
        dnorm(0, mean = mu.sf, sd = sigma, log = TRUE) -
            pnorm(
                0, mean = mu.sf, sd = sigma,
                lower.tail = FALSE, log.p = TRUE))
    
    o <- !is.na(signal)
    if (any(signal[o] < 0)) {
        warning("Limit of numerical accuracy reached with very
low intensity or very high background:\nsetting adjusted intensities
to small value")
        signal[o] <- pmax(signal[o], 1e-06)
    }
    signal
}

## Noob background correction for one channel
## pp is the parameters
## ib array of in-band signal
## oob array of out-of-band-signal
## ctl control probe signals
## offset padding for normalized signal
## return normalized in-band signal
.backgroundCorrCh1 <- function(x, pp, alpha, offset=15) {

    mu.bg <- pp$mu
    sigma <- pp$sigma
    sigma2 <- sigma * sigma

    if (alpha <= 0)
        stop("alpha must be positive")
    if (any(na.omit(sigma) <= 0))
        stop("sigma must be positive")

    mu.sf <- x - mu.bg - sigma2/alpha
    signal <- mu.sf + sigma2 * exp(
        dnorm(0, mean = mu.sf, sd = sigma, log = TRUE) -
            pnorm(
                0, mean = mu.sf, sd = sigma,
                lower.tail = FALSE, log.p = TRUE))
    
    o <- !is.na(signal)
    if (any(signal[o] < 0)) {
        warning("Limit of numerical accuracy reached with
very low intensity or very high background:
setting adjusted intensities to small value")
        signal[o] <- pmax(signal[o], 1e-06)
    }
    
    signal.min <- min(signal, na.rm = TRUE)
    signal <- signal - signal.min
    offset + signal
}

## train model using linear model with log transformed response
train.model.lm <- function(input, output) {

    fitdata <- data.frame(IB=as.vector(input)+1, LOB=log(as.vector(output)+1))
    m <- lm(LOB~IB, data=fitdata)
    ## m <- MASS::rlm(LOB~IB, data=fitdata)

    function(d) {
        force(m)
        pp <- predict(
            m, newdata=data.frame(IB=as.vector(d)),
            interval='prediction', level=0.8)
        
        list(mu=exp(pp[,'fit']), sigma=(exp(pp[,'upr'])-exp(pp[,'lwr']))/10.13)
        ## use upper bound for mu since true signal
        ## is often much higher than noise
        ## list(mu=exp(pp[,'upr']),
        ## sigma=(exp(pp[,'upr'])-exp(pp[,'lwr']))/10.13)
        ## list(mu=exp(pp[,'upr']),
        ## sigma=log(exp(pp[,'upr'])-exp(pp[,'lwr'])))
    }
}
