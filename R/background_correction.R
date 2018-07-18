#' Noob background correction
#'
#' The function takes a \code{SigSet} and returns a modified \code{SigSet}
#' with background subtracted. Background was modelled in a normal distribution
#' and true signal in an exponential distribution. The Norm-Exp deconvolution
#' is parameterized using Out-Of-Band (oob) probes
#' 
#' @param sset a \code{SigSet}
#' @param offset offset
#' @return a new \code{SigSet} with noob background correction
#' @examples
#' sset <- makeExampleTinyEPICDataSet()
#' sset.nb <- noob(sset)
#' @export
noob <- function(sset, offset=15) {

    ## sort signal based on channel
    ibR <- c(IR(sset), II(sset)[,'U'])    # in-band red signal
    ibG <- c(IG(sset), II(sset)[,'M'])    # in-band green signal

    ## set signal to 1 if 0
    ibR[ibR==0] <- 1
    ibG[ibG==0] <- 1

    ## oobG and oobR are untouched besides the 0>1 switch
    oobR(sset)[oobR(sset)==0] <- 1
    oobG(sset)[oobG(sset)==0] <- 1
    
    ## do background correction in each channel
    ibR.nl <- .backgroundCorrectionNoobCh1(
        ibR, oobR(sset), ctl(sset)$R, offset=offset)
    ibG.nl <- .backgroundCorrectionNoobCh1(
        ibG, oobG(sset), ctl(sset)$G, offset=offset)

    ## build back the list
    ## type IG
    if (length(IG(sset))>0)
        IG(sset) <- matrix(
            ibG.nl$i[seq_along(IG(sset))],
            nrow=nrow(IG(sset)), dimnames=dimnames(IG(sset)))
    else
        IG(sset) <- matrix(ncol=2, nrow=0, dimnames=list(NULL,c('M','U')))

    ## type IR
    if (length(IR(sset))>0)
        IR(sset) <- matrix(
            ibR.nl$i[seq_along(IR(sset))],
            nrow=nrow(IR(sset)), dimnames=dimnames(IR(sset)))
    else
        IR(sset) <- matrix(ncol=2, nrow=0, dimnames=list(NULL,c('M','U')))

    ## type II
    if (nrow(II(sset)) > 0)
        II(sset) <- as.matrix(data.frame(
            M=ibG.nl$i[(length(IG(sset))+1):length(ibG)],
            U=ibR.nl$i[(length(IR(sset))+1):length(ibR)],
            row.names=rownames(II(sset))))
    else
        II(sset) <- matrix(ncol=2, nrow=0, dimnames=list(NULL,c('M','U')))

    ## controls
    ctl(sset)$G <- ibG.nl$c
    ctl(sset)$R <- ibR.nl$c

    ## out-of-band
    oobR(sset) <- ibR.nl$o
    oobG(sset) <- ibG.nl$o

    sset
}

#' The openSeSAMe pipeline
#'
#' This function takes a \code{SigSet} and output beta values.It is
#' a simple wrapper of noob + nonlinear dye bias correction + pOOBAH masking.
#' 
#' @param sset a \code{SigSet}
#' @return a numeric vector for processed beta values
#' sset <- make
#' @examples
#' sset <- sesameDataGet('HM450.1.TCGA.PAAD')$sset
#' betas <- openSeSAMe(sset)
#' @export
openSeSAMe <- function(sset) {
    getBetas(dyeBiasCorrTypeINorm(noob(sset)))
}

## Noob background correction for one channel
## ib array of in-band signal
## oob array of out-of-band-signal
## ctl control probe signals
## offset padding for normalized signal
## return normalized in-band signal
.backgroundCorrectionNoobCh1 <- function(ib, oob, ctl, offset=15) {

    e <- MASS::huber(oob)
    mu <- e$mu
    sigma <- e$s
    alpha <- pmax(MASS::huber(ib)$mu-mu, 10)
    list(
        i=offset+.normExpSignal(mu, sigma, alpha, ib),
        c=offset+.normExpSignal(mu, sigma, alpha, ctl),
        o=offset+.normExpSignal(mu, sigma, alpha, oob))
}

## the following is adapted from Limma
## normal-exponential deconvolution (conditional expectation of
## xs|xf; WEHI code)
.normExpSignal <- function (mu, sigma, alpha, x)  {

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

#' Background subtraction with bleeding-through subtraction
#'
#' The function takes a \code{SigSet} and returns a modified \code{SigSet}
#' with background subtracted. Signal bleed-through was modelled using a
#' linear model with error estimated from cross-channel regression.
#' Norm-Exp deconvolution using Out-Of-Band (oob) probes.
#' 
#' @param sset a \code{SigSet}
#' @param offset offset
#' @param detailed if TRUE, return a list of \code{SigSet}
#' and regression function
#' @return a modified \code{SigSet} with background correction
#' @examples
#' sset <- makeExampleSeSAMeDataSet('HM450')
#' sset.nb <- noobsb(sset)
#' @export
noobsb <- function(sset, offset=15, detailed=FALSE) {

    ## sanitize
    ## sort signal based on channel
    ibR <- c(IR(sset), II(sset)[,'U'])    # in-band red signal
    ibR.other.channel <- c(oobG(sset), II(sset)[,'M'])
    ibG <- c(IG(sset), II(sset)[,'M'])    # in-band green signal
    ibG.other.channel <- c(oobR(sset), II(sset)[,'U'])
    ## set signal to 1 if 0
    ibR[ibR==0] <- 1
    ibG[ibG==0] <- 1
    
    ## oobG and oobR are untouched besides the 0>1 switch
    oobR(sset)[oobR(sset)==0] <- 1
    oobG(sset)[oobG(sset)==0] <- 1

    ## background correction Red
    ## use the other channel to predict background mean
    bg.GpredictR <- train.model.lm(IG(sset), oobR(sset))
    pp.bg.ibR <- bg.GpredictR(ibR.other.channel)
    pp.bg.oobR <- bg.GpredictR(IG(sset))
    pp.bg.ctlR <- bg.GpredictR(ctl(sset)$G)
    ## parameter estimation
    ## sigma.bgR <- MASS::huber(oobR(sset))$s # for now, use the global variance
    alphaR <- pmax(MASS::huber(ibR - pp.bg.ibR$mu)$mu, 10) # in-band variance
    ## correction
    ibR <- .backgroundCorrCh1(ibR, pp.bg.ibR, alphaR, offset=offset)
    oobR <- .backgroundCorrCh1(oobR(sset), pp.bg.oobR, alphaR, offset=offset)
    ctlR <- .backgroundCorrCh1(ctl(sset)$R, pp.bg.ctlR, alphaR, offset=offset)
    
    ## background correction Grn
    ## use the other channel to predict background mean
    bg.RpredictG <- train.model.lm(IR(sset), oobG(sset))
    pp.bg.ibG <- bg.RpredictG(ibG.other.channel)
    pp.bg.oobG <- bg.RpredictG(IR(sset))
    pp.bg.ctlG <- bg.RpredictG(ctl(sset)$R)
    ## parameter estimation
    ## sigma.bgG <- MASS::huber(oobG(sset))$s # for now, use the global variance
    alphaG <- pmax(MASS::huber(ibG - pp.bg.ibG$mu)$mu, 10) # in-band variance
    ## correction
    ibG <- .backgroundCorrCh1(ibG, pp.bg.ibG, alphaG, offset=offset)
    oobG <- .backgroundCorrCh1(oobG(sset), pp.bg.oobG, alphaG, offset=offset)
    ctlG <- .backgroundCorrCh1(ctl(sset)$G, pp.bg.ctlG, alphaG, offset=offset)

    ## build back the list
    ## type IG
    if (length(IG(sset))>0)
        IG(sset) <- matrix(
            ibG[seq_along(IG(sset))],
            nrow=nrow(IG(sset)), dimnames=dimnames(IG(sset)))
    else
        IG(sset) <- matrix(ncol=2, nrow=0, dimnames=list(NULL,c('M','U')))

    ## type IR
    if (length(IR(sset))>0)
        IR(sset) <- matrix(
            ibR[seq_along(IR(sset))],
            nrow=nrow(IR(sset)), dimnames=dimnames(IR(sset)))
    else
        IR(sset) <- matrix(ncol=2, nrow=0, dimnames=list(NULL,c('M','U')))

    ## type II
    if (nrow(II(sset)) > 0)
        II(sset) <- as.matrix(data.frame(
            M=ibG[(length(IG(sset))+1):length(ibG)],
            U=ibR[(length(IR(sset))+1):length(ibR)],
            row.names=rownames(II(sset))))
    else
        II(sset) <- matrix(ncol=2, nrow=0, dimnames=list(NULL,c('M','U')))

    ## controls
    ctl(sset)$G <- ctlG
    ctl(sset)$R <- ctlR

    ## out-of-band
    oobR(sset) <- oobR
    oobG(sset) <- oobG

    if (detailed)
        list(sset=sset, bg.RpredictG=bg.RpredictG, bg.GpredictR=bg.GpredictR)
    else
        sset
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
