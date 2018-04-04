#' Noob background correction
#' 
#' Norm-Exp deconvolution using Out-Of-Band (oob) probes
#' Note p-values are unchanged (based on the raw signal intensities).
#' @param sset a \code{SignalSet}
#' @param offset offset
#' @return a new \code{SignalSet} with noob background correction
#' @examples
#' sset <- makeExampleTinyEPICDataSet()
#' sset.nb <- noob(sset)
#' @export
noob <- function(sset, offset=15) {

    ## sort signal based on channel
    ibR <- c(sset@IR, sset@II[,'U'])    # in-band red signal
    ibG <- c(sset@IG, sset@II[,'M'])    # in-band green signal

    ## set signal to 1 if 0
    ibR[ibR==0] <- 1
    ibG[ibG==0] <- 1

    ## oobG and oobR are untouched besides the 0>1 switch
    sset@oobR[sset@oobR==0] <- 1
    sset@oobG[sset@oobG==0] <- 1
    
    ## do background correction in each channel
    ibR.nl <- .backgroundCorrectionNoobCh1(
        ibR, sset@oobR, sset@ctl$R, offset=offset)
    ibG.nl <- .backgroundCorrectionNoobCh1(
        ibG, sset@oobG, sset@ctl$G, offset=offset)

    ## build back the list
    ## type IG
    if (length(sset@IG)>0)
        sset@IG <- matrix(
            ibG.nl$i[seq_along(sset@IG)],
            nrow=nrow(sset@IG), dimnames=dimnames(sset@IG))
    else
        sset@IG <- matrix(ncol=2, nrow=0, dimnames=list(NULL,c('M','U')))

    ## type IR
    if (length(sset@IR)>0)
        sset@IR <- matrix(
            ibR.nl$i[seq_along(sset@IR)],
            nrow=nrow(sset@IR), dimnames=dimnames(sset@IR))
    else
        sset@IR <- matrix(ncol=2, nrow=0, dimnames=list(NULL,c('M','U')))

    ## type II
    if (nrow(sset@II) > 0)
        sset@II <- as.matrix(data.frame(
            M=ibG.nl$i[(length(sset@IG)+1):length(ibG)],
            U=ibR.nl$i[(length(sset@IR)+1):length(ibR)],
            row.names=rownames(sset@II)))
    else
        sset@II <- matrix(ncol=2, nrow=0, dimnames=list(NULL,c('M','U')))

    ## controls
    sset@ctl$G <- ibG.nl$c
    sset@ctl$R <- ibR.nl$c

    ## out-of-band
    sset@oobR <- ibR.nl$o
    sset@oobG <- ibG.nl$o

    sset
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

.getNormCtls <- function(sset) {
    
    if (sset@platform == 'hm27') {
        ## controversial to use, maybe mean of all signals in each channel?
        normctl.G <- sset@ctl[grep('norm.green', tolower(rownames(sset@ctl))),]
        normctl.R <- sset@ctl[grep('norm.red', tolower(rownames(sset@ctl))),]
    } else {                              # HM450 and EPIC
        normctl.G <- sset@ctl[grep('norm_(c|g)',tolower(rownames(sset@ctl))),]
        normctl.R <- sset@ctl[grep('norm_(a|t)',tolower(rownames(sset@ctl))),]
    }
    c(G=mean(normctl.G$G, na.rm=TRUE), R=mean(normctl.R$R, na.rm=TRUE))
    
}

#' background subtraction with bleeding-through subtraction
#'
#' Norm-Exp deconvolution using Out-Of-Band (oob) probes.
#' Mean and standard deviations are estimated from cross-channel regression
#' 
#' @param sset a \code{SignalSet}
#' @param offset offset
#' @param detailed if TRUE, return a list of \code{SignalSet}
#' and regression function
#' @return a modified \code{SignalSet} with background correction
#' @examples
#' sset <- makeExampleSeSAMeDataSet('HM450')
#' sset.nb <- noobsb(sset)
#' @export
noobsb <- function(sset, offset=15, detailed=FALSE) {

    ## sanitize
    ## sort signal based on channel
    ibR <- c(sset@IR, sset@II[,'U'])    # in-band red signal
    ibR.other.channel <- c(sset@oobG, sset@II[,'M'])
    ibG <- c(sset@IG, sset@II[,'M'])    # in-band green signal
    ibG.other.channel <- c(sset@oobR, sset@II[,'U'])
    ## set signal to 1 if 0
    ibR[ibR==0] <- 1
    ibG[ibG==0] <- 1
    
    ## oobG and oobR are untouched besides the 0>1 switch
    sset@oobR[sset@oobR==0] <- 1
    sset@oobG[sset@oobG==0] <- 1

    ## background correction Red
    ## use the other channel to predict background mean
    bg.GpredictR <- train.model.lm(sset@IG, sset@oobR)
    pp.bg.ibR <- bg.GpredictR(ibR.other.channel)
    pp.bg.oobR <- bg.GpredictR(sset@IG)
    pp.bg.ctlR <- bg.GpredictR(sset@ctl$G)
    ## parameter estimation
    ## sigma.bgR <- MASS::huber(sset@oobR)$s # for now, use the global variance
    alphaR <- pmax(MASS::huber(ibR - pp.bg.ibR$mu)$mu, 10) # in-band variance
    ## correction
    ibR <- .backgroundCorrCh1(ibR, pp.bg.ibR, alphaR, offset=offset)
    oobR <- .backgroundCorrCh1(sset@oobR, pp.bg.oobR, alphaR, offset=offset)
    ctlR <- .backgroundCorrCh1(sset@ctl$R, pp.bg.ctlR, alphaR, offset=offset)
    
    ## background correction Grn
    ## use the other channel to predict background mean
    bg.RpredictG <- train.model.lm(sset@IR, sset@oobG)
    pp.bg.ibG <- bg.RpredictG(ibG.other.channel)
    pp.bg.oobG <- bg.RpredictG(sset@IR)
    pp.bg.ctlG <- bg.RpredictG(sset@ctl$R)
    ## parameter estimation
    ## sigma.bgG <- MASS::huber(sset@oobG)$s # for now, use the global variance
    alphaG <- pmax(MASS::huber(ibG - pp.bg.ibG$mu)$mu, 10) # in-band variance
    ## correction
    ibG <- .backgroundCorrCh1(ibG, pp.bg.ibG, alphaG, offset=offset)
    oobG <- .backgroundCorrCh1(sset@oobG, pp.bg.oobG, alphaG, offset=offset)
    ctlG <- .backgroundCorrCh1(sset@ctl$G, pp.bg.ctlG, alphaG, offset=offset)

    ## build back the list
    ## type IG
    if (length(sset@IG)>0)
        sset@IG <- matrix(
            ibG[seq_along(sset@IG)],
            nrow=nrow(sset@IG), dimnames=dimnames(sset@IG))
    else
        sset@IG <- matrix(ncol=2, nrow=0, dimnames=list(NULL,c('M','U')))

    ## type IR
    if (length(sset@IR)>0)
        sset@IR <- matrix(
            ibR[seq_along(sset@IR)],
            nrow=nrow(sset@IR), dimnames=dimnames(sset@IR))
    else
        sset@IR <- matrix(ncol=2, nrow=0, dimnames=list(NULL,c('M','U')))

    ## type II
    if (nrow(sset@II) > 0)
        sset@II <- as.matrix(data.frame(
            M=ibG[(length(sset@IG)+1):length(ibG)],
            U=ibR[(length(sset@IR)+1):length(ibR)],
            row.names=rownames(sset@II)))
    else
        sset@II <- matrix(ncol=2, nrow=0, dimnames=list(NULL,c('M','U')))

    ## controls
    sset@ctl$G <- ctlG
    sset@ctl$R <- ctlR

    ## out-of-band
    sset@oobR <- oobR
    sset@oobG <- oobG

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
