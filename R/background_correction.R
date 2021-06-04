

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
#' sdf <- makeExampleTinyEPICDataSet()
#' sdf.nb <- noob(sdf)
#' sdf.nb.scrub <- scrub(sdf.nb)
#' @export
scrub <- function(sdf) {
    bG <- median(oobGpass(sdf))
    bR <- median(oobRpass(sdf))
    sdf@IR <- pmax(sdf@IR - bR,1)
    sdf@IG <- pmax(sdf@IG - bG,1)
    sdf@II <- cbind(M=pmax(sdf@II[,'M'] - bG,1), U=pmax(sdf@II[,'U'] - bR,1))
    sdf@oobR <- pmax(sdf@oobR - bR,1) # subtract oobR itself
    sdf@oobG <- pmax(sdf@oobG - bG,1) # subtract oobG itself
    sdf
}

noobSub <- function(sig, bg) {
    e <- MASS::huber(bg)
    mu <- e$mu
    sigma <- e$s
    alpha <- pmax(MASS::huber(sig)$mu-mu, 10)
    .normExpSignal(mu, sigma, alpha, sig)
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
#' sdf <- makeExampleTinyEPICDataSet()
#' sdf.nb <- noob(sdf)
#' sdf.nb.scrubSoft <- scrubSoft(sdf.nb)
#' @export
scrubSoft <- function(sdf) {
    oobR1 <- oobRpass(sdf)
    oobG1 <- oobGpass(sdf)
    sdf@IR <- noobSub(sdf@IR, oobR1)
    sdf@IG <- noobSub(sdf@IG, oobG1)
    sdf@II <- cbind(
        M=noobSub(sdf@II[,'M'], oobG1),
        U=noobSub(sdf@II[,'U'], oobR1))
    sdf@oobR <- noobSub(oobR(sdf), oobR1) # subtract oobR itself
    sdf@oobG <- noobSub(oobG(sdf), oobG1) # subtract oobG itself
    sdf
}

getBackgroundR <- function(sdf, bgR = NULL) {
    if (!is.null(bgR)) {
        oobR(sdf)[bgR,]
    } else {
        oobRpass(sdf)
    }
}

getBackgroundG <- function(sdf, bgG = NULL) {
    if (!is.null(bgG)) {
        oobG(sdf)[bgG,]
    } else {
        oobGpass(sdf)
    }
}

#' Noob background correction
#'
#' The function takes a \code{SigDF} and returns a modified \code{SigDF}
#' with background subtracted. Background was modelled in a normal distribution
#' and true signal in an exponential distribution. The Norm-Exp deconvolution
#' is parameterized using Out-Of-Band (oob) probes
#' 
#' @param sdf a \code{SigDF}
#' @param offset offset
#' @param bgR background red probes, if not given use all oobR
#' @param bgG background grn probes, if not given use all oobG
#' @return a new \code{SigDF} with noob background correction
#' @import stats
#' @examples
#' sdf <- makeExampleTinyEPICDataSet()
#' sdf.nb <- noob(sdf)
#' @export
noob <- function(sdf, bgR = NULL, bgG = NULL, offset=15) {

    ## if no Infinium-I probes
    if (nrow(IG(sdf)) == 0 && nrow(IR(sdf)) == 0) { return(sdf) }
    ## if not enough out-of-band signal
    oobG = oobG(noMask(sdf))
    oobR = oobR(noMask(sdf))
    if (sum(oobG > 0, na.rm=TRUE) < 100 ||
            sum(oobR > 0, na.rm=TRUE) < 100) { return(sdf) }
    oobR[oobR == 0] = 1
    oobG[oobG == 0] = 1

    ## grn channel: IG and oobG
    idx_mg = sdf$col == "G"
    idx_ug = sdf$col == "G" | sdf$col == "2"
    ibG = c(sdf$MG[idx_mg], sdf$UG[idx_ug])
    ibG[ibG == 0] = 1 # set signal to 1 if 0
    idx_r = sdf$col == "R"
    oob_g = c(sdf$MG[idx_r], sdf$UG[idx_r])
    oob_g[oob_g == 0] = 1
    out = .backgroundCorrectionNoobCh1(ibG, oob_g, controls(sdf)$G, oobG, offset=offset)
    sdf$MG[idx_mg] = out$i[1:sum(idx_mg)]
    sdf$UG[idx_ug] = out$i[(sum(idx_mg)+1):length(out$i)]
    sdf$MG[idx_r] = out$o[1:sum(idx_r)]
    sdf$UG[idx_r] = out$o[(sum(idx_r)+1):length(out$o)]

    ## red channel: IR and oobR
    idx_mr = sdf$col == "R"
    idx_ur = sdf$col == "R" | sdf$col == "2"
    ibR = c(sdf$MR[idx_mr], sdf$UR[idx_ur])
    ibR[ibR == 0] = 1 # set signal to 1 if 0
    idx_g = sdf$col == "G"
    oob_r = c(sdf$MR[idx_g], sdf$UR[idx_g])
    oob_r[oob_r == 0] = 1
    out = .backgroundCorrectionNoobCh1(ibR, oob_r, controls(sdf)$R, oobR, offset=offset)
    sdf$MR[idx_mr] = out$i[1:sum(idx_mr)]
    sdf$UR[idx_ur] = out$i[(sum(idx_mr)+1):length(out$i)]
    sdf$MR[idx_g] = out$o[1:sum(idx_g)]
    sdf$UR[idx_g] = out$o[(sum(idx_g)+1):length(out$o)]

    ## x2 = .backgroundCorrectionNoobCh1(x, ibR, oobR, offset=offset)
    ## sdf$MR[idx_g] = x2[1:length(idx_g)]
    ## sdf$UR[idx_g] = x2[(length(idx_g)+1):2*length(idx_g)]
    ## x2 = .backgroundCorrectionNoobCh1(x, ibG, oobG, offset=offset)
    ## sdf$MG[idx_r] = x2[1:length(idx_r)]
    ## sdf$UG[idx_r] = x2[(length(idx_r)+1):2*length(idx_r)]
        
    ## ibR.nl <- .backgroundCorrectionNoobCh1(
    ##     ibR, oobR(sdf), controls(sdf)$R, getBackgroundR(sdf, bgR), offset=offset)
    ## ibG.nl <- .backgroundCorrectionNoobCh1(
    ##     ibG, oobG(sdf), controls(sdf)$G, getBackgroundG(sdf, bgG), offset=offset)

    ## ## build back the list
    ## ## type IG
    ## if (length(IG(sdf))>0)
    ##     IG(sdf) <- matrix(
    ##         ibG.nl$i[seq_along(IG(sdf))],
    ##         nrow=nrow(IG(sdf)), dimnames=dimnames(IG(sdf)))
    ## else
    ##     IG(sdf) <- matrix(ncol=2, nrow=0, dimnames=list(NULL,c('M','U')))

    ## ## type IR
    ## if (length(IR(sdf))>0) {
    ##     IR(sdf) <- matrix(
    ##         ibR.nl$i[seq_along(IR(sdf))],
    ##         nrow=nrow(IR(sdf)), dimnames=dimnames(IR(sdf)))
    ## } else {
    ##     IR(sdf) <- matrix(ncol=2, nrow=0, dimnames=list(NULL,c('M','U')))
    ## }

    ## ## type II
    ## if (nrow(II(sdf)) > 0) {
    ##     II(sdf) <- as.matrix(data.frame(
    ##         M=ibG.nl$i[(length(IG(sdf))+1):length(ibG)],
    ##         U=ibR.nl$i[(length(IR(sdf))+1):length(ibR)],
    ##         row.names=rownames(II(sdf))))
    ## } else {
    ##     II(sdf) <- matrix(ncol=2, nrow=0, dimnames=list(NULL,c('M','U')))
    ## }

    ## ## controls
    ## ctl(sdf)$G <- ibG.nl$c
    ## ctl(sdf)$R <- ibR.nl$c

    ## ## out-of-band
    ## oobR(sdf) <- ibR.nl$o
    ## oobG(sdf) <- ibG.nl$o

    sdf
}

## Noob background correction for one channel
## ib array of in-band signal
## oob array of out-of-band-signal
## ctl control probe signals
## offset padding for normalized signal
## return normalized in-band signal
.backgroundCorrectionNoobCh1 <- function(ib, oob, ctl, bg, offset=15) {

    e <- MASS::huber(bg)
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
