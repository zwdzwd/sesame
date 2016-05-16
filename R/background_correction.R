
#' Noob background correction
#'
#' Norm-Exp deconvolution using Out-Of-Band (oob) probes
#'
#' Note p-values are unchanged (based on the raw signal intensities).
#' @param offset 
#' @import MASS
.backgroundCorrectionNoob <- function(offset=15) {

  ## sort signal based on channel
  ibR <- c(self$IR, self$II[,'U'])              # in-band red signal
  ibG <- c(self$IG, self$II[,'M'])              # in-band green signal

  ## set signal to 1 if 0
  ibR[ibR==0] <- 1
  ibG[ibG==0] <- 1

  ## oobG and oobR are untouched besides the 0>1 switch
  self$oobR[self$oobR==0] <- 1
  self$oobG[self$oobG==0] <- 1
  
  ## do background correction in each channel
  ibR.nl <- .backgroundCorrectionNoobCh1(ibR, self$oobR, self$ctl$R, offset=offset)
  ibG.nl <- .backgroundCorrectionNoobCh1(ibG, self$oobG, self$ctl$G, offset=offset)

  ## build back the list
  ## type IG
  if (length(self$IG)>0)
    self$IG <- matrix(ibG.nl$i[1:length(self$IG)],
                        nrow=nrow(self$IG), dimnames=dimnames(self$IG))
  else
    self$IG <- matrix(ncol=2, nrow=0, dimnames=list(NULL,c('M','U')))

  ## type IR
  if (length(self$IR)>0)
    self$IR <- matrix(ibR.nl$i[1:length(self$IR)],
                      nrow=nrow(self$IR), dimnames=dimnames(self$IR))
  else
    self$IR <- matrix(ncol=2, nrow=0, dimnames=list(NULL,c('M','U')))

  ## type II
  if (nrow(self$II) > 0)
    self$II <- as.matrix(data.frame(
      M=ibG.nl$i[(length(self$IG)+1):length(ibG)],
      U=ibR.nl$i[(length(self$IR)+1):length(ibR)],
      row.names=rownames(self$II)))
  else
    self$II <- matrix(ncol=2, nrow=0, dimnames=list(NULL,c('M','U')))

  ## controls
  self$ctl$G <- ibG.nl$c
  self$ctl$R <- ibR.nl$c

  invisible()
}

## Noob background correction for one channel
.backgroundCorrectionNoobCh1 <- function(ib, oob, ctl, offset=15) {
  ## @param ib array of in-band signal
  ## @param oob array of out-of-band-signal
  ## @param ctl control probe signals
  ## @param offset padding for normalized signal
  ## @return normalized in-band signal

  suppressPackageStartupMessages(library(MASS))
  e <- huber(oob)
  mu <- e$mu
  sigma <- e$s
  alpha <- pmax(huber(ib)$mu-mu, 10)
  return(list(i=offset+.normExpSignal(mu, sigma, alpha, ib),
              c=offset+.normExpSignal(mu, sigma, alpha, ctl)))
}

## the following is adapted from Limma
## normal-exponential deconvolution (conditional expectation of xs|xf; WEHI code)
.normExpSignal <- function (mu, sigma, alpha, x)  {
  sigma2 <- sigma * sigma
  if (alpha <= 0)
    stop("alpha must be positive")
  if (sigma <= 0)
    stop("sigma must be positive")
  mu.sf <- x - mu - sigma2/alpha
  signal <- mu.sf + sigma2 * exp(
    dnorm(0, mean = mu.sf, sd = sigma, log = TRUE) -
      pnorm(0, mean = mu.sf, sd = sigma, lower.tail = FALSE, log = TRUE))
  o <- !is.na(signal)
  if (any(signal[o] < 0)) {
    warning("Limit of numerical accuracy reached with very low intensity or very high background:\nsetting adjusted intensities to small value")
    signal[o] <- pmax(signal[o], 1e-06)
  }
  signal
}

.getNormCtls <- function(sset) {
  if (sset$platform == 'hm27') {
    ## controversial to use, maybe the mean of all signals in each channel?
    normctl.G <- sset$ctl[grep('norm.green', tolower(rownames(sset$ctl))),]
    normctl.R <- sset$ctl[grep('norm.red', tolower(rownames(sset$ctl))),]
  } else {                              # hm450 and EPIC
    normctl.G <- sset$ctl[grep('norm_(c|g)',tolower(rownames(sset$ctl))),]
    normctl.R <- sset$ctl[grep('norm_(a|t)',tolower(rownames(sset$ctl))),]
  }
  c(G=mean(normctl.G$G, na.rm=TRUE), R=mean(normctl.R$R, na.rm=TRUE))
}

## #' Correct dye bias using most balanced sample
## #'
## #' In practice, it doesn't matter as long as the reference level does not deviate much
## #' 
## #' @param ssets a list of normalized \code{SignalSet}s
## #' @return a list of normalized \code{SignalSet}s
## #' @export
## DyeBiasCorrectionMostBalanced <- function(ssets) {
##   normctls <- vapply(ssets, .getNormCtls, numeric(2))
##   most.balanced <- which.min(abs(normctls['G',] / normctls['R',] - 1))
##   ref <- mean(normctls[,most.balanced], na.rm=TRUE)
##   lapply(ssets, function(sset) DyeBiasCorrection(sset, ref))
## }

#' Correct dye bias
#'
#' @param sset a \code{SignalSet}s
#' @param ref reference signal level
#' @return a normalized \code{SignalSet}s
.dyeBiasCorrection <- function(ref=7000) {
  normctl <- .getNormCtls(self)
  fR <- ref/normctl['R']
  fG <- ref/normctl['G']

  self$IG <- matrix(c(fG*self$IG[,'M'], fG*self$IG[,'U']),
                    nrow=nrow(self$IG), ncol=ncol(self$IG), dimnames=dimnames(self$IG))
  self$IR <- matrix(c(fR*self$IR[,'M'], fR*self$IR[,'U']),
                    nrow=nrow(self$IR), ncol=ncol(self$IR), dimnames=dimnames(self$IR))
  self$II <- matrix(c(fG*self$II[,'M'], fR*self$II[,'U']),
                    nrow=nrow(self$II), ncol=ncol(self$II), dimnames=dimnames(self$II))
  self$ctl <- transform(self$ctl, G=fG*G, R=fR*R)
  self$oobG <- fG*self$oobG
  self$oobR <- fR*self$oobR
  invisible()
}
