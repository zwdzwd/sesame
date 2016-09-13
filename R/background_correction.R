#' Noob background correction
#' 
#' Norm-Exp deconvolution using Out-Of-Band (oob) probes
#' Note p-values are unchanged (based on the raw signal intensities).
#' @param sset a \code{SignalSet}
#' @param in.place modify \code{SignalSet} in place, faster
#' @param offset offset
#' @return a new \code{SignalSet} with noob background correction
#' @export
noob <- function(sset, in.place=FALSE, offset=15) {

  if (!in.place)
    sset <- sset$clone()
  ## sort signal based on channel
  ibR <- c(sset$IR, sset$II[,'U'])              # in-band red signal
  ibG <- c(sset$IG, sset$II[,'M'])              # in-band green signal

  ## set signal to 1 if 0
  ibR[ibR==0] <- 1
  ibG[ibG==0] <- 1

  ## oobG and oobR are untouched besides the 0>1 switch
  sset$oobR[sset$oobR==0] <- 1
  sset$oobG[sset$oobG==0] <- 1
  
  ## do background correction in each channel
  ibR.nl <- .backgroundCorrectionNoobCh1(ibR, sset$oobR, sset$ctl$R, offset=offset)
  ibG.nl <- .backgroundCorrectionNoobCh1(ibG, sset$oobG, sset$ctl$G, offset=offset)

  ## build back the list
  ## type IG
  if (length(sset$IG)>0)
    sset$IG <- matrix(ibG.nl$i[1:length(sset$IG)],
                      nrow=nrow(sset$IG), dimnames=dimnames(sset$IG))
  else
    sset$IG <- matrix(ncol=2, nrow=0, dimnames=list(NULL,c('M','U')))

  ## type IR
  if (length(sset$IR)>0)
    sset$IR <- matrix(ibR.nl$i[1:length(sset$IR)],
                      nrow=nrow(sset$IR), dimnames=dimnames(sset$IR))
  else
    sset$IR <- matrix(ncol=2, nrow=0, dimnames=list(NULL,c('M','U')))

  ## type II
  if (nrow(sset$II) > 0)
    sset$II <- as.matrix(data.frame(
      M=ibG.nl$i[(length(sset$IG)+1):length(ibG)],
      U=ibR.nl$i[(length(sset$IR)+1):length(ibR)],
      row.names=rownames(sset$II)))
  else
    sset$II <- matrix(ncol=2, nrow=0, dimnames=list(NULL,c('M','U')))

  ## controls
  sset$ctl$G <- ibG.nl$c
  sset$ctl$R <- ibR.nl$c

  ## out-of-band
  sset$oobR <- ibR.nl$o
  sset$oobG <- ibG.nl$o

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
  return(list(i=offset+.normExpSignal(mu, sigma, alpha, ib),
              c=offset+.normExpSignal(mu, sigma, alpha, ctl),
              o=offset+.normExpSignal(mu, sigma, alpha, oob)))
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
      pnorm(0, mean = mu.sf, sd = sigma, lower.tail = FALSE, log.p = TRUE))
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

#' Correct dye bias
#'
#' @param sset a \code{SignalSet}
#' @param ref reference signal level
#' @param in.place modify \code{SignalSet} in place, faster
#' @return a normalized \code{SignalSet}
#' @export
dyeBiasCorr <- function(sset, ref=5000, in.place=FALSE) {

  if (!in.place)
    sset <- sset$clone()
  normctl <- .getNormCtls(sset)
  fR <- ref/normctl['R']
  fG <- ref/normctl['G']

  sset$IG <- matrix(c(fG*sset$IG[,'M'], fG*sset$IG[,'U']),
                    nrow=nrow(sset$IG), ncol=ncol(sset$IG), dimnames=dimnames(sset$IG))
  sset$IR <- matrix(c(fR*sset$IR[,'M'], fR*sset$IR[,'U']),
                    nrow=nrow(sset$IR), ncol=ncol(sset$IR), dimnames=dimnames(sset$IR))
  sset$II <- matrix(c(fG*sset$II[,'M'], fR*sset$II[,'U']),
                    nrow=nrow(sset$II), ncol=ncol(sset$II), dimnames=dimnames(sset$II))
  sset$ctl$G <- fG*sset$ctl$G
  sset$ctl$R <- fR*sset$ctl$R
  sset$oobG <- fG*sset$oobG
  sset$oobR <- fR*sset$oobR
  sset
}


#' Correct dye bias using most balanced sample
#'
#' In practice, it doesn't matter as long as the
#' reference level does not deviate much.
#' 
#' @param ssets a list of normalized \code{SignalSet}s
#' @return a list of normalized \code{SignalSet}s
#' @export
dyeBiasCorrMostBalanced <- function(ssets) {
  normctls <- vapply(ssets, .getNormCtls, numeric(2))
  most.balanced <- which.min(abs(normctls['G',] / normctls['R',] - 1))
  ref <- mean(normctls[,most.balanced], na.rm=TRUE)
  lapply(ssets, function(sset) dyeBiasCorr(sset, ref))
}
