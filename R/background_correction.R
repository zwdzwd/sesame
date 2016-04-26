
#' Noob background correction
#'
#' Norm-Exp deconvolution using Out-Of-Band (oob) probes
#' @param sset a \code{SignalSet}
#' @import MASS
#' @return the normalized \code{SignalSet}
#' @export
BackgroundCorrectionNoob <- function(sset, offset=15) {

  ## sort signal based on channel
  ibR <- c(sset$IR, sset$II[,'U'])              # in-band red signal
  ibG <- c(sset$IG, sset$II[,'M'])              # in-band green signal

  ## set signal to 1 if 0
  ibR[ibR==0] <- 1
  ibG[ibG==0] <- 1
  sset$oobR[sset$oobR==0] <- 1
  sset$oobG[sset$oobG==0] <- 1
  
  ## do background correction in each channel
  ibR.nl <- .BackgroundCorrectionNoobCh1(ibR, sset$oobR, sset$ctl$R, offset=offset)
  ibG.nl <- .BackgroundCorrectionNoobCh1(ibG, sset$oobG, sset$ctl$G, offset=offset)

  ## build back the list
  IR.n <- matrix(ibR.nl$i[1:length(sset$IR)],
                 nrow=nrow(sset$IR), dimnames=dimnames(sset$IR))

  IG.n <- matrix(ibG.nl$i[1:length(sset$IG)],
                 nrow=nrow(sset$IG), dimnames=dimnames(sset$IG))

  II <- as.matrix(data.frame(
    M=ibG.nl$i[(length(sset$IG)+1):length(ibG)],
    U=ibR.nl$i[(length(sset$IR)+1):length(ibR)],
    row.names=rownames(sset$II)))

  ctl <- sset$ctl
  ctl$G <- ibG.nl$c
  ctl$R <- ibR.nl$c

  SignalSet(sset$platform,
            IG=IG.n, IR=IR.n,
            oobG=sset$oobG, oobR=sset$oobR, # oobG and oobR are untouched
            II=II, ctl=ctl)
}

## Noob background correction for one channel
.BackgroundCorrectionNoobCh1 <- function(ib, oob, ctl, offset=15) {
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
  return(list(i=offset+.NormExpSignal(mu, sigma, alpha, ib),
              c=offset+.NormExpSignal(mu, sigma, alpha, ctl)))
}

## the following is adapted from Limma
## normal-exponential deconvolution (conditional expectation of xs|xf; WEHI code)
.NormExpSignal <- function (mu, sigma, alpha, x)  {
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

.GetNormCtls <- function(sset) {
  normctl.G <- sset$ctl[grep('norm_(c|g)',tolower(rownames(sset$ctl))),]
  normctl.R <- sset$ctl[grep('norm_(a|t)',tolower(rownames(sset$ctl))),]
  c(G=mean(normctl.G$G, na.rm=TRUE), R=mean(normctl.R$R, na.rm=TRUE))
}

#' Correct dye bias using most balanced sample
#'
#' Correct dye bias using most balanced sample
#' 
#' @param ssets a list of normalized \code{SignalSet}s
#' @return a list of normalized \code{SignalSet}s
#' @export
DyeBiasCorrectionMostBalanced <- function(ssets) {
  normctls <- vapply(ssets, .GetNormCtls, numeric(2))
  most.balanced <- which.min(abs(normctls['G',] / normctls['R',] - 1))
  ref <- mean(normctls[,most.balanced], na.rm=TRUE)
  DyeBiasCorrection(ssets, ref, normctls=normctls)
}

#' Correct dye bias
#'
#' Correct dye bias
#'
#' @param ssets a list of \code{SignalSet}s
#' @param ref reference signal level
#' @return a list of normalized \code{SignalSet}s
#' @export
DyeBiasCorrection <- function(ssets, ref, normctls=NULL) {
  if (is.null(normctls))       # recompute normal control if not given
    normctls <- vapply(ssets, .GetNormCtls, numeric(2))
  factor <- ref / normctls
  ssets.n <- lapply(names(ssets), function(sample.name) {
    n <- list()
    sset <- ssets[[sample.name]]
    fR <- factor['R', sample.name]
    fG <- factor['G', sample.name]
    n$IR <- matrix(c(fR*sset$IR[,'M'], fR*sset$IR[,'U']),
                   nrow=nrow(sset$IR), dimnames=dimnames(sset$IR))
    n$IG <- matrix(c(fG*sset$IG[,'M'], fG*sset$IG[,'U']),
                   nrow=nrow(sset$IG), dimnames=dimnames(sset$IG))
    n$II <- matrix(c(fG*sset$II[,'M'], fR*sset$II[,'U']),
                   nrow=nrow(sset$II), dimnames=dimnames(sset$II))
    n$ctl <- transform(sset$ctl, G=fG*G, R=fR*R)

    ## not sure if we need to normalize oobs, it is dye-bias corrected anyway
    n$oobG <- fG*sset$oobG
    n$oobR <- fR*sset$oobR
    n
  })
  names(ssets.n) <- names(ssets)
  ssets.n
}
