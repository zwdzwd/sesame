#' background subtraction with bleeding-through subtraction
#'
#' Norm-Exp deconvolution using Out-Of-Band (oob) probes.
#' Mean and standard deviations are estimated from cross-channel regression
#' 
#' @param sset a \code{SignalSet}
#' @param in.place modify in place
#' @param offset offset
#' @param detailed if TRUE, return a list of \code{SignalSet} and regression function
#' @return a modified \code{SignalSet} with background correction
#' @export
noobsb <- function(sset, in.place=FALSE, offset=15, detailed=FALSE) {

  if (!in.place)
    sset <- sset$clone()
  
  ## sanitize
  ## sort signal based on channel
  ibR <- c(sset$IR, sset$II[,'U'])              # in-band red signal
  ibR.other.channel <- c(sset$oobG, sset$II[,'M'])
  ibG <- c(sset$IG, sset$II[,'M'])              # in-band green signal
  ibG.other.channel <- c(sset$oobR, sset$II[,'U'])
  ## set signal to 1 if 0
  ibR[ibR==0] <- 1
  ibG[ibG==0] <- 1
  
  ## oobG and oobR are untouched besides the 0>1 switch
  sset$oobR[sset$oobR==0] <- 1
  sset$oobG[sset$oobG==0] <- 1

  ## background correction Red
  # use the other channel to predict background mean
  bg.GpredictR <- train.model.glm(sset$IG, sset$oobR)
  mu.bg.ibR <- bg.GpredictR(ibR.other.channel)
  mu.bg.oobR <- bg.GpredictR(sset$IG)
  mu.bg.ctlR <- bg.GpredictR(sset$ctl$G)
  # parameter estimation
  sigma.bgR <- MASS::huber(sset$oobR)$s # for now, use the global variance
  alphaR <- pmax(MASS::huber(ibR - mu.bg.ibR)$mu, 10) # in-band variance
  # correction
  ibR <- .backgroundCorrCh1(ibR, mu.bg.ibR, sigma.bgR, alphaR, offset=offset)
  oobR <- .backgroundCorrCh1(sset$oobR, mu.bg.oobR, sigma.bgR, alphaR, offset=offset)
  ctlR <- .backgroundCorrCh1(sset$ctl$R, mu.bg.ctlR, sigma.bgR, alphaR, offset=offset)
  
  ## background correction Grn
  # use the other channel to predict background mean
  bg.RpredictG <- train.model.glm(sset$IR, sset$oobG)
  mu.bg.ibG <- bg.RpredictG(ibG.other.channel)
  mu.bg.oobG <- bg.RpredictG(sset$IR)
  mu.bg.ctlG <- bg.RpredictG(sset$ctl$R)
  # parameter estimation
  sigma.bgG <- MASS::huber(sset$oobG)$s # for now, use the global variance
  alphaG <- pmax(MASS::huber(ibG - mu.bg.ibG)$mu, 10) # in-band variance
  # correction
  ibG <- .backgroundCorrCh1(ibG, mu.bg.ibG, sigma.bgG, alphaG, offset=offset)
  oobG <- .backgroundCorrCh1(sset$oobG, mu.bg.oobG, sigma.bgG, alphaG, offset=offset)
  ctlG <- .backgroundCorrCh1(sset$ctl$G, mu.bg.ctlG, sigma.bgG, alphaG, offset=offset)

  ## build back the list
  ## type IG
  if (length(sset$IG)>0)
    sset$IG <- matrix(ibG[1:length(sset$IG)],
                      nrow=nrow(sset$IG), dimnames=dimnames(sset$IG))
  else
    sset$IG <- matrix(ncol=2, nrow=0, dimnames=list(NULL,c('M','U')))

  ## type IR
  if (length(sset$IR)>0)
    sset$IR <- matrix(ibR[1:length(sset$IR)],
                      nrow=nrow(sset$IR), dimnames=dimnames(sset$IR))
  else
    sset$IR <- matrix(ncol=2, nrow=0, dimnames=list(NULL,c('M','U')))

  ## type II
  if (nrow(sset$II) > 0)
    sset$II <- as.matrix(data.frame(
      M=ibG[(length(sset$IG)+1):length(ibG)],
      U=ibR[(length(sset$IR)+1):length(ibR)],
      row.names=rownames(sset$II)))
  else
    sset$II <- matrix(ncol=2, nrow=0, dimnames=list(NULL,c('M','U')))

  ## controls
  sset$ctl$G <- ctlG
  sset$ctl$R <- ctlR

  ## out-of-band
  sset$oobR <- oobR
  sset$oobG <- oobG

  if (detailed)
    list(sset=sset, bg.RpredictG=bg.RpredictG, bg.GpredictR=bg.GpredictR)
  else
    sset
}

## Noob background correction for one channel
## ib array of in-band signal
## oob array of out-of-band-signal
## ctl control probe signals
## offset padding for normalized signal
## return normalized in-band signal
.backgroundCorrCh1 <- function(x, mu.bg, sigma, alpha, offset=15) {
  sigma2 <- sigma * sigma
  if (alpha <= 0)
    stop("alpha must be positive")
  if (sigma <= 0)
    stop("sigma must be positive")
  mu.sf <- x - mu.bg - sigma2/alpha
  signal <- mu.sf + sigma2 * exp(
    dnorm(0, mean = mu.sf, sd = sigma, log = TRUE) -
      pnorm(0, mean = mu.sf, sd = sigma, lower.tail = FALSE, log.p = TRUE))
  o <- !is.na(signal)
  if (any(signal[o] < 0)) {
    warning("Limit of numerical accuracy reached with very low intensity or very high background:\nsetting adjusted intensities to small value")
    signal[o] <- pmax(signal[o], 1e-06)
  }
  signal.min <- min(signal, na.rm = T)
  signal <- signal - signal.min
  offset + signal
}

train.model.glm <- function(input, output) {
  fitdata <- data.frame(IB=as.vector(input)+1,OB=as.vector(output)+1)
  rg <- seq(min(fitdata$IB, na.rm=T), max(fitdata$IB, na.rm=T), length.out=1000)
  fitdata <- do.call(rbind, lapply(split(fitdata, cut(fitdata$IB, rg)), function(x) {
    lims <- quantile(x$OB, c(0.2,0.8))
    x[x$OB>lims[1] & x$OB<lims[2],]
  }))
  m <- glm(OB~IB, data=fitdata, family=gaussian(link='log'))
  function(d) {
    force(m)
    pdt <- predict(m, newdata=data.frame(IB=as.vector(d)), se.fit=T, type='response')
    pdt$fit # + 1000 #pdt$residual.scale * 10
  }
}

subtractBleeding <- function(sset, in.place=FALSE) {
  
  ## this subtract bleeding for type II probes using the same predict function as type I
  library(mgcv)
  if (!in.place)
    sset <- sset$clone()
  
  bg.GpredictR <- train.model.glm(sset$IG, sset$oobR)
  bg.RpredictG <- train.model.glm(sset$IR, sset$oobG)

  # type II
  if (dim(sset$II)[1] > 0) {
    newU <- sset$II[,'U'] - bg.GpredictR(sset$II[,'M'])
    newM <- sset$II[,'M'] - bg.RpredictG(sset$II[,'U'])
    sset$II <- cbind(M=ifelse(newM>0,newM,1), U=ifelse(newU>0,newU,1))
  }
  
  # type I
  bgR <- bg.GpredictR(as.vector(sset$oobG))
  bgG <- bg.RpredictG(as.vector(sset$oobR))
  bg.oobR <- bg.GpredictR(as.vector(sset$IG))
  bg.oobG <- bg.RpredictG(as.vector(sset$IR))
  dim(bgR) <- dim(sset$IR)
  dim(bgG) <- dim(sset$IG)
  dim(bg.oobR) <- dim(sset$oobR)
  dim(bg.oobG) <- dim(sset$oobG)
  sset$IR <- sset$IR - bgR
  sset$IG <- sset$IG - bgG
  sset$oobR <- sset$oobR - bg.oobR
  sset$oobG <- sset$oobG - bg.oobG
  
  sset$IR <- ifelse(sset$IR>0, sset$IR, 1)
  sset$IG <- ifelse(sset$IG>0, sset$IG, 1)
  sset$oobR <- ifelse(sset$oobR>0, sset$oobR, 1)
  sset$oobG <- ifelse(sset$oobG>0, sset$oobG, 1)
  
  sset
}

## subtract bleeding-through separately for type II
subtractBleedingTypeIIsep <- function(sset, in.place=FALSE) {

  library(mgcv)
  if (!in.place)
    sset <- sset$clone()
  
  bg.GpredictR <- train.model.glm(sset$IG, sset$oobR)
  bg.RpredictG <- train.model.glm(sset$IR, sset$oobG)

  ## type II
  ## separate out methylated and unmethylated signal, may fail on extremes like TGCT.
  if (dim(sset$II)[1] > 0) {
    II.meth <- sset$II[sset$II[,'M'] / (sset$II[,'M']+sset$II[,'U']) > 0.8,]
    II.unme <- sset$II[sset$II[,'M'] / (sset$II[,'M']+sset$II[,'U']) < 0.2,]
    bg.GpredictR.II <- train.model.glm(II.meth[,'M'], II.meth[,'U'])
    bg.RpredictG.II <- train.model.glm(II.unme[,'U'], II.unme[,'M'])

    newU <- sset$II[,'U'] - bg.GpredictR.II(sset$II[,'M'])
    newM <- sset$II[,'M'] - bg.RpredictG.II(sset$II[,'U'])
    sset$II <- cbind(M=ifelse(newM>0,newM,1), U=ifelse(newU>0,newU,1))
  } else {
    bg.GpredictR.II <- function(x) x
    bg.RpredictG.II <- function(x) x
  }
  
  # type I
  bgR <- bg.GpredictR(as.vector(sset$oobG))
  bgG <- bg.RpredictG(as.vector(sset$oobR))
  bg.oobR <- bg.GpredictR(as.vector(sset$IG))
  bg.oobG <- bg.RpredictG(as.vector(sset$IR))
  dim(bgR) <- dim(sset$IR)
  dim(bgG) <- dim(sset$IG)
  dim(bg.oobR) <- dim(sset$oobR)
  dim(bg.oobG) <- dim(sset$oobG)
  sset$IR <- sset$IR - bgR
  sset$IG <- sset$IG - bgG
  sset$oobR <- sset$oobR - bg.oobR
  sset$oobG <- sset$oobG - bg.oobG
  
  sset$IR <- ifelse(sset$IR>0, sset$IR, 0)
  sset$IG <- ifelse(sset$IG>0, sset$IG, 0)
  sset$oobR <- ifelse(sset$oobR>0, sset$oobR, 0)
  sset$oobG <- ifelse(sset$oobG>0, sset$oobG, 0)
  
  sset=sset
}


