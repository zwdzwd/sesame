#' Analyze DNA methylation data
#' 
#' biscuitR allows analysing DNA methylation data in R.
#' 
#' This package complements the biscuit software.
#' @aliases biscuitr
#' @author Wanding Zhou
#' @references To appear
#' @seealso To appear
#' @examples To appear
"_PACKAGE"
#> [1] "_PACKAGE"


#samples <- read.csv("samples1k.csv", stringsAsFactors=FALSE)
#samples <- read.csv("samples.csv", stringsAsFactors=FALSE)
## options("mc.cores"=20)
## options("mc.preschedule"=FALSE)

#' Import 1 IDAT file
#'
#' Import 1 IDAT file
#' 
#' @param fn IDAT file name
#' @import illuminaio
#' @return a data frame with 2 columns, corresponding to cy3 and cy5 color channel signal
readIDAT1 <- function(idat.name) {
  library(illuminaio)
  ida.grn <- readIDAT(paste0(idat.name,"_Grn.idat"));
  ida.red <- readIDAT(paste0(idat.name,"_Red.idat"));
  d <- cbind(cy3=ida.grn$Quants[,"Mean"], cy5=ida.red$Quants[,"Mean"])
  colnames(d) <- c(paste0(idat.name,'.cy3'), paste0(idat.name, '.cy5'))
  d
}

#' Import all IDATs
#' 
#' Import all IDATs from a directory
#' 
#' @param dir directory name.
#' @author Wanding Zhou
#' @return a data frame with signal intensity at each address
readAllIDATs <- function(dir.name) {
  fns <- list.files()
  sample.names <- unique(sub("_(Grn|Red).idat", "", fns[grep(".idat$", fns)]))

  do.call('cbind', lapply(sample.names, readIDAT1))
}

batchReadIDATs <- function(dir.name) {

  fns <- list.files()
  for (j in 0:10) {
    i <- 1
    dms <- do.call('cbind', lapply(samples$barcode[(j*1000+1):min((j+1)*1000,length(samples$barcode))], function(x) {
      i <<- i+1;
      if (i%%100==0) {
        cat(i,"\n");
        gc();
      }
      ida.grn <- readIDAT(paste0(x,"_Grn.idat"))
      ida.red <- readIDAT(paste0(x,"_Red.idat"))
      dm1 <- cbind(cy3=ida.grn$Quants[,"Mean"], cy5=ida.red$Quants[,"Mean"])
      colnames(dm1) <- c(paste0(x,'.cy3'), paste0(x, '.cy5'))
      dm1
    }))
    save(dms, file=paste0("datam.",j,".rda"))
  }
}

## data(hm450.controls)
#' Lookup address
#'
#' Lookup address and transform address to probe
#'
#' Translate data in chip address to probe address
#'
#' @param dm data frame in chip address, 2 columns: cy3/Grn and cy5/Red
#' @import FDb.InfiniumMethylation.hg19
#' @return a list with I, II, oob
#' I and II are matrices with columns M and U. oob is an array of integer
chipAddressToProbe <- function(dm) {
  library("FDb.InfiniumMethylation.hg19")
  data(hm450.ordering)

  Iord <- hm450.ordering[hm450.ordering$DESIGN=="I",]
  Iu2ch <- dm[match(Iord$U, rownames(dm)),] # typeI 2-channel
  Iu <- ifelse(Iord$col=='R', Iu2ch[,2], Iu2ch[,1])
  oob.u <- ifelse(Iord$col=='R', Iu2ch[,1], Iu2ch[,2])
  Im2ch <- dm[match(Iord$M, rownames(dm)),] # typeII 2-channel
  Im <- ifelse(Iord$col=='R', Im2ch[,2], Im2ch[,1])
  oob.m <- ifelse(Iord$col=='R', Im2ch[,1], Im2ch[,2])
  signal.I <- as.matrix(data.frame(M=Im, U=Iu, row.names=Iord$Probe_ID))

  IIord <- hm450.ordering[hm450.ordering$DESIGN=="II",]
  signal.II <- dm[match(IIord$U, rownames(dm)),]
  colnames(signal.II) <- c('M', 'U')
  rownames(signal.II) <- IIord$Probe_ID

  signal.oob <- c(oob.u, oob.m)

  return(list(I=signal.I, II=signal.II, oob=signal.oob))
}

#' Noob background correction in each channel
#'
#' background correction using Norm-Exp deconvolution
#' @param dmp a list with I, II, oob
#' @import MASS
backgroundCorrectionNoob <- function(dmp) {
  ib <- c(dmp$I, dmp$II)                # in-band signal
  library(MASS)
  ests <- huber(dmp$oob)
  mu <- ests$mu
  sigma <- ests$s
  alpha <- pmax(huber(ib)$mu-mu, 10)
  normexp.signal(mu, sigma, alpha, ib)
}

## the following from Limma
## normal-exponential deconvolution (conditional expectation of xs|xf; WEHI code)
normexp.signal <- function (mu, sigma, alpha, x)  {
  par = unlist(par)
  mu <- par[1]
  ## sigma <- exp(par[2])
  sigma2 <- sigma * sigma
  ## alpha <- exp(par[3])
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
