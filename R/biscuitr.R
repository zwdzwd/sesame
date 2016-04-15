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
  colnames(d) <- c('G', 'R')
  d
}

#' Import IDATs
#'
#' Import IDATs from a sample list
#'
#' @param samples sample list
#' @return a list of signal intensities
#'
#' Each element of the returned list contains a matrix
#' having signal intensity addressed by chip address
readIDATs <- function(sample.names, base.dir=NULL) {
  if (!is.null(base.dir))
    sample.paths <- paste0(base.dir,'/',sample.names)
  else
    sample.paths <- sample.names
  dm <- lapply(sample.paths, readIDAT1)
  names(dm) <- sample.names
  dm
}

#' Import IDATs
#' 
#' Import IDATs from a directory
#' 
#' @param dir directory name.
#' @author Wanding Zhou
#' @return a list of signal intensities
#'
#' Each element of the returned list contains a matrix
#' having signal intensity addressed by chip address
readIDATs.fromdir <- function(dir.name) {
  fns <- list.files(dir.name)
  sample.names <- unique(sub("_(Grn|Red).idat", "", fns[grep(".idat$", fns)]))

  readIDATs(paste0(dir.name,'/',sample.names))
}

#' Import IDATs
#' 
#' Import IDATs from a sample sheet
#' 
#' @param dir directory name.
#' @author Wanding Zhou
#' @return a list of signal intensities
#'
#' Each element of the returned list contains a matrix
#' having signal intensity addressed by chip address
readIDATs.from.samplesheet <- function(sample.sheet, base.dir=NULL) {
  sample.names <- read.csv(sample.sheet, stringsAsFactors=F)
  readIDATs(sample.names$barcode, base.dir=base.dir)
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

#' Lookup address in one sample
#'
#' Lookup address and transform address to probe
#'
#' Translate data in chip address to probe address.
#' Type I probes can be separated into Red and Grn channels. The
#' methylated allele and unmethylated allele are at different
#' addresses. For type II probes methylation allele and unmethylated allele are
#' at the same address. Grn channel is for methylated allele and Red channel is
#' for unmethylated allele. The out-of-band signals are type I probes measured
#' using the "wrong" channel.
#'
#' @param dm data frame in chip address, 2 columns: cy3/Grn and cy5/Red
#' @import FDb.InfiniumMethylation.hg19
#' @return a list with: IR, IG, oobR, oobG, II, ctrl
#' IR, IG and II are matrices with columns M and U.
#' oobR and oobG are arrays of integer.
#' ctrl is a data frame with columns G, R, col, type
chipAddressToProbe <- function(dm) {
  suppressPackageStartupMessages(library("FDb.InfiniumMethylation.hg19"))

  data(hm450.ordering)
  ## type I red channel / oob green channel
  IordR <- subset(hm450.ordering, DESIGN=='I' & col=='R')
  ## 2-channel for red probes' m allele
  IuR2ch <- dm[match(IordR$U, rownames(dm)),]
  IuR <- IuR2ch[,2]
  ## 2-channel for red probes' u allele
  ImR2ch <- dm[match(IordR$M, rownames(dm)),]
  ImR <- ImR2ch[,2]
  oob.G <- c(IuR2ch[,1], ImR2ch[,1])
  signal.I.R <- as.matrix(data.frame(M=ImR, U=IuR, row.names=IordR$Probe_ID))

  ## type I green channel / oob red channel
  IordG <- subset(hm450.ordering, DESIGN=='I' & col=='G')
  ## 2-channel for green probes' M allele
  IuG2ch <- dm[match(IordG$U, rownames(dm)),]
  IuG <- IuG2ch[,1]
  ## 2-channel for green probes' U allele
  ImG2ch <- dm[match(IordG$M, rownames(dm)),]
  ImG <- ImG2ch[,1]
  oob.R <- c(IuG2ch[,2], ImG2ch[,2])
  signal.I.G <- as.matrix(data.frame(M=ImG, U=IuG, row.names=IordG$Probe_ID))

  IIord <- hm450.ordering[hm450.ordering$DESIGN=="II",]
  signal.II <- dm[match(IIord$U, rownames(dm)),]
  colnames(signal.II) <- c('M', 'U')
  rownames(signal.II) <- IIord$Probe_ID

  ## control probes
  data(hm450.controls)
  ctrl <- as.data.frame(dm[match(hm450.controls$Address, rownames(dm)),])
  rownames(ctrl) <- make.names(hm450.controls$Name,unique=T)
  ctrl <- cbind(ctrl,hm450.controls[, c("Color_Channel","Type")])
  colnames(ctrl) <- c('G','R','col','type')

  return(list(IR=signal.I.R, IG=signal.I.G,
              oobR=oob.R, oobG=oob.G, II=signal.II, ctrl=ctrl))
}

#' Compute detection p-value
#'
#' Compute detection p-value for in-band signals
#' @param dmp a list with I, II, oob, ctrl
#' @import stats
detectionPval <- function(dmp) {
  negctrls <- dmp$ctrl[grep('negative', tolower(rownames(dmp$ctrl))),]
  negctrls <- subset(negctrls, col!=-99)

  library(stats)
  funcG <- ecdf(negctrls$G)
  funcR <- ecdf(negctrls$R)

  ## p-value is the minimium detection p-value of the 2 alleles
  IR <- 1-rowMax(cbind(funcR(dmp$IR[,'M']), funcR(dmp$IR[,'U'])))
  IG <- 1-rowMax(cbind(funcG(dmp$IG[,'M']), funcG(dmp$IG[,'U'])))
  II <- 1-rowMax(cbind(funcG(dmp$II[,'M']), funcR(dmp$II[,'U'])))

  names(IR) <- rownames(dmp$IR)
  names(IG) <- rownames(dmp$IG)
  names(II) <- rownames(dmp$II)
  ## return(list(IR=IR, IG=IG, II=II))
  pval <- c(IR,IG,II)
  pval[order(names(pval))]
}

#' Noob background correction
#'
#' background correction using Norm-Exp deconvolution
#' @param dmp a list with I, II, oob, ctrl
#' @import MASS
backgroundCorrNoob <- function(dmp, offset=15) {

  ## sort signal based on channel
  ibR <- c(dmp$IR, dmp$II[,'U'])              # in-band red signal
  ibG <- c(dmp$IG, dmp$II[,'M'])              # in-band green signal

  ## set signal to 1 if 0
  ibR[ibR==0] <- 1
  ibG[ibG==0] <- 1
  dmp$oobR[dmp$oobR==0] <- 1
  dmp$oobG[dmp$oobG==0] <- 1
  
  ## do background correction in each channel
  ibR.nl <- backgroundCorrNoobCh1(ibR, dmp$oobR, dmp$ctrl$R, offset=offset)
  ibG.nl <- backgroundCorrNoobCh1(ibG, dmp$oobG, dmp$ctrl$G, offset=offset)

  ## build back the list
  IR.n <- matrix(ibR.nl$i[1:length(dmp$IR)],
                 nrow=nrow(dmp$IR), dimnames=dimnames(dmp$IR))

  IG.n <- matrix(ibG.nl$i[1:length(dmp$IG)],
                 nrow=nrow(dmp$IG), dimnames=dimnames(dmp$IG))

  II <- as.matrix(data.frame(
    M=ibG.nl$i[(length(dmp$IG)+1):length(ibG)],
    U=ibR.nl$i[(length(dmp$IR)+1):length(ibR)],
    row.names=rownames(dmp$II)))

  ctrl <- dmp$ctrl
  ctrl$G <- ibG.nl$c
  ctrl$R <- ibR.nl$c

  return(list(IR=IR.n, IG=IG.n,
              oobR=dmp$oobR, oobG=dmp$oobG, II=II, ctrl=ctrl))
}

#' Noob background correction for 1 channel
#'
#' Noob background correction for 1 channel
#' @param ib array of in-band signal
#' @param oob array of out-of-band-signal
#' @return normalized in-band signal
backgroundCorrNoobCh1 <- function(ib, oob, ctrl, offset) {
  suppressPackageStartupMessages(library(MASS))
  e <- huber(oob)
  mu <- e$mu
  ## sigma <- log(e$s)
  ## alpha <- log(pmax(huber(ib)$mu-mu, 10))
  sigma <- e$s
  alpha <- pmax(huber(ib)$mu-mu, 10)
  ## cat(mu, sigma, alpha)
  return(list(i=offset+normexp.signal(mu, sigma, alpha, ib),
              c=offset+normexp.signal(mu, sigma, alpha, ctrl)))
}

## the following from Limma
## normal-exponential deconvolution (conditional expectation of xs|xf; WEHI code)
normexp.signal <- function (mu, sigma, alpha, x)  {
  ## par = unlist(par)
  ## mu <- par[1]
  ## sigma <- exp(par[2])
  ## sigma <- exp(lsigma)
  sigma2 <- sigma * sigma
  ## alpha <- exp(lalpha)
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

get.normctls <- function(dmp) {
  normctl.G <- dmp$ctrl[grep('norm_(c|g)',tolower(rownames(dmp$ctrl))),]
  normctl.R <- dmp$ctrl[grep('norm_(a|t)',tolower(rownames(dmp$ctrl))),]
  c(G=mean(normctl.G$G), R=mean(normctl.R$R))
}

#' Correct dye bias with most balanced sample
#'
#' Correct dye bias using normalization controls
#' @param dmps list of probe-addressed signals
#' @return normalized probe-addressed signals
dyeBiasCorr.most.balanced <- function(dmps) {
  normctls <- sapply(dmps, get.normctls)
  most.balanced <- which.min(abs(normctls['G',] / normctls['R',] - 1))
  ref <- mean(normctls[,most.balanced])
  dyeBiasCorr(dmps, ref, normctls=normctls)
}

#' Correct dye bias
#'
#' Correct dye bias
#'
#' @param dmps list of probe-addressed signals
#' @param ref reference signal level
#' @return normalized probe-addressed signals
dyeBiasCorr <- function(dmps, ref, normctls=NULL) {
  if (is.null(normctls))
    normctls <- sapply(dmps, get.normctls)
  factor <- ref / normctls
  dmps.n <- lapply(names(dmps), function(sample.name) {
    n <- list()
    dmp <- dmps[[sample.name]]
    fR <- factor['R', sample.name]
    fG <- factor['G', sample.name]
    n$IR <- matrix(c(fR*dmp$IR[,'M'], fR*dmp$IR[,'U']),
                   nrow=nrow(dmp$IR), dimnames=dimnames(dmp$IR))
    n$IG <- matrix(c(fG*dmp$IG[,'M'], fG*dmp$IG[,'U']),
                   nrow=nrow(dmp$IG), dimnames=dimnames(dmp$IG))
    n$II <- matrix(c(fG*dmp$II[,'M'], fR*dmp$II[,'U']),
                   nrow=nrow(dmp$II), dimnames=dimnames(dmp$II))
    n
  })
  names(dmps.n) <- names(dmps)
  dmps.n
}

#' Convert signal to beta
#'
#' Convert signal to beta
#'
#' @param dmp a list which corresponds to the probe signal
#' @return a list which corresponds to the beta value
probeSignalToBeta <- function(dmp, pval) {
  betas <- pmax(dmp$IR[,'M'],1) / pmax(dmp$IR[,'M']+dmp$IR[,'U'],2)
  betas <- c(betas,
             pmax(dmp$IG[,'M'],1) / pmax(dmp$IG[,'M']+dmp$IG[,'U'],2))
  betas <- c(betas,
             pmax(dmp$II[,'M'],1) / pmax(dmp$II[,'M']+dmp$II[,'U'],2))
  ## betas[c(pval$IR, pval$IG, pval$II)>0.05] <- NA
  betas[pval[names(betas)]>0.05] <- NA
  betas[order(names(betas))]
}

mprint <- function(...) {
  cat('[', as.character(Sys.time()),'] ', ..., '\n', sep='')
}
