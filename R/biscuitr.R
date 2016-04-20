#' @title
#' Analyze DNA methylation data
#' 
#' @description
#' biscuitR allows analysing DNA methylation data in R.
#' 
#' @details
#' This package complements array functionalities in complement to the biscuit software.
#' This implementation allow processing >10,000 samples in parallel on clusters.
#' @aliases biscuitr
#' @author
#' Wanding Zhou \email{Wanding.Zhou@vai.org},
#' Hui Shen \email{Hui.Shen@vai.org}
#' 
## @references To appear
## @seealso To appear
#' @examples
#' library(biscuitr)
#' 
#' ## read IDATs
#' dms <- ReadIDATsFromDir(sample.dir)
#'
#' ## translate chip address to probe address
#' dmps <- lapply(dmps, ChipAddressToProbe)
#'
#' ## detect p-value
#' pvals <- lapply(dmps, DetectPValue)
#'
#' ## normalization
#' dmps <- lapply(dmps, BackgroundCorrectionNoob)
#' dmps <- DyeBiasCorrectionMostBalanced(dmps)
#' dmps <- Funnorm(dmps)
#'
#' ## convert signal to beta values
#' betas <- Map(SignalToBeta, dmps, pvals)
#'
#' ## mask repeat and snp
#' betas <- MaskRepeatSnpHM450(betas)
#' 
#' @keywords DNAMethylation Microarray QualityControl
#' 
"_PACKAGE"

#' SignalSet
#' 
#' Construct a SignalSet object
#' @param IG intensity table for type I probes in green channel
#' @param IR intensity table for type I probes in red channel
#' @param II intensity table for type II probes
#' @param oobG out-of-band probes in green channel
#' @param oobR out-of-band probes in red channel
#' @param ctl all the control probe intensities
#' 
#' @return a SignalSet with IG, IR, II, oobG, oobR and ctl
#'   \item{IR, IG and II}{matrices with columns M and U}
#'   \item{oobR and oobG}{arrays of integer}
#'   \item{ctl}{a data frame with columns G, R, col, type}
#' @export
SignalSet <- function(IG=NULL, IR=NULL, II=NULL, oobG=NULL, oobR=NULL, ctl=NULL) {
  structure(list(IG=IG, IR=IR, II=II, oobG=oobG, oobR=oobR, ctl=ctl), class="SignalSet")
}

#' Select a chromosome
#'
#' Select a subset of signal set from a given chromosome
#'
#' @param dmp an object of class \code{SignalSet}
#' @param chrm target chromosome
#' @return an object of class \code{SignalSet}
#' @export
SelectChromosome <- function(dmp, chrm) {
  data(hm450.hg19.probe2chr)
  SignalSet(IG = dmp$IG[hm450.hg19.probe2chr[rownames(dmp$IG)] == chrm,],
            IR = dmp$IR[hm450.hg19.probe2chr[rownames(dmp$IR)] == chrm,],
            II = dmp$II[hm450.hg19.probe2chr[rownames(dmp$II)] == chrm,])
}

#' Merge two signal sets
#'
#' Merge two signal sets
#'
#' @param dmp1 an object of class \code{SignalSet}
#' @param dmp2 an object of class \code{SignalSet}
#' @return an object of class \code{SignalSet} after merging
MergeSignals <- function(dmp1, dmp2) {
  SignalSet(IG = rbind(dmp1$IG, dmp2$IG),
            IR = rbind(dmp1$IR, dmp2$IR),
            II = rbind(dmp1$II, dmp2$II))
}

#' Import one IDAT file
#'
#' Import one IDAT file
#' 
#' @param fn IDAT file name
#' @import illuminaio
#' @return a data frame with 2 columns, corresponding to cy3 and cy5 color channel signal
#' @export
ReadIDAT1 <- function(idat.name) {
  library(illuminaio)
  ida.grn <- readIDAT(paste0(idat.name,"_Grn.idat"));
  ida.red <- readIDAT(paste0(idat.name,"_Red.idat"));
  d <- cbind(cy3=ida.grn$Quants[,"Mean"], cy5=ida.red$Quants[,"Mean"])
  colnames(d) <- c('G', 'R')
  d
}

#' Import IDATs from a list of samples
#'
#' Import IDATs from a list of samples
#' 
#' Each element of the returned list contains a matrix
#' having signal intensity addressed by chip address
#' 
#' @param samples sample list
#' @param base.dir base directory
#' @return a list of data frames. Each data frame corresponds to a sample.
#' Data are signal intensities indexed by chip address.
#' @export
ReadIDATs <- function(sample.names, base.dir=NULL) {
  if (!is.null(base.dir))
    sample.paths <- paste0(base.dir,'/',sample.names)
  else
    sample.paths <- sample.names
  dm <- lapply(sample.paths, ReadIDAT1)
  names(dm) <- basename(sample.names)
  dm
}

#' Import IDATs from a directory
#' 
#' Import IDATs from a directory
#' 
#' Each element of the returned list contains a matrix
#' having signal intensity addressed by chip address
#' 
#' @param dir directory name.
#' @return a list of data frames. Each data frame corresponds to a sample.
#' Data are signal intensities indexed by chip address.
#' @export
ReadIDATsFromDir <- function(dir.name) {
  fns <- list.files(dir.name)
  sample.names <- unique(sub("_(Grn|Red).idat", "", fns[grep(".idat$", fns)]))

  ReadIDATs(paste0(dir.name,'/',sample.names))
}

#' Import IDATs from a sample sheet
#' 
#' Import IDATs from a sample sheet
#' 
#' Each element of the returned list contains a matrix
#' having signal intensity addressed by chip address
#' 
#' @param sample.sheet.path path to sample sheet
#' @param base.dir directory on which the \code{sample.sheet.path} is based
#' @return a list of data frames. Each data frame corresponds to a sample. 
#' Data are signal intensities indexed by chip address.
#' @export
ReadIDATsFromSampleSheet <- function(sample.sheet, base.dir=NULL) {
  sample.names <- read.csv(sample.sheet, stringsAsFactors=F)
  ReadIDATs(sample.names$barcode, base.dir=base.dir)
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
#' using the other channel.
#'
#' @param dm data frame in chip address, 2 columns: cy3/Grn and cy5/Red
#' @return a SignalSet, indexed by probe ID address
#' @import FDb.InfiniumMethylation.hg19
#' @export
ChipAddressToProbe <- function(dm) {

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
  ctl <- as.data.frame(dm[match(hm450.controls$Address, rownames(dm)),])
  rownames(ctl) <- make.names(hm450.controls$Name,unique=T)
  ctl <- cbind(ctl,hm450.controls[, c("Color_Channel","Type")])
  colnames(ctl) <- c('G','R','col','type')

  SignalSet(IR=signal.I.R, IG=signal.I.G, oobR=oob.R, oobG=oob.G, II=signal.II, ctl=ctl)
}

#' Compute detection p-value
#'
#' Compute detection p-value for in-band signals
#' @param dmp a SignalSet
#' @import stats
#' @return array of p-values for each probe
#' @export
DetectPValue <- function(dmp) {
  negctls <- dmp$ctl[grep('negative', tolower(rownames(dmp$ctl))),]
  negctls <- subset(negctls, col!=-99)

  library(stats)
  funcG <- ecdf(negctls$G)
  funcR <- ecdf(negctls$R)

  ## p-value is the minimium detection p-value of the 2 alleles
  IR <- 1-apply(cbind(funcR(dmp$IR[,'M']), funcR(dmp$IR[,'U'])),1,max)
  IG <- 1-apply(cbind(funcG(dmp$IG[,'M']), funcG(dmp$IG[,'U'])),1,max)
  II <- 1-apply(cbind(funcG(dmp$II[,'M']), funcR(dmp$II[,'U'])),1,max)

  names(IR) <- rownames(dmp$IR)
  names(IG) <- rownames(dmp$IG)
  names(II) <- rownames(dmp$II)
  
  pval <- c(IR,IG,II)
  pval[order(names(pval))]
}

#' Noob background correction
#'
#' Norm-Exp deconvolution using Out-Of-Band (oob) probes
#' @param dmp a \code{SignalSet}
#' @import MASS
#' @return the normalized \code{SignalSet}
#' @export
BackgroundCorrectionNoob <- function(dmp, offset=15) {

  ## sort signal based on channel
  ibR <- c(dmp$IR, dmp$II[,'U'])              # in-band red signal
  ibG <- c(dmp$IG, dmp$II[,'M'])              # in-band green signal

  ## set signal to 1 if 0
  ibR[ibR==0] <- 1
  ibG[ibG==0] <- 1
  dmp$oobR[dmp$oobR==0] <- 1
  dmp$oobG[dmp$oobG==0] <- 1
  
  ## do background correction in each channel
  ibR.nl <- .BackgroundCorrectionNoobCh1(ibR, dmp$oobR, dmp$ctl$R, offset=offset)
  ibG.nl <- .BackgroundCorrectionNoobCh1(ibG, dmp$oobG, dmp$ctl$G, offset=offset)

  ## build back the list
  IR.n <- matrix(ibR.nl$i[1:length(dmp$IR)],
                 nrow=nrow(dmp$IR), dimnames=dimnames(dmp$IR))

  IG.n <- matrix(ibG.nl$i[1:length(dmp$IG)],
                 nrow=nrow(dmp$IG), dimnames=dimnames(dmp$IG))

  II <- as.matrix(data.frame(
    M=ibG.nl$i[(length(dmp$IG)+1):length(ibG)],
    U=ibR.nl$i[(length(dmp$IR)+1):length(ibR)],
    row.names=rownames(dmp$II)))

  ctl <- dmp$ctl
  ctl$G <- ibG.nl$c
  ctl$R <- ibR.nl$c

  SignalSet(IR=IR.n, IG=IG.n, oobR=dmp$oobR, oobG=dmp$oobG, II=II, ctl=ctl)
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

.GetNormCtls <- function(dmp) {
  normctl.G <- dmp$ctl[grep('norm_(c|g)',tolower(rownames(dmp$ctl))),]
  normctl.R <- dmp$ctl[grep('norm_(a|t)',tolower(rownames(dmp$ctl))),]
  c(G=mean(normctl.G$G, na.rm=TRUE), R=mean(normctl.R$R, na.rm=TRUE))
}

#' Correct dye bias using most balanced sample
#'
#' Correct dye bias using most balanced sample
#' 
#' @param dmps a list of normalized \code{SignalSet}s
#' @return a list of normalized \code{SignalSet}s
#' @export
DyeBiasCorrectionMostBalanced <- function(dmps) {
  normctls <- vapply(dmps, .GetNormCtls, numeric(2))
  most.balanced <- which.min(abs(normctls['G',] / normctls['R',] - 1))
  ref <- mean(normctls[,most.balanced], na.rm=TRUE)
  DyeBiasCorrection(dmps, ref, normctls=normctls)
}

#' Correct dye bias
#'
#' Correct dye bias
#'
#' @param dmps a list of \code{SignalSet}s
#' @param ref reference signal level
#' @return a list of normalized \code{SignalSet}s
#' @export
DyeBiasCorrection <- function(dmps, ref, normctls=NULL) {
  if (is.null(normctls))
    normctls <- vapply(dmps, .GetNormCtls, numeric(2))
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
#' Convert signal to beta and mask low p-value probes
#'
#' @param dmp a \code{SignalSet}
#' @return beta values
#' @export
SignalToBeta <- function(dmp, pval) {
  betas <- pmax(dmp$IR[,'M'],1) / pmax(dmp$IR[,'M']+dmp$IR[,'U'],2)
  betas <- c(betas,
             pmax(dmp$IG[,'M'],1) / pmax(dmp$IG[,'M']+dmp$IG[,'U'],2))
  betas <- c(betas,
             pmax(dmp$II[,'M'],1) / pmax(dmp$II[,'M']+dmp$II[,'U'],2))
  ## betas[c(pval$IR, pval$IG, pval$II)>0.05] <- NA
  betas[pval[names(betas)]>0.05] <- NA
  betas[order(names(betas))]
}

#' @export
MPrint <- function(...) {
  cat('[', as.character(Sys.time()),'] ', ..., '\n', sep='')
}

#' Mask SNP and repeat HM450
#'
#' Mask SNP and repeat for HM450
#' @param betas a matrix of betas
#' @return masked matrix of beta values
#' @export
MaskRepeatSnpHM450 <- function(betas) {
  data(hm450.mask)
  betas[names(hm450.mask),] <- NA
  MPrint('Masked ',length(hm450.mask), ' probes.')
  betas
}
