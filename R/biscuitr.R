#' @title
#' Analyze DNA methylation data
#' 
#' @description
#' biscuitR allows analysing DNA methylation data in R.
#' 
#' @details
#' This package complements array functionalities in complement to the biscuit software.
#' This implementation allow processing >10,000 samples in parallel on clusters.
#' It supports both HM450 and EPIC platform.
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
#' ssets <- lapply(ssets, ChipAddressToSignal)
#'
#' ## detect p-value
#' pvals <- lapply(ssets, DetectPValue)
#'
#' ## normalization
#' ssets <- lapply(ssets, BackgroundCorrectionNoob)
#' ssets <- DyeBiasCorrectionMostBalanced(ssets)
#' ssets <- Funnorm(ssets)
#'
#' ## convert signal to beta values
#' betas <- mapply(SignalToBeta, ssets, pvals)
#'
#' ## mask repeats and SNPs
#' betas <- MaskRepeatSNPs(betas, 'hm450')
#' 
#' @keywords DNAMethylation Microarray QualityControl
#' 
"_PACKAGE"

#' SignalSet
#' 
#' Construct a SignalSet object
#' @param platform character() for platform of the array
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
SignalSet <- function(platform,
                      IG=NULL, IR=NULL,
                      II=NULL,
                      oobG=NULL, oobR=NULL,
                      ctl=NULL) {
  sset <- list(IG=IG, IR=IR,
               II=II,
               oobG=oobG, oobR=oobR,
               ctl=ctl,
               platform=platform)
  class(sset) <- "SignalSet"
  sset
}

#' Select a chromosome
#'
#' Select a subset of signal set from a given chromosome
#'
#' @param sset an object of class \code{SignalSet}
#' @param chrm target chromosome
#' @return an object of class \code{SignalSet}
#' @export
SelectChromosome <- function(sset, chrm) {
  probe2chr <- GetBuiltInData('hg19.probe2chr', sset$platform)
  SignalSet(sset$platform,
            IG = sset$IG[probe2chr[rownames(sset$IG)] == chrm,],
            IR = sset$IR[probe2chr[rownames(sset$IR)] == chrm,],
            II = sset$II[probe2chr[rownames(sset$II)] == chrm,])
}

#' Merge two signal sets
#'
#' Merge two signal sets
#'
#' @param sset1 an object of class \code{SignalSet}
#' @param sset2 an object of class \code{SignalSet}
#' @return an object of class \code{SignalSet} after merging
#' @export
MergeSignals <- function(sset1, sset2) {
  SignalSet(sset1$platform,
            IG = rbind(sset1$IG, sset2$IG),
            IR = rbind(sset1$IR, sset2$IR),
            II = rbind(sset1$II, sset2$II))
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
  chip.type <- switch(
    ida.red$ChipType,
    'BeadChip 8x5'='EPIC',
    'BeadChip 12x8'='hm450',
    'BeadChip 12x1'='hm27')
  attr(d, 'platform') <- chip.type
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
ChipAddressToSignal <- function(dm) {

  platform <- attr(dm, 'platform')
  dm.ordering <- GetBuiltInData('ordering', platform)

  ## type I red channel / oob green channel
  IordR <- subset(dm.ordering, DESIGN=='I' & col=='R')
  ## 2-channel for red probes' m allele
  IuR2ch <- dm[match(IordR$U, rownames(dm)),]
  IuR <- IuR2ch[,2]
  ## 2-channel for red probes' u allele
  ImR2ch <- dm[match(IordR$M, rownames(dm)),]
  ImR <- ImR2ch[,2]
  oob.G <- c(IuR2ch[,1], ImR2ch[,1])
  signal.I.R <- as.matrix(data.frame(M=ImR, U=IuR, row.names=IordR$Probe_ID))

  ## type I green channel / oob red channel
  IordG <- subset(dm.ordering, DESIGN=='I' & col=='G')
  ## 2-channel for green probes' M allele
  IuG2ch <- dm[match(IordG$U, rownames(dm)),]
  IuG <- IuG2ch[,1]
  ## 2-channel for green probes' U allele
  ImG2ch <- dm[match(IordG$M, rownames(dm)),]
  ImG <- ImG2ch[,1]
  oob.R <- c(IuG2ch[,2], ImG2ch[,2])
  signal.I.G <- as.matrix(data.frame(M=ImG, U=IuG, row.names=IordG$Probe_ID))

  IIord <- dm.ordering[dm.ordering$DESIGN=="II",]
  signal.II <- dm[match(IIord$U, rownames(dm)),]
  colnames(signal.II) <- c('M', 'U')
  rownames(signal.II) <- IIord$Probe_ID

  ## control probes
  dm.controls <- GetBuiltInData('controls', platform)

  ctl <- as.data.frame(dm[match(dm.controls$Address, rownames(dm)),])
  rownames(ctl) <- make.names(dm.controls$Name,unique=T)
  ctl <- cbind(ctl, dm.controls[, c("Color_Channel","Type")])
  colnames(ctl) <- c('G','R','col','type')

  SignalSet(platform,
            IG = signal.I.G, IR=signal.I.R,
            oobG = oob.G, oobR=oob.R,
            II = signal.II, ctl=ctl)
}

#' Compute detection p-value
#'
#' Compute detection p-value for in-band signals
#' @param sset a SignalSet
#' @import stats
#' @return array of p-values for each probe
#' @export
DetectPValue <- function(sset) {
  negctls <- sset$ctl[grep('negative', tolower(rownames(sset$ctl))),]
  negctls <- subset(negctls, col!=-99)

  library(stats)
  funcG <- ecdf(negctls$G)
  funcR <- ecdf(negctls$R)

  ## p-value is the minimium detection p-value of the 2 alleles
  IR <- 1-apply(cbind(funcR(sset$IR[,'M']), funcR(sset$IR[,'U'])),1,max)
  IG <- 1-apply(cbind(funcG(sset$IG[,'M']), funcG(sset$IG[,'U'])),1,max)
  II <- 1-apply(cbind(funcG(sset$II[,'M']), funcR(sset$II[,'U'])),1,max)

  names(IR) <- rownames(sset$IR)
  names(IG) <- rownames(sset$IG)
  names(II) <- rownames(sset$II)
  
  pval <- c(IR,IG,II)
  pval[order(names(pval))]
}

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
            oobG=sset$oobG, oobR=sset$oobR,
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
  if (is.null(normctls))
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
    n
  })
  names(ssets.n) <- names(ssets)
  ssets.n
}

#' Convert signal to beta
#'
#' Convert signal to beta
#'
#' Convert signal to beta and mask low p-value probes
#'
#' @param sset a \code{SignalSet}
#' @return beta values
#' @export
SignalToBeta <- function(sset, pval) {
  betas1 <- pmax(sset$IR[,'M'],1) / pmax(sset$IR[,'M']+sset$IR[,'U'],2)
  betas2 <- pmax(sset$IG[,'M'],1) / pmax(sset$IG[,'M']+sset$IG[,'U'],2)
  betas3 <- pmax(sset$II[,'M'],1) / pmax(sset$II[,'M']+sset$II[,'U'],2)
  betas <- c(betas1, betas2, betas3)
  ## betas[c(pval$IR, pval$IG, pval$II)>0.05] <- NA
  betas[pval[names(betas)]>0.05] <- NA
  betas[order(names(betas))]
}

#' @export
MPrint <- function(...) {
  cat('[', as.character(Sys.time()),'] ', ..., '\n', sep='')
}

#' Mask SNPs and repeats
#'
#' Mask SNPs and repeats
#' @param betas a matrix of betas
#' @return masked matrix of beta values
#' @export
MaskRepeatSNPs <- function(betas, platform) {
  dm.mask <- GetBuiltInData('mask', platform)
  betas[names(dm.mask),] <- NA
  MPrint('Masked ',length(dm.mask), ' probes.')
  betas
}
