#' @title
#' Analyze DNA methylation data
#' 
#' @description
#' SEnsible and step-wise analysis of DNA methylation data
#' 
#' @details
#' This package complements array functionalities that allow
#' processing >10,000 samples in parallel on clusters.
#' It supports both HM450 and EPIC platform.
#' @aliases sesame
#' @author
#' Wanding Zhou \email{Wanding.Zhou@vai.org},
#' Hui Shen \email{Hui.Shen@vai.org}
#' 
## @references To appear
## @seealso To appear
#' @examples
#' library(sesame)
#' 
#' ## read IDATs
#' dms <- ReadIDATsFromDir(sample.dir)
#'
#' ## translate chip address to probe address
#' ssets <- lapply(dms, ChipAddressToSignal)
#'
#' ## detect p-value
#' ssets <- lapply(ssets, DetectPValue)
#'
#' ## normalization
#' ssets <- lapply(ssets, BackgroundCorrectionNoob)
#' ssets <- DyeBiasCorrectionMostBalanced(ssets)
#' ssets <- Funnorm(ssets)
#'
#' ## convert signal to beta values
#' betas <- sapply(ssets, SignalToBeta)
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
#' @param pval named numeric vector of p-values
#' @return a SignalSet with IG, IR, II, oobG, oobR and ctl
#'   \item{IR, IG and II}{matrices with columns M and U}
#'   \item{oobR and oobG}{arrays of integer}
#'   \item{ctl}{a data frame with columns G, R, col, type}
#' @export
SignalSet <- function(platform,
                      IG=NULL, IR=NULL,
                      II=NULL,
                      oobG=NULL, oobR=NULL,
                      ctl=NULL, pval=NULL) {
  sset <- list(IG=IG, IR=IR,
               II=II,
               oobG=oobG, oobR=oobR,
               ctl=ctl, pval=pval,
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
  dms <- lapply(sample.paths, ReadIDAT1)
  names(dms) <- basename(sample.names)
  dms
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
#'
#' Compute detection p-value. This is typically done before any normalization.
#' @param sset a SignalSet
#' @import stats
#' @return a SignalSet with pval
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

  sset$pval <- c(IR,IG,II)
  sset$pval[order(names(sset$pval))]
  sset
}

#' Convert signal to beta
#' 
#' Convert signal to beta and mask low p-value probes
#' 
#' @param sset a \code{SignalSet}
#' @return a vector of beta values
#' @export
SignalToBeta <- function(sset) {
  betas1 <- pmax(sset$IG[,'M'],1) / pmax(sset$IG[,'M']+sset$IG[,'U'],2)
  betas2 <- pmax(sset$IR[,'M'],1) / pmax(sset$IR[,'M']+sset$IR[,'U'],2)
  betas3 <- pmax(sset$II[,'M'],1) / pmax(sset$II[,'M']+sset$II[,'U'],2)
  betas <- c(betas1, betas2, betas3)
  betas[sset$pval[names(betas)]>0.05] <- NA
  betas[order(names(betas))]
}

#' Convert signal to M-value
#'
#' Convert signal to M-value and mask low p-value probes
#'
#' @param sset a \code{SignalSet}
#' @return a vector of M values
#' @export
SignalToMValue <- function(sset) {
  m1 <- log2(pmax(sset$IG[,'M'],1) / pmax(sset$IG[,'U']))
  m2 <- log2(pmax(sset$IR[,'M'],1) / pmax(sset$IR[,'U']))
  m3 <- log2(pmax(sset$II[,'M'],1) / pmax(sset$II[,'U']))
  m <- c(m1, m2, m3)
  m[sset$pval[names(m)]>0.05] <- NA
  m[order(names(m))]
}

#' Convert beta-value to M-value
#'
#' Convert beta-value to M-value (aka logit transform)
#' @param a vector of beta values
#' @return a vector of M values
#' @export
BetaValueToMValue <- function(b) {
  log2(b/(1-b))
}

#' Convert M-value to beta-value
#'
#' Convert M-value to beta-value (aka inverse logit transform)
#' @param a vector of M values
#' @param a vector of beta values
#' @export
MValueToBetaValue <- function(m) {
  2^m/(1+2^m)
}

#' @export
Message <- function(...) {
  cat('[', as.character(Sys.time()),'] ', ..., '\n', sep='')
}

#' @export
MaskRepeatSNPs <- function(x, platform) {
  UseMethod('MaskRepeatSNPs', x)
}

#' @export
MaskRepeatSNPs.numeric <- function(betas, platform) {
  dm.mask <- GetBuildInData('mask', platform)
  n.na <- sum(is.na(betas))
  betas[names(betas) %in% dm.mask] <- NA
  Message('Masked probes: ', n.na, ' (before) ', sum(is.na(betas)), ' (after).')
  betas
}

#' Mask SNPs and repeats
#'
#' Mask SNPs and repeats
#' @param betas a matrix (multiple samples) or a vector of betas (one sample)
#' @return masked matrix of beta values
#' @export
MaskRepeatSNPs.matrix <- function(betas, platform) {
  dm.mask <- GetBuiltInData('mask', platform)
  n.na <- sum(is.na(betas))
  betas[rownames(betas) %in% dm.mask,] <- NA
  Message('Masked probes: ', n.na, ' (before) ', sum(is.na(betas)), ' (after).')
  betas
}

#' @export
MaskRepeatSNPs.SignalSet <- function(sset, platform) {
  dm.mask <- GetBuiltInData('mask', platform)
  lapply(c('IG','IR','II'), function(nm.cat) 
    sset[[nm.cat]][rownames(sset[[nm.cat]]) %in% dm.mask,] <<- NA)
  sset
}

#' Infer gender from signals
#'
#' Infer gender from signals
#' @param ssets a list of \code{SignalSet}s
#' @return named vector of 0 (for female) and 1 (for male)
#' @export
InferGenders <- function(ssets) {
  probe2chr <- GetBuiltInData('hg19.probe2chr', ssets[[1]]$platform)
  g <- t(vapply(ssets, function(sset) {
    all.signals <- c(apply(sset$IR,1,max), apply(sset$IG,1,max), apply(sset$II,1,max))
    all <- rbind(sset$IG, sset$IR, sset$II)
    all.X <- all[(probe2chr[rownames(all)] == 'chrX'),]
    all.X.betas <- all.X[,'M']/(all.X[,'M']+all.X[,'U'])
    c(median(all.signals[probe2chr[names(all.signals)] == 'chrX']),
      median(all.signals[probe2chr[names(all.signals)] == 'chrY']),
      median(all.X.betas, na.rm=TRUE),
      sum(all.X.betas>0.3 & all.X.betas<0.7, na.rm=TRUE) /
        sum(!is.na(all.X.betas)))
  }, numeric(4)))
  colnames(g) <- c('x.median','y.median','x.beta.median','x.intermed.frac')
  g <- as.data.frame(g)

  ## the simpler classification
  ## g$genders <- ifelse(g$y.median < 500, 0, 1)
  
  ## more accurate but "might" overfit
  g$genders <- ifelse(
    g$y.median < 500,
    ifelse(g$y.median > 300 & g$x.intermed.frac<0.2,1,0),
    ifelse(g$y.median<2000 & g$x.intermed.frac>0.5, 0, 1))
  g
}

#' Filter signal set by beta
#'
#' Filter signal set by beta value
#' @param a SignalSet
#' @return a SignalSet
#' @export
FilterByBeta <- function(sset, max=1.1, min=-0.1) {
  lapply(c('IG','IR','II'), function(nm.cat) {
    b <- sset[[nm.cat]][,'M']/(sset[[nm.cat]][,'M']+sset[[nm.cat]][,'U'])
    sset[[nm.cat]] <<- sset[[nm.cat]][(b>min & b<max),]
  })
  sset
}

#' Measure input
#'
#' log2 mean of the M and U signals in all the probes
#' 
#' This measure roughly quantifies the input DNA quantities to each sample.
#'
#' @param a SignalSet
#' @return input measure
#' @export
MeasureInput <- function(sset) {
  log2(mean(c(sset$IG, sset$IR, sset$II)))
}
