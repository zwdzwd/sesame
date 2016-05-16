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
#' \dontrun{
#' library(sesame)
#' 
#' ## read IDATs
#' dms <- ReadIDATsFromDir(sample.dir)
#'
#' ## translate chip address to probe address
#' ssets <- mclapply(dms, ChipAddressToSignal)
#'
#' ## detect p-value
#' ssets <- mclapply(ssets, DetectPValue)
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
#' }
#' @keywords DNAMethylation Microarray QualityControl
#' 
"_PACKAGE"

#' SignalSet class
#' 
#' @docType class
#' @importFrom R6 R6Class
#' @importFrom stats ecdf
#' @export
#' @return Object of class \code{SignalSet}
#' @format An \code{\link{R6Class}} object.
#' @examples
#' SignalSet$new("EPIC")
#' @field platform platform name, currently supports "EPIC", "hm450" and "hm27"
#' @field IG intensity table for type I probes in green channel
#' @field IR intensity table for type I probes in red channel
#' @field II intensity table for type II probes
#' @field oobG out-of-band probes in green channel
#' @field oobR out-of-band probes in red channel
#' @field ctl all the control probe intensities
#' @field pval named numeric vector of p-values
#' \describe{
#'   \item{Documentation}{For full documentation of each method go to }
#'   \item{\code{new(platform)}}{Create a SignalSet in the specified platform}
#'   \item{\code{detectPValue()}}{Detect P-value for each probe}
#' }
SignalSet <- R6Class(
  'SignalSet',
  portable = FALSE,
  public = list(
    platform = 'EPIC',
    IG = NULL,
    IR = NULL,
    II = NULL,
    oobG = NULL,
    oobR = NULL,
    ctl = NULL,
    pval = NULL,
    
    initialize = function(x) self$platform <- x,
    
    detectPValue = function() {
      negctls <- ctl[grep('negative', tolower(rownames(ctl))),]
      negctls <- subset(negctls, col!=-99)

      funcG <- ecdf(negctls$G)
      funcR <- ecdf(negctls$R)

      ## p-value is the minimium detection p-value of the 2 alleles
      pIR <- 1-apply(cbind(funcR(IR[,'M']), funcR(IR[,'U'])),1,max)
      pIG <- 1-apply(cbind(funcG(IG[,'M']), funcG(IG[,'U'])),1,max)
      pII <- 1-apply(cbind(funcG(II[,'M']), funcR(II[,'U'])),1,max)

      names(pIR) <- rownames(IR)
      names(pIG) <- rownames(IG)
      names(pII) <- rownames(II)

      pval <- c(pIR,pIG,pII)
      pval <<- pval[order(names(pval))]
      invisible()
    },
    
    toBeta = function() {
      betas1 <- pmax(IG[,'M'],1) / pmax(IG[,'M']+IG[,'U'],2)
      betas2 <- pmax(IR[,'M'],1) / pmax(IR[,'M']+IR[,'U'],2)
      betas3 <- pmax(II[,'M'],1) / pmax(II[,'M']+II[,'U'],2)
      betas <- c(betas1, betas2, betas3)
      betas[self$pval[names(betas)]>0.05] <- NA
      betas[order(names(betas))]
    },
    
    toM = function() {
      m1 <- log2(pmax(IG[,'M'],1) / pmax(IG[,'U']))
      m2 <- log2(pmax(IR[,'M'],1) / pmax(IR[,'U']))
      m3 <- log2(pmax(II[,'M'],1) / pmax(II[,'U']))
      m <- c(m1, m2, m3)
      m[pval[names(m)]>0.05] <- NA
      m[order(names(m))]
    },

    maskRepeatSNPs = function() {
      dm.mask <- getBuiltInData('mask', platform)
      lapply(c('IG','IR','II'), function(nm.cat) 
        get(nm.cat)[rownames(sset[[nm.cat]]) %in% dm.mask,] <<- NA)
    },
    
    inferGender = function() {
      probe2chr <- getBuiltInData('hg19.probe2chr', ssets[[1]]$platform)
      all.signals <- c(apply(IR,1,max), apply(IG,1,max), apply(II,1,max))
      all <- rbind(IG, IR, II)
      all.X <- all[(probe2chr[rownames(all)] == 'chrX'),]
      all.X.betas <- all.X[,'M']/(all.X[,'M']+all.X[,'U'])
      x.median <- median(all.signals[probe2chr[names(all.signals)] == 'chrX'])
      y.median <- median(all.signals[probe2chr[names(all.signals)] == 'chrY'])
      x.beta.median <- median(all.X.betas, na.rm=TRUE)
      x.intermed.frac <- sum(all.X.betas>0.3 & all.X.betas<0.7, na.rm=TRUE) /
        sum(!is.na(all.X.betas))
      
      ## more accurate but "might" overfit
      if (y.median < 500) {
        if (y.median > 300 & x.intermed.frac<0.2) 1
        else 0
      } else {
        if (y.median < 2000 & x.intermed.frac>0.5) 0
        else 1
      }
    },

    measureInput = function() {
      log2(mean(c(IG, IR, II)))
    },

    noob = .backgroundCorrectionNoob,
    dyebias = .dyeBiasCorrection
  )
)


#' Import one IDAT file
#'
#' @param idat.name IDAT file name
#' @importFrom illuminaio readIDAT
#' @return a data frame with 2 columns, corresponding to
#' cy3 (Grn) and cy5 (Red) color channel signal
#' @export
readIDAT1 <- function(idat.name) {
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
#' Each element of the returned list contains a matrix
#' having signal intensity addressed by chip address
#' 
#' @param sample.names a sample list
#' @param base.dir base directory
#' @param mc use multiple cores
#' @param mc.cores number of cores to use
#' @import parallel
#' @return a list of data frames. Each data frame corresponds to a sample.
#' Data are signal intensities indexed by chip address.
#' @export
readIDATs <- function(sample.names, base.dir=NULL, mc=FALSE, mc.cores=8) {
  if (!is.null(base.dir))
    sample.paths <- paste0(base.dir,'/',sample.names)
  else
    sample.paths <- sample.names
  if (mc)
    dms <- mclapply(sample.paths, readIDAT1, mc.cores=mc.cores)
  else
    dms <- lapply(sample.paths, readIDAT1)
  names(dms) <- basename(sample.names)
  dms
}

#' Import IDATs from a directory
#' 
#' Each element of the returned list contains a matrix
#' having signal intensity addressed by chip address
#' 
#' @param dir.name directory name.
#' @param ... multiple core parameters: mc and mc.cores see \code{readIDATs}
#' @return a list of data frames. Each data frame corresponds to a sample.
#' Data are signal intensities indexed by chip address.
#' @export
readIDATsFromDir <- function(dir.name, ...) {
  fns <- list.files(dir.name)
  sample.names <- unique(sub("_(Grn|Red).idat", "", fns[grep(".idat$", fns)]))
  readIDATs(paste0(dir.name,'/',sample.names), ...)
}

#' Import IDATs from a sample sheet
#' 
#' Each element of the returned list contains a matrix
#' having signal intensity addressed by chip address
#' 
#' @param sample.sheet path to sample sheet
#' @param base.dir directory on which the \code{sample.sheet.path} is based
#' @param ... multiple core parameters: mc and mc.cores see \code{readIDATs}
#' @return a list of data frames. Each data frame corresponds to a sample. 
#' Data are signal intensities indexed by chip address.
#' @export
readIDATsFromSampleSheet <- function(sample.sheet, base.dir=NULL, ...) {
  sample.names <- read.csv(sample.sheet, stringsAsFactors=F)
  readIDATs(sample.names$barcode, base.dir=base.dir, ...)
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
chipAddressToSignal <- function(dm) {

  platform <- attr(dm, 'platform')
  dm.ordering <- getBuiltInData('ordering', platform)

  sset <- SignalSet$new(platform)

  ## type I green channel / oob red channel
  IordG <- dm.ordering[((dm.ordering$DESIGN=='I')&(dm.ordering$col=='G')),]
  ## 2-channel for green probes' M allele
  IuG2ch <- dm[match(IordG$U, rownames(dm)),]
  IuG <- IuG2ch[,1]
  ## 2-channel for green probes' U allele
  ImG2ch <- dm[match(IordG$M, rownames(dm)),]
  ImG <- ImG2ch[,1]
  sset$oobR <- c(IuG2ch[,2], ImG2ch[,2])
  sset$IG <- as.matrix(data.frame(M=ImG, U=IuG, row.names=IordG$Probe_ID))

  ## type I red channel / oob green channel
  IordR <- dm.ordering[((dm.ordering$DESIGN=='I')&(dm.ordering$col='R')),]
  ## 2-channel for red probes' m allele
  IuR2ch <- dm[match(IordR$U, rownames(dm)),]
  IuR <- IuR2ch[,2]
  ## 2-channel for red probes' u allele
  ImR2ch <- dm[match(IordR$M, rownames(dm)),]
  ImR <- ImR2ch[,2]
  sset$oobG <- c(IuR2ch[,1], ImR2ch[,1])
  sset$IR <- as.matrix(data.frame(M=ImR, U=IuR, row.names=IordR$Probe_ID))

  ## type II
  IIord <- dm.ordering[dm.ordering$DESIGN=="II",]
  signal.II <- dm[match(IIord$U, rownames(dm)),]
  colnames(signal.II) <- c('M', 'U')
  rownames(signal.II) <- IIord$Probe_ID
  sset$II <- signal.II

  ## control probes
  dm.controls <- getBuiltInData('controls', platform)
  ctl <- as.data.frame(dm[match(dm.controls$Address, rownames(dm)),])
  rownames(ctl) <- make.names(dm.controls$Name,unique=T)
  ctl <- cbind(ctl, dm.controls[, c("Color_Channel","Type")])
  colnames(ctl) <- c('G','R','col','type')
  sset$ctl <- ctl

  sset
}

subsetBeta <- function(sset, max=1.1, min=-0.1) {
  sset <- sset$clone()
  lapply(c('IG','IR','II'), function(nm.cat) {
    s <- sset[[nm.cat]]
    b <- s[,'M']/(s[,'M']+s[,'U'])
    sset[[nm.cat]] <<- s[(b>min & b<max),]
    invisible()
  })
  invisible()
}

subsetChromosome <- function(sset, chrm) {
  sset <- sset$clone()
  probe2chr <- getBuiltInData('hg19.probe2chr', sset$platform)
  sset$IG <- sset$IG[probe2chr[rownames(sset$IG)] == chrm,]
  sset$IR <- sset$IR[probe2chr[rownames(sset$IR)] == chrm,]
  sset$II <- sset$II[probe2chr[rownames(sset$II)] == chrm,]
  sset
}
