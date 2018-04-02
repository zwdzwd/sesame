#' @title
#' Analyze DNA methylation data
#' 
#' @description
#' SEnsible and step-wise analysis of DNA methylation data
#' 
#' @details
#' This package complements array functionalities that allow
#' processing >10,000 samples in parallel on clusters.
#' @aliases sesame
#' @author
#' Wanding Zhou \email{Wanding.Zhou@vai.org},
#' Hui Shen \email{Hui.Shen@vai.org}
#' Timothy J Triche Jr \email{Tim.Triche@vai.org}
#' 
## @references To appear
## @seealso To appear
#' @examples
#' 
#' library(sesame)
#' library(BiocParallel)
#' 
#' ## read IDATs
#' ssets <- readIDATsFromDir(system.file('extdata', '', package='sesameData'))
#'
#' ## normalization
#' ssets <- bplapply(ssets, noob)
#' ssets <- bplapply(ssets, dyeBiasCorrTypeINorm)
#' 
#' ## convert signal to beta values
#' betas <- do.call(cbind, bplapply(ssets, getBetas))
#' 
#' @keywords DNAMethylation Microarray QualityControl
#' 
"_PACKAGE"

#' SignalSet class
#' 
#' @docType class
#' @importFrom R6 R6Class
#' @importFrom stats ecdf
#' @importFrom utils data
#' @importFrom utils download.file
#' @export
#' @return Object of class \code{SignalSet}
#' @format An \code{\link{R6Class}} object.
#' @examples
#' SignalSet$new("EPIC")
#' @field platform platform name, supports "EPIC", "HM450" and "HM27"
#' @field IG intensity table for type I probes in green channel
#' @field IR intensity table for type I probes in red channel
#' @field II intensity table for type II probes
#' @field oobG out-of-band probes in green channel
#' @field oobR out-of-band probes in red channel
#' @field ctl all the control probe intensities
#' @field pval named numeric vector of p-values
#' @field mask probe mask
#' \describe{
#'   \item{Documentation}{For full documentation of each method go to }
#'   \item{\code{new(platform)}}{Create a SignalSet in the specified platform}
#'   \item{\code{detectPValue()}}{Detect P-value for each probe}
#'   \item{\code{toM()}}{Convert to M values}
#'   \item{\code{totalIntensities()}}{Total intensity on each probe}
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
    
        toM = function() {
            m1 <- log2(pmax(IG[,'M'],1) / pmax(IG[,'U']))
            m2 <- log2(pmax(IR[,'M'],1) / pmax(IR[,'U']))
            m3 <- log2(pmax(II[,'M'],1) / pmax(II[,'U']))
            m <- c(m1, m2, m3)
            m[pval[names(m)]>0.05] <- NA
            m[order(names(m))]
        },
    
        totalIntensities = function() {
            rowSums(rbind(IG, IR, II))
        }
    )
)

## get negative control probes
negControls <- function(sset) {
    
    negctls <- sset$ctl[grep('negative', tolower(rownames(sset$ctl))),]
    negctls <- subset(negctls, col!=-99)
    negctls
}

#' mean intensity
#'
#' mean intensity
#'
#' @param sset a \code{SignalSet}
#' @return mean of all intensities
#' @examples
#' sset <- makeExampleSeSAMeDataSet()
#' intensities <- meanIntensity(sset)
#' @export
meanIntensity <- function(sset) {
    mean(c(sset$IG, sset$IR, sset$II), na.rm=TRUE)
}

#' get sex information
#'
#' @param sset a \code{SignalSet}
#' @return medianY and fracXlinked
#' @examples
#' sset <- makeExampleSeSAMeDataSet()
#' getSexInfo(sset)
#' @export
getSexInfo <- function(sset) {
    cleanY <- get(paste0(sset$platform, '.female.clean.chrY.probes'))
    xLinked <- get(paste0(sset$platform, '.female.xlinked.chrX.probes'))
    xLinkedBeta <- getBetas(sset[xLinked], quality.mask=FALSE)
    c(
        medianY=median(sset[cleanY]$totalIntensities()),
        medianX=median(sset[xLinked]$totalIntensities()),
        fracXlinked=(sum(
            xLinkedBeta>0.3 & xLinkedBeta<0.7, na.rm = TRUE) /
                sum(!(is.na(xLinkedBeta)))))
}

#' infer sex
#'
#' @param sset a \code{SignalSet}
#' @return 'F' or 'M'
#' We established our sex calling based on the median intensity of
#' chromosome X, Y and the fraction of intermediately methylated probes
#' among the identified X-linked probes. This is similar to the
#' approach by Minfi (Aryee et al., 2014) but also different in that
#' we used the fraction of intermediate beta value rather than
#' median intensity for all chromosome X probes. Instead of using
#' all probes from the sex chromosomes, we used our curated set of Y
#' chromosome probes and X-linked probes which exclude potential
#' cross-hybridization and pseudo-autosomal effect.
#'
#' XXY male (Klinefelter's), 45,X female (Turner's) can confuse the
#' model sometimes.
#' Our function works on a single sample.
#' @importFrom randomForest randomForest
#' @import sesameData
#' @examples
#' sset <- makeExampleSeSAMeDataSet()
#' sex <- inferSex(sset)
#' @export
inferSex <- function(sset) {
    sex.info <- getSexInfo(sset)
    as.character(predict(sesameData::sex.inference.model, sex.info))
}

#' infer ethnicity
#'
#' this uses both the built-in rsprobes as well as the type I
#' Color-Channel-Switching probes
#'
#' sset better be background subtracted and dyebias corrected for
#' best accuracy
#'
#' @param sset a \code{SignalSet}
#' @return string of ethnicity
#' @importFrom randomForest randomForest
#' @import sesameData
#' @examples
#' sset <- makeExampleSeSAMeDataSet("HM450")
#' inferEthnicity(sset)
#' @export
inferEthnicity <- function(sset) {
    ccsprobes <- sesameData::ethnicity.ccs.probes
    rsprobes <- sesameData::ethnicity.rs.probes
    ethnicity.model <- sesameData::ethnicity.model
    af <- c(
        getBetas(
            sset[rsprobes], quality.mask = FALSE,
            nondetection.mask=FALSE),
        getAFTypeIbySumAlleles(
            sset[ccsprobes],
            quality.mask = FALSE,
            nondetection.mask = FALSE))
    
    as.character(predict(ethnicity.model, af))
}

#' subset a SignalSet
#'
#' @param sset a \code{SignalSet}
#' @param probes probe names
#' @return a \code{SignalSet} with only probes
#' @export
`[.SignalSet` <- function(sset, probes) {
    sset <- sset$clone()
    sset$IR <- sset$IR[rownames(sset$IR) %in% probes,]
    sset$IG <- sset$IG[rownames(sset$IG) %in% probes,]
    sset$II <- sset$II[rownames(sset$II) %in% probes,]
    sset$oobR <- sset$oobR[rownames(sset$oobR) %in% probes,]
    sset$oobG <- sset$oobG[rownames(sset$oobG) %in% probes,]
    sset
}

#' get beta values
#'
#' @param sset \code{SignalSet}
#' @param quality.mask whether to mask low quality probes
#' @param nondetection.mask whether to mask nondetection
#' @param mask.use.tcga whether to use TCGA masking, only applies to HM450
#' @param pval.threshold p-value threshold for nondetection mask
#' @return a numeric vector, beta values
#' @examples
#' sset <- makeExampleSeSAMeDataSet('HM450')
#' betas <- getBetas(sset)
#' @export
getBetas <- function(
    sset, quality.mask=TRUE, nondetection.mask=TRUE,
    mask.use.tcga=FALSE, pval.threshold=0.05) {
    
    betas1 <- pmax(sset$IG[,'M'],1) / pmax(sset$IG[,'M']+sset$IG[,'U'],2)
    betas2 <- pmax(sset$IR[,'M'],1) / pmax(sset$IR[,'M']+sset$IR[,'U'],2)
    betas3 <- pmax(sset$II[,'M'],1) / pmax(sset$II[,'M']+sset$II[,'U'],2)
    betas <- c(betas1, betas2, betas3)
    if (nondetection.mask)
        betas[sset$pval[names(betas)] > pval.threshold] <- NA
    if (quality.mask) {
        if (mask.use.tcga) {
            mask <- get(paste0(sset$platform, '.mask.tcga'))
        } else {
            mask <- get(paste0(sset$platform, '.mask'))
        }
        betas[names(betas) %in% mask] <- NA
    }
    betas[order(names(betas))]
}

#' get beta values treating type I by summing channels
#'
#' used for rescuing beta values on Color-Channel-Switching CCS probes
#'
#' @param sset \code{SignalSet}
#' @param quality.mask whether to mask low quality probes
#' @param nondetection.mask whether to mask nondetection
#' @param pval.threshold p-value threshold for nondetection mask
#' @return beta values
#' @examples
#' sset <- makeExampleSeSAMeDataSet()
#' betas <- getBetasTypeIbySumChannels(sset)
#' @export
getBetasTypeIbySumChannels <- function(
    sset, quality.mask=TRUE, nondetection.mask=TRUE, pval.threshold=0.05) {
    
    ## .oobR <- oobR[rownames(IG),]
    ## .oobG <- oobG[rownames(IR),]
    betas1 <- pmax(sset$IG[,'M']+sset$oobR[,'M'],1) /
        pmax(sset$IG[,'M']+sset$oobR[,'M']+sset$IG[,'U']+sset$oobR[,'U'],2)
    betas2 <- pmax(sset$IR[,'M']+sset$oobG[,'M'],1) /
        pmax(sset$IR[,'M']+sset$oobG[,'M']+sset$IR[,'U']+sset$oobG[,'U'],2)
    betas3 <- pmax(sset$II[,'M'],1) / pmax(sset$II[,'M']+sset$II[,'U'],2)
    betas <- c(betas1, betas2, betas3)
    if (nondetection.mask)
        betas[sset$pval[names(betas)]>pval.threshold] <- NA
    if (quality.mask) {
        mask <- get(paste0(sset$platform, '.mask'))
        betas[names(betas) %in% mask] <- NA
    }
    betas[order(names(betas))]
}

#' get allele frequency treating type I by summing alleles
#'
#' used for extra allele frequency on Color-Channel-Switching CCS probes
#'
#' @param sset \code{SignalSet}
#' @param quality.mask whether to mask low quality probes
#' @param nondetection.mask whether to mask nondetection
#' @param pval.threshold p-value threshold for nondetection mask
#' @return beta values
#' @examples
#' sset <- readRDS(system.file(
#'     'extdata','EPIC.sset.LNCaP.Rep1.rds',package='sesameData'))
#' betas <- getAFTypeIbySumAlleles(sset)
#' @export
getAFTypeIbySumAlleles <- function(
    sset, quality.mask=TRUE, nondetection.mask=TRUE, pval.threshold=0.05) {

    ## .oobR <- oobR[rownames(IG),]
    ## .oobG <- oobG[rownames(IR),]
    betas1 <- pmax(rowSums(sset$oobR),1) /
        pmax(rowSums(sset$oobR) + rowSums(sset$IG), 2)
    betas2 <- pmax(rowSums(sset$oobG),1) /
        pmax(rowSums(sset$oobG) + rowSums(sset$IR), 2)
    betas3 <- pmax(sset$II[,'M'],1) / pmax(sset$II[,'M']+sset$II[,'U'],2)
    betas <- c(betas1, betas2, betas3)
    if (nondetection.mask)
        betas[sset$pval[names(betas)]>pval.threshold] <- NA
    if (quality.mask) {
        mask <- get(paste0(sset$platform, '.mask'))
        betas[names(betas) %in% mask] <- NA
    }
    betas <- betas[
        intersect(names(betas), sesameData::ethnicity.ccs.probes)]
    betas[order(names(betas))]
}

## Import one IDAT file
## return a data frame with 2 columns, corresponding to
## cy3 (Grn) and cy5 (Red) color channel signal
readIDAT1 <- function(idat.name) {
    ida.grn <- illuminaio::readIDAT(paste0(idat.name,"_Grn.idat"));
    ida.red <- illuminaio::readIDAT(paste0(idat.name,"_Red.idat"));
    d <- cbind(cy3=ida.grn$Quants[,"Mean"], cy5=ida.red$Quants[,"Mean"])
    colnames(d) <- c('G', 'R')
    chip.type <- switch(
        ida.red$ChipType,
        'BeadChip 8x5'='EPIC',
        'BeadChip 12x8'='HM450',
        'BeadChip 12x1'='HM27')
    attr(d, 'platform') <- chip.type
    d
}

#' Import IDATs from a list of samples
#'
#' Each element of the returned list contains a matrix
#' having signal intensity addressed by chip address
#'
#' sample.names is a vector of common prefixes between
#' the _Grn.idat and _Red.idat.
#' 
#' @param sample.names a sample list
#' @param base.dir base directory
#' @param raw to return raw data without mapping to signal
#' @param mc use multiple cores
#' @param mc.cores number of cores to use
#' @return a list of \code{SignalSet}s or a list of matrices if `raw=TRUE`
#' @importFrom BiocParallel bplapply
#' @importFrom BiocParallel MulticoreParam
#' @examples
#' ssets <- readIDATs(sub('_Grn.idat','',system.file(
#'  "extdata", "4207113116_A_Grn.idat", package = "sesameData")))
#' @export
readIDATs <- function(
    sample.names, base.dir=NULL, raw=FALSE, mc=FALSE, mc.cores=NULL) {

    if (is.null(mc.cores)) {
        if (is.null(getOption('mc.cores')))
            mc.cores <- 4
        else
            mc.cores <- getOption('mc.cores')
    }
    
    if (!is.null(base.dir))
        sample.paths <- paste0(base.dir,'/',sample.names)
    else
        sample.paths <- sample.names
    
    if (mc || !is.null(mc.cores))
        dms <- bplapply(
            sample.paths, readIDAT1,
            BPPARAM = MulticoreParam(workers=mc.cores))
    else
        dms <- lapply(sample.paths, readIDAT1)

    names(dms) <- basename(sample.names)
    if (!raw) {
        if (mc) {
            bplapply(dms, chipAddressToSignal, mc.cores=mc.cores)
        } else {
            lapply(dms, chipAddressToSignal)
        }
    } else {
        dms
    }
}

#' Import IDATs from a directory
#' 
#' Each element of the returned list contains a matrix
#' having signal intensity addressed by chip address
#'
#' dir.name is the folder containing the IDAT files.
#' 
#' @param dir.name directory name.
#' @param ... multiple core parameters: mc and mc.cores see \code{readIDATs}
#' @return a list of \code{SignalSet}s
#' @examples
#' ssets <- readIDATsFromDir(
#'   system.file("extdata", "", package = "sesameData"))
#' @export
readIDATsFromDir <- function(dir.name, ...) {
    fns <- list.files(dir.name)
    sample.names <- unique(sub(
        "_(Grn|Red).idat", "", fns[grep(".idat$", fns)]))
    readIDATs(paste0(dir.name,'/',sample.names), ...)
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
chipAddressToSignal <- function(dm) {

    platform <- attr(dm, 'platform')
    dm.ordering <- get(paste0(platform, '.ordering'))

    sset <- SignalSet$new(platform)

    ## type I green channel / oob red channel
    IordG <- dm.ordering[((dm.ordering$DESIGN=='I')&(dm.ordering$col=='G')),]
    ## 2-channel for green probes' M allele
    IuG2ch <- dm[match(IordG$U, rownames(dm)),]
    IuG <- IuG2ch[,1]
    ## 2-channel for green probes' U allele
    ImG2ch <- dm[match(IordG$M, rownames(dm)),]
    ImG <- ImG2ch[,1]
    sset$oobR <- as.matrix(
        data.frame(M=ImG2ch[,2], U=IuG2ch[,2], row.names=IordG$Probe_ID))
    sset$IG <- as.matrix(data.frame(M=ImG, U=IuG, row.names=IordG$Probe_ID))

    ## type I red channel / oob green channel
    IordR <- dm.ordering[((dm.ordering$DESIGN=='I')&(dm.ordering$col=='R')),]
    ## 2-channel for red probes' m allele
    IuR2ch <- dm[match(IordR$U, rownames(dm)),]
    IuR <- IuR2ch[,2]
    ## 2-channel for red probes' u allele
    ImR2ch <- dm[match(IordR$M, rownames(dm)),]
    ImR <- ImR2ch[,2]
    sset$oobG <- as.matrix(
        data.frame(M=ImR2ch[,1], U=IuR2ch[,1], row.names=IordR$Probe_ID))
    sset$IR <- as.matrix(data.frame(M=ImR, U=IuR, row.names=IordR$Probe_ID))

    ## type II
    IIord <- dm.ordering[dm.ordering$DESIGN=="II",]
    signal.II <- dm[match(IIord$U, rownames(dm)),]
    colnames(signal.II) <- c('M', 'U')
    rownames(signal.II) <- IIord$Probe_ID
    sset$II <- signal.II

    ## control probes
    dm.controls <- get(paste0(platform, '.controls'))
    ctl <- as.data.frame(dm[match(dm.controls$Address, rownames(dm)),])
    rownames(ctl) <- make.names(dm.controls$Name,unique=TRUE)
    ctl <- cbind(ctl, dm.controls[, c("Color_Channel","Type")])
    colnames(ctl) <- c('G','R','col','type')
    sset$ctl <- ctl

    sset$pval <- detectionPoobEcdf(sset)
    sset
}

#' compute internal bisulfite conversion control
#'
#' compute GCT score for internal bisulfite conversion control
#'
#' The higher the GCT score, the more likely the incomplete conversion.
#' The lower the GCT score, the more likely over-conversion.
#' 
#' @param sset signal set
#' @param use.median use median to compute GCT
#' @return GCT score (the higher, the more incomplete conversion)
#' @examples
#' sset <- makeExampleSeSAMeDataSet('HM450')
#' bisConversionControl(sset)
#' 
#' @export
bisConversionControl <- function(sset, use.median=FALSE) {
    extC <- get(paste0(sset$platform, '.typeI.extC'))
    extT <- get(paste0(sset$platform, '.typeI.extT'))
    if (use.median) {
        median(sset$oobG[extC,], na.rm=TRUE) /
            median(sset$oobG[extT,], na.rm=TRUE)
    } else {
        mean(sset$oobG[extC,], na.rm=TRUE) /
            mean(sset$oobG[extT,], na.rm=TRUE)
    }
}
