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
#' sset <- readIDATs(sub('_Grn.idat','',system.file(
#'     'extdata','4207113116_A_Grn.idat',package='sesameData')))[[1]]
#' 
#' ## The TCGA standard pipeline
#' betas <- getBetas(dyeBiasCorrTypeINorm(noob(sset)))
#' 
#' @keywords DNAMethylation Microarray QualityControl
#' 
"_PACKAGE"

#' SigSet class
#'
#' This is the main data class for SeSAMe. The class holds different
#' classes of signal intensities.
#'
#' @slot IG intensity table for type I probes in green channel
#' @slot IR intensity table for type I probes in red channel
#' @slot II intensity table for type II probes
#' @slot oobG out-of-band probes in green channel
#' @slot oobR out-of-band probes in red channel
#' @slot ctl all the control probe intensities
#' @slot pval named numeric vector of p-values
#' @slot platform "EPIC", "HM450" or "HM27"
#' @return a \code{SigSet} object
#'
#' @name SigSet-class
#' @rdname SigSet-class
#' @examples
#' ## Create an empty EPIC object.
#' SigSet("EPIC")
#' @exportClass SigSet
setClass(
    "SigSet",
    representation(
        IG = 'matrix',
        IR = 'matrix',
        II = 'matrix',
        oobG = 'matrix',
        oobR = 'matrix',
        ctl = 'data.frame',
        pval = 'numeric',
        platform = 'character'))

#' Constructor method of SigSet Class.
#'
#' The function takes a string describing the platform of the data. It can be
#' one of "HM27", "HM450" or "EPIC".
#' 
#' @name SigSet
#' @param .Object target object
#' @param platform "EPIC", "HM450" or "HM27"
#' @rdname SigSet-class
#' @return a \code{SigSet} object
#' @aliases initialize,SigSet-method
setMethod(
    "initialize", "SigSet",
    function(.Object, platform = c('EPIC', 'HM450', 'HM27'), ...)  {
        .Object <- callNextMethod()
        .Object@platform <- match.arg(platform)
        .Object
    })

#' Wrapper function for building a new \code{SigSet}
#'
#' The function takes a string describing the platform of the data. It can be
#' one of "HM27", "HM450" or "EPIC".
#'
#' @param ... additional arguments
#' @name SigSet
#' @rdname SigSet-class
#' @examples
#' SigSet('EPIC')
#' @export
SigSet <- function(...) new("SigSet", ...)

#' The display method for SigSet
#'
#' The function outputs the number of probes in each category and the first
#' few signal measurements.
#'
#' @param object displayed object
#' @return message of number of probes in each category.
#' @rdname show-methods
#' @aliases show,SigSet-method
#' @examples
#' sset <- readRDS(system.file(
#'     'extdata','EPIC.sset.LNCaP.Rep1.rds',package='sesameData'))
#' print(sset)
setMethod(
    "show", "SigSet",
    function(object) {
        cat("SigSet", object@platform, "\n - @IG probes:",
            nrow(object@IG), '-', as.numeric(head(object@IG, n=3)),
            "...\n  - @IR probes:", nrow(object@IR),
            '-', as.numeric(head(object@IR, n=3)),
            "...\n  - @II probes:", nrow(object@II),
            '-', as.numeric(head(object@II, n=3)),
            "...\n  - @oobG probes:", nrow(object@oobG),
            '-', as.numeric(head(object@oobG, n=3)),
            "...\n  - @oobR probes:", nrow(object@oobR),
            '-', as.numeric(head(object@oobR, n=3)),
            "...\n  - @ctl probes:", nrow(object@ctl),
            "...\n  - @pval:", length(object@pval),
            "-", as.numeric(head(object@pval, n=3)), "...\n")
    })

#' Select a subset of probes
#'
#' The function takes a \code{SigSet} as input and output another
#' \code{SigSet} with probes from the given probe selection.
#'
#' @param sset a \code{SigSet}
#' @param probes target probes
#' @return another sset with probes specified
#' @examples
#' sset <- readRDS(system.file(
#'     'extdata','EPIC.sset.LNCaP.Rep1.rds',package='sesameData'))
#' subsetSignal(sset, rownames(slot(sset, 'IR')))
#' @export
subsetSignal <- function(sset, probes) {
    stopifnot(is(sset, "SigSet"))
    sset@IR <- sset@IR[rownames(sset@IR) %in% probes,]
    sset@IG <- sset@IG[rownames(sset@IG) %in% probes,]
    sset@II <- sset@II[rownames(sset@II) %in% probes,]
    sset@oobR <- sset@oobR[rownames(sset@oobR) %in% probes,]
    sset@oobG <- sset@oobG[rownames(sset@oobG) %in% probes,]
    sset
}

## get negative control probes
negControls <- function(sset) {
    stopifnot(is(sset, "SigSet"))
    negctls <- sset@ctl[grep('negative', tolower(rownames(sset@ctl))),]
    negctls <- subset(negctls, col!=-99)
    negctls
}

#' Mean Intensity
#'
#' The function takes one single \code{SigSet} and computes mean
#' intensity of all the in-band measurements. This includes all Type-I
#' in-band measurements and all Type-II probe measurements. Both methylated
#' and unmethylated alleles are considered. This function outputs a single
#' numeric for the mean.
#'
#' @param sset a \code{SigSet}
#' @return mean of all intensities
#' @examples
#' sset <- makeExampleSeSAMeDataSet()
#' meanIntensity(sset)
#' @export
meanIntensity <- function(sset) {
    stopifnot(is(sset, "SigSet"))
    mean(c(sset@IG, sset@IR, sset@II), na.rm=TRUE)
}

#' M+U Intensities for All Probes
#'
#' The function takes one single \code{SigSet} and computes total
#' intensity of all the in-band measurements by summing methylated and
#' unmethylated alleles. This function outputs a single numeric for the mean.
#' @param sset a \code{SigSet}
#' @return a vector of M+U signal for each probe
#' @examples
#' sset <- makeExampleSeSAMeDataSet()
#' totalIntensities(sset)
#' @export
totalIntensities <- function(sset) {
    rowSums(rbind(sset@IG, sset@IR, sset@II))
}

#' Get sex-related information
#'
#' The function takes a \code{SigSet} and returns a vector of three
#' numerics: the median intensity of chrY probes; the median intensity of
#' chrX probes; and fraction of intermediate chrX probes. chrX and chrY
#' probes excludes pseudo-autosomal probes.
#'
#' @param sset a \code{SigSet}
#' @return medianY and medianX, fraction of XCI, methylated and unmethylated X
#' probes, median intensities of auto-chromosomes.
#' @examples
#' sset <- makeExampleSeSAMeDataSet()
#' getSexInfo(sset)
#' @export
getSexInfo <- function(sset) {
    stopifnot(is(sset, "SigSet"))
    cleanY <- get(paste0(sset@platform, '.female.clean.chrY.probes'))
    xLinked <- get(paste0(sset@platform, '.female.xlinked.chrX.probes'))
    probe2chr <- get(paste0(sset@platform, '.hg19.probe2chr'))
    xLinkedBeta <- getBetas(subsetSignal(sset, xLinked), quality.mask=FALSE)
    intens <- totalIntensities(sset)
    probes <- intersect(names(intens), names(probe2chr))
    intens <- intens[probes]
    probe2chr <- probe2chr[probes]
    c(
        medianY=median(totalIntensities(subsetSignal(sset, cleanY))),
        medianX=median(totalIntensities(subsetSignal(sset, xLinked))),
        fracXlinked=(sum(
            xLinkedBeta>0.3 & xLinkedBeta<0.7, na.rm = TRUE) /
                sum(!(is.na(xLinkedBeta)))),
        fracXmeth=(
            sum(xLinkedBeta > 0.7, na.rm = TRUE) / sum(!(is.na(xLinkedBeta)))),
        fracXunmeth=(
            sum(xLinkedBeta < 0.3, na.rm = TRUE) / sum(!(is.na(xLinkedBeta)))),
        tapply(intens, probe2chr, median))
}


#' Infer Sex Karyotype
#'
#' The function takes a \code{SigSet} and infers the sex chromosome Karyotype
#' and presence/absence of X-chromosome inactivation (XCI). chrX, chrY and XCI
#' are inferred relatively independently. This function gives a more detailed
#' look of potential sex chromosome aberrations.
#'
#' @param sset a \code{SigSet}
#' @return Karyotype string, with XCI
#' @examples
#' sset <- readRDS(system.file(
#'     'extdata','EPIC.sset.LNCaP.Rep1.rds',package='sesameData'))
#' inferSexKaryotypes(sset)
#' @export
inferSexKaryotypes <- function(sset) {
    stopifnot(is(sset, "SigSet"))
    sex.info <- getSexInfo(sset)
    auto.median <- median(sex.info[paste0('chr',seq_len(22))], na.rm=TRUE)
    XdivAuto <- sex.info['medianX'] / auto.median
    YdivAuto <- sex.info['medianY'] / auto.median
    if (XdivAuto > 1.2) {
        if (sex.info['fracXlinked'] >= 0.5) 
            sexX <- 'XaXi'
        else if (sex.info['fracXmeth'] > sex.info['fracXunmeth'])
            sexX <- 'XiXi'
        else
            sexX <- 'XaXa'
    } else {
        if (sex.info['fracXmeth'] > sex.info['fracXunmeth'])
            sexX <- 'Xi'
        else
            sexX <- 'Xa'
    }

    ## adjust X copy number by fraction of XCI
    if ((sexX == 'Xi' || sexX == 'Xa') && XdivAuto >= 1.0 &&
        sex.info['fracXlinked'] >= 0.5)
        sexX <- 'XaXi'
    
    if (YdivAuto > 0.3 || sex.info['medianY'] > 2000)
        sexY <- 'Y'
    else
        sexY <- ''

    karyotype <- paste0(sexX, sexY)
    karyotype
}

#' Infer Sex
#'
#' @param sset a \code{SigSet}
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
#' sset <- readRDS(system.file(
#'     'extdata','EPIC.sset.LNCaP.Rep1.rds',package='sesameData'))
#' inferSex(sset)
#' @export
inferSex <- function(sset) {
    stopifnot(is(sset, "SigSet"))
    sex.info <- getSexInfo(sset)[seq_len(3)]
    as.character(predict(sesameData::sex.inference.model, sex.info))
}

#' infer ethnicity
#'
#' This function uses both the built-in rsprobes as well as the type I
#' Color-Channel-Switching probes to infer ethnicity.
#'
#' sset better be background subtracted and dyebias corrected for
#' best accuracy
#'
#' @param sset a \code{SigSet}
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
            subsetSignal(sset, rsprobes), quality.mask = FALSE,
            nondetection.mask=FALSE),
        getAFTypeIbySumAlleles(
            subsetSignal(sset, ccsprobes)))
    
    as.character(predict(ethnicity.model, af))
}

#' get beta values
#'
#' @param sset \code{SigSet}
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

    stopifnot(is(sset, "SigSet"))
    
    betas1 <- pmax(sset@IG[,'M'],1) / pmax(sset@IG[,'M']+sset@IG[,'U'],2)
    betas2 <- pmax(sset@IR[,'M'],1) / pmax(sset@IR[,'M']+sset@IR[,'U'],2)
    betas3 <- pmax(sset@II[,'M'],1) / pmax(sset@II[,'M']+sset@II[,'U'],2)
    betas <- c(betas1, betas2, betas3)
    if (nondetection.mask) {
        sset@pval <- sset@pval[match(names(betas), names(sset@pval))]
        betas[sset@pval > pval.threshold] <- NA
    }

    if (quality.mask) {
        if (mask.use.tcga) {
            mask <- get(paste0(sset@platform, '.mask.tcga'))
        } else {
            mask <- get(paste0(sset@platform, '.mask'))
        }
        betas[names(betas) %in% mask] <- NA
    }
    betas[order(names(betas))]
}

#' Get beta values treating type I by summing channels
#'
#' This function is used for rescuing beta values on
#' Color-Channel-Switching CCS probes. The function takes a \code{SigSet}
#' and returns beta value except that Type-I in-band signal and out-of-band
#' signal are combined. This prevents color-channel switching due to SNPs.
#'
#' @param sset \code{SigSet}
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

    stopifnot(is(sset, "SigSet"))
    
    ## .oobR <- oobR[rownames(IG),]
    ## .oobG <- oobG[rownames(IR),]
    betas1 <- pmax(sset@IG[,'M']+sset@oobR[,'M'],1) /
        pmax(sset@IG[,'M']+sset@oobR[,'M']+sset@IG[,'U']+sset@oobR[,'U'],2)
    betas2 <- pmax(sset@IR[,'M']+sset@oobG[,'M'],1) /
        pmax(sset@IR[,'M']+sset@oobG[,'M']+sset@IR[,'U']+sset@oobG[,'U'],2)
    betas3 <- pmax(sset@II[,'M'],1) / pmax(sset@II[,'M']+sset@II[,'U'],2)
    betas <- c(betas1, betas2, betas3)
    if (nondetection.mask) {
        sset@pval <- sset@pval[match(names(betas), names(sset@pval))]
        betas[sset@pval > pval.threshold] <- NA
    }
    if (quality.mask) {
        mask <- get(paste0(sset@platform, '.mask'))
        betas[names(betas) %in% mask] <- NA
    }
    betas[order(names(betas))]
}

#' Get allele frequency treating type I by summing alleles
#'
#' Takes a \code{SigSet} as input and returns a numeric vector containing
#' extra allele frequencies based on Color-Channel-Switching (CCS) probes.
#' If no CCS probes exist in the \code{SigSet}, then an numeric(0) is
#' returned.
#'
#' @param sset \code{SigSet}
#' @return beta values
#' @examples
#' sset <- readRDS(system.file(
#'     'extdata','EPIC.sset.LNCaP.Rep1.rds',package='sesameData'))
#' betas <- getAFTypeIbySumAlleles(sset)
#' @export
getAFTypeIbySumAlleles <- function(sset) {

    stopifnot(is(sset, "SigSet"))
    
    if (any(rownames(sset@oobR) != rownames(sset@IG)))
        stop("oobR-IG not matched. Likely a malformed sset.");
    if (any(rownames(sset@oobG) != rownames(sset@IR)))
        stop("oobG-IR not matched. Likely a malformed sset.");
    
    ## .oobG <- oobG[rownames(IR),]
    af <- c(
        pmax(rowSums(sset@oobR),1)/pmax(rowSums(sset@oobR)+rowSums(sset@IG),2),
        pmax(rowSums(sset@oobG),1)/pmax(rowSums(sset@oobG)+rowSums(sset@IR),2))
    af <- af[
        intersect(names(af), sesameData::ethnicity.ccs.probes)]
    af[order(names(af))]
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
#' If `raw = FALSE` (default), a list of \code{SigSet}s are returned. Each
#' \code{SigSet} contains all the signal measurement on the platform. If
#' `raw = TRUE`, a list is returned with each element of the returned list
#' containing a matrix having signal intensity indexed by chip address.
#'
#' Sample.names is a vector of common prefixes between the *_Grn.idat and
#' *_Red.idat. `num.processes` controls the number of parallel workers. It
#' is default to 1 which means serial.
#' 
#' @param sample.names a sample list
#' @param base.dir base directory
#' @param raw to return raw data without mapping to signal
#' @param num.processes number of parallele processes, serial if 1
#' @return a list of \code{SigSet}s or a list of matrices if `raw=TRUE`
#' @importFrom BiocParallel bplapply
#' @importFrom BiocParallel MulticoreParam
#' @examples
#' ssets <- readIDATs(sub('_Grn.idat','',system.file(
#'     "extdata", "4207113116_A_Grn.idat", package = "sesameData")))
#' @export
readIDATs <- function(
    sample.names, base.dir=NULL, raw=FALSE, num.processes=1) {

    if (!is.null(base.dir))
        sample.paths <- paste0(base.dir,'/',sample.names)
    else
        sample.paths <- sample.names
    
    if (num.processes > 1)
        dms <- bplapply(
            sample.paths, readIDAT1,
            BPPARAM = MulticoreParam(workers = num.processes))
    else
        dms <- lapply(sample.paths, readIDAT1)

    names(dms) <- make.names(basename(sample.names), unique=TRUE)
    if (!raw) {
        if (num.processes > 1) {
            bplapply(
                dms, chipAddressToSignal,
                BPPARAM = MulticoreParam(workers = num.processes))
        } else {
            lapply(dms, chipAddressToSignal)
        }
    } else {
        dms
    }
}

#' Import IDATs from a directory
#' 
#' The input is the directory name as a string. The function identifies all
#' the IDAT files under the directory. The function returns a list of
#' \code{SigSet}s each processed from a pair of IDAT files under the
#' directory.
#' 
#' @param dir.name the directory containing the IDAT files.
#' @param max.num.samples maximum number of samples to process, to protect
#' from memory crash.
#' @param recursive search IDAT files recursively
#' @param ... see \code{readIDATs}
#' @return a list of \code{SigSet}s, if more than `max.num.samples` samples
#' found, then only return the prefixes (a vector of character strings).
#' Consider using `readIDATs` to process each chunk separately.
#' 
#' @examples
#' ## only search what are directly under
#' ssets <- readIDATsFromDir(
#'     system.file("extdata", "", package = "sesameData"))
#' 
#' ## search files recursively
#' ssets <- readIDATsFromDir(
#'     system.file(package = "sesameData"), recursive=TRUE)
#' @export
readIDATsFromDir <- function(
    dir.name, max.num.samples = 100, recursive = FALSE, ...) {

    stopifnot(file.exists(dir.name))
    
    prefixes <- unique(sub(
        '_(Grn|Red).idat$', '',
        list.files(dir.name, '\\.idat$', recursive = recursive)))

    if (length(prefixes) == 0)
        stop("No IDAT file found.")
    if (!all(file.exists(file.path(dir.name, paste0(prefixes, '_Grn.idat')))))
        stop("IDAT names unmatched.")
    if (!all(file.exists(file.path(dir.name, paste0(prefixes, '_Red.idat')))))
        stop("IDAT names unmatched.")

    if (length(prefixes) > max.num.samples) {
        warning(sprintf(
            "%d (>%d) samples found. Returning the prefixes.",
            length(prefixes), max.num.samples))
        file.path(dir.name, prefixes)
    } else {
        readIDATs(file.path(dir.name, prefixes), ...)
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
#' using the other channel.
#'
#' @param dm data frame in chip address, 2 columns: cy3/Grn and cy5/Red
#' @return a SigSet, indexed by probe ID address
chipAddressToSignal <- function(dm) {

    platform <- attr(dm, 'platform')
    dm.ordering <- get(paste0(platform, '.ordering'))

    sset <- SigSet(platform)

    ## type I green channel / oob red channel
    IordG <- dm.ordering[((dm.ordering$DESIGN=='I')&(dm.ordering$col=='G')),]
    ## 2-channel for green probes' M allele
    IuG2ch <- dm[match(IordG$U, rownames(dm)),]
    IuG <- IuG2ch[,1]
    ## 2-channel for green probes' U allele
    ImG2ch <- dm[match(IordG$M, rownames(dm)),]
    ImG <- ImG2ch[,1]
    sset@oobR <- as.matrix(
        data.frame(M=ImG2ch[,2], U=IuG2ch[,2], row.names=IordG$Probe_ID))
    sset@IG <- as.matrix(data.frame(M=ImG, U=IuG, row.names=IordG$Probe_ID))

    ## type I red channel / oob green channel
    IordR <- dm.ordering[((dm.ordering$DESIGN=='I')&(dm.ordering$col=='R')),]
    ## 2-channel for red probes' m allele
    IuR2ch <- dm[match(IordR$U, rownames(dm)),]
    IuR <- IuR2ch[,2]
    ## 2-channel for red probes' u allele
    ImR2ch <- dm[match(IordR$M, rownames(dm)),]
    ImR <- ImR2ch[,2]
    sset@oobG <- as.matrix(
        data.frame(M=ImR2ch[,1], U=IuR2ch[,1], row.names=IordR$Probe_ID))
    sset@IR <- as.matrix(data.frame(M=ImR, U=IuR, row.names=IordR$Probe_ID))

    ## type II
    IIord <- dm.ordering[dm.ordering$DESIGN=="II",]
    signal.II <- dm[match(IIord$U, rownames(dm)),]
    colnames(signal.II) <- c('M', 'U')
    rownames(signal.II) <- IIord$Probe_ID
    sset@II <- signal.II

    ## control probes
    dm.controls <- get(paste0(platform, '.controls'))
    ctl <- as.data.frame(dm[match(dm.controls$Address, rownames(dm)),])
    rownames(ctl) <- make.names(dm.controls$Name, unique=TRUE)
    ctl <- cbind(ctl, dm.controls[, c("Color_Channel","Type")])
    colnames(ctl) <- c('G','R','col','type')
    sset@ctl <- ctl

    sset <- detectionPoobEcdf(sset)
    sset
}

#' Compute internal bisulfite conversion control
#'
#' Compute GCT score for internal bisulfite conversion control. The function
#' takes a \code{SigSet} as input. The higher the GCT score, the more likely
#' the incomplete conversion. The lower the GCT score, the more likely
#' over-conversion.
#' 
#' @param sset signal set
#' @param use.median use median to compute GCT instead of mean
#' @return GCT score (the higher, the more incomplete conversion)
#' @examples
#' sset <- makeExampleSeSAMeDataSet('HM450')
#' bisConversionControl(sset)
#' 
#' @export
bisConversionControl <- function(sset, use.median=FALSE) {
    extC <- get(paste0(sset@platform, '.typeI.extC'))
    extT <- get(paste0(sset@platform, '.typeI.extT'))
    if (use.median) {
        median(sset@oobG[extC,], na.rm=TRUE) /
            median(sset@oobG[extT,], na.rm=TRUE)
    } else {
        mean(sset@oobG[extC,], na.rm=TRUE) /
            mean(sset@oobG[extT,], na.rm=TRUE)
    }
}

#### Following are obsolete and exist for backward-compatibility #####

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

#' sesame R6 to S4
#'
#' SignalSet-SigSet conversion
#'
#' @param sset signalset in R6
#' @return signalset in S4
#' @import methods
#' @examples
#' 
#' sset <- SignalSet$new('EPIC')
#' sset$IG <- matrix(1:4, nrow=2)
#' sset$IR <- matrix(1:4, nrow=2)
#' sset$II <- matrix(1:4, nrow=2)
#' sset$oobG <- matrix(1:4, nrow=2)
#' sset$oobR <- matrix(1:4, nrow=2)
#' sset$ctl <- data.frame(G=1:2,R=3:4)
#' sset$pval <- rep(0,2)
#'
#' signalR6toS4(sset)
#' @export
signalR6toS4 <- function(sset) {
    new('SigSet', IG=sset$IG, IR=sset$IR, II=sset$II,
        oobG=sset$oobG, oobR=sset$oobR, ctl=sset$ctl,
        pval=sset$pval, platform=toupper(sset$platform))
}
