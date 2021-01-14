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
#' sset <- readIDATpair(sub('_Grn.idat','',system.file(
#'     'extdata','4207113116_A_Grn.idat',package='sesameData')))
#'
#' ## The OpenSesame pipeline
#' betas <- openSesame(sset)
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
#' @slot extra extra data, currently,
#' IGG => Type-I green that is inferred to be green
#' IRR => Type-I red that is inferred to be red
#' pvals => list of other pvals
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
        extra = 'list', # extra data, extended to allow additional data
        platform = 'character'))

## update old SigSet without slot extra
updateSigSet <- function(sset) {
    sset2 <- SigSet(sset@platform)
    for(sname in c("IG", "IR", "II", "oobG", "oobR", "ctl", "pval", "extra", "platform")) {
        if (sname == 'extra')
            slot(sset2, sname) <- list()
        else
            slot(sset2, sname) <- slot(sset, sname)
    }
    sset2
}

#' Constructor method of SigSet class.
#'
#' The function takes a string describing the platform of the data. It can be
#' one of "HM27", "HM450" or "EPIC".
#'
#' @name SigSet
#' @param .Object target object
#' @param platform "EPIC", "HM450", "HM27" or other strings for custom arrays
#' @rdname SigSet-class
#' @return a \code{SigSet} object
#' @aliases initialize,SigSet-method
setMethod(
    "initialize", "SigSet",
    function(.Object, platform, ...)  {
        .Object <- callNextMethod()
        .Object@extra <- list()
        .Object@platform <- platform # match.arg(platform)
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
#' @import methods
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
#' sset <- sesameDataGet('EPIC.1.LNCaP')$sset
#' print(sset)
setMethod(
    "show", "SigSet",
    function(object) {
        cat("SigSet", object@platform, "\n - @IG probes:",
            nrow(object@IG), '-', as.numeric(head(object@IG, n=3)),
            "...\n - @IR probes:", nrow(object@IR),
            '-', as.numeric(head(object@IR, n=3)),
            "...\n - @II probes:", nrow(object@II),
            '-', as.numeric(head(object@II, n=3)),
            "...\n - @oobG probes:", nrow(object@oobG),
            '-', as.numeric(head(object@oobG, n=3)),
            "...\n - @oobR probes:", nrow(object@oobR),
            '-', as.numeric(head(object@oobR, n=3)),
            "...\n - @ctl probes:", nrow(object@ctl),
            "...\n - @pval:", length(object@pval),
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
#' sset <- sesameDataGet('EPIC.1.LNCaP')$sset
#' subsetSignal(sset, rownames(slot(sset, 'IR')))
#' @export
subsetSignal <- function(sset, probes) {
    stopifnot(is(sset, "SigSet"))
    IR(sset) <- IR(sset)[rownames(IR(sset)) %in% probes,,drop=FALSE]
    IG(sset) <- IG(sset)[rownames(IG(sset)) %in% probes,,drop=FALSE]
    II(sset) <- II(sset)[rownames(II(sset)) %in% probes,,drop=FALSE]
    oobR(sset) <- oobR(sset)[rownames(oobR(sset)) %in% probes,,drop=FALSE]
    oobG(sset) <- oobG(sset)[rownames(oobG(sset)) %in% probes,,drop=FALSE]
    sset
}

## get negative control probes
negControls <- function(sset) {
    stopifnot(is(sset, "SigSet"))
    negctls <- ctl(sset)[grep('negative', tolower(rownames(ctl(sset)))),]
    negctls <- subset(negctls, col!=-99)
    negctls
}

#' report M and U for regular probes
#'
#' @param sset a \code{SigSet}
#' @return a data frame of M and U columns
#' @examples
#' sset <- sesameDataGet('EPIC.1.LNCaP')$sset
#' signalMU(sset)
#' @export
signalMU <- function(sset) {
    rbind(IR(sset), IG(sset), II(sset))
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
#' @param mask.use.manifest use mask column in the manifest to filter probes
#' attributes set in extra(sset)
#' @return mean of all intensities
#' @examples
#' sset <- makeExampleSeSAMeDataSet()
#' meanIntensity(sset)
#' @export
meanIntensity <- function(sset, mask.use.manifest = TRUE) {
    stopifnot(is(sset, "SigSet"))
    IGset <- IG(sset)
    IRset <- IR(sset)
    IIset <- II(sset)
    if (mask.use.manifest && extraHas(sset, 'mask')) {
        pmask <- extraGet(sset, 'mask')
        IGset <- IGset[!pmask[rownames(IGset)],]
        IRset <- IRset[!pmask[rownames(IRset)],]
        IIset <- IIset[!pmask[rownames(IIset)],]
    }
    
    mean(c(IG(sset), IR(sset), II(sset)), na.rm=TRUE)
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
#' intensities <- totalIntensities(sset)
#' @export
totalIntensities <- function(sset) {
    rowSums(rbind(IG(sset), IR(sset), II(sset)))
}

#' Calculate intensity Z-score
#'
#' This function compute intensity Z-score with respect to the mean.
#' Log10 transformation is done first. Probes of each design type are
#' grouped before Z-scores are computed.
#'
#' @param sset a \code{SigSet}
#' @return a vector of Z-score for each probe
#' @examples
#' sset <- sesameDataGet('HM450.1.TCGA.PAAD')$sset
#' head(totalIntensityZscore(sset))
#' @export
totalIntensityZscore <- function(sset) {
    Zscore <- rbind(
        scale(log10(1+rowSums(IG(sset)))),
        scale(log10(1+rowSums(IR(sset)))),
        scale(log10(1+rowSums(II(sset)))))[,1]
    Zscore[sort(names(Zscore))]
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

    if (is(sset, "SigSetList"))
        return(do.call(cbind, lapply(sset, getSexInfo)))
    
    stopifnot(is(sset, "SigSet"))

    cleanY <- sesameDataGet(paste0(
        sset@platform,'.probeInfo'))$chrY.clean

    xLinked <- sesameDataGet(paste0(
        sset@platform,'.probeInfo'))$chrX.xlinked

    probe2chr <- sesameDataGet(paste0(
        sset@platform,'.probeInfo'))$probe2chr.hg19

    xLinkedBeta <- getBetas(subsetSignal(sset, xLinked), mask=FALSE)
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
#' sset <- sesameDataGet('EPIC.1.LNCaP')$sset
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
#' sset <- sesameDataGet('EPIC.1.LNCaP')$sset
#' inferSex(sset)
#' @export
inferSex <- function(sset) {
    stopifnot(is(sset, "SigSet"))
    stopifnot(sset@platform %in% c('EPIC','HM450'))
    sex.info <- getSexInfo(sset)[seq_len(3)]
    as.character(predict(
        sesameDataGet('sex.inference'), sex.info))
}

#' Infer Ethnicity
#'
#' This function uses both the built-in rsprobes as well as the type I
#' Color-Channel-Switching probes to infer ethnicity.
#'
#' sset better be background subtracted and dyebias corrected for
#' best accuracy
#'
#' Please note: the betas should come from sigset *without*
#' channel inference.
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

    if (is(sset, "SigSetList"))
        return(vapply(sset, inferEthnicity, character(1)))
    
    stopifnot(is(sset, 'SigSet'))
    stopifnot(sset@platform %in% c('EPIC','HM450'))

    ethnicity.inference <- sesameDataGet('ethnicity.inference')
    ccsprobes <- ethnicity.inference$ccs.probes
    rsprobes <- ethnicity.inference$rs.probes
    ethnicity.model <- ethnicity.inference$model
    af <- c(
        getBetas(sset, mask=FALSE)[rsprobes],
        getAFTypeIbySumAlleles(
            subsetSignal(sset, ccsprobes)))

    as.character(predict(ethnicity.model, af))
}

#' Mask beta values by design quality
#' 
#' Currently quality masking only supports three platforms
#' 
#' @param sset a \code{SigSet} object
#' @param mask.use.tcga whether to use TCGA masking, only applies to HM450
#' @return a filtered \code{SigSet}
#' @examples
#' sset <- sesameDataGet('EPIC.1.LNCaP')$sset
#' sset.masked <- qualityMask(sset)
#' @export 
qualityMask <- function(
    sset,
    mask.use.tcga = FALSE) {

    if(!(sset@platform %in% c('HM27','HM450','EPIC'))) {
        message(sprintf(
            "Quality masking is not supported for %s.", sset@platform))
        return(sset)
    }
        
    
    if (mask.use.tcga) {
        stopifnot(sset@platform == 'HM450')
        masked <- sesameDataGet('HM450.probeInfo')$mask.tcga
    } else {
        masked <- sesameDataGet(paste0(sset@platform, '.probeInfo'))$mask
    }

    if (!extraHas(sset, 'mask')) {
        resetMask(sset);
    }
    sset@extra$mask[masked] <- TRUE

    return(sset)
}

#' Reset Masking
#'
#' @param sset a \code{SigSet}
#' @return a new \code{SigSet} with mask reset to empty
#' @examples
#' sset <- sesameDataGet('EPIC.1.LNCaP')$sset
#' sset.no.mask <- resetMask(sset)
#' @export
resetMask <- function(sset) {
    probes <- probeIDs(sset)
    sset@extra$mask <- setNames(rep(FALSE, length(probes)), probes)
}

#' Mask Sigset by detection p-value
#'
#' @param sset a \code{SigSet}
#' @param pval.method which method to use in calculating p-values
#' @param pval.threshold the p-value threshold
#' @return a filtered \code{SigSet}
#' @examples
#' sset <- sesameDataGet('EPIC.1.LNCaP')$sset
#' sset.masked <- detectionMask(sset)
#' @export
detectionMask <- function(
    sset, pval.method=NULL, pval.threshold=0.05) {
    if (is.null(pval.method)) {
        pv <- pval(sset)
    } else {
        stopifnot(
            extraHas(sset, 'pvals') &&
                pval.method %in% names(sset@extra$pvals))
        pv <- sset@extra$pvals[[pval.method]]
    }

    if (!extraHas(sset, 'mask')) {
        resetMask(sset);
    }
    sset@extra$mask[pv[names(sset@extra$mask)] > pval.threshold] <- TRUE

    sset
}

#' Get beta Values
#'
#' sum.typeI is used for rescuing beta values on
#' Color-Channel-Switching CCS probes. The function takes a \code{SigSet}
#' and returns beta value except that Type-I in-band signal and out-of-band
#' signal are combined. This prevents color-channel switching due to SNPs.
#' 
#' @param sset \code{SigSet}
#' @param sum.TypeI whether to sum type I channels
#' @param mask whether to use mask
#' @return a numeric vector, beta values
#' @examples
#' sset <- sesameDataGet('EPIC.1.LNCaP')$sset
#' betas <- getBetas(sset)
#' @export
getBetas <- function(sset, mask=TRUE, sum.TypeI = FALSE) {

    if (is(sset, "SigSetList")) {
        return(do.call(cbind, lapply(
            sset, getBetas, mask=mask, sum.TypeI=sum.TypeI)))
    }
    
    stopifnot(is(sset, "SigSet"))

    IGs <- IG(sset)
    IRs <- IR(sset)

    ## optionally summing channels protects
    ## against channel misspecification
    if (sum.TypeI) {
        IGs <- IGs + oobR2(sset)
        IRs <- IRs + oobG2(sset)
    } else if (!is.null(sset@extra$IGG) && !is.null(sset@extra$IRR)) {
        IGs[!sset@extra$IGG,] <- sset@oobR[!sset@extra$IGG,]
        IRs[!sset@extra$IRR,] <- sset@oobG[!sset@extra$IRR,]
    }

    betas <- c(
        pmax(IGs[,'M'],1) / pmax(IGs[,'M']+IGs[,'U'],2),
        pmax(IRs[,'M'],1) / pmax(IRs[,'M']+IRs[,'U'],2),
        pmax(II(sset)[,'M'],1) / pmax(II(sset)[,'M']+II(sset)[,'U'],2))

    if (mask) {
        betas[sset@extra$mask[names(betas)]] <- NA
    }

    betas
}

#' Get allele frequency treating type I by summing alleles
#'
#' Takes a \code{SigSet} as input and returns a numeric vector containing
#' extra allele frequencies based on Color-Channel-Switching (CCS) probes.
#' If no CCS probes exist in the \code{SigSet}, then an numeric(0) is
#' returned.
#'
#' @param sset \code{SigSet}
#' @param known.ccs.only consider only known CCS probes
#' @return beta values
#' @examples
#' sset <- sesameDataGet('EPIC.1.LNCaP')$sset
#' betas <- getAFTypeIbySumAlleles(sset)
#' @export
getAFTypeIbySumAlleles <- function(sset, known.ccs.only = TRUE) {

    stopifnot(is(sset, "SigSet"))

    if (any(rownames(oobR(sset)) != rownames(IG(sset))))
        stop("oobR-IG not matched. Likely a malformed sset.");
    if (any(rownames(oobG(sset)) != rownames(IR(sset))))
        stop("oobG-IR not matched. Likely a malformed sset.");

    ## .oobG <- oobG[rownames(IR),]
    af <- c(
        pmax(rowSums(oobR(sset)),1)/(
            pmax(rowSums(oobR(sset))+rowSums(IG(sset)),2)),
        pmax(rowSums(oobG(sset)),1)/(
            pmax(rowSums(oobG(sset))+rowSums(IR(sset)),2)))

    if (known.ccs.only)
        af <- af[intersect(
            names(af),
            sesameDataGet('ethnicity.inference')$ccs.probes)]

    af[order(names(af))]
}

## res is the output of illuminaio::readIDAT
inferPlatform <- function(res) {
    switch(res$ChipType,
        'BeadChip 8x5'='EPIC',
        'BeadChip 12x8'='HM450',
        'BeadChip 12x1'='HM27',
        "Custom")
}

## Import one IDAT file
## return a data frame with 2 columns, corresponding to
## cy3 (Grn) and cy5 (Red) color channel signal
readIDAT1 <- function(grn.name, red.name, platform='') {
    ida.grn <- suppressWarnings(illuminaio::readIDAT(grn.name));
    ida.red <- suppressWarnings(illuminaio::readIDAT(red.name));
    d <- cbind(
        cy3=ida.grn$Quants[,"Mean"],
        cy5=ida.red$Quants[,"Mean"])
    colnames(d) <- c('G', 'R')

    if (platform != '') {
        attr(d, 'platform') <- platform
    } else {
        ## this is not always accurate
        ## TODO should identify unique tengo IDs.
        attr(d, 'platform') <- inferPlatform(ida.red)
    }
    d
}

#' Import a pair of IDATs from one sample
#'
#' The function takes a prefix string that are shared with _Grn.idat
#' and _Red.idat. The function returns a \code{SigSet}.
#'
#' @param prefix.path sample prefix without _Grn.idat and _Red.idat
#' @param manifest optional design manifest file
#' @param controls optional control probe manifest file
#' @param verbose     be verbose?  (FALSE)
#' @param platform EPIC, HM450 and HM27 etc.
#' 
#' @return a \code{SigSet}
#' 
#' @examples
#' sset <- readIDATpair(sub('_Grn.idat','',system.file(
#'     "extdata", "4207113116_A_Grn.idat", package = "sesameData")))
#' @export
readIDATpair <- function(
    prefix.path, platform = '', manifest = NULL,
    controls = NULL, verbose=FALSE) {

    if (file.exists(paste0(prefix.path, '_Grn.idat'))) {
        grn.name <- paste0(prefix.path, '_Grn.idat')
    } else if (file.exists(paste0(prefix.path, '_Grn.idat.gz'))) {
        grn.name <- paste0(prefix.path, '_Grn.idat.gz')
    } else {
        stop('Grn IDAT does not exist')
    }
    
    if (file.exists(paste0(prefix.path, '_Red.idat'))) {
        red.name <- paste0(prefix.path, '_Red.idat')
    } else if (file.exists(paste0(prefix.path, '_Red.idat.gz'))) {
        red.name <- paste0(prefix.path, '_Red.idat.gz')
    } else {
        stop('Red IDAT does not exist')
    }

    if (verbose == TRUE) {
        message("Reading IDATs for ", basename(prefix.path), "...")
    }

    dm <- readIDAT1(grn.name, red.name, platform=platform)

    if (is.null(manifest)) { # pre-built platforms, EPIC, HM450, HM27 etc
        df_address <- sesameDataGet(paste0(
            attr(dm, 'platform'), '.address'))
        manifest <- df_address$ordering
        controls <- df_address$controls
    }

    ## this is critical, sset must have default one p-value
    sset <- pOOBAH(chipAddressToSignal(dm, manifest, controls))
    pval(sset) <- extra(sset)[['pvals']][['pOOBAH']]
    sset
}

#' Identify IDATs from a directory
#'
#' The input is the directory name as a string. The function identifies all
#' the IDAT files under the directory. The function returns a vector of such
#' IDAT prefixes under the directory.
#'
#' @param dir.name the directory containing the IDAT files.
#' @param recursive search IDAT files recursively
#' @param use.basename basename of each IDAT path is used as sample name
#' This won't work in rare situation where there are duplicate IDAT files.
#' @return the IDAT prefixes (a vector of character strings).
#'
#' @examples
#' ## only search what are directly under
#' IDATprefixes <- searchIDATprefixes(
#'     system.file("extdata", "", package = "sesameData"))
#'
#' ## search files recursively is by default
#' IDATprefixes <- searchIDATprefixes(
#'     system.file(package = "sesameData"), recursive=TRUE)
#' @export
searchIDATprefixes <- function(dir.name,
    recursive = TRUE, use.basename = TRUE) {

    stopifnot(file.exists(dir.name))

    paths <- list.files(dir.name, '\\.idat(.gz)?$', recursive = recursive)
    prefixes <- unique(sub('_(Grn|Red).idat(.gz)?', '', paths))

    df <- data.frame(
        paths=paths, 
        prefix=sub('_(Grn|Red).idat.*', '', paths),
        channel=sub('.*_(Grn|Red).idat.*', '\\1', paths))
    
    byprefix <- split(df, df$prefix)
    ## valid IDAT should has both Grn and Red
    is.valid <- vapply(byprefix, function(x) all(
        sort(x[,'channel']) == c('Grn','Red')), logical(1))
    
    prefixes <- names(is.valid)[is.valid]
    if (length(prefixes) == 0)
        stop("No IDAT file found.")

    prefixes <- file.path(dir.name, prefixes)
    
    ## set name attributes so that names are auto-assigned for
    ## lapply and mclapply
    if (use.basename) {
        names(prefixes) <- basename(prefixes)
    } else {
        names(prefixes) <- prefixes
    }
    prefixes
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
#' @param manifest a data frame with columns Probe_ID, M, U and col
#' @param controls a data frame with columns Address and Name. This is optional
#' but might be necessary for some preprocessing methods that depends on these
#' control probes. This is left for backward compatibility. Updated version
#' should have controls consolidated into manifest.
#' @return a SigSet, indexed by probe ID address
chipAddressToSignal <- function(
    dm, manifest, controls = NULL) {

    sset <- SigSet(attr(dm, 'platform'))

    ## Infinium I green channel / oob red channel
    IordG <- manifest[(!is.na(manifest$col))&(manifest$col=='G'),]
    ## 2-channel for green probes' M allele
    IuG2ch <- dm[match(IordG$U, rownames(dm)),]
    ## 2-channel for green probes' U allele
    ImG2ch <- dm[match(IordG$M, rownames(dm)),]
    oobR(sset) <- as.matrix(
        data.frame(M=ImG2ch[,'R'], U=IuG2ch[,'R'], row.names=IordG$Probe_ID))
    IG(sset) <- as.matrix(data.frame(
        M = ImG2ch[,'G'], U = IuG2ch[,'G'],
        row.names = IordG$Probe_ID))

    ## Infinium I red channel / oob green channel
    IordR <- manifest[(!is.na(manifest$col))&(manifest$col=='R'),]
    ## 2-channel for red probes' m allele
    IuR2ch <- dm[match(IordR$U, rownames(dm)),]
    ## 2-channel for red probes' u allele
    ImR2ch <- dm[match(IordR$M, rownames(dm)),]
    oobG(sset) <- as.matrix(
        data.frame(M=ImR2ch[,'G'], U=IuR2ch[,'G'], row.names=IordR$Probe_ID))
    IR(sset) <- as.matrix(data.frame(
        M = ImR2ch[,'R'], U = IuR2ch[,'R'],
        row.names = IordR$Probe_ID))
    IR(sset)[rowSums(is.na(IR(sset))) > 0] <- 0
    
    ## Infinium II
    if (isUniqProbeID(manifest$Probe_ID)) { # mouse array probe ID
        is.InfII <- probeID_designType(manifest$Probe_ID) == '2'
    } else { # traditional probe ID with just the cg number
        is.InfII <- is.na(manifest$col)
    }
    IIord <- manifest[is.InfII,]
    signal.II <- dm[match(IIord$U, rownames(dm)),c('G','R')]
    colnames(signal.II) <- c('M', 'U')
    rownames(signal.II) <- IIord$Probe_ID
    II(sset) <- signal.II

    ## control probes
    ctl_idx <- grep('^ctl',manifest$Probe_ID)
    if (length(ctl_idx) > 0) { # control probes are included in manifest
        ctl_ord <- manifest[ctl_idx,]
        ctl <- as.data.frame(dm[match(ctl_ord$U, rownames(dm)),])
        ctl <- cbind(ctl, ctl_ord$col, ctl_ord$Probe_ID)
        rownames(ctl) <- ctl_ord$Probe_ID
        colnames(ctl) <- c('G','R','col','type')
        ctl(sset) <- ctl
    } else if (!is.null(controls)) {
        ctl <- as.data.frame(dm[match(controls$Address, rownames(dm)),])
        rownames(ctl) <- make.names(controls$Name, unique=TRUE)
        ctl <- cbind(ctl, controls[, c("Color_Channel","Type")])
        colnames(ctl) <- c('G','R','col','type')
        ctl(sset) <- ctl
    }

    ## additional annotation in manifest
    if ('mask' %in% colnames(manifest)) {
        sset <- extraSet(
            sset, 'mask', setNames(manifest$mask, manifest$Probe_ID))
    }

    ## backward support for mouse array, to remove in the future
    if ('mapUniq' %in% colnames(manifest)) {
        sset <- extraSet(
            sset, 'mask', setNames(!manifest$mapUniq, manifest$Probe_ID))
    }
    
    sset
}

#' Compute internal bisulfite conversion control
#'
#' Compute GCT score for internal bisulfite conversion control. The function
#' takes a \code{SigSet} as input. The higher the GCT score, the more likely
#' the incomplete conversion.
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

    stopifnot(sset@platform %in% c('EPIC','HM450'))
    extC <- sesameDataGet(paste0(sset@platform, '.probeInfo'))$typeI.extC
    extT <- sesameDataGet(paste0(sset@platform, '.probeInfo'))$typeI.extT
    prbs <- rownames(oobG(sset))
    extC <- intersect(prbs, extC)
    extT <- intersect(prbs, extT)
    if (use.median) {
        median(oobG(sset)[extC,], na.rm=TRUE) /
            median(oobG(sset)[extT,], na.rm=TRUE)
    } else {
        mean(oobG(sset)[extC,], na.rm=TRUE) /
            mean(oobG(sset)[extT,], na.rm=TRUE)
    }
}

## R6 utility functions deleted after 1.5.0

