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

#' Infer and reset color channel for Type-I probes instead of
#' using what is specified in manifest
#' 
#' @param sset a \code{SigSet}
#' @param verbose whether to print correction summary
#' @param switch_failed whether to switch failed probes
#' @import matrixStats
#' @return a \code{SigSet}
#' @examples
#'
#' sset <- sesameDataGet('EPIC.1.LNCaP')$sset
#' inferTypeIChannel(sset)
#' 
#' @export
inferTypeIChannel <- function(sset, switch_failed = TRUE, verbose = TRUE) {
    red_channel <- rbind(IR(sset), oobR(sset))
    grn_channel <- rbind(oobG(sset), IG(sset))
    red_idx0 <- seq_len(nrow(red_channel)) <= nrow(IR(sset)) # old red index
    red_max <- rowMaxs(red_channel)
    grn_max <- rowMaxs(grn_channel)
    red_idx <- red_max > grn_max # new red index

    ## stop inference when in-band signal is lower than a minimum
    min_ib <- quantile(pmin(rowMins(red_channel), rowMins(grn_channel)), 0.95)
    big_idx <- pmax(red_max, grn_max) > min_ib # in-band is big enough? 

    if (verbose) {
        message('Type-I color channel reset:')
        message('R>R: ', sum(red_idx0 & red_idx & big_idx))
        message('G>G: ', sum(!red_idx0 & !red_idx & big_idx))
        message('R>G: ', sum(red_idx0 & !red_idx & big_idx))
        message('G>R: ', sum(!red_idx0 & red_idx & big_idx))
        message('Failed: ', sum(!big_idx))
    }

    if (switch_failed) {
        IR(sset) <- red_channel[red_idx,]
        oobG(sset) <- grn_channel[red_idx,]
        IG(sset) <- grn_channel[!red_idx,]
        oobR(sset) <- red_channel[!red_idx,]
    } else {
        IR(sset) <- rbind(
            red_channel[red_idx & big_idx,],
            red_channel[red_idx0 & !big_idx,])
        oobG(sset) <- rbind(
            grn_channel[red_idx & big_idx,],
            grn_channel[red_idx0 & !big_idx,])
        IG(sset) <- rbind(
            grn_channel[!red_idx & big_idx,],
            grn_channel[!red_idx0 & !big_idx,])
        oobR(sset) <- rbind(
            red_channel[!red_idx & big_idx,],
            red_channel[!red_idx0 & !big_idx,])
    }
    sset
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

    ethnicity.inference <- sesameDataGet('ethnicity.inference')
    ccsprobes <- ethnicity.inference$ccs.probes
    rsprobes <- ethnicity.inference$rs.probes
    ethnicity.model <- ethnicity.inference$model
    af <- c(
        getBetas(
            subsetSignal(sset, rsprobes), quality.mask = FALSE,
            nondetection.mask=FALSE),
        getAFTypeIbySumAlleles(
            subsetSignal(sset, ccsprobes)))

    as.character(predict(ethnicity.model, af))
}

#' Get beta Values
#'
#' sum.typeI is used for rescuing beta values on
#' Color-Channel-Switching CCS probes. The function takes a \code{SigSet}
#' and returns beta value except that Type-I in-band signal and out-of-band
#' signal are combined. This prevents color-channel switching due to SNPs.
#' 
#' @param sset \code{SigSet}
#' @param quality.mask whether to mask low quality probes
#' @param nondetection.mask whether to mask nondetection
#' @param mask.use.tcga whether to use TCGA masking, only applies to HM450
#' @param pval.threshold p-value threshold for nondetection mask
#' @param sum.TypeI whether to sum type I channels
#' @return a numeric vector, beta values
#' @examples
#' sset <- sesameDataGet('EPIC.1.LNCaP')$sset
#' betas <- getBetas(sset)
#' @export
getBetas <- function(
    sset,
    quality.mask = TRUE,
    nondetection.mask = TRUE,
    mask.use.tcga = FALSE,
    pval.threshold = 0.05,
    sum.TypeI = FALSE) {
    
    if (is(sset, "SigSetList")) {
        return(do.call(cbind, lapply(
            sset, getBetas, quality.mask = quality.mask, 
            nondetection.mask = nondetection.mask,
            mask.use.tcga = mask.use.tcga, pval.threshold = pval.threshold)))
    }
        
    stopifnot(is(sset, "SigSet"))

    ## optionally summing channels protects
    ## against channel misspecification
    IGs <- IG(sset)
    IRs <- IR(sset)
    if (sum.TypeI) {
        IGs <- IGs + oobR(sset)
        IRs <- IRs + oobG(sset)
    }
    
    betas1 <- pmax(IGs[,'M'],1) / pmax(IGs[,'M']+IGs[,'U'],2)
    betas2 <- pmax(IRs[,'M'],1) / pmax(IRs[,'M']+IRs[,'U'],2)
    betas3 <- pmax(II(sset)[,'M'],1) / pmax(II(sset)[,'M']+II(sset)[,'U'],2)
    betas <- c(betas1, betas2, betas3)
    if (nondetection.mask) {
        pval(sset) <- pval(sset)[match(names(betas), names(pval(sset)))]
        betas[pval(sset) > pval.threshold] <- NA
    }

    ## currently quality masking only supports three platforms
    if (quality.mask && sset@platform %in% c('HM27','HM450','EPIC')) {
        if (mask.use.tcga) {
            stopifnot(sset@platform == 'HM450')
            mask <- sesameDataGet('HM450.probeInfo')$mask.tcga
        } else {
            mask <- sesameDataGet(paste0(sset@platform, '.probeInfo'))$mask
        }
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

## Import one IDAT file
## return a data frame with 2 columns, corresponding to
## cy3 (Grn) and cy5 (Red) color channel signal
readIDAT1 <- function(grn.name, red.name, platform='') {
    ida.grn <- illuminaio::readIDAT(grn.name);
    ida.red <- illuminaio::readIDAT(red.name);
    d <- cbind(cy3=ida.grn$Quants[,"Mean"], cy5=ida.red$Quants[,"Mean"])
    colnames(d) <- c('G', 'R')

    if (platform != '') {
        attr(d, 'platform') <- platform
    } else {
        ## this is not always accurate
        ## TODO should identify unique tengo IDs.
        attr(d, 'platform') <- switch(
            ida.red$ChipType,
            'BeadChip 8x5'='EPIC',
            'BeadChip 12x8'='HM450',
            'BeadChip 12x1'='HM27')
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
    prefix.path, platform = '',
    manifest = NULL, verbose=FALSE) {

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

    chipAddressToSignal(dm, manifest, controls)
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
chipAddressToSignal <- function(dm, manifest, controls = NULL) {

    platform <- attr(dm, 'platform')

    sset <- SigSet(platform)

    ## type I green channel / oob red channel
    ## IordG <- manifest[((manifest$DESIGN=='I')&(manifest$col=='G')),]
    IordG <- manifest[(!is.na(manifest$col))&(manifest$col=='G'),]
    ## 2-channel for green probes' M allele
    IuG2ch <- dm[match(IordG$U, rownames(dm)),]
    IuG <- IuG2ch[,1]
    ## 2-channel for green probes' U allele
    ImG2ch <- dm[match(IordG$M, rownames(dm)),]
    ImG <- ImG2ch[,1]
    oobR(sset) <- as.matrix(
        data.frame(M=ImG2ch[,2], U=IuG2ch[,2], row.names=IordG$Probe_ID))
    IG(sset) <- as.matrix(data.frame(M=ImG, U=IuG, row.names=IordG$Probe_ID))

    ## type I red channel / oob green channel
    ## IordR <- manifest[((manifest$DESIGN=='I')&(manifest$col=='R')),]
    IordR <- manifest[(!is.na(manifest$col))&(manifest$col=='R'),]
    ## 2-channel for red probes' m allele
    IuR2ch <- dm[match(IordR$U, rownames(dm)),]
    IuR <- IuR2ch[,2]
    ## 2-channel for red probes' u allele
    ImR2ch <- dm[match(IordR$M, rownames(dm)),]
    ImR <- ImR2ch[,2]
    oobG(sset) <- as.matrix(
        data.frame(M=ImR2ch[,1], U=IuR2ch[,1], row.names=IordR$Probe_ID))
    IR(sset) <- as.matrix(data.frame(M=ImR, U=IuR, row.names=IordR$Probe_ID))

    ## type II
    ## IIord <- manifest[manifest$DESIGN=="II",]
    IIord <- manifest[is.na(manifest$col),]
    signal.II <- dm[match(IIord$U, rownames(dm)),]
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
    extC <- sesameDataGet(paste0(sset@platform, '.probeInfo'))$typeI.extC
    extT <- sesameDataGet(paste0(sset@platform, '.probeInfo'))$typeI.extT
    if (use.median) {
        median(oobG(sset)[extC,], na.rm=TRUE) /
            median(oobG(sset)[extT,], na.rm=TRUE)
    } else {
        mean(oobG(sset)[extC,], na.rm=TRUE) /
            mean(oobG(sset)[extT,], na.rm=TRUE)
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
