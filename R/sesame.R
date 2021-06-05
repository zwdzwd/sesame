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
#' sdf <- readIDATpair(sub('_Grn.idat','',system.file(
#'     'extdata','4207113116_A_Grn.idat',package='sesameData')))
#'
#' ## The OpenSesame pipeline
#' betas <- openSesame(sdf)
#'
#' @keywords DNAMethylation Microarray QualityControl
#'
"_PACKAGE"

#' report M and U for regular probes
#'
#' @param sdf a \code{SigDF}
#' @param mask whether to apply mask
#' @return a data frame of M and U columns
#' @examples
#' sesameDataCache("EPIC") # if not done yet
#' sdf <- sesameDataGet('EPIC.1.SigDF')
#' head(signalMU(sdf))
#' @export
signalMU <- function(sdf, mask = TRUE) {
    s2 = rbind(
        with(sdf[sdf$col=="G",], data.frame(M=MG, U=UG, Probe_ID=Probe_ID)),
        with(sdf[sdf$col=="R",], data.frame(M=MR, U=UR, Probe_ID=Probe_ID)),
        with(sdf[sdf$col=="2",], data.frame(M=UG, U=UR, Probe_ID=Probe_ID)))
    s2 = s2[match(sdf$Probe_ID, s2$Probe_ID),]
    if (mask) { s2 = s2[!sdf$mask,] }
    rownames(s2) = NULL
    s2
}

#' Mean Intensity
#'
#' The function takes one single \code{SigDF} and computes mean
#' intensity of all the in-band measurements. This includes all Type-I
#' in-band measurements and all Type-II probe measurements. Both methylated
#' and unmethylated alleles are considered. This function outputs a single
#' numeric for the mean.
#'
#' @param sdf a \code{SigDF}
#' @param mask whether to mask probes using mask column
#' 
#' @return mean of all intensities
#' @examples
#' sesameDataCache("EPIC") # if not done yet
#' sdf <- sesameDataGet('EPIC.1.SigDF')
#' meanIntensity(sdf)
#' @export
meanIntensity <- function(sdf, mask = TRUE) {
    stopifnot(is(sdf, "SigDF"))
    with(signalMU(sdf, mask = mask), mean(c(M,U), na.rm=TRUE))
}

#' M+U Intensities for All Probes
#'
#' The function takes one single \code{SigDF} and computes total
#' intensity of all the in-band measurements by summing methylated and
#' unmethylated alleles. This function outputs a single numeric for the mean.
#' @param sdf a \code{SigDF}
#' @return a vector of M+U signal for each probe
#' @examples
#' sesameDataCache("EPIC") # if not done yet
#' sdf <- sesameDataGet('EPIC.1.SigDF')
#' intensities <- totalIntensities(sdf)
#' @export
totalIntensities <- function(sdf, mask = TRUE) {
    with(signalMU(sdf), setNames(M+U, Probe_ID))
}

subsetvec <- function(vec, vecnames) {
    vec[names(vec) %in% vecnames]
}

#' Get sex-related information
#'
#' The function takes a \code{SigDF} and returns a vector of three
#' numerics: the median intensity of chrY probes; the median intensity of
#' chrX probes; and fraction of intermediate chrX probes. chrX and chrY
#' probes excludes pseudo-autosomal probes.
#'
#' @param sdf a \code{SigDF}
#' @return medianY and medianX, fraction of XCI, methylated and unmethylated X
#' probes, median intensities of auto-chromosomes.
#' @examples
#' sesameDataCache("EPIC") # if not done yet
#' sdf <- sesameDataGet('EPIC.1.SigDF')
#' getSexInfo(sdf)
#' @export
getSexInfo <- function(sdf) {

    stopifnot(is(sdf, "SigDF"))

    cleanY = sesameDataGet(paste0(
        platform(sdf),'.probeInfo'))$chrY.clean

    xLinked = sesameDataGet(paste0(
        platform(sdf),'.probeInfo'))$chrX.xlinked

    probe2chr = sesameDataGet(paste0(
        platform(sdf),'.probeInfo'))$probe2chr.hg19

    xLinkedBeta = getBetas(sdf, mask=FALSE)[xLinked]
    intens = totalIntensities(sdf)
    probes = intersect(names(intens), names(probe2chr))
    intens = intens[probes]
    probe2chr = probe2chr[probes]

    c(
        medianY=median(subsetvec(intens, cleanY)),
        medianX=median(subsetvec(intens, xLinked)),
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
#' The function takes a \code{SigDF} and infers the sex chromosome Karyotype
#' and presence/absence of X-chromosome inactivation (XCI). chrX, chrY and XCI
#' are inferred relatively independently. This function gives a more detailed
#' look of potential sex chromosome aberrations.
#'
#' @param sdf a \code{SigDF}
#' @return Karyotype string, with XCI
#' @examples
#' sesameDataCache("EPIC") # if not done yet
#' sdf <- sesameDataGet('EPIC.1.SigDF')
#' inferSexKaryotypes(sdf)
#' @export
inferSexKaryotypes <- function(sdf) {
    stopifnot(is(sdf, "SigDF"))
    sex.info <- getSexInfo(sdf)
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
#' @param sdf a \code{SigDF}
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
#' sesameDataCache("EPIC") # if not done yet
#' sdf <- sesameDataGet('EPIC.1.SigDF')
#' inferSex(sdf)
#' @export
inferSex <- function(sdf) {
    stopifnot(is(sdf, "SigDF"))
    stopifnot(platform(sdf) %in% c('EPIC','HM450'))
    sex.info <- getSexInfo(sdf)[seq_len(3)]
    as.character(predict(
        sesameDataGet('sex.inference'), sex.info))
}

#' Infer Ethnicity
#'
#' This function uses both the built-in rsprobes as well as the type I
#' Color-Channel-Switching probes to infer ethnicity.
#'
#' s better be background subtracted and dyebias corrected for
#' best accuracy
#'
#' Please note: the betas should come from SigDF *without*
#' channel inference.
#'
#' @param sdf a \code{SigDF}
#' @return string of ethnicity
#' @importFrom randomForest randomForest
#' @import sesameData
#' @examples
#' sdf <- sesameDataGet('EPIC.1.SigDF')
#' inferEthnicity(sdf)
#' @export
inferEthnicity <- function(sdf) {

    stopifnot(is(sdf, 'SigDF'))
    stopifnot(platform(sdf) %in% c('EPIC','HM450'))

    ethnicity.inference <- sesameDataGet('ethnicity.inference')
    ccsprobes <- ethnicity.inference$ccs.probes
    rsprobes <- ethnicity.inference$rs.probes
    ethnicity.model <- ethnicity.inference$model
    af <- c(
        getBetas(sdf, mask=FALSE)[rsprobes],
        getAFTypeIbySumAlleles(sdf, known.ccs.only = FALSE)[ccsprobes])

    as.character(predict(ethnicity.model, af))
}

#' Get beta Values
#'
#' sum.typeI is used for rescuing beta values on
#' Color-Channel-Switching CCS probes. The function takes a \code{SigDF}
#' and returns beta value except that Type-I in-band signal and out-of-band
#' signal are combined. This prevents color-channel switching due to SNPs.
#' 
#' @param sdf \code{SigDF}
#' @param sum.TypeI whether to sum type I channels
#' @param mask whether to use mask
#' @return a numeric vector, beta values
#' @examples
#' sesameDataCache("EPIC") # if not done yet
#' sdf <- sesameDataGet('EPIC.1.SigDF')
#' betas <- getBetas(sdf)
#' @export
getBetas <- function(sdf, mask=TRUE, sum.TypeI = FALSE) {

    stopifnot(is(sdf, "SigDF"))

    if (sum.TypeI) {
        betas = c(
            with(sdf[sdf$col != "2",],
                setNames(pmax(MG+MR,1) / pmax(MG+MR+UG+UR,2), Probe_ID)),
            with(sdf[sdf$col == "2",],
                setNames(pmax(UG,1) / pmax(UG+UR,2), Probe_ID)))
    } else {
        betas = c(
            with(sdf[sdf$col == "G",],
                setNames(pmax(MG,1) / pmax(MG+UG,2), Probe_ID)),
            with(sdf[sdf$col == "R",],
                setNames(pmax(MR,1) / pmax(MR+UR,2), Probe_ID)),
            with(sdf[sdf$col == "2",],
                setNames(pmax(UG,1) / pmax(UG+UR,2), Probe_ID)))
    }

    betas = betas[match(sdf$Probe_ID, names(betas))] # always use the original order
    if (mask) {
        betas[sdf$mask] = NA
    }
    betas
}

#' Get allele frequency treating type I by summing alleles
#'
#' Takes a \code{SigDF} as input and returns a numeric vector containing
#' extra allele frequencies based on Color-Channel-Switching (CCS) probes.
#' If no CCS probes exist in the \code{SigDF}, then an numeric(0) is
#' returned.
#'
#' @param sdf \code{SigDF}
#' @param known.ccs.only consider only known CCS probes
#' @return beta values
#' @examples
#' sesameDataCache("EPIC") # if not done yet
#' sdf <- sesameDataGet('EPIC.1.SigDF')
#' af <- getAFTypeIbySumAlleles(sdf)
#' @export
getAFTypeIbySumAlleles <- function(sdf, known.ccs.only = TRUE) {

    stopifnot(is(sdf, "SigDF"))

    af = c(
        with(sdf[sdf$col == "G",],
            setNames(pmax(MR+UR,1) / pmax(MR+UR+MG+UG,2), Probe_ID)),
        with(sdf[sdf$col == "R",],
            setNames(pmax(MG+UG,1) / pmax(MR+UR+MG+UG,2), Probe_ID)))

    if (known.ccs.only) {
        af = af[intersect(
            names(af), sesameDataGet('ethnicity.inference')$ccs.probes)]
    }
    
    af[order(names(af))] # TODO: avoid sorting
}

## res is the output of illuminaio::readIDAT
## Infer platform from IDATs
inferPlatform <- function(res) {
    sig <- sesameDataGet('idatSignature')
    names(which.max(vapply(
        sig, function(x) sum(x %in% rownames(res$Quants)), integer(1))))
}

inferPlatformFromProbeIDs <- function(probeIDs) {
    sig = sesameDataGet("probeIDSignature")
    names(which.max(vapply(
        sig, function(x) sum(probeIDs %in% x), integer(1))))
}

defaultAssembly <- function(platform) {
    platform2build = c(
        "HM27"="hg38",
        "HM450"="hg38",
        "EPIC"="hg38",
        "MM285"="mm10",
        "Mammal40"="hg38"
    )
    if (!(platform %in% names(platform2build))) {
        stop(sprintf(
            "Platform %s not supported. Try custom manifest?", platform))
    }
    platform2build[platform]
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
#' and _Red.idat. The function returns a \code{SigDF}.
#'
#' @param prefix.path sample prefix without _Grn.idat and _Red.idat
#' @param manifest optional design manifest file
#' @param controls optional control probe manifest file
#' @param verbose     be verbose?  (FALSE)
#' @param platform EPIC, HM450 and HM27 etc.
#' 
#' @return a \code{SigDF}
#' 
#' @examples
#' sdf <- readIDATpair(sub('_Grn.idat','',system.file(
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

    sdf = chipAddressToSignal(dm, manifest)
    if (!is.null(controls)) {
        attr(sdf, "controls") = readControls(dm, controls)
    }
    sdf
}

readControls <- function(dm, controls) {
    ctl <- as.data.frame(dm[match(controls$Address, rownames(dm)),])
    rownames(ctl) <- make.names(controls$Name, unique=TRUE)
    ctl <- cbind(ctl, controls[, c("Color_Channel","Type")])
    colnames(ctl) <- c('G','R','col','type')
    ctl
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
#' @param mft a data frame with columns Probe_ID, M, U and col
#' @return a SigDF, indexed by probe ID address
chipAddressToSignal <- function(dm, mft) {

    mft1 = mft[!is.na(mft$col),]
    tmpM = dm[match(mft1$M, rownames(dm)),]
    tmpU = dm[match(mft1$U, rownames(dm)),]
    sdf = data.frame(
        Probe_ID=mft1$Probe_ID,
        MG=tmpM[,"G"], MR=tmpM[,"R"], UG=tmpU[,"G"], UR=tmpU[,"R"],
        col=mft1$col, mask=FALSE)
    if ("mask" %in% colnames(mft1)) { sdf$mask = mft1$mask; }

    mft2 = mft[is.na(mft$col),]
    if (nrow(mft2) > 0) {
        tmp = dm[match(mft2$U, rownames(dm)),]
        s2 = data.frame(
            Probe_ID=mft2$Probe_ID,
            MG=NA, MR=NA, UG=tmp[,"G"], UR=tmp[,"R"], col="2", mask=FALSE)
        if ("mask" %in% colnames(mft2)) { s2$mask = mft2$mask; }
        sdf = rbind(sdf, s2)
    }
    sdf$col = factor(sdf$col, levels=c("G","R","2"))
    sdf = sdf[match(mft$Probe_ID, sdf$Probe_ID),] # always the mft order
    sdf = structure(sdf, class=c("SigDF", "data.frame"))
    attr(sdf, "platform") = attr(dm, 'platform')
    rownames(sdf) = NULL
    sdf
}

SigSetToSigDF = function(sset) {
    df = rbind(
        data.frame(
            Probe_ID = rownames(sset@IG),
            MG = sset@IG[,"M"], MR = sset@oobR[,"M"],
            UG = sset@IG[,"U"], UR = sset@oobR[,"U"], col="G", mask=FALSE),
        data.frame(
            Probe_ID = rownames(sset@IR),
            MG = sset@oobG[,"M"], MR = sset@IR[,"M"],
            UG = sset@oobG[,"U"], UR = sset@IR[,"U"], col="R", mask=FALSE),
        data.frame(
            Probe_ID = rownames(sset@II),
            MG = NA, MR = NA,
            UG = sset@II[,"M"], UR = sset@II[,"U"], col="2", mask=FALSE))
    sdf = structure(df, class=c("SigDF", "data.frame"))
    sdf$col = factor(sdf$col, levels=c("G","R","2"))
    attr(sdf, "platform") = sset@platform
    attr(sdf, "controls") = sset@ctl
    rownames(sdf) = NULL
    sdf
}

## retired functions:
## after 1.10.5
## bisConversionControl, SigSet class, updateSigSet, subsetSignal,
## makeExampleSeSAMeDataSet, makeExampleTinyEPICDataSet, saveMask, 
## restoreMask, buildControlMatrix450k, detectionPoobEcdf2,
## detectionPnegNormTotal, detectionPnegNormGS, detectionPfixedNorm,
## detectionPnegNorm, noobsb, IG2,IGR2,oobG2, oobR2，SigSetMethod,
## SigSetList, SigSetsToRGChannelSet, RGChannelSetToSigSets,
## SigSetToRatioSet, parseGEOSignalABFile
## after 1.5.0
## R6 utility functions
