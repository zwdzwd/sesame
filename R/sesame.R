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
    stopifnot(is(sdf, "SigDF"))
    dG = InfIG(sdf); dR = InfIR(sdf); d2 = InfII(sdf)
    sdf2 = rbind(
        data.frame(M = dG$MG, U = dG$UG, Probe_ID = dG$Probe_ID),
        data.frame(M = dR$MR, U = dR$UR, Probe_ID = dR$Probe_ID),
        data.frame(M = d2$UG, U = d2$UR, Probe_ID = d2$Probe_ID))
    sdf2 = sdf2[match(sdf$Probe_ID, sdf2$Probe_ID),]
    if (mask) { sdf2 = sdf2[!sdf$mask,] }
    rownames(sdf2) = NULL
    sdf2
}

#' Whole-dataset-wide Mean Intensity
#'
#' The function takes one single \code{SigDF} and computes mean
#' intensity of all the in-band measurements. This includes all Type-I
#' in-band measurements and all Type-II probe measurements. Both methylated
#' and unmethylated alleles are considered. This function outputs a single
#' numeric for the mean.
#'
#' Note: mean in this case is more informative than median because
#' methylation level is mostly bimodal.
#'
#' @param sdf a \code{SigDF}
#' @param mask whether to mask probes using mask column
#' @return mean of all intensities
#' @examples
#' sesameDataCache("EPIC") # if not done yet
#' sdf <- sesameDataGet('EPIC.1.SigDF')
#' meanIntensity(sdf)
#' @export
meanIntensity <- function(sdf, mask = TRUE) {
    stopifnot(is(sdf, "SigDF"))
    s = signalMU(sdf, mask = mask)
    mean(c(s$M,s$U), na.rm=TRUE)
}

#' Whole-dataset-wide Median Total Intensity (M+U)
#'
#' The function takes one single \code{SigDF} and computes median
#' intensity of M+U for each probe. This function outputs a single
#' numeric for the median.
#'
#' @param sdf a \code{SigDF}
#' @param mask whether to mask probes using mask column
#' @return median of all intensities
#' @examples
#' sesameDataCache("EPIC") # if not done yet
#' sdf <- sesameDataGet('EPIC.1.SigDF')
#' medianTotalIntensity(sdf)
#' @export
medianTotalIntensity = function(sdf, mask = TRUE) {
    stopifnot(is(sdf, "SigDF"))
    s = signalMU(sdf, mask = mask)
    median(c(s$M + s$U), na.rm=TRUE)
}

#' Whole-dataset-wide Probe Success Rate
#'
#' This function calculates the probe success rate using
#' pOOBAH detection p-values. Probes that has a detection p-value
#' higher than a specific threshold are considered failed probes.
#' @param sdf a \code{SigDF}
#' @param max_pval the maximum p-value to consider detection success
#' @param mask whether or not we count the masked probes in SigDF
#' @return a fraction number as probe success rate
#' @examples
#' sesameDataCache("EPIC") # if not done yet
#' sdf = sesameDataGet('EPIC.1.SigDF')
#' probeSuccessRate(sdf)
#' @export
probeSuccessRate = function(sdf, mask = TRUE, max_pval = 0.05) {
    pval = pOOBAH(sdf, return.pval = TRUE)
    if (mask) { pval = pval[!sdf$mask] }
    pval = na.omit(pval)
    stopifnot(length(pval) > 100)
    sum(pval < max_pval) / length(pval)
}

#' M+U Intensities Array
#'
#' The function takes one single \code{SigDF} and computes total
#' intensity of all the in-band measurements by summing methylated and
#' unmethylated alleles. This function outputs a single numeric for the mean.
#' @param sdf a \code{SigDF}
#' @param mask whether to mask probes using mask column
#' @return a vector of M+U signal for each probe
#' @examples
#' sesameDataCache("EPIC") # if not done yet
#' sdf <- sesameDataGet('EPIC.1.SigDF')
#' intensities <- totalIntensities(sdf)
#' @export
totalIntensities <- function(sdf, mask = FALSE) {
    stopifnot(is(sdf, "SigDF"))
    s = signalMU(sdf, mask = mask)
    setNames(s$M+s$U, s$Probe_ID)
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
        d1 = InfI(sdf); d2 = InfII(sdf)
        betas = c(setNames(
            pmax(d1$MG+d1$MR,1)/pmax(d1$MG+d1$MR+d1$UG+d1$UR,2), d1$Probe_ID),
            setNames(pmax(d2$UG,1) / pmax(d2$UG+d2$UR,2), d2$Probe_ID))
    } else {
        dG = InfIG(sdf); dR = InfIR(sdf); d2 = InfII(sdf)
        betas = c(
            setNames(pmax(dG$MG,1) / pmax(dG$MG+dG$UG,2), dG$Probe_ID),
            setNames(pmax(dR$MR,1) / pmax(dR$MR+dR$UR,2), dR$Probe_ID),
            setNames(pmax(d2$UG,1) / pmax(d2$UG+d2$UR,2), d2$Probe_ID))
    }
    
    ## always use the original order
    betas = setNames(betas[match(sdf$Probe_ID, names(betas))], sdf$Probe_ID)
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

    dG = InfIG(sdf); dR = InfIR(sdf)
    af = c(setNames(
        pmax(dG$MR+dG$UR,1)/pmax(dG$MR+dG$UR+dG$MG+dG$UG,2), dG$Probe_ID),
        setNames(
            pmax(dR$MG+dR$UG,1)/pmax(dR$MR+dR$UR+dR$MG+dR$UG,2), dR$Probe_ID))

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
    ida.grn <- suppressWarnings(illuminaio::readIDAT(grn.name))
    ida.red <- suppressWarnings(illuminaio::readIDAT(red.name))
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
#' @param verbose be verbose?  (FALSE)
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
    controls = NULL, verbose = FALSE) {
    
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
        df_address = sesameDataGet(paste0(
            attr(dm, 'platform'), '.address'))
        manifest <- df_address$ordering
        controls <- df_address$controls
    }

    sdf = chipAddressToSignal(dm, manifest)
    if (!is.null(controls) && nrow(controls) > 0) {
        attr(sdf, "controls") = readControls(dm, controls)
    }
    sdf
}

readControls <- function(dm, controls) {
    if ("Color_Channel" %in% colnames(controls)) { # legacy control data
        ctl <- as.data.frame(dm[match(controls$Address, rownames(dm)),])
        rownames(ctl) <- make.names(controls$Name, unique=TRUE)
        ctl <- cbind(ctl, controls[, c("Color_Channel","Type")])
        colnames(ctl) <- c('G','R','col','type')
        ctl = ctl[!(is.na(ctl$G)|is.na(ctl$R)),] # no NA in controls
    } else {
        ctl = as.data.frame(chipAddressToSignal(dm, controls))
    }
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

#' Compute internal bisulfite conversion control
#'
#' Compute GCT score for internal bisulfite conversion control. The function
#' takes a \code{SigSet} as input. The higher the GCT score, the more likely
#' the incomplete conversion.
#'
#' @param sdf a SigDF
#' @return GCT score (the higher, the more incomplete conversion)
#' @examples
#' sesameDataCache("EPIC") # if not done yet
#' sdf <- sesameDataGet('EPIC.1.SigDF')
#' bisConversionControl(sdf)
#'
#' @export
bisConversionControl <- function(sdf) {

    stopifnot(sdfPlatform(sdf) %in% c('EPICplus','EPIC','HM450'))
    extC <- sesameDataGet(paste0(sdfPlatform(sdf), '.probeInfo'))$typeI.extC
    extT <- sesameDataGet(paste0(sdfPlatform(sdf), '.probeInfo'))$typeI.extT
    ## prbs <- rownames(oobG(sset))
    df = InfIR(sdf)
    extC = intersect(df$Probe_ID, extC)
    extT = intersect(df$Probe_ID, extT)
    dC = df[match(extC, df$Probe_ID),]
    dT = df[match(extT, df$Probe_ID),]
    mean(c(dC$MG, dC$UG), na.rm=TRUE) / mean(c(dT$MG, dT$UG), na.rm=TRUE)
}

## retired functions:
## after 1.11
## bisConversionControl, SigSet class, updateSigSet, subsetSignal,
## makeExampleSeSAMeDataSet, makeExampleTinyEPICDataSet, saveMask, 
## restoreMask, buildControlMatrix450k, detectionPoobEcdf2,
## detectionPnegNormTotal, detectionPnegNormGS, detectionPfixedNorm,
## detectionPnegNorm, noobsb, IG2,IGR2,oobG2, oobR2，SigSetMethod,
## SigSetList, SigSetsToRGChannelSet, RGChannelSetToSigSets,
## SigSetToRatioSet, parseGEOSignalABFile
## after 1.5.0
## R6 utility functions
