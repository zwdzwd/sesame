getIntensityRatioYvsAuto = function(sdf) {
    intens = totalIntensities(sdf)
    prbA = getAutosomeProbes(sdfPlatform(sdf))
    prbY = getProbesByChromosome("chrY", sdfPlatform(sdf))
    median(intens[prbY], na.rm=TRUE) / median(intens[prbA], na.rm=TRUE)
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
        sdfPlatform(sdf),'.probeInfo'))$chrY.clean

    xLinked = sesameDataGet(paste0(
        sdfPlatform(sdf),'.probeInfo'))$chrX.xlinked

    probe2chr = sesameDataGet(paste0(
        sdfPlatform(sdf),'.probeInfo'))$probe2chr.hg19

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
#' @param x either a raw \code{SigDF} or a beta value vector named by probe ID
#' SigDF is preferred over beta values.
#' @param pfm platform Only MM285, EPIC and HM450 are supported.
#' @return 'F' or 'M'
#' We established our sex calling based on the CpGs hypermethylated in
#' inactive X (XiH), CpGs hypomethylated in inactive X (XiL) and signal
#' intensity ratio of Y-chromosome over autosomes. Currently human inference
#' uses a random forest and mouse inference uses a support vector machine.
#'
#' The function checks the sample quality. If the sample is of poor quality
#' the inference return NA.
#' 
#' Note many factors such as Dnmt genotype, XXY male (Klinefelter's),
#' 45,X female (Turner's) can confuse the model sometimes.
#' This function works on a single sample.
#' @importFrom randomForest randomForest
#' @import sesameData
#' @import e1071
#' @examples
#' sesameDataCache("EPIC") # if not done yet
#' sdf = sesameDataGet('EPIC.1.SigDF')
#' inferSex(sdf)
#' @export
inferSex = function(x, pfm = NULL) {

    if (is.null(pfm)) {
        if (is.numeric(x)) {
            pfm = inferPlatformFromProbeIDs(names(x))
        } else if (is(x, "SigDF")) {
            pfm = sdfPlatform(x)
        }
    }
    
    stopifnot(pfm %in% c('EPIC','HM450','MM285'))
    if (pfm == 'MM285'){
        pkgTest("e1071")
        inf = sesameDataGet("MM285.inferences")$sex
        if (is(x, "SigDF")) { # if possible use signal intensity
            intensYvsAuto = getIntensityRatioYvsAuto(x)

            ## assuming unnormalized data
            betas = getBetas(pOOBAH(resetMask(dyeBiasNL(noob(x)))))

            ## if probe success rate is low, give up
            if (sum(is.na(betas)) / length(betas) > 0.3) { return(NA); }

            betasXiH = median(betas[inf$XiH], na.rm=TRUE)
            return(predict(inf$XY,
                data.frame(
                    betasXiH=betasXiH,
                    intensYvsAuto=intensYvsAuto))[[1]])
            
        } else if (is.numeric(x)) { # use beta value
            
            ## if probe success rate is low, give up
            if (sum(is.na(x)) / length(x) > 0.3) { return(NA); }
            betasXiH = median(x[inf$XiH], na.rm=TRUE)
            betasXiL = median(x[inf$XiL], na.rm=TRUE)
            
            return(predict(inf$X2,
                data.frame(betasXiL=betasXiL, betasXiH=betasXiH))[[1]])
            
        } else { return(NA); }
    } else { # need to update human sex inference to match the mouse format
        sex.info <- getSexInfo(x)[seq_len(3)]
        as.character(predict(
            sesameDataGet('sex.inference'), sex.info))
    }
}
