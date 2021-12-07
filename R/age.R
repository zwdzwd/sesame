
Hv.age2response <- function(x, adult.age=20) {
    ## trafo
    x=(x+1)/(adult.age+1)
    ifelse(x<=1,log(x),x-1)
}

Hv.response2age <- function(x, adult.age=20) {
    ## anti.trafo
    ifelse(
        x<0,
        (1+adult.age)*exp(x)-1,
        (1+adult.age)*x+adult.age)
}

#' Horvath 353 age predictor
#'
#' The function takes a named numeric vector of beta values. The name attribute
#' contains the probe ID (cg, ch or rs IDs). The function looks for overlapping
#' probes and estimate age using Horvath aging model (Horvath 2013
#' Genome Biology). The function outputs a single numeric of age in years.
#'
#' @param betas a probeID-named vector of beta values
#' @return age in years
#' @export
#' @examples
#'
#' betas <- sesameDataGet('HM450.1.TCGA.PAAD')$betas
#' predictAgeHorvath353(betas)
#' sesameDataClearCache()
#' 
predictAgeHorvath353 <- function(betas) {
    predictAge(betas, sesameDataGet('age.inference')$Horvath353)
}

#' Horvath Skin and Blood age predictor
#'
#' The function takes a named numeric vector of beta values. The name attribute
#' contains the probe ID (cg, ch or rs IDs). The function looks for overlapping
#' probes and estimate age using Horvath aging model (Horvath et al. 2018
#' Aging, 391 probes). The function outputs a single numeric of age in years.
#'
#' @param betas a probeID-named vector of beta values
#' @return age in years
#' @export
#' @examples
#'
#' betas <- sesameDataGet('HM450.1.TCGA.PAAD')$betas
#' predictAgeSkinBlood(betas)
#' sesameDataClearCache()
#' 
predictAgeSkinBlood <- function(betas) {
    predictAge(betas, sesameDataGet('age.inference')$SkinBlood)
}

predictAge <- function(betas, cf) {
    probes <- intersect(names(na.omit(betas)), cf$CpGmarker[-1])

    if (length(probes) < 10) {
        stop('Fewer than 10 matching probes left. Age prediction abort.')
    }

    drop(Hv.response2age(
        cf$CoefficientTraining[1] +
            cf$CoefficientTraining[
                match(probes, cf$CpGmarker)] %*% betas[probes]))
}

#' Mouse age predictor
#'
#' The function takes a named numeric vector of beta values. The name attribute
#' contains the probe ID. The function looks for overlapping
#' probes and estimate age using an aging model built from 321 MM285 probes.
#' The function outputs a single numeric of age in months. The clock is most
#' accurate with the sesame preprocessing.
#'
#' @param betas a probeID-named vector of beta values
#' @param na_fallback use the fallback default for NAs.
#' @return age in month
#' @examples
#'
#' betas = sesameDataGet('MM285.10.tissue')$betas
#' predictMouseAgeInMonth(betas[,1])
#' sesameDataClearCache()
#' 
#' @export
predictMouseAgeInMonth = function(betas, na_fallback=TRUE) {
    coefs = sesameDataGet("MM285.clock347")
    dat = betas[names(coefs$slopes)]
    if (sum(is.na(dat)) > 0 && na_fallback) {
        k = is.na(dat)
        dat[k] = coefs$na_fallback[names(k[k])]
    }
    sum(dat * coefs$slopes) + coefs$intercept
}
