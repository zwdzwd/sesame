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
#' @param betas a probeID-named vector of beta values
#' @return age in years
#' @export
#' @examples
#' betas <- readRDS(system.file(
#'     'extdata','HM450.betas.TCGA-2L-AAQA-01A-21D-A38H-05.rds',
#'     package='sesameData'))
#' predictAgeHorvath353(betas)
predictAgeHorvath353 <- function(betas) {

    Hv <- sesameData::agePredHorvath353
    probes <- intersect(names(na.omit(betas)), Hv$CpGmarker[-1])

    if (length(probes) < 10) {
        stop('Fewer than 10 matching probes left. Age prediction abort.')
    }

    drop(Hv.response2age(
        Hv$CoefficientTraining[1] + 
            Hv$CoefficientTraining[match(probes, Hv$CpG)] %*% betas[probes]))
}
