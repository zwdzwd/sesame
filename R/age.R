
#' Predict age using linear models
#'
#' The function takes a named numeric vector of beta values. The name attribute
#' contains the probe ID (cg, ch or rs IDs). The function looks for overlapping
#' probes and estimate age using different models.
#'
#' You can get the models such as the Horvath aging model (Horvath 2013
#' Genome Biology) from sesameDataGet. The function outputs a single numeric
#' of age in years.
#'
#' Here are some built-in age models:
#' Anno/HM450/Clock_Horvath353.rds
#' Anno/HM450/Clock_Hannum.rds
#' Anno/HM450/Clock_SkinBlood.rds
#' Anno/EPIC/Clock_PhenoAge.rds
#' Anno/MM285/Clock_Zhou347.rds
#' see vignette inferences.html#Age__Epigenetic_Clock for details
#' 
#' @param betas a probeID-named vector of beta values
#' @param model a model object from sesameDataGet. should contain
#' param, intercept, response2age. default to the Horvath353 model.
#' @param na_fallback use fall back values if na
#' @param min_nonna the minimum number of non-NA values.
#' @return age in the unit specified in the model (usually in year, but
#' sometimes can be month, like in the mouse clocks).
#' @examples
#' betas <- sesameDataGet('HM450.1.TCGA.PAAD')$betas
#' \dontrun{
#' ## download age models from
#' ## https://github.com/zhou-lab/InfiniumAnnotationV1/tree/main/Anno
#' ## e.g., Anno/HM450/Clock_Horvath353.rds
#' predictAge(betas, model)
#' }
#' @export
predictAge <- function(betas, model, na_fallback=FALSE, min_nonna = 10) {

    betas <- betas[model$param$Probe_ID]
    if (sum(!is.na(betas)) < min_nonna) {
        stop("Fewer than 10 matching probes left. Age prediction abort.")
    }
    if (sum(is.na(betas)) > 0) {
        if (na_fallback) {
            k <- is.na(betas)
            betas[k] <- model$param$na_fallback[k]
        } else {
            probes <- intersect(names(na.omit(betas)), model$param$Probe_ID)
            betas <- betas[probes]
            model$param <- model$param[match(probes, model$param$Probe_ID),]
        }
    }
    drop(model$response2age(betas %*% model$param$slope + model$intercept))
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
#' cat("Deprecated. See predictAge")
#' @export
predictMouseAgeInMonth <- function(betas, na_fallback=TRUE) {
    .Deprecated("predictAge")
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
#' @examples
#' cat("Deprecated. See predictAge")
#' @export
predictAgeHorvath353 <- function(betas) {
    .Deprecated("predictAge")
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
#' @examples
#' cat("Deprecated. See predictAge")
#' @export
predictAgeSkinBlood <- function(betas) {
    .Deprecated("predictAge")
}


## Hv.age2response <- function(x, adult.age=20) {
##     ## trafo
##     x <- (x+1)/(adult.age+1)
##     ifelse(x<=1,log(x),x-1)
## }

## Hv.response2age <- function(x, adult.age=20) {
##     ## anti.trafo
##     ifelse(
##         x<0,
##         (1+adult.age)*exp(x)-1,
##         (1+adult.age)*x+adult.age)
## }
