
#' Extract the probe type field from probe ID
#' This only works with the new probe ID system.
#' See https://github.com/zhou-lab/InfiniumAnnotation for illustration
#'
#' @param Probe_ID Probe ID
#' @return a vector of '1' and '2' suggesting Infinium-I and Infinium-II
#' @import stringr
#' @examples
#' probeID_designType("cg36609548_TC21")
#' @export
probeID_designType <- function(Probe_ID) {
    stopifnot(all(grepl('_', Probe_ID))) # make sure it's the new ID system
    vapply(Probe_ID, function(x) substr(
        strsplit(x,'_')[[1]][2],3,3), character(1))
}

#' Convert beta-value to M-value
#'
#' Logit transform a beta value vector to M-value vector.
#'
#' Convert beta-value to M-value (aka logit transform)
#' @param b vector of beta values
#' @return a vector of M values
#' @examples
#' BetaValueToMValue(c(0.1, 0.5, 0.9))
#' @export
BetaValueToMValue <- function(b) {
    log2(b/(1-b))
}

#' Convert M-value to beta-value
#'
#' Convert M-value to beta-value (aka inverse logit transform)
#' 
#' @param m a vector of M values
#' @return a vector of beta values
#' @examples
#' MValueToBetaValue(c(-3, 0, 3))
#' @export
MValueToBetaValue <- function(m) {
    2^m/(1+2^m)
}

#' Check SeSAMe versions
#'
#' print package verison of sesame and depended packages to help troubleshoot
#' installation issues.
#'
#' @return print the version of sesame, sesameData, biocondcutor and R
#' @importFrom utils packageVersion
#' @export
#' @examples
#' sesame_checkVersion()
sesame_checkVersion <- function() {
    rv <- R.Version()
    msg <- paste0(
        "SeSAMe requires matched versions of ",
        "R, sesame, sesameData and ExperimentHub.\n",
        "Here is the current versions installed:\n",
        sprintf("R: %s.%s\n", rv$major, rv$minor),
        sprintf("Bioconductor: %s\n", BiocManager::version()),
        sprintf("sesame: %s\n", packageVersion("sesame")),
        sprintf("sesameData: %s\n", packageVersion("sesameData")),
        sprintf("ExperimentHub: %s\n", packageVersion("ExperimentHub")))
    message(msg)
}


#' sesamize function is deprecated.
#' Please check https://github.com/zwdzwd/sesamize for previous scripts
#'
#' @param ... arguments for sesamize
#' @return a message text for deprecated function
#' @export
#' @examples
#' cat("Deprecated. see https://github.com/zwdzwd/sesamize")
sesamize <- function(...) {
    .Deprecated("https://github.com/zwdzwd/sesamize")
}
