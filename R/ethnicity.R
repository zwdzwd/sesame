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
#' @param verbose print more messages
#' @return string of ethnicity
#' @import sesameData
#' @examples
#' sdf <- sesameDataGet('EPIC.1.SigDF')
#' ## inferEthnicity(sdf)
#' @export
inferEthnicity <- function(sdf, verbose = FALSE) {
    .Deprecated("Please use CytoMethIC::cmi_classify.")
}
