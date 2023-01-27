
#' Generate some additional color palettes
#'
#' @param pal a string for adhoc pals
#' @param n the number of colors for interpolation
#' @param space rgb or Lab
#' @return a palette-generating function
#' @examples
#' library(pals)
#' pal.bands(palgen("whiteturbo"))
#' @export
palgen <- function(pal, n=150, space = "Lab") {

    requireNamespace("pals")
    adhoc_pals <- list(
        whiteturbo = c("white","white",pals::turbo(10)[seq(2,10)]),
        whitejet = c("white","white","lightblue",
            "blue","green","yellow","orange","red","darkred"),
        whiteblack = c("white", "black"))

    if (length(pal) == 1 && is.character(pal) &&
            (pal %in% names(adhoc_pals))) {
        pal <- adhoc_pals[[pal]]
    }

    if (is.character(pal)) {
        requireNamespace("grDevices")
        grDevices::colorRampPalette(pal, space = space)
    } else if (is.function(pal)) {
        pal
    } else {
        stop("Please provide the right pal format.")
    }
}
