#' Customize an existing plot
#'
#' @param mar.left left margin [0.03]
#' @param mar.right right margin [0.03]
#' @param mar.bottom bottom margin [0.03]
#' @param mar.top top margin [0.03]
#' @param mar margin in all directions [0.03]
#' @return an object of class WCustomize
#' @examples
#' WHeatmap(matrix(c('fred','frank','brad',
#'                 'frank','fred','frank'), ncol=2)) +
#'   WLegendV(NULL, RightOf(), label.fontsize = 20) +
#'   WCustomize(mar.right=0.1)
#' @export
WCustomize <- function(mar.left=NULL, mar.right=NULL, mar.top=NULL, mar.bottom=NULL, mar=NULL) {
  wc <- list()
  invisible(lapply(names(as.list(match.call()))[-1], function (nm) {
    wc[[nm]] <<- get(nm)
  }))
  force(wc)
  structure(function(group) {
    if (!is.null(wc[['mar']])) {
      group$mar$left   <- wc[['mar']]
      group$mar$right  <- wc[['mar']]
      group$mar$top    <- wc[['mar']]
      group$mar$bottom <- wc[['mar']]
    }
    if (!is.null(wc$mar.left))
      group$mar$left <- wc$mar.left
    if (!is.null(wc$mar.right))
      group$mar$right <- wc$mar.right
    if (!is.null(wc$mar.top))
      group$mar$top <- wc$mar.top
    if (!is.null(wc$mar.bottom))
      group$mar$bottom <- wc$mar.bottom

    group
  }, class = c('WCustomize', 'WObject'))
}
