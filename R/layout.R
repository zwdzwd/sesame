
#' plot multiple figures in a matrix
#'
#' This function can take WObject, or gg (from ggplot)
#' since the coordinates are not set, gg can be converted to WGG
#' @param objs a list of plotting objects either WObject or gg
#' @param ncols number of columns
#' @return WGroup
#' @export
WMatrix <- function(objs, ncols=1) {
  objs <- lapply(objs, function(x) {
    if ('gg' %in% class(x)) {
      x <- WGG(x);
    }
    x <- Resolve(x);
    x
  })

  num.plots <- length(objs)
  nrows <- ceiling(num.plots / ncols)
  
  width <- 1.0 / ncols
  height <- 1.0 / nrows

  mtx <- NULL
  for (i in 1:num.plots) {
    bottom <- 1.0 / nrows * (nrows-ceiling(i / ncols))
    left <- ((i-1) %% ncols) * width
    objs[[i]]$dm <- WDim(left, bottom, width, height, nr=1, nc=1)
    if (i==1) mtx <- ResolvedWGroup(objs[[i]])
    else mtx <- mtx + objs[[i]]
  }
  mtx
}
