
#' print a WGenerator
#'
#' This calls WGenerator and creates a WGroup to enclose the produced object.
#'
#' @param x a WGenerator object
#' @param ... additional options
#' @return the WGroup containing the plotting object
#' @export
print.WGenerator <- function(x, ...) {
  x <- x(NULL)
  group <- ResolvedWGroup(x)
  print(group)
  return(group)
}
