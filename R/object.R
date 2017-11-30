#' Construct a WObject
#'
#' @param dm position
#' @param name name
#' @return a WObject
#' @export
WObject <- function(dm=NULL, name='') {
  structure(list(dm=dm, name=name), class='WObject')
}


#' merge plotting objects
#'
#' @param group a WGroup or a plotting object
#' @param p a new plotting object
#' @return a WGroup
#' @export
`+.WObject` <- function(group, p) {

  ## if first plotting object, create a WGroup to enclose it
  if (!('WGroup' %in% class(group))) {
    group <- Resolve(group, NULL)
    group <- ResolvedWGroup(group)
  }

  if ('WCustomize' %in% class(p)) {
    group <- p(group)
    return (group)
  }

  p <- Resolve(p, group)
  group <- AddWGroup(group, p)
  group
}

