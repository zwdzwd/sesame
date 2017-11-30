#' WGG object
#' form ggplot with coordinates
#' 
#' @param ggobj ggplot plotting object
#' @param dm dimension
#' @param name name
#' @return WGG object
#' @export
WGG <- function(ggobj, dm=NULL, name='') {
  if (is.null(dm)) {
    dm <- WDim(0,0,1,1,nr=1, nc=1)
  }
  ggp <- list(ggobj=ggobj, dm=dm, name=name)
  class(ggp) <- c('WGG', 'WObject')

  force(ggp);
  structure(function(group) {
    ggp$dm <- Resolve(ggp$dm, group, nr=1, nc=1)
    ggp
  }, class=c('WGenerator', 'WObject'))
}

#' plot WGG object
#' 
#' @param x WGG
#' @param cex scaling factor for text
#' @param layout.only plot layout
#' @param stand.alone produce a stand.alone plot
#' @param ... extra options
#' @export
print.WGG <- function(x, cex=1, layout.only=FALSE, stand.alone=TRUE, ...) {

  print(x$ggobj, vp=viewport(x=unit(x$dm$left,'npc'), y=unit(x$dm$bottom,'npc'),
                             width=unit(x$dm$width,'npc'), height=unit(x$dm$height,'npc'),
                             just=c('left','bottom')))
}

CalcTextBounding.WGG <- function(x, group) {
  dm <- DimToTop(x, group)
  dm <- DimNPCToPoints(dm)
  dm
}
