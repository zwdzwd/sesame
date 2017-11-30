
#' WGrob object
#' plot from a gList of grob objects
#'
#' @param glist gList object
#' @param dm dimension
#' @param name name
#' @return WGrob object
#' @export
WGrob <- function(glist, dm=NULL, name='') {
  if (is.null(dm)) {
    dm <- WDim(0,0,1,1,nr=1,nc=1)
  }

  wg <- structure(list(
    glist = glist,
    dm=dm, name=name), class=c("WGrob", "WObject"))
  force(wg);
  structure(function(group) {
    wg$dm <- Resolve(wg$dm, group, nr=1, nc=1)
    wg
  }, class=c("WGenerator", "WObject"))
}

#' plot WGrob object
#'
#' @param x WGrob
#' @param cex scaling factor for text
#' @param layout.only plot layout
#' @param stand.alone produce a stand.alone plot
#' @param ... extra options
#' @export
print.WGrob <- function(x, cex=1, layout.only=FALSE, stand.alone=TRUE, ...) {

  if (stand.alone) {
    group <- ResolvedWGroup(x)
    print(group)
    return(group)
  }

  if (layout.only)
    return(.print.layout(x))
  
  pushViewport(viewport(x=unit(x$dm$left,'npc'), y=unit(x$dm$bottom,'npc'),
                        width=unit(x$dm$width,'npc'), height=unit(x$dm$height,'npc'),
                        just=c('left','bottom')))
  grid.draw(x$glist)
  upViewport()
}

CalcTextBounding.WGrob <- function(x, group) {
  dm <- DimToTop(x, group)
  dm <- DimNPCToPoints(dm)
  dm
}
