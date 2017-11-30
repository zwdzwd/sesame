
#' WLegendV
#'
#' a vertical legend
#'
#' @param x a name or a plotting object, if NULL use the last plotting object
#' @param dm position
#' @param name name of the plotted legend
#' @param label.fontsize label font size
#' @param n.stops number of stops in computing continuous legend
#' @param n.text number of text labels in continuous legend
#' @param width width of each unit in plotted legend
#' @param height height of each unit in plotted legend
#' @param ... additional options to WHeatmap
#' @return an object of class WLegendV
#' @examples
#' WHeatmap(matrix(1:4,nrow=2))+WLegendV(NULL, RightOf())
#' @export
WLegendV <- function(x=NULL, dm=NULL, name='',
                     n.stops=20, n.text=5, label.fontsize=12,
                     width=0.05, height=0.02, ...) {

  kargs <- list(...)
  kargs$dm <- dm
  kargs$name <- name
  force(x); force(kargs);
  force(n.stops); force(n.text); force(label.fontsize);
  structure(function(group) {
    if (is.null(x))
      x <- group$children[[length(group$children)]]$name
    
    x <- Resolve(x, group)
    if (x$continuous) {
      d <- seq(from=x$cm$dmin, to=x$cm$dmax, length.out=n.stops)
      kargs$data <- matrix(
        d, dimnames=list(format(d, digits=2, trim=TRUE)))
    } else {
      d <- x$cm$mapper
      d <- d[order(names(d))]
      kargs$data <- matrix(names(d), dimnames=list(names(d), NULL))
      kargs$continuous <- FALSE
    }

    kargs$cm <- x$cm
    legend <- do.call(WHeatmap, kargs)(group)
    nr <- nrow(kargs$data)
    nc <- ncol(kargs$data)
    ## when dm is from TopOf etc use nr and nc
    ## when dm is from TopLeftOf etc use hard.dm
    legend$dm <- Resolve(dm, group, nr=nr, nc=nc,
                         hard.dm=WDim(0,0,width*nc,height*nr,nr=nr,nc=nc))
    legend$yticklabels <- TRUE
    if (x$continuous)
      legend$yticklabels.n <- n.text
    else
      legend$yticklabels.n <- nr
    legend$yticklabel.fontsize <- label.fontsize
    if (is.null(kargs$yticklabel.side))
      legend$yticklabel.side <- 'r'
    if (is.null(kargs$orientation))
      legend$orientation <- 'v'
    class(legend) <- c('WLegendV', class(legend))
    legend
  }, class=c('WGenerator','WObject'))
}

#' WLegendH
#'
#' a horizontal legend
#'
#' @param x a name or a plotting object, if NULL use the last plotting object
#' @param dm position
#' @param name name of the plotted legend
#' @param label.fontsize label font size
#' @param n.stops number of stops in computing continuous legend
#' @param n.text number of text labels in continuous legend
#' @param width width of each unit in plotted legend
#' @param height height of each unit in plotted legend
#' @param ... additional options to WHeatmap
#' @return an object of class WLegendH
#' @examples
#' WHeatmap(matrix(1:4,nrow=2))+WLegendH(NULL, Beneath())
#' @export
WLegendH <- function(x=NULL, dm=NULL, name='',
                     n.stops=20, n.text=5, label.fontsize=12,
                     width=0.02, height=0.05, ...) {

  kargs <- list(...)
  kargs$dm <- dm
  kargs$name <- name
  force(x); force(kargs);
  force(n.stops); force(n.text); force(label.fontsize);
  structure(function(group) {
    if (is.null(x))
      x <- group$children[[length(group$children)]]$name
    x <- Resolve(x, group)
    if (x$continuous) {
      d <- seq(from=x$cm$dmin, to=x$cm$dmax, length.out=n.stops)
      kargs$data <- matrix(
        d, nrow=1, dimnames=list(NULL, format(d, digits=2, trim=TRUE)))
    } else {
      d <- x$cm$mapper
      d <- d[order(names(d))]
      kargs$data <- matrix(names(d), dimnames=list(NULL, names(d)), nrow=1)
      kargs$continuous <- FALSE
    }
    kargs$cm <- x$cm
    legend <- do.call(WHeatmap, kargs)(group)
    nr <- nrow(kargs$data)
    nc <- ncol(kargs$data)
    ## when dm is from TopOf etc use nr and nc
    ## when dm is from TopLeftOf etc use hard.dm
    legend$dm <- Resolve(dm, group, nr=nr, nc=nc,
                         hard.dm=WDim(0,0,width*nc,height*nr, nr=nr, nc=nc))
    legend$xticklabels <- TRUE
    if (x$continuous)
      legend$xticklabels.n <- n.text
    else
      legend$xticklabels.n <- nc
    legend$xticklabel.fontsize <- label.fontsize
    if (is.null(kargs$xticklabel.side))
      legend$xticklabel.side <- 'b'
    if (is.null(kargs$orientation))
      legend$orientation <- 'h'
    class(legend) <- c('WLegendH', class(legend))
    legend
  }, class=c('WGenerator','WObject'))
}
