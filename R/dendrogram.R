
#' WDendrogram class
#'
#' WDendrogram class
#'
#' @param clust hclust object
#' @param dm plotting dimension
#' @param facing direction of the dendrogram plot
#' @param name name of the dendrogram plot
#' @return an object of class WDendrogram
#' @export
WDendrogram <- function(clust=NULL, dm=WDim(0,0,1,1), name='', facing=c("bottom", "top", "left", "right")) {

  stopifnot('hclust' %in% class(clust))
  dd <- lapply(formals(), eval)
  invisible(lapply(names(as.list(match.call()))[-1], function (nm) {
    dd[[nm]] <<- get(nm)
  }))
  class(dd) <- c('WDendrogram', 'WObject')
  force(dd)
  structure(function(group) {
    dd$dm <- Resolve(dm,group) # nr=1 nc=1
    dd
  }, class=c('WGenerator', 'WObject'))
}

#' print a dendrogram
#'
#' @param x a dendrogram
#' @param stand.alone plot is stand alone
#' @param layout.only plot layout only
#' @param cex factor to scaling texts
#' @param ... additional options (ignored)
#' @return NULL
#' @export
print.WDendrogram <- function(x, stand.alone=TRUE, layout.only=FALSE, cex=1, ...) {
  if (layout.only)
    return(.print.layout(x))
  pushViewport(viewport(x=unit(x$dm$left,'npc'), y=unit(x$dm$bottom,'npc'),
                        width=unit(x$dm$width,'npc'), height=unit(x$dm$height,'npc'),
                        just=c('left','bottom'), gp=gpar(col='black')))
  grid.dendrogram(as.dendrogram(x$clust), facing=x$facing)
  upViewport()
}

plot.WDendrogram <- print.WDendrogram

CalcTextBounding.WDendrogram <- function(dd, group) {
  dm <- DimToTop(dd, group)
  dm <- DimNPCToPoints(dm)
  dm
}

#' Draw dendrogram under grid system
#'
#' The dendrogram can be renderred. A viewport is created which contains the dendrogram.
#'
#' -order should leaves of dendrogram be put in the normal order (1, ..., n) or reverse order (n, ..., 1)?
#' -... pass to `grid::viewport` which contains the dendrogram.
#'
#'
#' This function only plots the dendrogram without adding labels. The leaves of the dendrogram
#' locates at unit(c(0.5, 1.5, ...(n-0.5))/n, "npc").
#'
#' @param dend a stats::dendrogram object.
#' @param facing facing of the dendrogram.
#' @param max_height maximum height of the dendrogram.
#' @param order order
#' @param ... additional options
#' @import grid
#' @source adapted from the ComplexHeatmap package authored by Zuguang Gu <z.gu@dkfz.de>
#' @export
grid.dendrogram = function(dend, facing = c("bottom", "top", "left", "right"),
                           max_height = NULL, order = c("normal", "reverse"), ...) {

  facing = match.arg(facing)[1]

  if(is.null(max_height)) {
    max_height = attr(dend, "height")
  }

  is.leaf = function(object) {
    leaf = attr(object, "leaf")
    if(is.null(leaf)) {
      FALSE
    } else {
      leaf
    }
  }

  draw.d = function(dend, max_height, facing = "bottom", order = "normal", max_width = 0, env = NULL) {
    leaf = attr(dend, "leaf")
    d1 = dend[[1]]  # child tree 1
    d2 = dend[[2]]  # child tree 2
    height = attr(dend, "height")
    midpoint = attr(dend, "midpoint")

    if(is.leaf(d1)) {
      x1 = x[as.character(attr(d1, "label"))]
    } else {
      x1 = attr(d1, "midpoint") + x[as.character(labels(d1))[1]]
    }
    y1 = attr(d1, "height")

    if(is.leaf(d2)) {
      x2 = x[as.character(attr(d2, "label"))]
    } else {
      x2 = attr(d2, "midpoint") + x[as.character(labels(d2))[1]]
    }
    y2 = attr(d2, "height")

    # graphic parameters for current branch
    edge_gp1 = as.list(attr(d1, "edgePar"))
    edge_gp2 = as.list(attr(d2, "edgePar"))

    if(is.null(env)) {
      begin = TRUE
      env = new.env()
      n = nobs(dend)
      env$x0 = NULL
      env$y0 = NULL
      env$x1 = NULL
      env$y1 = NULL
      env$col = NULL
      env$lty = NULL
      env$lwd = NULL
    } else {
      begin = FALSE
    }

    for(gp_name in c("col", "lwd", "lty")) {
      if(is.null(edge_gp1[[gp_name]])) {
        env[[gp_name]] = c(env[[gp_name]], rep(get.gpar(gp_name)[[gp_name]], 2))
      } else {
        env[[gp_name]] = c(env[[gp_name]], rep(edge_gp1[[gp_name]], 2))
      }
      if(is.null(edge_gp2[[gp_name]])) {
        env[[gp_name]] = c(env[[gp_name]], rep(get.gpar(gp_name)[[gp_name]], 2))
      } else {
        env[[gp_name]] = c(env[[gp_name]], rep(edge_gp2[[gp_name]], 2))
      }
    }


    ## plot the connection line
    if(order == "normal") {
      if(facing == "bottom") {
        env$x0 = c(env$x0, c(x1, x1, x2, x2))
        env$y0 = c(env$y0, c(y1, height, y2, height))
        env$x1 = c(env$x1, c(x1, (x1+x2)/2, x2, (x1+x2)/2))
        env$y1 = c(env$y1, c(height, height, height, height))
      } else if(facing == "top") {
        env$x0 = c(env$x0, c(x1, x1, x2, x2))
        env$y0 = c(env$y0, max_height - c(y1, height, y2, height))
        env$x1 = c(env$x1, c(x1, (x1+x2)/2, x2, (x1+x2)/2))
        env$y1 = c(env$y1, max_height - c(height, height, height, height))
      } else if(facing == "right") {
        env$x0 = c(env$x0, max_height - c(y1, height, y2, height))
        env$y0 = c(env$y0, max_width - c(x1, x1, x2, x2))
        env$x1 = c(env$x1, max_height - c(height, height, height, height))
        env$y1 = c(env$y1, max_width - c(x1, (x1+x2)/2, x2, (x1+x2)/2))
      } else if(facing == "left") {
        env$x0 = c(env$x0, c(y1, height, y2, height))
        env$y0 = c(env$y0, max_width - c(x1, x1, x2, x2))
        env$x1 = c(env$x1, c(height, height, height, height))
        env$y1 = c(env$y1, max_width - c(x1, (x1+x2)/2, x2, (x1+x2)/2))
      }
    } else {
      if(facing == "bottom") {
        env$x0 = c(env$x0, max_width - c(x1, x1, x2, x2))
        env$y0 = c(env$y0, c(y1, height, y2, height))
        env$x1 = c(env$x1, max_width - c(x1, (x1+x2)/2, x2, (x1+x2)/2))
        env$y1 = c(env$y1, c(height, height, height, height))
      } else if(facing == "top") {
        env$x0 = c(env$x0, max_width - c(x1, x1, x2, x2))
        env$y0 = c(env$y0, max_height - c(y1, height, y2, height))
        env$x1 = c(env$x1, max_width - c(x1, (x1+x2)/2, x2, (x1+x2)/2))
        env$y1 = c(env$y1, max_height - c(height, height, height, height))
      } else if(facing == "right") {
        env$x0 = c(env$x0, max_height - c(y1, height, y2, height))
        env$y0 = c(env$y0, c(x1, x1, x2, x2))
        env$x1 = c(env$x1, max_height - c(height, height, height, height))
        env$y1 = c(env$y1, c(x1, (x1+x2)/2, x2, (x1+x2)/2))
      } else if(facing == "left") {
        env$x0 = c(env$x0, c(y1, height, y2, height))
        env$y0 = c(env$y0, c(x1, x1, x2, x2))
        env$x1 = c(env$x1, c(height, height, height, height))
        env$y1 = c(env$y1, c(x1, (x1+x2)/2, x2, (x1+x2)/2))
      }
    }
    ## do it recursively
    if(!is.leaf(d1)) {
      draw.d(d1, max_height, facing, order, max_width, env = env)
    }
    if(!is.leaf(d2)) {
      draw.d(d2, max_height, facing, order, max_width, env = env)
    }

    if(begin) {
      grid.segments(env$x0, env$y0, env$x1, env$y1, default.units = "native", gp = gpar(col = env$col, lty = env$lty, lwd = env$lwd))
    }
  }

  labels = as.character(labels(dend))
  x = seq_along(labels) - 0.5 # leaves are places at x = 0.5, 1.5, ..., n - 0.5

  names(x) = labels
  n = length(labels)

  order = match.arg(order)[1]

  if(facing %in% c("top", "bottom")) {
    pushViewport(viewport(xscale = c(0, n), yscale = c(0, max_height), ...))
    draw.d(dend, max_height, facing, order, max_width = n)
    upViewport()
  } else if(facing %in% c("right", "left")) {
    pushViewport(viewport(yscale = c(0, n), xscale = c(0, max_height), ...))
    draw.d(dend, max_height, facing, order, max_width = n)
    upViewport()
  }
}
