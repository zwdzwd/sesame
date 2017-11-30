## colorspace.name:
## diverge_hcl, diverge_hsv, terrain_hcl,
## heat_hcl, sequential_hcl and rainbow_hcl

## brewer.name
## see display.brewer.all() for the brewer colors

#' Color Map Parameters
#'
#' Create color map parameters
#'
#' @param dmin minimum for continuous color map
#' @param dmax maximum for continuous color map
#' @param use.data use data as color, data must be either common color names
#' or hexdecimal color names
#' @param label2color a named vector or list that defines label to color mapping
#' explicitly for discrete color mapping
#' @param brewer.name palette name for RColorbrewer
#' @param brewer.n number of stop points in RColorBrewer for continuous color map
#' @param colorspace.name colorspace name
#' @param colorspace.n number of stops in colorspace palettes
#' @param cmap customized colormap name
#' @param stop.points custome stop points
#' @param grey.scale whether to use grey scale
#' @return an object of class CMPar
#' @export
CMPar <- function(dmin = NULL, dmax = NULL, # color scale max and min
                  brewer.name=NULL, brewer.n=3,
                  colorspace.name=NULL, colorspace.n=2,
                  cmap=NULL, label2color=NULL, use.data=FALSE,
                  stop.points=NULL, # color names at stop points
                  grey.scale=FALSE) {
  cmp <- lapply(formals(), eval)
  invisible(lapply(names(as.list(match.call()))[-1], function (nm) {
    cmp[[nm]] <<- get(nm)
  }))
  cmp
}

#' Constructor for ColoMap object
#'
#' Create color maps
#'
#' @param continuous whether colormap is continuous
#' @param colors colors for each data point
#' @param dmin miminum in continuous color map
#' @param dmax maximum in continuous color map
#' @param scaler scaler function from data range to 0-1
#' @param mapper function that maps data to color
#' @return an object of class ColorMap
#' @export
ColorMap <- function(continuous=TRUE,
                     colors=NULL,
                     dmin=NULL, dmax=NULL,
                     scaler=NULL, mapper=NULL) {
  cm <- lapply(formals(), eval)
  invisible(lapply(names(as.list(match.call()))[-1], function (nm) {
    cm[[nm]] <<- get(nm)
  }))
  class(cm) <- 'ColorMap'
  cm
}

#' map data to continuous color
#'
#' @param data numeric vector
#' @param cmp an color map parameter object of class CMPar
#' @param given.cm given colormap
#' @import colorspace
#' @import RColorBrewer
#' @importFrom grDevices rgb colorRamp
#' @importFrom grDevices col2rgb
#' @return an object of ColorMap
#' @export
MapToContinuousColors <- function(data, cmp=CMPar(), given.cm=NULL) {

  if (!is.null(given.cm)) {
    given.cm$colors <- apply(
      given.cm$mapper(given.cm$scaler(data)), 1,
      function(x) do.call(rgb, c(as.list(x), maxColorValue=255)))
    return(given.cm)
  }

  if (is.null(cmp$stop.points)) {
    if (is.null(cmp$cmap) &&
        is.null(cmp$brewer.name) &&
        is.null(cmp$colorspace.name)) {
      cmp$cmap <- 'jet'
    }
    if (!is.null(cmp$cmap)) {
      ## get(cmp$cmap)
      ## data(cmp$cmap)
      cmp$stop.points <- get(paste0(cmp$cmap,'.stops'))
    } else if (!is.null(cmp$brewer.name)) {
      ## use display.brewer.all for the brewer colors
      ## note that brewer.n cannot be >8 typically
      if (cmp$brewer.n < 3)
        cmp$brewer.n <- 3
      cmp$stop.points <- brewer.pal(cmp$brewer.n, cmp$brewer.name)
    } else {
      ## colorspace.name can be
      ## diverge_hcl, diverge_hsv, terrain_hcl, heat_hcl, sequential_hcl and rainbow_hcl
      ## colorspace.n can be very large
      cmp$stop.points <- get(cmp$colorspace.name)(cmp$colorspace.n)
    }
  }

  ## cap data
  if (!is.null(cmp$dmax))
    data[data>=cmp$dmax] <- cmp$dmax
  if (!is.null(cmp$dmin))
    data[data<=cmp$dmin] <- cmp$dmin

  .dmax <- max(cmp$dmax, data, na.rm=TRUE)
  .dmin <- min(cmp$dmin, data, na.rm=TRUE)
  if (.dmax==.dmin) # when range==0
    .dmax <- .dmax+1
  data <- (data - .dmin) / (.dmax-.dmin)

  cm <- ColorMap(
    dmin = .dmin, dmax = .dmax,
    scaler = function(x) {(x-.dmin)/(.dmax-.dmin)},
    mapper = colorRamp(cmp$stop.points, alpha=TRUE))
  cm$colors = apply(cm$mapper(data), 1, function(x) {
    if (any(is.na(x)))
      x <- col2rgb('#C0C0C0', alpha=TRUE)
    do.call(rgb, c(as.list(x), maxColorValue=255))
  })
  cm
}

#' map data to discrete color
#'
#' @param data numeric vector
#' @param cmp an color map parameter object of class CMPar
#' @param given.cm given color map
#' @return an object of ColorMap
#' @import RColorBrewer
#' @import colorspace
#' @importFrom grDevices rgb colorRamp
#' @export
MapToDiscreteColors <- function(data, cmp=CMPar(), given.cm=NULL) {

  if (!is.null(given.cm)) {
    given.cm$colors <- given.cm$mapper[as.character(data)]
    return(given.cm)
  }

  alphabet <- as.character(unique(as.vector(data)))
  if (is.null(cmp$cmap) &&
      is.null(cmp$brewer.name) &&
      is.null(cmp$colorspace.name)) {
    cmp$brewer.name <- 'Accent'
  }

  ## when label2color is explicitly set
  if (cmp$use.data) {
    mapped.colors <- alphabet
  } else if (!is.null(cmp$label2color)) {
    mapped.colors <- cmp$label2color[alphabet]
  } else if (!is.null(cmp$brewer.name) &&
      length(alphabet)<=brewer.pal.info[cmp$brewer.name,'maxcolors']) {
    ## use grey scale for binary and unary data
    if (length(alphabet)<3)
      mapped.colors <- c('#C0C0C0','#808080')[1:length(alphabet)]
    else
      mapped.colors <- brewer.pal(length(alphabet), cmp$brewer.name)
  } else if (!is.null(cmp$cmap)) {
    ## get(cmp$cmap)
    ## data(cmp$cmap)
    stop.points <- get(paste0(cmp$cmap,'.stops'))
    mapped.colors <- colorRamp(stop.points, alpha=TRUE)(length(alphabet))
  } else {
    if (is.null(cmp$colorspace.name))
      cmp$colorspace.name <- 'rainbow_hcl'
    mapped.colors <- get(cmp$colorspace.name)(length(alphabet))
  }

  cm <- ColorMap(
    continuous=FALSE,
    mapper=setNames(mapped.colors, alphabet))
  cm$colors=cm$mapper[as.character(data)]
  cm
}
