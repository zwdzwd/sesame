#' WHeatmap object
#'
#' Create a heatmap
#'
#' @param data data matrix
#' @param dm plotting dimension (a WDim or a WDimGenerator object)
#' @param name name of the plot
#' @param continuous whether the data should be treated as continuous or discrete
#' @param cmp a CMPar object, for tunning color mapping parameters
#' @param cm a given color map
#' @param xticklabels to plot xtick labels, one may supply characters to plot just a subset of xtick labels
#' @param xticklabels.n number of xtick labels to plot (resample for aethetics by default)
#' @param xticklabel.side xticklabel side (t or b)
#' @param xticklabel.fontsize xticklabel font size
#' @param xticklabel.pad padding between xticklabel and x-axis
#' @param xticklabel.rotat xticklabel rotation
#' @param yticklabels to plot ytick labels, one may supply characters to plot just a subset of ytick labels
#' @param yticklabels.n number of ytick labels to plot (resample for aethetics by default)
#' @param yticklabel.side yticklabel side (l or r)
#' @param yticklabel.fontsize yticklabel font size
#' @param yticklabel.pad padding between yticklabel and y-axis
#' @param yticklabel.rotat yticklabel rotation
#' @param gp a list of graphical parameters
#' @param sub.name subclass name
#' @return one or a list of heatmaps (depends on whether dimension is split)
#' @export
WHeatmap <- function(data=NULL, dm=NULL, name='', continuous=NULL,
                     cmp = CMPar(), # colormapping parameters
                     cm = NULL,

                     ## tick label on x-axis
                     xticklabels = NULL,
                     xticklabels.n = NULL,
                     xticklabel.side = 'b',
                     xticklabel.fontsize = 12,
                     xticklabel.rotat = 90,
                     xticklabel.pad = 0.005,

                     ## tick label on y-axis
                     yticklabels = NULL,
                     yticklabels.n = NULL,
                     yticklabel.side = 'l',
                     yticklabel.fontsize = 12,
                     yticklabel.rotat = 0,
                     yticklabel.pad = 0.005,

                     ## subclass name
                     sub.name = NULL,

                     ## graph parameters
                     gp = NULL) {

  if(!('matrix' %in% class(data))) {
    data <- tryCatch({
      as.matrix(data)
    }, error = function(e) {
      message('data argument must be matrix-like. Abort.')
      stop()
    })
  }

  hm <- lapply(formals(), eval)
  invisible(lapply(names(as.list(match.call()))[-1], function (nm) {
    hm[[nm]] <<- get(nm)
  }))

  if (is.null(hm$dm))
    hm$dm <- WDim(0,0,1,1,nr=nrow(data), nc=ncol(data))

  ## auto-infer continuous/discrete
  if (is.null(continuous)) {
    if (!is.null(cm))
      hm$continuous <- cm$continuous
    else if (!is.numeric(data) || length(unique(data)) < 5)
      hm$continuous <- FALSE
    else
      hm$continuous <- TRUE
  }

  ## graph parameters
  hm$gp <- list()
  hm$gp$col <- 'white'
  hm$gp$lty <- 'blank'
  lapply(names(gp), function(x) {hm$gp[[x]] <<- gp[[x]]})

  ## map to colors
  if (hm$continuous)
    hm$cm <- MapToContinuousColors(hm$data, cmp=hm$cmp, given.cm=cm)
  else
    hm$cm <- MapToDiscreteColors(hm$data, cmp=hm$cmp, given.cm=cm)

  class(hm) <- c('WHeatmap', 'WObject')
  if (!is.null(sub.name))
    class(hm) <- c(sub.name, class(hm))

  force(hm);
  structure(function(group) {
    hm$dm <- Resolve(hm$dm, group, nr=nrow(hm$data), nc=ncol(hm$data))
    ## split if dimension indicates so
    if (!is.null(hm$dm$column.split) || !is.null(hm$dm$row.split)) {
      return(SplitWHeatmap(hm, hm$dm, cm, group))
    }
    hm
  }, class=c('WGenerator', 'WObject'))
}

SplitWHeatmap <- function(hm, dm, cm, group) {

  if (is.null(dm$column.split)) {
    column.split <- list(dm)
    column.split[[1]]$left <- 0
    column.split[[1]]$width <- 1
  } else {
    column.split <- dm$column.split
  }
  if (is.null(dm$row.split)) {
    row.split <- list(dm)
    row.split[[1]]$bottom <- 0
    row.split[[1]]$height <- 1
  } else {
    row.split <- dm$row.split
  }

  all.nc <- sapply(column.split, function(dm) dm$nc)
  all.nr <- sapply(row.split, function(dm) dm$nr)
  sum.nc <- sum(all.nc)
  sum.nr <- sum(all.nr)
  nc.data <- ncol(hm$data)
  nr.data <- nrow(hm$data)
  col.inds <- c(0,round(cumsum(all.nc) * nc.data / sum.nc))
  row.inds <- c(0,round(cumsum(all.nr) * nr.data / sum.nr))

  sub.dms.col <- column.split[order(sapply(column.split, function(dm) dm$left))]
  sub.dms.row <- rev(row.split[order(sapply(row.split, function(dm) dm$bottom))])
  sub.dms <- expand.grid(seq_along(sub.dms.row), seq_along(sub.dms.col))
  k <- apply(sub.dms, 1, function(dm.i) {
    ir <- dm.i[1]
    ic <- dm.i[2]
    sub.dm.row <- sub.dms.row[[ir]]
    sub.dm.col <- sub.dms.col[[ic]]
    sub.hm <- hm
    sub.hm$dm <- WDim(sub.dm.col$left, sub.dm.row$bottom, sub.dm.col$width, sub.dm.row$height,
                      row.split=sub.dm.row$row.split, column.split=sub.dm.col$column.split,
                      nr=sub.dm.row$nr, nc=sub.dm.col$nc)
    sub.hm$data <- hm$data[(row.inds[ir]+1):row.inds[ir+1],
                           (col.inds[ic]+1):col.inds[ic+1], drop=FALSE]
    sub.hm$cmp$cm <- cm
    sub.hm$name <- paste0(hm$name, '.', ir, '.', ic)
    do.call(WHeatmap, sub.hm)(group)
  })
  k$name <- hm$name
  k$group.dm <- dm
  k$affine <- TRUE
  w.group <- do.call(ResolvedWGroup, k)
  # w.group$dm$row.split <- sub.dms.row
  # w.group$dm$column.split <- sub.dms.col

  return(w.group)
}

#' Calculate Texting Bounding for WHeatmap
#' 
#' @param hm an object of class WHeatmap
#' @param group an object of class WGroup
#' @return an object of class WDim in coordinate points
#' @export
CalcTextBounding.WHeatmap <- function(hm, group) {

  ## this needs be called at the ROOT view port
  dm <- DimToTop(hm, group)
  ## bottom, left, top, right
  left <- NPCToPoints(dm$left)
  bottom <- NPCToPoints(dm$bottom)
  top <- bottom + NPCToPoints(dm$height)
  right <- left + NPCToPoints(dm$width)

  if (!(is.null(hm[['yticklabels']]) ||
          (length(hm[['yticklabels']])==1 &&
             hm[['yticklabels']] == FALSE))) {
    if (hm$yticklabel.side=='l') {
      if (is.null(rownames(hm$data))) {
        .text.margin <- 0
      } else {
        .text.margin <- max(sapply(
          rownames(hm$data), function(t) text.width(t, hm$yticklabel.fontsize))) -
            NPCToPoints(LengthToTop(hm, group, hm$yticklabel.pad))
      }
      left <- left - .text.margin
    } else {
      if (is.null(rownames(hm$data))) {
        .text.margin <- 0
      } else {
        .text.margin <- max(sapply(
          rownames(hm$data), function(t) text.width(t, hm$yticklabel.fontsize))) +
            NPCToPoints(LengthToTop(hm, group, hm$yticklabel.pad))
      }
      right <- right + .text.margin
    }
  }

  #  if (!is.null(hm$title)) {
  #    if (hm$title.side=='l') {
  #      rg$left <- rg$left - text.width(hm$title, hm$title.fontsize) - hm$title.pad
  #    } else {
  #      rg$right <- rg$right + text.width(hm$title, hm$title.fontsize) + hm$title.pad
  #    }
  #  }

  if (!(is.null(hm[['xticklabels']]) ||
          (length(hm[['xticklabels']])==1 &&
             hm[['xticklabels']] == FALSE))) {
    if (hm$xticklabel.side=='b') {
      if (is.null(colnames(hm$data))) {
        .text.margin <- 0
      } else {
        .text.margin <- max(sapply(
          colnames(hm$data), function(t) text.width(t, hm$xticklabel.fontsize))) -
            NPCToPoints(LengthToTop(hm, group, hm$xticklabel.pad))
      }
      bottom <- bottom - .text.margin
    } else {
      if (is.null(colnames(hm$data))) {
        .text.margin <- 0
      } else {
        .text.margin <- max(sapply(
          colnames(hm$data), function(t) text.width(t, hm$xticklabel.fontsize))) +
            NPCToPoints(LengthToTop(hm, group, hm$xticklabel.pad))
      }
      top <- top + .text.margin
    }
  }
  dm$left <- left
  dm$bottom <- bottom
  dm$width <- right-left
  dm$height <- top-bottom
  dm
}

#' plot WHeatmap
#'
#' @param x a WHeatmap
#' @param stand.alone plot is stand alone
#' @param layout.only plot layout only
#' @param cex factor to scaling texts
#' @param ... additional options
#' @return \code{NULL}
#' @import grid
#' @export
print.WHeatmap <- function(x, cex=1, layout.only=FALSE, stand.alone=TRUE, ...) {

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

  nc = ncol(x$data)
  nr = nrow(x$data)
  xc = (seq_len(nc)-1)/nc
  yc = (rev(seq_len(nr))-1)/nr
  expand.index <- expand.grid(seq_len(nr), seq_len(nc))
  grid.rect(xc[expand.index[[2]]], yc[expand.index[[1]]],
            width=unit(1/nc, 'npc'), height=unit(1/nr, 'npc'),
            gp=do.call('gpar', c(list(fill=x$cm$colors), x$gp)), just=c('left','bottom'))

  ## x tick labels
  if (!is.null(x[['xticklabels']])) {
    .WPrintXTickLabels(x, x[['xticklabels']], cex=max(cex))
  }

  ## y tick labels
  if (!is.null(x[['yticklabels']])) {
    .WPrintYTickLabels(x, x[['yticklabels']], cex=max(cex))
  }

  upViewport()
}

plot.WHeatmap <- print.WHeatmap

