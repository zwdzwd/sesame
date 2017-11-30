#' Get dimensions
#' 
#' @param x WDim object or a plotting object
#' @export
getdim <- function(x) {
  if('WDim' %in% class(x))
    dm <- x
  else
    dm <- x$dm
  cat(dm$left, ' ', dm$bottom, ' ', dm$width, ' ', dm$height,'\n')
}

## #' format WDim
## #'
## #' @param object dimension
## #' @param ... additional parameters
## #' @return NULL, output WDim representation
## #' @export
## str.WDim <- function(object, ...) {
##   browser()
##   message("WDim object")
##   message(sprintf(" - l: %1.2f | b: %1.2f | w: %1.2f | h: %1.2f",
##                   object$left, object$bottom, object$width, object$height))
##   message(sprintf(' - nr: %d | nc: %d', object$nr, object$nc))
##   invisible(NULL)
## }

Resolve.WGenerator <- function(x, group) {
  x(group)
}

Resolve.WDim <- function(x, ...) x

Resolve.WHeatmap <- function(x, ...) x

Resolve.WDendrogram <- function(x, ...) x

Resolve.WGroup <- function(x, ...) x

Resolve.WObject <- function(x, ...) x

Resolve.WDimGenerator <- function(x, group, ...) {
  x(group, ...)
  ## if (is.null(x$dm)) {
  ##   message('Error dm is NULL. Abort.')
  ##   stop()
  ## }
  ## if ('WDim' %in% class(x$dm))
  ##   return(x)
  ## if ('function' %in% class(x$dm)) {
  ##   x$dm <- x$dm(x, nr, nc, group)
  ##   return(x)
  ## }
  ## message('Unknown dimension')
  ## stop()
}

Resolve.character <- function(x, group) {
  GroupNameGet(group, x)
}

## float length to top level
LengthToTop <- function(obj, root, .length) {
  parent <- GetParentIn(obj, root)
  if (is.null(parent)) {
    return(.length)
  }
  .length <- .length * max(obj$dm$width, obj$dm$height) # take max here, may be optimized.
  return(LengthToTop(parent, root, .length))
}

## x is inside obj's dimension
XToTop <- function(obj, root, x) {
  parent <- GetParentIn(obj, root)
  x <- XFromAffine(x, obj$dm)
  if (is.null(parent)) {
    return(x)
  }
  return(XToTop(parent, root, x))
}

## y is inside obj's dimension
YToTop <- function(obj, root, y) {
  parent <- GetParentIn(obj, root)
  y <- YFromAffine(y, obj$dm)
  if (is.null(parent)) {
    return(y)
  }
  return(YToTop(parent, root, y))
}

## float dimension to top level
DimToTop <- function(obj, root, dm=NULL) {
  if (is.null(dm)) {
    dm <- obj$dm
  }
  parent <- GetParentIn(obj, root)
  if (is.null(parent)) {
    return(dm)
  }
  dm <- FromAffine(dm, parent$dm)
  return(DimToTop(parent, root, dm))
}

## describe a dimension in points
DimInPoints <- function(dm) {
  dm.new <- dm
  dm.new$left <- NPCToPoints(dm$left)
  dm.new$bottom <- NPCToPoints(dm$bottom)
  dm.new$width <- NPCToPoints(dm$width)
  dm.new$height <- NPCToPoints(dm$height)
  dm.new
}

#' class WDim
#'
#' class WDim
#'
#' @param left left coordinate
#' @param bottom bottom coordinate
#' @param width width
#' @param height height
#' @param nr number of row
#' @param nc number of column
#' @param text.x x anchor for text
#' @param text.y y anchor for text
#' @param text.just just for text
#' @param column.split a list of WDim objects for column split
#' @param row.split a list of WDim objects for row split
#' @export
WDim <- function(left=0, bottom=0, width=1, height=1, nr=1, nc=1,
                 text.x=0, text.y=0, text.just=c('center','center'),
                 column.split=NULL, row.split=NULL) {
  dm <- list(left=left, bottom=bottom, width=width, height=height, nr=nr, nc=nc,
             column.split=column.split, row.split=row.split,
             text.x=0, text.y=0, text.just=text.just)
  class(dm) <- 'WDim'
  dm
}

.DimRight <- function(dm) {
  dm$left+dm$width
}

.DimTop <- function(dm) {
  dm$bottom+dm$height
}

## Draw height
##
## Height without inter-subdimension space
.DimDrawHeight <- function(dm) {
  if (is.null(dm$row.split)) return(dm$height)
  else return(sum(sapply(dm$row.split, .DimDrawHeight)))
}

## Draw width
##
## Width without inter-subdimension space
.DimDrawWidth <- function(dm) {
  if (is.null(dm$column.split)) return(dm$width)
  else return(sum(sapply(dm$column.split, .DimDrawWidth)))
}

## group dimensions
.DimGroup <- function(...) {
  dms <- list(...)
  left = min(sapply(dms, function(dm) dm$left))
  bottom = min(sapply(dms, function(dm) dm$bottom))
  width = max(sapply(dms, function(dm) dm$left+dm$width)) - left
  height = max(sapply(dms, function(dm) dm$bottom+dm$height)) - bottom
  WDim(left=left, bottom=bottom, width=width, height=height)
}


ResolveToTopDim <- function(x, group) {
  if (is.null(x))
    return(NULL)
  else
    return(DimToTop(Resolve(x, group), group))
}

#' Top left of
#'
#' Place a new object to the top left corner of another.
#' @param x target object, either a name, a object or NULL which refers to the last plotting object
#' @param just the part from the new object that should be attached to
#' @param v.pad vertical translational padding [0.0]
#' @param h.pad horizontal translational padding [0.0]
#' @return a WDimGenerator
#' @export
TopLeftOf <- function(x=NULL,
                      just=c('right','bottom'),
                      v.pad=0.0, h.pad=0.0) {

  stopifnot(just[1] %in% c('left','center','right'))
  stopifnot(just[2] %in% c('bottom','center','top'))
  force(x); force(just); force(v.pad); force(h.pad);
  structure(function(group, nr=1, nc=1, hard.dm=NULL) {
    if (is.null(x))
      x <- group$children[[length(group$children)]]$name
    x <- ResolveToTopDim(x, group)

    if (is.null(hard.dm)) {
      dm <- WDim()
      dm$width <- x$width/x$nc*nc
      dm$height <- x$height/x$nr*nr
      dm$nc <- nc
      dm$nr <- nr
    } else {
      dm <- hard.dm
    }

    if (just[1]=='right') {
      dm$left <- x$left - dm$width
    } else if (just[1]=='center') {
      dm$left <- x$left - dm$width/2
    } else {
      dm$left <- x$left
    }
    dm$left <- dm$left + h.pad

    if (just[2]=='bottom') {
      dm$bottom <- x$bottom + x$height
    } else if (just[2]=='center') {
      dm$bottom <- x$bottom + x$height - dm$height/2
    } else {
      dm$bottom <- x$bottom + x$height - dm$height
    }
    dm$bottom <- dm$bottom + v.pad

    dm$text.x <- x$left + h.pad
    dm$text.y <- x$bottom + x$height + v.pad
    dm$text.just <- just

    dm
  }, class='WDimGenerator')
}

#' Top right of
#'
#' Place a new object to the top right corner of another.
#' @param x target object, either a name, a object or NULL which refers to the last plotting object
#' @param just the part from the new object that should be attached to
#' @param v.pad vertical translational padding [0.0]
#' @param h.pad horizontal translational padding [0.0]
#' @return a WDimGenerator
#' @export
TopRightOf <- function(x=NULL,
                       just=c('left','bottom'),
                       v.pad=0.0, h.pad=0.0) {

  stopifnot(just[1] %in% c('left','center','right'))
  stopifnot(just[2] %in% c('bottom','center','top'))
  force(x); force(just); force(v.pad); force(h.pad);
  structure(function(group, nr=1, nc=1, hard.dm=NULL) {
    if (is.null(x))
      x <- group$children[[length(group$children)]]$name
    x <- ResolveToTopDim(x, group)

    if (is.null(hard.dm)) {
      dm <- WDim()
      dm$width <- x$width/x$nc*nc
      dm$height <- x$height/x$nr*nr
      dm$nc <- nc
      dm$nr <- nr
    } else {
      dm <- hard.dm
    }

    if (just[1]=='right') {
      dm$left <- x$left + x$width - dm$width
    } else if (just[1]=='center') {
      dm$left <- x$left + x$width - dm$width/2
    } else {
      dm$left <- x$left + x$width
    }
    dm$left <- dm$left + h.pad

    if (just[2]=='bottom') {
      dm$bottom <- x$bottom + x$height
    } else if (just[2]=='center') {
      dm$bottom <- x$bottom + x$height - dm$height/2
    } else {
      dm$bottom <- x$bottom + x$height - dm$height
    }
    dm$bottom <- dm$bottom + v.pad

    dm$text.x <- x$left + x$width + h.pad
    dm$text.y <- x$bottom + x$height + v.pad
    dm$text.just <- just

    dm
  }, class='WDimGenerator')
}

#' Bottom left of
#'
#' Place a new object to the bottom left corner of another.
#' @param x target object, either a name, a object or NULL which refers to the last plotting object
#' @param just the part from the new object that should be attached to
#' @param v.pad vertical translational padding [0.0]
#' @param h.pad horizontal translational padding [0.0]
#' @return a WDimGenerator
#' @export
BottomLeftOf <- function(x=NULL,
                         just=c('right', 'bottom'),
                         v.pad=0.0, h.pad=0.0) {

  stopifnot(just[1] %in% c('left','center','right'))
  stopifnot(just[2] %in% c('bottom','center','top'))
  force(x); force(just); force(v.pad); force(h.pad);
  structure(function(group, nr=1, nc=1, hard.dm=NULL) {
    if (is.null(x))
      x <- group$children[[length(group$children)]]$name
    x <- ResolveToTopDim(x, group)

    if (is.null(hard.dm)) {
      dm <- WDim()
      dm$width <- x$width/x$nc*nc
      dm$height <- x$height/x$nr*nr
      dm$nc <- nc
      dm$nr <- nr
    } else {
      dm <- hard.dm
    }

    if (just[1]=='right') {
      dm$left <- x$left - dm$width
    } else if (just[1]=='center') {
      dm$left <- x$left - dm$width/2
    } else {
      dm$left <- x$left
    }
    dm$left <- dm$left + h.pad

    if (just[2]=='bottom') {
      dm$bottom <- x$bottom
    } else if (just[2]=='center') {
      dm$bottom <- x$bottom - dm$height/2
    } else {
      dm$bottom <- x$bottom - dm$height
    }
    dm$bottom <- dm$bottom + v.pad

    dm$text.x <- x$left + h.pad
    dm$text.y <- x$bottom + v.pad
    dm$text.just <- just

    dm
  }, class='WDimGenerator')
}

#' Bottom right of
#'
#' Place a new object to the bottom right corner of another.
#' @param x target object, either a name, a object or NULL which refers to the last plotting object
#' @param just the part from the new object that should be attached to
#' @param v.pad vertical translational padding [0.0]
#' @param h.pad horizontal translational padding [0.0]
#' @return a WDimGenerator
#' @export
BottomRightOf <- function(x=NULL,
                          just=c('left','bottom'),
                          v.pad=0.0, h.pad=0.0) {

  stopifnot(just[1] %in% c('left','center','right'))
  stopifnot(just[2] %in% c('bottom','center','top'))
  force(x); force(just); force(v.pad); force(h.pad);
  structure(function(group, nr=1, nc=1, hard.dm=NULL) {
    if (is.null(x))
      x <- group$children[[length(group$children)]]$name
    x <- ResolveToTopDim(x, group)

    if (is.null(hard.dm)) {
      dm <- WDim()
      dm$width <- x$width/x$nc*nc
      dm$height <- x$height/x$nr*nr
      dm$nc <- nc
      dm$nr <- nr
    } else {
      dm <- hard.dm
    }

    if (just[1]=='right') {
      dm$left <- x$left + x$width - dm$width
    } else if (just[1]=='center') {
      dm$left <- x$left + x$width - dm$width/2
    } else {
      dm$left <- x$left + x$width
    }
    dm$left <- dm$left + h.pad

    if (just[2]=='bottom') {
      dm$bottom <- x$bottom
    } else if (just[2]=='center') {
      dm$bottom <- x$bottom - dm$height/2
    } else {
      dm$bottom <- x$bottom - dm$height
    }
    dm$bottom <- dm$bottom + v.pad

    dm$text.x <- x$left + x$width + h.pad
    dm$text.y <- x$bottom + v.pad
    dm$text.just <- just

    dm
  }, class='WDimGenerator')
}


#' Top of
#'
#' Generate dimension top of another object
#'
#' @param x an object with dimension
#' @param height the height of the new object (when NULL, set to proportional to data)
#' @param pad padding between the target and current
#' @param min.ratio minimum ratio of dimensions when auto-scale
#' @param h.aln object for horizontal alignment (when NULL, set to x)
#' @param v.scale object for vertical scaling (when NULL, set to x)
#' @param v.scale.proportional when v.scale is provided, whether to make proportional to data
#' @return a dimension generator on top of x
#' @export
TopOf <- function(x=NULL, height=NULL, pad=0.01, min.ratio=0.02,
                  h.aln=NULL, v.scale=NULL, v.scale.proportional=FALSE) {

  force(x); force(h.aln); force(v.scale);
  force(v.scale.proportional)
  force(height); force(pad); force(min.ratio);
  structure(function(group, nr=1, nc=1, hard.dm=NULL) {
    if (is.null(x))
      x <- group$children[[length(group$children)]]$name
    x <- ResolveToTopDim(x, group)
    h.aln <- ResolveToTopDim(h.aln, group)
    v.scale <- ResolveToTopDim(v.scale, group)

    dm <- x
    dm$nr <- nr
    dm$nc <- nc

    ## vertical alignment
    if (is.null(height)) {
      if (is.null(v.scale)) {
        v.scale <- x
        v.scale.proportional <- TRUE
      } else {
        dm$row.split <- v.scale$row.split
      }
      if (v.scale.proportional) {
        .height <- 1 / v.scale$nr * nr
        .height <- max(min.ratio, .height)
        .height <- min(1/min.ratio, .height)
        dm$height <- .height * .DimDrawHeight(v.scale)
      } else {
        dm$height <- v.scale$height
      }
    } else {
      dm$height <- height
    }
    dm$bottom <- x$bottom+pad+x$height

    ## horizontal alignment
    if (is.null(h.aln)) h.aln <- x
    dm$left <- h.aln$left
    dm$width <- h.aln$width
    dm$column.split <- h.aln$column.split

    dm$text.x <- x$left + x$width/2
    dm$text.y <- x$bottom + x$height + pad
    dm$text.just <- c('center','bottom')

    dm
  }, class='WDimGenerator')
}

#' Beneath
#'
#' Generate dimension beneath another object
#'
#' @param x an object with dimension
#' @param height the height of the new object (when NULL set proportional to the data)
#' @param pad padding between the target and current
#' @param min.ratio minimum ratio of dimensions when auto-scale
#' @param h.aln object for horizontal alignment (when NULL, set to x)
#' @param v.scale object for vertical scaling (when NULL, set to x)
#' @param v.scale.proportional when v.scale is provided, whether to make proportional to data
#' @return a dimension generator beneath x
#' @export
Beneath <- function(x=NULL, height=NULL, pad=0.01, min.ratio=0.02,
                    h.aln=NULL, v.scale=NULL, v.scale.proportional=FALSE) {

  force(x); force(h.aln); force(v.scale);
  force(v.scale.proportional)
  force(height); force(pad); force(min.ratio);
  structure(function(group, nr=1, nc=1, hard.dm=NULL) {
    if (is.null(x))
      x <- group$children[[length(group$children)]]$name
    x <- ResolveToTopDim(x, group)
    h.aln <- ResolveToTopDim(h.aln, group)
    v.scale <- ResolveToTopDim(v.scale, group)

    dm <- x
    dm$nr <- nr
    dm$nc <- nc

    ## vertical alignment
    if (is.null(height)) {
      if (is.null(v.scale)) {
        v.scale <- x
        v.scale.proportional <- TRUE
      } else {
        dm$row.split <- v.scale$row.split
      }
      if (v.scale.proportional) {
        .height <- 1 / v.scale$nr * nr
        .height <- max(min.ratio, .height)
        .height <- min(1/min.ratio, .height)
        dm$height <- .height * .DimDrawHeight(v.scale)
      } else {
        dm$height <- v.scale$height
      }
    } else {
      dm$height <- height
    }
    dm$bottom <- x$bottom-pad-dm$height

    ## horizontal alignment
    if (is.null(h.aln)) h.aln <- x
    dm$left <- h.aln$left
    dm$width <- h.aln$width
    dm$column.split <- h.aln$column.split

    dm$text.x <- x$left + x$width/2
    dm$text.y <- x$bottom - pad
    dm$text.just <- c('center','top')

    dm
  }, class='WDimGenerator')
}

#' LeftOf
#'
#' Generate dimension to the left of another object
#'
#' @param x an object with dimension
#' @param width the width of the new object (when NULL, set proportional to data)
#' @param pad padding between the target and current
#' @param min.ratio minimum ratio of dimensions when auto-scale
#' @param v.aln object for vertical alignment (when NULL, set to x)
#' @param h.scale object for horizontal scaling (when NULL, set to x)
#' @param h.scale.proportional when h.scale is provided, whether to make proportional to data
#' @return a dimension to the left of x
#' @export
LeftOf <- function(x=NULL, width=NULL, pad=0.01, min.ratio=0.02,
                   v.aln=NULL, h.scale=NULL, h.scale.proportional=FALSE) {

  force(x); force(v.aln); force(h.scale);
  force(h.scale.proportional)
  force(width); force(pad); force(min.ratio);
  structure(function(group, nr=1, nc=1, hard.dm=NULL) {
    if (is.null(x))
      x <- group$children[[length(group$children)]]$name
    x <- ResolveToTopDim(x, group)
    v.aln <- ResolveToTopDim(v.aln, group)
    h.scale <- ResolveToTopDim(h.scale, group)

    dm <- x
    dm$nr <- nr
    dm$nc <- nc

    ## horizontal alignment
    if (is.null(width)) {
      if (is.null(h.scale)) {
        h.scale <- x
        h.scale.proportional <- TRUE
      } else {
        dm$column.split <- h.scale$column.split
      }
      if (h.scale.proportional) {
        .width <- 1 / h.scale$nc * nc
        .width <- max(min.ratio, .width)
        .width <- min(1/min.ratio, .width)
        dm$width <- .width * .DimDrawWidth(h.scale)
      } else {
        dm$wdith <- h.scale$width
      }
    } else {
      dm$width <- width
    }
    dm$left <- x$left-pad-dm$width

    ## vertical alignment
    if (is.null(v.aln)) v.aln <- x
    dm$bottom <- v.aln$bottom
    dm$height <- v.aln$height
    dm$row.split <- v.aln$row.split

    dm$text.x <- x$left - pad
    dm$text.y <- x$bottom + x$height/2
    dm$text.just <- c('right','center')

    dm
  }, class='WDimGenerator')
}

#' RightOf
#'
#' Generate dimension to the right of another object
#'
#' @param x an object with dimension
#' @param width the width of the new object (when NULL, set proportional to data)
#' @param pad padding between the target and current
#' @param min.ratio minimum ratio of dimensions when auto-scale
#' @param v.aln object for vertical alignment (when NULL, set to x)
#' @param h.scale object for horizontal scaling (when NULL, set to x)
#' @param h.scale.proportional when h.scale is provided, whether to make proportional to data
#' @return a dimension to the right of x
#' @export
RightOf <- function(x=NULL, width=NULL, pad=0.01, min.ratio=0.02, 
                    v.aln=NULL, h.scale=NULL, h.scale.proportional=FALSE) {

  force(x); force(v.aln); force(h.scale);
  force(h.scale.proportional);
  force(width); force(pad); force(min.ratio);
  structure(function(group, nr=1, nc=1, hard.dm=NULL) {
    if (is.null(x))
      x <- group$children[[length(group$children)]]$name
    x <- ResolveToTopDim(x, group)
    v.aln <- ResolveToTopDim(v.aln, group)
    h.scale <- ResolveToTopDim(h.scale, group)

    dm <- x
    dm$nr <- nr
    dm$nc <- nc

    ## horizontal alignment
    if (is.null(width)) {
      if (is.null(h.scale)) {
        h.scale <- x
        h.scale.proportional <- TRUE
      } else {
        dm$column.split <- h.scale$column.split
      }
      if (h.scale.proportional) {
        .width <- 1 / h.scale$nc * nc
        .width <- max(min.ratio, .width)
        .width <- min(1/min.ratio, .width)
        dm$width <- .width * .DimDrawWidth(h.scale)
      } else {
        dm$width <- h.scale$width
      }
    } else {
      dm$width <- width
    }
    dm$left <- dm$left+pad+x$width

    ## vertical alignment
    if (is.null(v.aln)) v.aln <- x
    dm$bottom <- v.aln$bottom
    dm$height <- v.aln$height
    dm$row.split <- v.aln$row.split

    dm$text.x <- x$left + x$width + pad
    dm$text.y <- x$bottom + x$height/2
    dm$text.just <- c('left','center')

    dm
  }, class='WDimGenerator')
}

#' place an arbitrary position w.r.t a subplot
#'
#' @param anchor.x x coordinates
#' @param anchor.y y coordinates
#' @param x plotting object to anchor
#' @param just adjustment of new plot
#' @param data.coord whether the coordinates is in term of data
#' @return a WDimGenerator object
#' @export
WPosition <- function(anchor.x, anchor.y, x=NULL, just=c('left','bottom'), data.coord=FALSE) {
  stopifnot(just[1] %in% c('left','center','right'))
  stopifnot(just[2] %in% c('bottom','center','top'))
  force(anchor.x); force(anchor.y); force(x); force(just); force(data.coord);
  structure(function(group, nr=1, nc=1, hard.dm=NULL) {
    if (is.null(x))
      x <- group$children[[length(group$children)]]$name

    x <- Resolve(x, group)

    ## resolve anchor positions to top level
    if (data.coord) {
      anchor.x <- (anchor.x - 0.5) / x$dm$nc
      anchor.y <- (x$dm$nr - anchor.y + 0.5) / x$dm$nr
    }

    anchor.x <- XToTop(x, group, anchor.x)
    anchor.y <- YToTop(x, group, anchor.y)

    x <- DimToTop(x, group)

    if (is.null(hard.dm)) {
      dm <- WDim()
      dm$width <- x$width/x$nc*nc
      dm$height <- x$height/x$nr*nr
      dm$nc <- nc
      dm$nr <- nr
    } else {
      dm <- hard.dm
    }

    if (just[1]=='left') {
      dm$left <- anchor.x
    } else if (just[1]=='center') {
      dm$left <- anchor.x - dm$width/2
    } else {
      dm$left <- anchor.x - dm$width
    }

    if (just[2]=='bottom') {
      dm$bottom <- anchor.y
    } else if (just[2]=='center') {
      dm$bottom <- anchor.y - dm$height/2
    } else {
      dm$bottom <- anchor.y - dm$height
    }

    dm$text.x <- anchor.x
    dm$text.y <- anchor.y
    dm$text.just <- just

    dm
  }, class='WDimGenerator')
}
