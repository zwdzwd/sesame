#' Convert from absolute coordinates to affine coordinates
#'
#' @param dm dimension on the same coordinate system
#' as the affine system (absolute coordinates)
#' @param dm.sys dimension of the affine system
#' @return dimension on affine coordinates (relative coordinates)
ToAffine <- function(dm, dm.sys) {
  dm.affine <- dm
  dm.affine$left <- (dm$left - dm.sys$left) / dm.sys$width
  dm.affine$bottom <- (dm$bottom - dm.sys$bottom) / dm.sys$height
  dm.affine$width <- dm$width / dm.sys$width
  dm.affine$height <- dm$height / dm.sys$height
  dm.affine$text.x <- XToAffine(dm$text.x, dm.sys)
  dm.affine$text.y <- YToAffine(dm$text.y, dm.sys)
  dm.affine
}

#' Convert from affine coordinates to absolute coordinates
#'
#' @param dm.affine dimension on affine coordinates (relative coordinates)
#' @param dm.sys dimension of the affine system
#' @return dimension on the same coordinate system
FromAffine <- function(dm.affine, dm.sys) {
  dm <- dm.affine
  dm$left <- dm.sys$left + dm.affine$left * dm.sys$width
  dm$bottom <- dm.sys$bottom + dm.affine$bottom * dm.sys$height
  dm$width <- dm.sys$width * dm.affine$width
  dm$height <- dm.sys$height * dm.affine$height
  dm$text.x <- XFromAffine(dm$text.x, dm.sys)
  dm$text.y <- YFromAffine(dm$text.y, dm.sys)
  dm
}

XToAffine <- function(x, dm.sys) {
  (x - dm.sys$left) / dm.sys$width
}

YToAffine <- function(y, dm.sys) {
  (y - dm.sys$bottom) / dm.sys$height
}

XFromAffine <- function(x, dm.sys) {
  dm.sys$left + x * dm.sys$width
}

YFromAffine <- function(y, dm.sys) {
  dm.sys$bottom + y * dm.sys$height
}

#' Construct a WGroup
#'
#' @param ... plotting objects to be grouped
#' @param name name of the group
#' @param group.dm group dimension, by default use the dm of the merge of members
#' @param group.from.member group merged from member coordinates (require affine == FALSE),
#' the supplied group.dm is ignored
#' @param affine whether the group members are on affine coordinates already
#' @param nr number of rows
#' @param nc number of columns
#' @param mar a WMar object
#' @return a WGroup object
#' 
#' @export
WGroup <- function(..., name='', group.dm=NULL, group.from.member=FALSE, mar=WMar(), affine=FALSE, nr=NULL, nc=NULL) {
  ## row and column.split must be a set separately ??
  objs <- list(...)

  dms <- lapply(objs, function(o) o$dm)
  ## the grouped coordinates in the source system
  ## the bounding box coordinates before scaling
  group.dm.source <- do.call(.DimGroup, dms)
  if (is.null(nc))
    group.dm.source$nc <- max(sapply(objs, function(o) o$dm$nc))
  else
    group.dm.source$nc <- nc
    
  if (is.null(nr))
    group.dm.source$nr <- max(sapply(objs, function(o) o$dm$nr))
  else
    group.dm.source$nr <- nr

  ## by default, group.dm inherits nr and nc from group.dm.source
  if (is.null(group.dm)) {
    group.dm <- WDim(nr=group.dm.source$nr,nc=group.dm.source$nc)
  }
  
  ## convert to affine coordinates if not
  ## the group dm is the merge of the member
  ## dm before affine conversion
  if (!affine) {
    objs <- lapply(objs, function(obj) {
      obj$dm <- ToAffine(obj$dm, group.dm.source)
      obj
    })

    ## ignore supplied group.dm, set from the grouped members
    if (group.from.member)
      group.dm <- group.dm.source
  }

  force(name)
  force(group.dm)
  force(group.dm.source)
  force(mar)
  structure(function(group) {
    group.dm <- Resolve(group.dm, group, nr=group.dm.source$nr, nc=group.dm.source$nc)
    group.obj <- structure(list(
      children=objs,
      name=name,
      mar=mar,
      dm=group.dm), class=c('WGroup', 'WObject'))

    ## assign names if missing
    missing.inds <- which(sapply(objs, function(obj) obj$name==''))
    assigned.names <- GroupAssignNames(group.obj, length(missing.inds))
    lapply(seq_along(missing.inds), function(i) {
      group.obj$children[[missing.inds[i]]]$name <<- assigned.names[i]
    })

    stopifnot(GroupCheckNameUnique(group.obj))
    stopifnot(GroupCheckNameValid(group.obj))
    group.obj
  }, class=c('WGenerator', 'WObject'))
}

## a convenient wrapper for resolved WGroup
ResolvedWGroup <- function(...) {
  Resolve(WGroup(..., group.from.member=TRUE), NULL)
}

## calculate bounding box including texts
CalcTextBounding.WGroup <- function(group.obj, top.group=NULL) {
  if (is.null(top.group)) top.group <- group.obj
  group.dmb <- DimInPoints(group.obj$dm)
  dmb <- do.call(
    .DimGroup, lapply(group.obj$children,
                      function(obj) CalcTextBounding(obj, top.group)))
  ## dmb <- FromAffine(dmb, group.dmb)
  .DimGroup(dmb, group.dmb)
}

#' Add a plotting object to a group
#'
#' The object to be added are in the same coordinate system as the group.
#' @param group.obj WGroup object to be added to
#' @param new.obj plotting object to be added
#' @return a WGroup object where new.obj is added.
AddWGroup <- function(group.obj, new.obj) {

  if ('WAnnotate' %in% class(new.obj)) # WLabel does not contribute to dimensions
    dm <- group.obj$dm
  else {
    new.obj$dm <- Resolve(new.obj$dm, group.obj)
    dm <- .DimGroup(group.obj$dm, new.obj$dm)
  }
  dm$nc <- max(group.obj$dm$nc, new.obj$dm$nc)
  dm$nr <- max(group.obj$dm$nr, new.obj$dm$nr)

  ## put old and new children's dimensions
  ## to npc of the new dimension
  ## olds
  group.obj$children <- lapply(group.obj$children, function(obj) {
    obj$dm <- ToAffine(FromAffine(obj$dm, group.obj$dm), dm)
    obj
  })

  if (new.obj$name %in% GroupNames(group.obj)) {
    message('New object name ', new.obj$name, ' conflicts with existing names. Abort.')
    stop()
  }

  if (new.obj$name=='')
    new.obj$name <- GroupAssignNames(group.obj)

  ## new
  new.obj$dm <- ToAffine(new.obj$dm, dm)
  group.obj$dm <- dm
  group.obj$children[[length(group.obj$children)+1]] <- new.obj
  group.obj
}

#' Check whether group names are unique
#'
#' @param group.obj a WGroup
#' @return TRUE or FALSE
GroupCheckNameUnique <- function(group.obj) {
  if (!('WGroup' %in% class(group.obj)))
    return(TRUE)
  all.nms <- GroupNames(group.obj)
  if (length(all.nms)!=length(unique(all.nms)))
    return(FALSE)
  for(child in group.obj$children)
    if (!GroupCheckNameUnique(child))
      return(FALSE)
  return(TRUE)
}

GroupCheckNameValid <- function(group.obj) {
  all(GroupAllNames(group.obj)!='')
}

GroupNames <- function(group.obj) {
  sapply(group.obj$children, function(x) x$name)
}

GroupAllNames <- function(group.obj) {
  do.call(c, lapply(group.obj$children, function(x){
    if ('WGroup' %in% class(x))
      GroupAllNames(x)
    else
      x$name
  }))
}

GroupAssignNames <- function(group.obj, n=1) {
  i <- 0
  all.names <- GroupNames(group.obj)
  assigned <- NULL
  repeat{
    i <- i+1
    .name <- paste0('..internal.',i)
    if (!(.name %in% all.names) && n<=1) {
      assigned <- c(assigned, .name)
      n <- n-1
      break
    }
  }
  assigned
}

#' subset WGroup
#'
#' @param x a WGroup object
#' @param i integer indexing element
#' @export
`[.WGroup` <- function(x, i) {
  if(is.numeric(i))
    return(x$children[[i]])
  for (xx in x$children) {
    if (xx$name == i[1]) {
      if (length(i)>1 && 'WGroup' %in% class(xx))
        return(xx[i[2:length(i)]])
      return(xx)
    }
  }
  return(NULL)
}

#' Get an plotting object from a group's descendants
#'
#' @param x a WGroup object
#' @param nm name
#' @param force.unique assume the name is unique in the descendants and get one object instead of a list
#' @return if `force.unique==FALSE` return a list. Otherwise, one plotting object.
GroupDeepGet <- function(x, nm, force.unique=TRUE) {
  objs <- list()
  for (xx in x$children) {
    if (xx$name == nm)
      objs[[length(objs)+1]] <- xx
    if ('WGroup' %in% class(xx))
      objs <- c(objs, GroupDeepGet(xx, nm, force.unique=FALSE))
  }
  if (force.unique) {
    if (length(objs) > 1) {
      message('The name ',nm,' is ambiguous. Please provide full path.')
      stop()
    }
    return(objs[[1]])
  } else {
    return (objs)
  }
}

.GroupNameGet <- function(x, nm) {
  if (length(nm)==1) {
    if (!is.null(x[nm]))
      return(x[nm])
    else
      return(GroupDeepGet(x, nm))
  } else {
    return(x[nm])
  }
}

GroupNameGet <- function(x, nm) {
  obj <- .GroupNameGet(x, nm)
  if (is.null(obj)) {
    message('Object: ',nm,' unfound.')
    stop()
  }
  return(obj)
}

GetParentIn <- function(x, ancestor) {
  if (is.null(x$name)) # is root
    return(NULL)
  if ('WGroup' %in% class(ancestor)) {
    for (child in ancestor$children) {
      if (child$name==x$name)
        return(ancestor)
      if ('WGroup' %in% class(child)) {
        parent <- GetParentIn(x, child)
        if (!is.null(parent))
          return(parent)
      }
    }
  }
  return(NULL)
}

WFlatten <- function(.obs) {
  obs <- list()
  for(o in .obs){
    if ('WGroup' %in% class(o))
      obs <- c(obs, o$obs)
    else
      obs[[length(obs)+1]] <- o
  }
  obs
}

#' show layout
#'
#' @param x plot
#' @export
ly <- function(x) print(x, layout.only=TRUE)

#' Scale group
#'
#' Scale group to incorporate text on margins
#' @param group.obj group object that needs to be scaled
#' @return scaled group obj
ScaleGroup <- function(group.obj) {

  mar <- group.obj$mar

  dmb <- CalcTextBounding(group.obj)
  dmb$left <- dmb$left - dmb$width*mar$left
  dmb$bottom <- dmb$bottom - dmb$height*mar$bottom
  dmb$width <- dmb$width*(1+mar$left+mar$right)
  dmb$height <- dmb$height*(1+mar$bottom+mar$top)
  group.dmb <- DimInPoints(group.obj$dm)
  group.obj$dm <- ToAffine(group.dmb, dmb)
  cex <- c(group.dmb$width / dmb$width,
           group.dmb$height / dmb$height)

  list(cex=cex, group=group.obj)
}


.print.layout <- function(x) {
  pad <- 0.01
  pushViewport(viewport(x=unit(x$dm$left+x$dm$width*pad,'npc'),
                        y=unit(x$dm$bottom+x$dm$height*pad,'npc'),
                        width=unit(x$dm$width*(1-2*pad),'npc'),
                        height=unit(x$dm$height*(1-2*pad),'npc'),
                        just=c('left','bottom')))
  grid.rect(gp=gpar(col='red', lty='dashed', fill=NA))
  grid.text(x$name, gp=gpar(col='red'))
  return (upViewport())
}

#' Draw WGroup
#'
#' @param x a WGroup
#' @param cex factor for scaling fonts
#' @param layout.only to plot layout only
#' @param stand.alone to plot stand alone
#' @param ... additional options
#' @import grid
#' @export
print.WGroup <- function(x, stand.alone=TRUE, cex=1, layout.only=FALSE, ...) {

  if (stand.alone) {
    res <- ScaleGroup(x)
    cex <- res$cex
    x <- res$group
    grid.newpage()
  }

  pushViewport(viewport(x=unit(x$dm$left,'npc'),
                        y=unit(x$dm$bottom,'npc'),
                        width=unit(x$dm$width,'npc'),
                        height=unit(x$dm$height,'npc'),
                        just=c('left','bottom')))
  if (layout.only) {
    grid.rect(gp=gpar(col='green', lwd=5, alpha=0.6, fill=NULL))
    grid.text(x$name, gp=gpar(col='green', alpha=1))
  }

  for (child in x$children) {
    print(child, stand.alone=FALSE, cex=cex, layout.only=layout.only)
  }
  upViewport()
}

plot.WGroup <- print.WGroup

