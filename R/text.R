#' @import grid
NPCToPoints <- function(npc) {
  as.numeric(convertUnit(unit(npc,'npc'),'points'))
}

DimNPCToPoints <- function(dm) {
  dm$left <- NPCToPoints(dm$left)
  dm$bottom <- NPCToPoints(dm$bottom)
  dm$height <- NPCToPoints(dm$height)
  dm$width <- NPCToPoints(dm$width)
  dm
}

## TODO: the following doesn't work for polymath calls

## text width and scale to specified font size
##
## @param txt text string
## @param fontsize font size
## @import grid
text.width <- function(txt, fontsize=NULL) {
  w <- convertUnit(stringWidth(paste0(txt)),'points', valueOnly=TRUE)
  if (!is.null(fontsize)) ## scale to the given font size
    w <- w / get.gpar('fontsize')$fontsize * fontsize * 1.1 # a magic number that just works 
  w
}

## text height and scale to specified font size
## @import grid
text.height <- function(txt, fontsize=NULL) {
  w <- convertUnit(stringHeight(txt),'points', valueOnly=TRUE)
  if (!is.null(fontsize)) ## scale to the given font size
    w <- w / get.gpar('fontsize')$fontsize * fontsize * 1.1 # a magic number that just works
  w
}

rotate.just1 <- function(just) {
  if (just[1]=='left')
    y <- 'top'
  else if (just[1]=='right')
    y <- 'bottom'
  else
    y <- 'center'
  if (just[2]=='bottom')
    x <- 'left'
  else if (just[2]=='top')
    x <- 'right'
  else
    x <- 'center'
  c(x,y)
}

## change just by rotation
rotate.just <- function(just, rot) {
  n <- round((rot/90) %% 4)
  if (n>0) {
    for(i in 1:n) {
      just <- rotate.just1(just)
    }
  }
  just
}
