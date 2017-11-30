.TickLabelResample <- function(labels, ticklabels.n) {
  text.height1 <- convertUnit(stringHeight('a'),'npc',valueOnly=TRUE)
  total.height <- 1
  n.labels <- length(labels)
  if (!is.null(ticklabels.n))
    n.texts <- ticklabels.n
  else if (total.height*1.2 < text.height1*n.labels) {
    n.texts <- max(floor(total.height/text.height1*0.4),2)
  } else {
    n.texts <- n.labels
  }
  sample.inds <- round(seq(1, n.labels, length.out=n.texts))
}

.WPrintXTickLabels <- function(hm, labels=NULL, cex=1) {

  if (length(labels)==1 && is.logical(labels)) {
    if (labels) {
      labels <- colnames(hm$data)
    } else {
      labels <- NULL
    }
  }

  if (!is.null(labels)) {
    nc = ncol(hm$data)
    x.mid <- (seq_len(nc)-0.5)/nc
    if (hm$xticklabel.side == 'b') {
      .text.just = 'right'
      .text.y = - hm$xticklabel.pad
    } else {
      .text.just = 'left'
      .text.y = 1 + hm$xticklabel.pad
    }
    .text.rot = hm$xticklabel.rotat
    if (!is.logical(hm$xticklabels))
      sample.inds <- which(labels %in% hm$xticklabels)
    else
      sample.inds <- .TickLabelResample(labels, hm$xticklabels.n)
    grid.text(labels[sample.inds],
              x=x.mid[sample.inds], y=unit(.text.y,'npc'), rot=.text.rot,
              just=c(.text.just, 'center'), gp=gpar(fontsize=hm$xticklabel.fontsize*cex))
  }
}

.WPrintYTickLabels <- function(hm, labels=NULL, cex=1) {

  if (length(labels)==1 && is.logical(labels)) {
    if (labels) {
      labels <- rownames(hm$data)
    } else {
      labels <- NULL
    }
  }

  if (!is.null(labels)) {
    nr = nrow(hm$data)
    y.mid <- (rev(seq_len(nr))-0.5)/nr
    if (hm$yticklabel.side == 'l') {
      .text.just = 'right'
      .text.x = - hm$yticklabel.pad
    } else {
      .text.just = 'left'
      .text.x = 1 + hm$yticklabel.pad
    }
    .text.rot = hm$yticklabel.rotat
    if (!is.logical(hm$yticklabels))
      sample.inds <- which(labels %in% hm$yticklabels)
    else
      sample.inds <- .TickLabelResample(labels, hm$yticklabels.n)
    grid.text(labels[sample.inds],
              x=unit(.text.x,'npc'), y=y.mid[sample.inds], rot=.text.rot,
              just=c(.text.just,'center'), gp=gpar(fontsize=hm$yticklabel.fontsize*min(cex)))
  }
}
