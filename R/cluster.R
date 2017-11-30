
#' row cluster a matrix
#'
#' row cluster a matrix
#'
#' @param mat input matrix
#' @param ... extra color bars or matrix that needs row reordered.
#' @param hc.method method to use in hclust
#' @param dist.method method to use in dist
#' @return a list of clustered row, column and matrix
#' @export
row.cluster <- function(mat, ..., hc.method='ward.D2', dist.method='euclidean') {
  d.row <- dist(mat, method=dist.method)
  r <- list()
  r$row.clust <- hclust(d.row, method=hc.method)
  r$column.clust <- NULL
  r$mat <- mat[r$row.clust$order,]
  r$extra <- lapply(list(...), function(x) {
    if ('matrix' %in% class(x)) {
      x[r$row.clust$order, ]
    } else {
      x[r$row.clust$order]
    }
  })
  r
}

#' column cluster a matrix
#'
#' column cluster a matrix
#'
#' @param mat input matrix
#' @param ... extra color bars or matrix that needs column reordered
#' @param hc.method method to use in hclust
#' @param dist.method method to use in dist
#' @return a list of clustered row, column and matrix
#' @export
column.cluster <- function(mat, ..., hc.method='ward.D2', dist.method='euclidean') {
  d.column <- dist(t(mat), method=dist.method)
  r <- list()
  r$row.clust <- NULL
  r$column.clust <- hclust(d.column, method=hc.method)
  r$mat <- mat[,r$column.clust$order]
  r$extra <- lapply(list(...), function(x) {
    if ('matrix' %in% class(x)) {
      x[,r$column.clust$order]
    } else {
      x[r$column.clust$order]
    }
  })
  r
}

#' row- and column-cluster a matrix
#'
#' row- and column-cluster a matrix
#'
#' @param mat input matrix
#' @param hc.method method to use in hclust
#' @param dist.method method to use in dist
#' @param extra.row extra row reordering
#' @param extra.column extra column reordering
#' @return a list of clustered row, column and matrix
#' @import stats
#' @export
both.cluster <- function(mat, extra.row=NULL, extra.column=NULL, hc.method='ward.D2', dist.method='euclidean') {
  d.row <- dist(mat, method=dist.method)
  d.column <- dist(t(mat), method=dist.method)
  r <- list()
  r$row.clust <- hclust(d.row, method=hc.method)
  r$column.clust <- hclust(d.column, method=hc.method)
  r$mat <- mat[r$row.clust$order, r$column.clust$order]
  r$extra <- c(
    lapply(extra.row, function(x) {
      if ('matrix' %in% class(x)) {
        x[r$row.clust$order, ]
      } else {
        x[r$row.clust$order]
      }
    }), lapply(extra.column, function(x) {
      if ('matrix' %in% class(x)) {
        x[,r$column.clust$order]
      } else {
        x[r$column.clust$order]
      }
    }))
  r
}



