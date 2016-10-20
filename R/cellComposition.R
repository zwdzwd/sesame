## clean reference set to non NA sites
cleanRefSet <- function(g, platform) {
  mapinfo <- getBuiltInData('mapped.probes.hg19', platform)
  g.clean <- g[apply(g, 1, function(x) !any(is.na(x))),]
  g.clean <- g.clean[rownames(g.clean) %in% names(mapinfo),]
  g.clean <- g.clean[!(as.vector(seqnames(mapinfo[rownames(g.clean)])) %in% c('chrX','chrY','chrM')),]
  g.clean <- g.clean[grep('cg',rownames(g.clean)),]
  g.clean
}

#' restrict refset to differentially methylated probes
#' use with care, might introduce bias
#' 
#' @param g reference set
#' @return g
#'
#' @export
diffRefSet <- function(g) {
  g <- g[apply(g, 1, function(x) min(x) != max(x)),]
  message('Reference set is based on ', dim(g)[1], ' probes from ', dim(g)[2], ' cell types.')
  g
}

#' retrieve reference set
#'
#' @param cells reference cell types
#' @param platform EPIC or hm450
#' @return g
#' @export
getRefSet <- function(cells=NULL, platform='EPIC') {
  if (is.null(cells)) {
    cells <- c('CD4T', 'CD19B','CD56NK','CD14Monocytes', 'granulocytes');
  }

  g <- sapply(cells, function(cell) getBuiltInData(cell, platform, 'cellref'));
  g <- cleanRefSet(g, platform);
  message('Reference set is based on ', dim(g)[1], ' probes from ', dim(g)[2], ' cell types.')
  g
}

errFunc <- function(f, g, q) {
  gamma <- q - g %*% f[2:length(f)]
  sum(ifelse(gamma < f[1] / 2, abs(gamma), abs(gamma - f[1])), na.rm=T)
}

double.transform.f <- function(f, nu1, nu2, step.size) {
  if (f[nu1] + step.size > 1) return(NULL);
  if (f[nu2] - step.size < 0) return(NULL);
  f[nu1] <- f[nu1] + step.size
  f[nu2] <- f[nu2] - step.size
  f[1] <- 1 - sum(f[2:length(f)]) # renormalize to avoid numerical drift
  f
}

dichotomize <- function(q, rmin=0.15, rmax=0.85) {
  q[q > rmax] <- 1
  q[q < rmin] <- 0
  middle <- na.omit(which(q >= rmin & q <= rmax))
  q[middle] <- (q[middle]-rmin) / (rmax-rmin)
  q
}

getg0 <- function(f, g, q) {
  gamma <- q - g %*% f[2:length(f)]
  ifelse(gamma < f[1] / 2, 0, 1)
}

#' Estimate cell composition using reference
#'
#' @param g reference methylation
#' @param q target measurement: length(q) == nrow(g)
#' @param temp annealing temperature
#' @param maxIter maximum iteration number after nothing changed
#' @param delta delta score to reset counter
#' @param verbose output debug info
#' @return a list of fraction, min error and unknown component methylation state
#' @export
estimateCellComposition <- function(g, q, temp=0.5, maxIter=1000, delta=0.0001, verbose=FALSE) {

  q <- dichotomize(q)
  
  ## initialize
  M <- ncol(g)
  frac <- c(1, rep(0, M))
  errcurrent <- errFunc(frac, g, q)
  errmin <- errcurrent
  frac.min <- frac
  niter <- 1
  repeat {
    nu <- sample(1:(M+1), 2)
    step.size <- runif(1)
    frac.test <- double.transform.f(frac, nu[1], nu[2], step.size);
    if (!is.null(frac.test)) {

      if (verbose) {
        message(errcurrent, '\t')
        message(do.call(paste,
                        lapply(frac, function(x) sprintf('%1.2f', x))),
                '\t', step.size, '\t', temp, '\t',
                do.call(paste, lapply(frac.min, function(x) sprintf('%1.2f', x))), '\t', errmin)
      }

      errtest <- errFunc(frac.test, g, q)
      
      ## update best
      if (errtest < errmin) {
        errmin <- errtest
        frac.min <- frac.test
        if ((errmin-errtest) > errmin*delta)
          niter <- 1;
      } else {
        niter <- niter + 1
        if (niter > maxIter) {
          break;
        }
      }
      
      ## rejection sampling
      if (runif(1) < exp(-(errtest-errcurrent)/temp)) {
        errcurrent <- errtest
        frac <- frac.test
      }
    }
  }

  list(frac = setNames(frac.min, c("unknown", colnames(g))),
       err = errmin,
       g0 = getg0(frac.min, g, q))
}
