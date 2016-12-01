#' test differential methylation on each locus
#'
#' @param betas beta values
#' @param sample.data data frame for sample information, column names
#' are predictor variables (e.g., sex, age, treatment, tumor/normal etc)
#' and are referenced in formula. Rows are samples.
#' @param formula formula
#' @param se.lb lower bound to standard error of slope, lower this to get
#' more difference of small effect size.
#' @param balanced whether design is balanced or not. default to TRUE, when
#' unbalanced will use Welch's method to estimate standard error. This is slower
#' balanced=TRUE still give a good result in most cases.
#' @param cf.test factors to test (default to all factors in formula except
#' intercept). Use "all" for all factors.
#' @return cf coefficient table for each factor
#' @export
DML <- function(betas, sample.data, formula, se.lb=0.06, balanced=TRUE, cf.test=NULL) {

  design <- model.matrix(formula, sample.data)
  n.cpg <- dim(betas)[1]
  n.cf <- dim(design)[2]

  ## cf.test specify the factors to be reported
  if (is.null(cf.test)) {
    cf.test <- colnames(design)
    cf.test <- cf.test[cf.test != '(Intercept)']
  } else if (cf.test == 'all') {
    cf.test <- colnames(design)
  }

  ## preprare output
  cf <- lapply(cf.test, function(cfi) matrix(data=NA, nrow=n.cpg, ncol=5, dimnames=list(rownames(betas), c('Estimate', 'Std. Error', 't-stat', 'Pr(>|t|)', 'Effect size'))))
  names(cf) <- cf.test

  message('Testing differential methylation on each locus:')
  n.skip <- 0
  group <- factor(apply(design, 1, paste, collapse="_"))
  for (i in 1:n.cpg) {

    if (i%%ceiling(n.cpg/80)==0) message('.', appendLF=FALSE);
    
    ## filter NA
    sample.is.na <- is.na(betas[i,])
    if (all(sample.is.na)) next;
    
    design1 <- design[!sample.is.na,,drop=FALSE]
    betas1 <- betas[i,!sample.is.na]
    if (sum(apply(design1[,2:n.cf,drop=FALSE], 2, function(x) length(unique(x)))==1) > 0) {
      n.skip <- n.skip + 1
      next;
    }

    ## sigma is padded, boundary adjusted, and log2 damped (should use beta distribution)
    ## wts <- rep(1, length(betas1))
    pad <- 0.01
    betas1.padded <- pmin(pmax(betas1, pad), 1-pad)
    wts <- 1/(betas1.padded*(1-betas1.padded)) # 1/var
    stopifnot(all(wts>0))

    ## QR-solve weighted least square
    z <- .lm.fit(design1*wts, betas1*wts)
    names(z$coefficients) <- colnames(design1)
    p1 <- 1:z$rank
    coefs <- z$coefficients[z$pivot[p1]]
    residuals <- z$residuals / wts

    if (balanced) {
      se <- sum(residuals^2)/(length(residuals)-1)
    } else {
      ## Welch-Satterthwaite se
      rs <- residuals^2
      se <- apply(design1, 2, function(group1) {
        sqrt(sum(vapply(split(rs, group1), mean, numeric(1)), na.rm=TRUE))
      })
    }
    se <- pmax(se, se.lb)      # lower bound coefficient se, this selects against high effect size differences
    rdf <- apply(design1, 2, function(x) max(min(tabulate(x)),1))

    ## t-statistics
    t.stat <- coefs / se
    pval <- 2*pt(abs(t.stat), rdf, lower.tail=FALSE)
    stopifnot(!any(is.na(pval)))
    
    ## output
    fitted.rg <- range(betas1 - residuals)
    eff <- fitted.rg[2] - fitted.rg[1]  # effect size
    ans <- cbind(coefs, se, t.stat, pval, eff)
    lapply(cf.test, function(cfi)
      if (cfi %in% rownames(ans))
        cf[[cfi]][i,] <<- ans[cfi,])
    ## betas.fitted[i,!sample.is.na] <- betas1 - residuals
  }
  message('.\n', appendLF=FALSE)
  ## CpGs are correlated, should apply p-value adjustment on segments
  ## cf <- lapply(cf, function(cf1) cbind(cf1, P.adjusted=p.adjust(cf1[,'Pr(>|t|)'],method='BH')))
  class(cf) <- 'diffMeth'

  message('Significant loci (p<0.05):')
  sigcnts <- lapply(cf, function(x) sum(x[,'Pr(>|t|)'] < 0.05, na.rm=TRUE))
  sigmsg <- lapply(seq_along(sigcnts), function(i) sprintf(' - %s: %d', names(sigcnts[i]), sigcnts[i][[1]]))
  message(do.call(paste0, list(sigmsg, collapse='\n')))

  cf
}

#' find DMR
#'
#' This subroutine uses Euclidean distance to group CpGs and
#' then combine p-values for each segment.
#' 
#' @param cf coefficient table from diffMeth, when NULL will be computed from beta
#' @param betas beta values for distance calculation
#' @param sample.data data frame for sample information, column names
#' are predictor variables (e.g., sex, age, treatment, tumor/normal etc)
#' and are referenced in formula. Rows are samples.
#' @param formula formula
#' @param dist.cutoff distance cutoff (default to use dist.cutoff.quantile)
#' @param seg.per.locus number of segments per locus
#' higher value leads to more segments
#' @param platform EPIC or hm450
#' @param refversion hg38 or hg19
#' @param ... additional parameters to DML
#' @return coefficient table with segment ID and segment P-value
#' @export
DMR <- function(betas, sample.data=NULL, formula=NULL, cf=NULL, dist.cutoff=NULL, seg.per.locus=0.5, platform='EPIC', refversion='hg38', ...) {

  pkgTest('GenomicRanges')

  if (is.null(cf)) {
    if (is.null(sample.data) || is.null(formula))
      stop('Need either cf or sample data and formula to run DML.')
    cf <- DML(betas, sample.data, formula, ...)
  }
  
  ## filter NA
  betas.noNA <- betas[!apply(betas, 1, function(x) all(is.na(x))),]

  ## sort by coordinates
  probe.coords <- getBuiltInData(paste0(platform,'.mapped.probes.',refversion))
  cpg.ids <- intersect(rownames(betas.noNA), names(probe.coords))
  probe.coords <- GenomicRanges::sort(probe.coords[cpg.ids])
  betas.coord.srt <- betas.noNA[names(probe.coords),]

  message("Merging correlated CpGs...", appendLF=FALSE)
  cpg.ids <- rownames(betas.coord.srt)
  cpg.coords <- probe.coords[cpg.ids]
  cpg.chrm <- as.vector(GenomicRanges::seqnames(cpg.coords))
  cpg.start <- GenomicRanges::start(cpg.coords)
  cpg.end <- GenomicRanges::end(cpg.coords)

  n.cpg <- length(cpg.ids)

  ## euclidean distance suffcies
  beta.dist <- sapply(1:(n.cpg-1), function(i) sqrt(sum((betas.coord.srt[i,] - betas.coord.srt[i+1,])^2, na.rm=TRUE)))

  ## 1-correlation coefficient
  ## beta.dist <- sapply(1:(n.cpg-1), function(i) {
  ##   x <- cor(betas.coord.srt[i,],betas.coord.srt[i+1,],use='na.or.complete',method='spearman')
  ##   if (is.na(x)) x <- 0
  ##   1-x
  ## })
  
  chrm.changed <- (cpg.chrm[-1] != cpg.chrm[-n.cpg])
  if (is.null(dist.cutoff))
    dist.cutoff <- quantile(beta.dist, 1-seg.per.locus) # empirical cutoff based on quantiles
  change.points <- (beta.dist > dist.cutoff | chrm.changed)
  seg.ids <- cumsum(c(TRUE,change.points))
  message("Done.")

  ## add back unmapped
  all.cpg.ids <- rownames(betas)
  unmapped <- all.cpg.ids[!(all.cpg.ids %in% cpg.ids)]
  cpg.ids <- c(cpg.ids, unmapped)
  cpg.chrm <- c(cpg.chrm, rep('*', length(unmapped)))
  cpg.start <- c(cpg.start, rep(NA, length(unmapped)))
  cpg.end <- c(cpg.end, rep(NA, length(unmapped)))

  ## segment coordinates
  seg.ids <- c(seg.ids, seq.int(from=seg.ids[length(seg.ids)]+1, length.out=length(unmapped)))
  seg.chrm <- as.vector(tapply(cpg.chrm, seg.ids, function(x) x[1]))
  seg.start <- as.vector(tapply(cpg.start, seg.ids, function(x) x[1]))
  seg.end <- as.vector(tapply(cpg.end, seg.ids, function(x) x[length(x)]))

  ## combine p-value
  message("Combine p-values... ", appendLF=FALSE)
  cf <- lapply(cf, function(cf1) {
    seg.pval <- as.vector(tapply(cf1[cpg.ids, 'Pr(>|t|)'], seg.ids, function(x) pnorm(sum(qnorm(x))/sqrt(length(x))))) #  Stouffer's Z-score method
    seg.pval.adj <- p.adjust(seg.pval, method="BH")
    seg.ids.cf <- match(rownames(cf1), cpg.ids)
    cf1 <- as.data.frame(cf1)
    cf1$Seg.ID <- seg.ids[seg.ids.cf]
    cf1$Seg.chrm <- seg.chrm[cf1$Seg.ID]
    cf1$Seg.start <- seg.start[cf1$Seg.ID]
    cf1$Seg.end <- seg.end[cf1$Seg.ID]
    cf1$Seg.Pval <- seg.pval[cf1$Seg.ID]
    cf1$Seg.Pval.adj <- seg.pval.adj[cf1$Seg.ID]
    cf1
  })
  message("Done.")

  cf
}

#' top segments in differential methylation
#' 
#' This is a convenience function to show top differential methylated segments.
#' 
#' @param cf1 coefficient table of one factor from segmentDMR
#' @return coefficient table ordered by adjusted p-value of segments
#' @export
topSegments <- function(cf1) {
  x <- unique(cf1[order(cf1[,'Seg.Pval']),c('Seg.ID','Seg.chrm','Seg.start','Seg.end','Seg.Pval.adj')])
  rownames(x) <- x[['Seg.ID']]
  x
}

#' top loci in differential methylation
#' 
#' This is a convenience function to show top differential methylated segments.
#' 
#' @param cf1 coefficient table of one factor from diffMeth
#' @return coefficient table ordered by p-value of each locus
#' @export
topLoci <- function(cf1) {
  cf1[order(cf1[,'Pr(>|t|)']),]
}

#' select segment from coefficient table
#'
#' @param cf1 coefficient table of one factor from segmentDMR
#' @param seg.id segment ID
#' @return coefficient table from given segment
#' @export
getSegment <- function(cf1, seg.id) {
  cf1[cf1[,'Seg.ID']==seg.id,]
}
