#' test differential methylation on each locus
#'
#' @param betas beta values
#' @param sample.data data frame for sample information, column names
#' are predictor variables (e.g., sex, age, treatment, tumor/normal etc)
#' and are referenced in formula. Rows are samples.
#' @param formula formula
#' @param se.lb lower bound to standard error of slope
#' @param cf.test factors to test (default to all factors in formula except
#' intercept). Use "all" for all factors.
#' @return cf coefficient table for each factor
#' @export
diffMeth <- function(betas, sample.data, formula, se.lb=0.01, cf.test=NULL) {

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
  cf <- lapply(cf.test, function(cfi) matrix(data=NA, nrow=n.cpg, ncol=6, dimnames=list(rownames(betas), c('Estimate', 'Std. Error', 't-stat', 'Pr(>|t|)', 'residual d.f.', 'Effect size'))))
  names(cf) <- cf.test

  cat('Testing differential methylation on each locus:\n')
  n.skip <- 0
  group <- factor(apply(design, 1, paste, collapse="_"))
  for (i in 1:n.cpg) {

    if (i%%5000==0) message('.', appendLF=FALSE);
    if (i%%300000==0) message('\n', appendLF=FALSE);
    
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

    rss <- sum(z$residuals^2) # residual sum of squares
    if (rss<=0) {
      rdf <- 1 # arbitrarily set degree of freedom
      se <- 1/wts # slope se
    } else {
      ## Welch-Satterthwaite correction for residual degree of freedom
      group1 <- group[!sample.is.na]
      group1N <- pmax(tabulate(group1),1)
      rss.group <- vapply(split(z$residuals, group1), function(x) sum(x^2), numeric(1)) / group1N
      rdf <- sum(rss.group)^2/sum(rss.group^2/pmax(group1N-(z$rank-1),1))
      ## rdf <- length(betas1) - z$rank # only works for balanced design

      ## slope se
      XTXinv <- chol2inv(z$qr[p1,p1,drop=FALSE]) # (t(X)*X)^-1
      resvar <- rss/rdf # residual variance
      se <- sqrt(diag(XTXinv)*resvar) # slope se
    }
    se <- pmax(se, se.lb)      # lower bound coefficient se
    
    ## t-statistics
    t.stat <- coefs / se
    pval <- 2*pt(abs(t.stat), rdf, lower.tail=FALSE)
    stopifnot(!any(is.na(pval)))
    
    ## output
    fitted.rg <- range(betas1 - z$residuals / wts)
    eff <- fitted.rg[2] - fitted.rg[1]  # effect size
    ans <- cbind(coefs, se, t.stat, pval, rdf, eff)
    lapply(cf.test, function(cfi)
      if (cfi %in% rownames(ans))
        cf[[cfi]][i,] <<- ans[cfi,])
    ## z$residuals/wts
    ## betas.fitted[i,!sample.is.na] <- betas1 - z$residuals/wts
  }
  message('.\n', appendLF=FALSE)
  cf <- lapply(cf, function(cf1) cbind(cf1, P.adjusted=p.adjust(cf1[,'Pr(>|t|)'],method='BH')))
  class(cf) <- 'diffMeth'

  cf
}

#' segment DMR
#'
#' This subroutine uses Euclidean distance to group CpGs and
#' then combine p-values for each segment.
#' 
#' @param cf coefficient table from diffMeth
#' @param betas beta values for distance calculation
#' @param dist.cutoff distance cutoff (default to use dist.cutoff.quantile)
#' @param dist.cutoff.quantile quantile to use in selecting cutoff distance
#' @param platform EPIC or hm450
#' @param refversion hg38 or hg19
#' @return coefficient table with segment ID and segment P-value
#' @export
segmentDMR <- function(betas, cf, dist.cutoff=NULL, dist.cutoff.quantile=0.5, platform='EPIC', refversion='hg38') {

  pkgTest('GenomicRanges')
  
  ## filter NA
  betas.noNA <- betas[!apply(betas, 1, function(x) all(is.na(x))),]

  ## sort by coordinates
  probe.coords <- getBuiltInData(paste0(platform,'.mapped.probes.',refversion))
  cpg.ids <- intersect(rownames(betas.noNA), names(probe.coords))
  probe.coords <- GenomicRanges::sort(probe.coords[cpg.ids])
  betas.coord.srt <- betas.noNA[names(probe.coords),]

  message("Merging correlated CpGs...", appendLF=FALSE)
  cpg.ids <- rownames(betas.coord.srt)
  n.cpg <- length(cpg.ids)
  euc.dist <- sapply(1:(n.cpg-1), function(i) sqrt(sum((betas.coord.srt[i,] - betas.coord.srt[i+1,])^2, na.rm=TRUE))) # euclidean distance
  chrm <- GenomicRanges::seqnames(probe.coords[cpg.ids])
  chrm.changed <- as.vector(chrm[-1]!=chrm[-n.cpg])
  dist.cutoff <- quantile(euc.dist, dist.cutoff.quantile) # empirical cutoff based on quantiles
  change.points <- (euc.dist > dist.cutoff | chrm.changed)
  seg.ids <- cumsum(c(TRUE,change.points))
  message("Done.")

  ## add back unmapped
  all.cpg.ids <- rownames(betas)
  unmapped <- all.cpg.ids[!(all.cpg.ids %in% cpg.ids)]
  cpg.ids <- c(cpg.ids, unmapped)
  seg.ids <- c(seg.ids, seq.int(from=seg.ids[length(seg.ids)]+1, length.out=length(unmapped)))

  ## combine p-value
  message("Combine p-values...", appendLF=FALSE)
  cf <- lapply(cf, function(cf1) {
    seg.pval <- tapply(cf1[cpg.ids, 'P.adjusted'], seg.ids, function(x) pnorm(sum(qnorm(x))/sqrt(length(x)))) #  Stouffer's Z-score method
    seg.ids.cf <- seg.ids[match(rownames(cf1), cpg.ids)]
    cf1 <- cbind(cf1, Seg.ID = seg.ids.cf, Seg.Pval=seg.pval[seg.ids.cf])
    cf1
  })
  message("Done.")

  cf
}
