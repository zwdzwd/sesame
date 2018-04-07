#' Test differential methylation on each locus
#'
#' The function takes a beta value matrix with probes on the rows and
#' samples on the columns. It also takes a sample information data frame
#' (sample.data) and formula for testing. The function outputs a list of
#' coefficient tables for each factor tested.
#'
#' @param betas beta values
#' @param sample.data data frame for sample information, column names
#' are predictor variables (e.g., sex, age, treatment, tumor/normal etc)
#' and are referenced in formula. Rows are samples.
#' @param formula formula
#' @param se.lb lower bound to standard error of slope, lower this to get
#' more difference of small effect size.
#' @param balanced whether design is balanced or not. default to FALSE, when
#' unbalanced will use Welch's method to estimate standard error.
#' balance=TRUE is faster.
#' @param cf.test factors to test (default to all factors in formula except
#' intercept). Use "all" for all factors.
#' @return cf - a list of coefficient tables for each factor
#' @import stats
#' @examples
#' betas <- readRDS(system.file(
#'     'extdata','HM450.betas.76matchedTCGAchr20.rds',
#'     package='sesameData'))
#' sample.info <- readRDS(system.file(
#'     'extdata','HM450.sampleinfo.76matchedTCGAchr20.rds',
#'     package='sesameData'))
#' cf <- DML(betas, sample.info, ~type)
#' @export
DML <- function(
    betas, sample.data, formula, se.lb=0.06, balanced=FALSE, cf.test=NULL) {

    design <- model.matrix(formula, sample.data)
    ## convert to factor for faster processing
    design.fac <- data.frame(lapply(as.data.frame(design), as.factor))
    rdf0 <- unlist(
        lapply(design.fac, function(x) max(min(tabulate(x)),1)),
        recursive = FALSE, use.names=FALSE)
    
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
    cf <- lapply(cf.test, function(cfi) matrix(
        data=NA, nrow=n.cpg, ncol=5,
        dimnames=list(rownames(betas),
        c('Estimate', 'Std. Error', 't-stat', 'Pr(>|t|)', 'Effect size'))))
    
    names(cf) <- cf.test

    message('Testing differential methylation on each locus:')
    for (i in seq_len(n.cpg)) {

        if (i%%ceiling(n.cpg/80)==0) message('.', appendLF=FALSE);
        
        ## filter NA
        sample.is.na <- is.na(betas[i,])
        if (all(sample.is.na)) next;

        if (any(sample.is.na)) {      # this "if" improves performance
            design1 <- design[!sample.is.na,,drop=FALSE]
            design1.fac <- design.fac[!sample.is.na,]
            betas1 <- betas[i,!sample.is.na]

            if (sum(apply(
                design1[,2:n.cf,drop=FALSE], 2,
                function(x) length(unique(x)))==1) > 0) next;
            
            rdf <- unlist(lapply(design1.fac, function(x) max(min(
                tabulate(x)),1)), recursive = FALSE, use.names=FALSE)
            
        } else {
            design1 <- design
            design1.fac <- design.fac
            betas1 <- betas[i,]
            rdf <- rdf0
        }

        ## sigma is padded, boundary adjusted, and log2
        ## damped (should use beta distribution)
        ## wts <- rep(1, length(betas1))
        pad <- 0.01
        betas1.padded <- pmin(pmax(betas1, pad), 1-pad)
        wts <- 1/(betas1.padded*(1-betas1.padded)) # 1/var

        ## QR-solve weighted least square
        z <- .lm.fit(design1*wts, betas1*wts)
        names(z$coefficients) <- colnames(design1)
        p1 <- seq_len(z$rank)
        coefs <- z$coefficients[z$pivot[p1]]
        residuals <- z$residuals / wts

        if (balanced) {
            se <- sum(residuals^2)/(length(residuals)-1)
        } else {
            ## Welch-Satterthwaite se
            rs <- residuals^2
            se <- unlist(
                lapply(design1.fac, function(group1) {sqrt(sum(
                    vapply(split(
                        rs, group1), mean, numeric(1)), na.rm=TRUE))}),
                recursive = FALSE, use.names=FALSE)
        }
        
        ## lower bound coefficient se, this selects against
        ## high effect size differences
        se <- pmax(se, se.lb)

        ## t-statistics
        t.stat <- coefs / se
        pval <- 2*pt(abs(t.stat), rdf, lower.tail=FALSE)
        stopifnot(!any(is.na(pval)))

        ## output
        fitted.rg <- range(betas1 - residuals)
        eff <- fitted.rg[2] - fitted.rg[1]  # effect size
        ans <- cbind(coefs, se, t.stat, pval, eff)
        for (cfi in cf.test) {
            if (cfi %in% rownames(ans))
                cf[[cfi]][i,] <- ans[cfi,]
        }
        ## betas.fitted[i,!sample.is.na] <- betas1 - residuals
    }
    message('.\n', appendLF=FALSE)
    ## CpGs are correlated, should apply p-value adjustment on segments
    ## cf <- lapply(cf, function(cf1) cbind(cf1,
    ## P.adjusted=p.adjust(cf1[,'Pr(>|t|)'],method='BH')))

    message('Significant loci (p<0.05):')
    sigcnts <- lapply(cf, function(x) sum(x[,'Pr(>|t|)'] < 0.05, na.rm=TRUE))
    sigmsg <- lapply(seq_along(sigcnts), function(i) sprintf(
        ' - %s: %d significant loci.', names(sigcnts[i]), sigcnts[i][[1]]))
    message(do.call(paste0, list(sigmsg, collapse='\n')))

    cf
}

#' Find Differentially Methylated Region (DMR)
#'
#' This subroutine uses Euclidean distance to group CpGs and
#' then combine p-values for each segment. The function performs DML test first
#' if cf is NULL. It groups the probe testing results into differential
#' methylated regions in a coefficient table with additional columns
#' designating the segment ID and statistical significance (P-value) testing
#' the segment.
#' 
#' @param betas beta values for distance calculation
#' @param sample.data data frame for sample information, column names
#' are predictor variables (e.g., sex, age, treatment, tumor/normal etc)
#' and are referenced in formula. Rows are samples.
#' @param formula formula
#' @param cf coefficient table from diffMeth, when NULL will be computed from
#' beta. If cf is given, sample.data and formula are ignored.
#' @param dist.cutoff distance cutoff (default to use dist.cutoff.quantile)
#' @param seg.per.locus number of segments per locus
#' higher value leads to more segments
#' @param platform EPIC or HM450
#' @param refversion hg38 or hg19
#' @param ... additional parameters to DML
#' @return coefficient table with segment ID and segment P-value
#' @examples
#' betas <- readRDS(system.file(
#'     'extdata','HM450.betas.76matchedTCGAchr20.rds',
#'     package='sesameData'))
#' sample.info <- readRDS(system.file(
#'     'extdata','HM450.sampleinfo.76matchedTCGAchr20.rds',
#'     package='sesameData'))
#' cf <- DMR(betas, sample.info, ~type)
#' @export
DMR <- function(
    betas, sample.data=NULL, formula=NULL, cf=NULL, dist.cutoff=NULL,
    seg.per.locus=0.5, platform='EPIC', refversion='hg38', ...) {

    pkgTest('GenomicRanges')
    
    if (is.null(cf)) {
        if (is.null(sample.data) || is.null(formula))
            stop('Need either cf or sample data and formula to run DML.')
        cf <- DML(betas, sample.data, formula, ...)
    }
    
    ## filter NA
    betas.noNA <- betas[!apply(betas, 1, function(x) all(is.na(x))),]

    ## sort by coordinates
    probe.coords <- get(paste0(platform,'.mapped.probes.',refversion))
    cpg.ids <- intersect(rownames(betas.noNA), names(probe.coords))
    probe.coords <- GenomicRanges::sort(probe.coords[cpg.ids])
    betas.coord.srt <- betas.noNA[names(probe.coords),]

    message("Merging correlated CpGs ... ", appendLF=FALSE)
    cpg.ids <- rownames(betas.coord.srt)
    cpg.coords <- probe.coords[cpg.ids]
    cpg.chrm <- as.vector(GenomicRanges::seqnames(cpg.coords))
    cpg.start <- GenomicRanges::start(cpg.coords)
    cpg.end <- GenomicRanges::end(cpg.coords)

    n.cpg <- length(cpg.ids)

    ## euclidean distance suffcies
    beta.dist <- vapply(seq_len(n.cpg-1), function(i) sum(
        (betas.coord.srt[i,] - betas.coord.srt[i+1,])^2, na.rm=TRUE), 1)

    ## 1-correlation coefficient
    ## beta.dist <- sapply(seq_len(n.cpg-1), function(i) {
    ##   x <- cor(betas.coord.srt[i,],betas.coord.srt[i+1,],
    ## use='na.or.complete',method='spearman')
    ##   if (is.na(x)) x <- 0
    ##   1-x
    ## })
    
    chrm.changed <- (cpg.chrm[-1] != cpg.chrm[-n.cpg])
    ## empirical cutoff based on quantiles
    if (is.null(dist.cutoff))
        dist.cutoff <- quantile(beta.dist, 1-seg.per.locus)
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
    seg.ids <- c(
        seg.ids, seq.int(
            from=seg.ids[length(seg.ids)]+1, length.out=length(unmapped)))
    
    seg.chrm <- as.vector(tapply(cpg.chrm, seg.ids, function(x) x[1]))
    seg.start <- as.vector(tapply(cpg.start, seg.ids, function(x) x[1]))
    seg.end <- as.vector(tapply(cpg.end, seg.ids, function(x) x[length(x)]))

    message(sprintf('Generated %d segments.', seg.ids[length(seg.ids)]))

    ## combine p-value
    message("Combine p-values ... ")
    cfnames <- names(cf)
    cf <- lapply(seq_along(cf), function(i) {
        cf1 <- cf[[i]]
        ## Stouffer's Z-score method
        seg.pval <- as.vector(tapply(
            cf1[cpg.ids, 'Pr(>|t|)'], seg.ids,
            function(x) pnorm(sum(qnorm(x))/sqrt(length(x)))))
        
        seg.pval.adj <- p.adjust(seg.pval, method="BH")
        seg.ids.cf <- match(rownames(cf1), cpg.ids)
        cf1 <- as.data.frame(cf1)
        cf1$Seg.ID <- seg.ids[seg.ids.cf]
        cf1$Seg.chrm <- seg.chrm[cf1$Seg.ID]
        cf1$Seg.start <- seg.start[cf1$Seg.ID]
        cf1$Seg.end <- seg.end[cf1$Seg.ID]
        cf1$Seg.Pval <- seg.pval[cf1$Seg.ID]
        cf1$Seg.Pval.adj <- seg.pval.adj[cf1$Seg.ID]
        message(sprintf(
            ' - %s: %d significant segments.',
            cfnames[i], sum(seg.pval<0.05, na.rm=TRUE)))
        message(sprintf(
            ' - %s: %d significant segments (after BH).',
            cfnames[i], sum(seg.pval.adj<0.05, na.rm=TRUE)))
        cf1
    })
    names(cf) <- cfnames
    message("Done.")

    cf
}

#' Top segments in differential methylation
#' 
#' This is a utility function to show top differential methylated segments.
#' The function takes a coefficient table as input and output the same
#' table ordered by the sigificance of the segments.
#' 
#' @param cf1 coefficient table of one factor from DMR
#' @return coefficient table ordered by adjusted p-value of segments
#' @examples
#' betas <- readRDS(system.file(
#'     'extdata','HM450.betas.76matchedTCGAchr20.rds',
#'     package='sesameData'))
#' sample.info <- readRDS(system.file(
#'     'extdata','HM450.sampleinfo.76matchedTCGAchr20.rds',
#'     package='sesameData'))
#' cf <- DMR(betas, sample.info, ~type)
#' topSegments(cf[[1]])
#' @export
topSegments <- function(cf1) {
    x <- unique(cf1[order(cf1[,'Seg.Pval']), c(
        'Seg.ID','Seg.chrm','Seg.start','Seg.end','Seg.Pval','Seg.Pval.adj')])
    rownames(x) <- x[['Seg.ID']]
    x
}

#' Top loci in differential methylation
#' 
#' This is a convenience function to show top differential methylated segments.
#' The function takes a coefficient table as input and output the same
#' table ordered by the sigificance of the locus.
#' 
#' @param cf1 coefficient table of one factor from diffMeth
#' @return coefficient table ordered by p-value of each locus
#' @examples
#' betas <- readRDS(system.file(
#'     'extdata','HM450.betas.76matchedTCGAchr20.rds',
#'     package='sesameData'))
#' sample.info <- readRDS(system.file(
#'     'extdata','HM450.sampleinfo.76matchedTCGAchr20.rds',
#'     package='sesameData'))
#' cf <- DMR(betas, sample.info, ~type)
#' topLoci(cf[[1]])
#' @export
topLoci <- function(cf1) {
    cf1[order(cf1[,'Pr(>|t|)']),]
}

#' Select segment from coefficient table
#'
#' This function takes a coefficient table and returns a subset of the table
#' targeting only the specified segment using segment ID.
#'
#' @param cf1 coefficient table of one factor from DMR
#' @param seg.id segment ID
#' @return coefficient table from given segment
#' @examples
#' betas <- readRDS(system.file(
#'     'extdata','HM450.betas.76matchedTCGAchr20.rds',
#'     package='sesameData'))
#' sample.info <- readRDS(system.file(
#'     'extdata','HM450.sampleinfo.76matchedTCGAchr20.rds',
#'     package='sesameData'))
#' cf <- DMR(betas, sample.info, ~type)
#' getSegment(cf[[1]], cf[[1]][['Seg.ID']][1])
#' @export
getSegment <- function(cf1, seg.id) {
    cf1[cf1[,'Seg.ID']==seg.id,]
}
