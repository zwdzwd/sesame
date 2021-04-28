## 1-correlation coefficient
## beta.dist <- sapply(seq_len(n.cpg-1), function(i) {
##   x <- cor(betas.coord.srt[i,],betas.coord.srt[i+1,],
## use='na.or.complete',method='spearman')
##   if (is.na(x)) x <- 0; 1-x; })

DMLShrinkage1 <- function(
    betas, design, i, rdf0, design.fac, se.lb=0.06, balanced=FALSE) {
    
    n.cf <- dim(design)[2]
    ## filter NA
    sample.is.na <- is.na(betas[i,])
    if (all(sample.is.na)) return(NULL);

    if (any(sample.is.na)) {      # this "if" improves performance
        design1 <- design[!sample.is.na,,drop=FALSE]
        design1.fac <- design.fac[!sample.is.na,]
        betas1 <- betas[i,!sample.is.na]
        if (sum(apply(
            design1[,2:n.cf,drop=FALSE], 2,
            function(x) length(unique(x)))==1) > 0) return(NULL);
        
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
    ## coefs <- z$coefficients[z$pivot[seq_len(z$rank)]]
    coefs <- z$coefficients
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
    t.stat <- coefs / se[seq_along(coefs)]
    pval <- 2*pt(abs(t.stat), rdf, lower.tail=FALSE)
    stopifnot(!any(is.na(pval)))
    ## output
    fitted.rg <- range(betas1 - residuals)
    eff <- fitted.rg[2] - fitted.rg[1]  # effect size
    cbind(coefs, se, t.stat, pval, eff)
}


#' Test differential methylation on each locus Using Shrinkage Estimator
#'
#' The function takes a beta value matrix with probes on the rows and
#' samples on the columns. It also takes a sample information data frame
#' (meta) and formula for testing. The function outputs a list of
#' coefficient tables for each factor tested.
#'
#' @param betas beta values, matrix or SummarizedExperiment
#' @param formula formula
#' @param meta data frame for sample information, column names
#' are predictor variables (e.g., sex, age, treatment, tumor/normal etc)
#' and are referenced in formula. Rows are samples.
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
#' sesameDataCache("HM450") # in case not done yet
#' data <- sesameDataGet('HM450.76.TCGA.matched')
#' cf_list <- DMLShrinkage(data$betas, ~type, meta=data$sampleInfo)
#' @export
DMLShrinkage <- function(
    betas, formula, meta=NULL, se.lb=0.06, balanced=FALSE, cf.test=NULL) {

    if(is(betas, "SummarizedExperiment")) {
        betas0 = betas
        betas = assay(betas0)
        meta = colData(betas0)
    }
    design <- model.matrix(formula, meta)
    ## convert to factor for faster processing
    design.fac <- data.frame(lapply(as.data.frame(design), as.factor))
    rdf0 <- unlist(
        lapply(design.fac, function(x) max(min(tabulate(x)),1)),
        recursive = FALSE, use.names=FALSE)
    n.cpg <- dim(betas)[1]

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
        ans <- DMLShrinkage1(betas, design, i, rdf0, design.fac,
            se.lb=se.lb, balanced=balanced)
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
