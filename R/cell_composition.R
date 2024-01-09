## clean reference set to non NA sites
cleanRefSet <- function(g, platform = c('EPIC','HM450','HM27')) {

    platform <- match.arg(platform)
    mapinfo <- sesameDataGet(paste0(
        platform, '.probeInfo'))[[paste0('mapped.probes.hg19')]]
    ## mapinfo <- sesameDataGet(paste0(
    ## platform, ".address"))[["hg19"]]
    g <- g[GenomicRanges::intersect(rownames(g), names(mapinfo)),,drop=FALSE]
    g.clean <- g[apply(g, 1, function(x) !any(is.na(x))),,drop=FALSE]
    g.clean <- g.clean[rownames(g.clean) %in% names(mapinfo),,drop=FALSE]
    g.clean <- g.clean[!(as.vector(GenomicRanges::seqnames(
        mapinfo[rownames(g.clean)])) %in% c('chrX','chrY','chrM')),,drop=FALSE]
    g.clean <- g.clean[grep('cg',rownames(g.clean)),,drop=FALSE]
    g.clean
}

#' Restrict refset to differentially methylated probes
#' use with care, might introduce bias
#'
#' The function takes a matrix with probes on the rows and cell types on
#' the columns and output a subset matrix and only probes that show
#' discordant methylation levels among the cell types.
#' 
#' @param g a matrix with probes on the rows and cell types on the columns
#' @return g a matrix with a subset of input probes (rows)
#' @examples
#'
#' g = diffRefSet(getRefSet(platform='HM450'))
#' sesameDataGet_resetEnv()
#' 
#' @export
diffRefSet <- function(g) {
    g <- g[apply(g, 1, function(x) min(x) != max(x)),]
    message(
        'Reference set is based on ', dim(g)[1], ' differential probes from ',
        dim(g)[2], ' cell types.')
    g
}

#' Retrieve reference set
#'
#' The function retrieves the curated reference DNA methylation status for
#' a set of cell type names under the Infinium platform. Supported cell types
#' include "CD4T", "CD19B", "CD56NK", "CD14Monocytes", "granulocytes", "scFat",
#' "skin" etc. See package sesameData for more details. The function output a
#' matrix with probes on the rows and specified cell types on the columns.
#' 0 suggests unmethylation and 1 suggests methylation. Intermediate
#' methylation and nonclusive calls are left with NA. 
#'
#' @param cells reference cell types
#' @param platform EPIC or HM450
#' @return g, a 0/1 matrix with probes on the rows and specified cell types
#' on the columns.
#' @examples
#'
#' betas = getRefSet('CD4T', platform='HM450')
#' sesameDataGet_resetEnv()
#' 
#' @export
getRefSet <- function(cells=NULL, platform = c('EPIC','HM450')) {
    
    platform <- match.arg(platform)
    if (is.null(cells)) {
        cells <- c('CD4T', 'CD19B','CD56NK','CD14Monocytes', 'granulocytes');
    }

    refdata <- sesameDataGet('ref.methylation')[cells]
    probes <- Reduce(intersect, lapply(refdata, names))
    g <- do.call(cbind, lapply(refdata, function(x) x[probes]))
    g <- cleanRefSet(g, platform);
    message(
        'Reference set is based on ', dim(g)[1], ' probes from ',
        dim(g)[2], ' cell types.')
    g
}

errFunc <- function(f, g, q) {
    gamma <- q - g %*% f[2:length(f)]
    sum(ifelse(gamma < f[1] / 2, abs(gamma), abs(gamma - f[1])), na.rm=TRUE)
}

## transform fraction (f) by altering 2 components (nu1 and nu2) by step.size
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

## @param frac0 initial fraction
.optimizeCellComposition <- function(
    g, q, frac0=NULL, temp=0.5, maxIter=1000,
    delta=0.0001, step.max=1.0, verbose=FALSE) {

    M <- ncol(g) # number of reference
    if (is.null(frac0)) { frac <- c(1, rep(0, M))
    } else { frac <- frac0 } # use given fraction estimate
    
    ## initialize
    errcurrent <- errFunc(frac, g, q); errmin <- errcurrent
    frac.min <- frac; niter <- 1
    repeat {
        nu <- sample(seq_len(M+1), 2)
        step.size <- runif(1) * step.max
        frac.test <- double.transform.f(frac, nu[1], nu[2], step.size);
        if (!is.null(frac.test)) {
            if (verbose) {
                message('errcurrent=', errcurrent, 'frac=',
                    paste(lapply(frac, function(x) sprintf('%1.2f', x)),
                        collapse='-'),
                    ';stepsize=', step.size, ';temp=', temp, ';best=',
                    paste(lapply(frac.min,
                        function(x) sprintf('%1.2f', x)), collapse='-'),
                    ';err=', errmin) }

            errtest <- errFunc(frac.test, g, q)
            if (errtest < errmin) { # update best
                errmin <- errtest
                frac.min <- frac.test
                if ((errmin-errtest) > errmin*delta)
                    niter <- 1;
            } else {
                niter <- niter + 1
                if (niter > maxIter) {
                    break;
                }}
            
            ## rejection sampling
            if (runif(1) < exp(-(errtest-errcurrent)/temp)) {
                errcurrent <- errtest
                frac <- frac.test
            }}}

    list(frac.min=frac.min, errmin=errmin)
}

#' Estimate cell composition using reference
#'
#' This is a reference-based cell composition estimation. The function takes a
#' reference methylation status matrix (rows for probes and columns for cell
#' types, can be obtained by getRefSet function) and a query beta value
#' measurement. The length of the target beta values should be the same as
#' the number of rows of the reference matrix. The method assumes one unknown
#' component. It outputs a list containing the estimated cell fraction, the
#' error of optimization and methylation status of the unknown component.
#'
#' @param g reference methylation
#' @param q target measurement: length(q) == nrow(g)
#' @param dichotomize to dichotomize query beta value before estimate,
#' this relieves unclean background subtraction
#' @param refine to refine estimate, takes longer
#' @param ... extra parameters for optimization, this includes
#' temp - annealing temperature (0.5)
#' maxIter - maximum iteration to stop after converge (1000)
#' delta - delta score to reset counter (0.0001)
#' verbose - output debug info (FALSE)
#' @return a list of fraction, min error and unknown component methylation state
## @examples
## g <- diffRefSet(getRefSet(platform='HM450'))
## M <- ncol(g)
## trueFrac <- runif(M+1)
## trueFrac <- trueFrac / sum(trueFrac)
## g0 <- sample(c(0,1), nrow(g), replace=TRUE)
## q <- cbind(g0, g) %*% trueFrac + rnorm(length(g0), mean=0, sd = 0.0)
## q[q<0] <- 0
## q[q>1] <- 1
## est <- estimateCellComposition(g, q)
## @export - TODO fix colinearity
estimateCellComposition <- function(
    g, q, refine=TRUE, dichotomize=FALSE, ...) {

    if (dichotomize) {
        q <- dichotomize(q);
    }

    ## raw
    res <- .optimizeCellComposition(g, q, step.max=1, ...);
    ## refine
    res <- .optimizeCellComposition(
        g, q, frac0=res$frac.min, step.max=0.05, ...)

    list(
        frac = setNames(res$frac.min, c("unknown", colnames(g))),
        err = res$errmin,
        g0 = getg0(res$frac.min, g, q))
}


#' Estimate leukocyte fraction using a two-component model
#'
#' The method assumes only two components in the mixture: the leukocyte
#' component and the target tissue component. The function takes the beta
#' values matrix of the target tissue and the beta value matrix of the
#' leukocyte. Both matrices have probes on the row and samples on the column.
#' Row names should have probe IDs from the platform. The function outputs
#' a single numeric describing the fraction of leukocyte.
#'
#' @param betas.tissue tissue beta value matrix (#probes X #samples)
#' @param betas.leuko leukocyte beta value matrix,
#' if missing, use the SeSAMe default by infinium platform
#' @param betas.tumor optional, tumor beta value matrix
#' @param platform "HM450", "HM27" or "EPIC"
#' @return leukocyte estimate, a numeric vector
#' @importFrom utils head
#' @importFrom utils tail
#' @examples
#'
#' betas.tissue <- sesameDataGet('HM450.1.TCGA.PAAD')$betas
#' estimateLeukocyte(betas.tissue)
#' sesameDataGet_resetEnv()
#' 
#' @export
estimateLeukocyte<-function(betas.tissue, betas.leuko = NULL,
    betas.tumor = NULL, platform = c('EPIC','HM450','HM27')){
    
    platform <- match.arg(platform)
    if (!is.matrix(betas.tissue)) { betas.tissue <- as.matrix(betas.tissue) }
    if (is.null(betas.leuko)) {
        betas.leuko <- sesameDataGet('leukocyte.betas')[[platform]] }
    if (!is.matrix(betas.leuko)) { betas.leuko <- as.matrix(betas.leuko) }

    ## choose the probes to work with
    ave.leuko <- rowMeans(betas.leuko, na.rm=TRUE)
    ave.tissue <- rowMeans(betas.tissue, na.rm=TRUE)
    probes <- intersect(names(ave.leuko), names(ave.tissue))
    ave.leuko <- ave.leuko[probes]; ave.tissue <- ave.tissue[probes]

    if (toupper(platform) %in% c("HM450", "EPIC")) { nprobes <- 1000
    } else if (toupper(platform) == "HM27") { nprobes <- 100; }
    
    tt <- sort(ave.leuko - ave.tissue)
    probes.leuko.lo <- names(head(tt, n=nprobes)) # leuk-specific hypo
    probes.leuko.hi <- names(tail(tt, n=nprobes)) # leuk-specific hyper

    if (!is.null(betas.tumor)) { # test tumor if specified
        if (!is.matrix(betas.tumor))
            betas.tumor <- as.matrix(betas.tumor)
        betas.tissue <- betas.tumor }
        
    ## more refined range if dataset is big
    if (dim(betas.tissue)[2] >= 10) { # exclude NAs
        t.hi <- apply(betas.tissue[probes.leuko.hi,],1,min, na.rm=TRUE)
        t.lo <- apply(betas.tissue[probes.leuko.lo,],1,max, na.rm=TRUE)
    } else {
        t.hi <- rep(0.0, length(probes.leuko.hi))
        t.lo <- rep(1.0, length(probes.leuko.lo)) }
    l.hi <- as.numeric(as.matrix(ave.leuko[probes.leuko.hi]))
    l.lo <- as.numeric(as.matrix(ave.leuko[probes.leuko.lo]))

    ## calculating Leukocyte Percentage
    leuko.estimate <- vapply(seq_len(ncol(betas.tissue)), function(i) {
        s.hi <- betas.tissue[probes.leuko.hi, i]
        s.lo <- betas.tissue[probes.leuko.lo, i]
        p <- c((s.hi - t.hi) / (l.hi - t.hi), (s.lo - t.lo) / (l.lo - t.lo))
        if (sum(!is.na(p)) < 10) return(NA); # not enough data
        dd <- density(na.omit(p))
        dd$x[which.max(dd$y)]
    }, numeric(1))
    names(leuko.estimate) <- colnames(betas.tissue)
    leuko.estimate
}

twoCompsDiff <- function(pop1, pop2) {
    pb <- intersect(rownames(pop1), rownames(pop2))
    pop1 <- pop1[pb,]
    pop2 <- pop2[pb,]
    tt <- sort(rowMeans(pop1) - rowMeans(pop2))
    res <- list(
        diff_1m2u = names(tail(tt, n=1000)),
        diff_1u2m = names(head(tt, n=1000)))
    res
}


#' Estimate the fraction of the 2nd component in a 2-component mixture
#' 
#' @param pop1 Reference methylation level matrix for population 1
#' @param pop2 Reference methylation level matrix for population 2
#' @param target Target methylation level matrix to be analyzed
#' @param diff_1m2u A vector of differentially methylated probes (methylated
#' in population 1 but unmethylated in population 2)
#' @param diff_1u2m A vector of differentially methylated probes (unmethylated
#' in population 1 but methylated in population 2)
#' @param use.ave use population average in selecting differentially
#' methylated probes
#' @return Estimate of the 2nd component in the 2-component mixture
twoCompsEst2 <- function(
    pop1, pop2, target, use.ave=TRUE, diff_1m2u=NULL, diff_1u2m=NULL) {

    pb <- intersect(
        intersect(rownames(pop1), rownames(pop2)),
        rownames(target))
    message(length(pb), " probes shared. Starting from there.\n")
    pop1 <- pop1[pb,]; pop2 <- pop2[pb,]; target <- target[pb,]

    if (is.null(diff_1m2u) || is.null(diff_1u2m)) {
        if (use.ave) {
            tt <- sort(rowMeans(pop1) - rowMeans(pop2))
            diff_1m2u <- names(tail(tt, n=1000))
            diff_1u2m <- names(head(tt, n=1000))
        } else {
            diff_1u2m <- names(which(apply(pop1,1,function(x) {
                all(x<0.3, na.rm=TRUE) && sum(is.na(x)) / length(x) < 0.5
            }) & apply(pop2,1,function(x) {
                all(x>0.7, na.rm=TRUE) && sum(is.na(x)) / length(x) < 0.5
            })))
            
            diff_1m2u <- names(which(apply(pop1,1,function(x) {
                all(x>0.7, na.rm=TRUE) && sum(is.na(x)) / length(x) < 0.5
            }) & apply(pop2,1,function(x) {
                all(x<0.3, na.rm=TRUE) && sum(is.na(x)) / length(x) < 0.5
            })))
        }}
    
    message(length(diff_1u2m), " probes meth. in 2 and unmeth. in 1.\n")
    message(length(diff_1m2u), " probes meth. in 1 and unmeth. in 2.\n")

    d1u2m_hi <- apply(pop2[diff_1u2m,], 1, max, na.rm=TRUE)
    d1u2m_lo <- apply(pop1[diff_1u2m,], 1, min, na.rm=TRUE)
    d1m2u_hi <- apply(pop1[diff_1m2u,], 1, max, na.rm=TRUE)
    d1m2u_lo <- apply(pop2[diff_1m2u,], 1, min, na.rm=TRUE)

    est <- vapply(seq_len(ncol(target)), function(i) {
        xx <- c((target[diff_1u2m,i] - d1u2m_lo) / (d1u2m_hi - d1u2m_lo),
            1-(target[diff_1m2u,i] - d1m2u_lo) / (d1m2u_hi - d1m2u_lo))
        dd <- density(na.omit(xx))
        dd$x[which.max(dd$y)]
    }, numeric(1))
    names(est) <- colnames(target)
    est
}

