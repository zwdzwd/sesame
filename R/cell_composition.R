## clean reference set to non NA sites
cleanRefSet <- function(g, platform = c('EPIC','HM450','HM27')) {

    platform <- match.arg(platform)
    pkgTest('GenomicRanges')
    mapinfo <- get(paste0(platform,'.mapped.probes.hg19'))
    g <- g[intersect(rownames(g), names(mapinfo)),,drop=FALSE]
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
#' g <- diffRefSet(getRefSet(platform='HM450'))
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
#' betas <- getRefSet('CD4T', platform='HM450')
#' @export
getRefSet <- function(cells=NULL, platform = c('EPIC','HM450')) {
    
    platform <- match.arg(platform)
    if (is.null(cells)) {
        cells <- c('CD4T', 'CD19B','CD56NK','CD14Monocytes', 'granulocytes');
    }

    refdata <- mget(paste0('cellref.', cells), inherits=TRUE)
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
    if (is.null(frac0)) {
        frac <- c(1, rep(0, M))
    } else { # use given fraction estimate
        frac <- frac0
    }
    
    ## initialize
    errcurrent <- errFunc(frac, g, q)
    errmin <- errcurrent
    frac.min <- frac
    niter <- 1

    ## initialize
    repeat {
        nu <- sample(seq_len(M+1), 2)
        step.size <- runif(1) * step.max
        frac.test <- double.transform.f(frac, nu[1], nu[2], step.size);
        if (!is.null(frac.test)) {

            if (verbose) {
                message(
                    'errcurrent=', errcurrent, 'frac=',
                    paste(lapply(
                        frac,
                        function(x) sprintf('%1.2f', x)), collapse='-'),
                    ';stepsize=', step.size, ';temp=', temp, ';best=',
                    paste(lapply(
                        frac.min,
                        function(x) sprintf('%1.2f', x)), collapse='-'),
                    ';err=', errmin)
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
#' @examples
#' g <- diffRefSet(getRefSet(platform='HM450'))
#' M <- ncol(g)
#' trueFrac <- runif(M+1)
#' trueFrac <- trueFrac / sum(trueFrac)
#' g0 <- sample(c(0,1), nrow(g), replace=TRUE)
#' q <- cbind(g0, g) %*% trueFrac + rnorm(length(g0), mean=0, sd = 0.0)
#' q[q<0] <- 0
#' q[q>1] <- 1
#' est <- estimateCellComposition(g, q)
#' @export
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
#' betas.tissue <- readRDS(system.file(
#'     'extdata','HM450.betas.TCGA-2L-AAQA-01A-21D-A38H-05.rds',
#'     package='sesameData'))
#' estimateLeukocyte(betas.tissue)
#' 
#' @export
estimateLeukocyte<-function(
    betas.tissue, betas.leuko = NULL,
    betas.tumor = NULL,
    platform = c('EPIC','HM450','HM27')){
    
    platform <- match.arg(platform)

    if (!is.matrix(betas.tissue))
        betas.tissue <- as.matrix(betas.tissue)

    if (is.null(betas.leuko)) {
        betas.leuko <- get(paste0(platform,'.betas.leuko.whole'))
    }

    if (!is.matrix(betas.leuko))
        betas.leuko <- as.matrix(betas.leuko)

    ## choose the probes to work with
    ave.leuko <- rowMeans(betas.leuko, na.rm=TRUE)
    ave.tissue <- rowMeans(betas.tissue, na.rm=TRUE)

    probes <- intersect(names(ave.leuko), names(ave.tissue))
    ave.leuko <- ave.leuko[probes]
    ave.tissue <- ave.tissue[probes]

    if (toupper(platform) %in% c("HM450", "EPIC"))
        nprobes <- 1000
    else if (toupper(platform) == "HM27")
        nprobes <- 100
    else
        stop("Platform not supported.")
    
    tt <- sort(ave.leuko - ave.tissue)
    ## leukocyte-specific hypo-methylation
    probes.leuko.low <- names(head(tt, n=nprobes))
    ## leukocyte-specific hyper-methylation
    probes.leuko.high <- names(tail(tt, n=nprobes))

    ## test tumor if specified
    if (!is.null(betas.tumor)) {
        if (!is.matrix(betas.tumor))
            betas.tumor <- as.matrix(betas.tumor)
        betas.tissue <- betas.tumor
    }
        
    ## more refined range if dataset is big
    if (dim(betas.tissue)[2] >= 10) {
        ## the old code doesn't exclude NAs, this should be better
        t.high <- apply(betas.tissue[probes.leuko.high,],1,min, na.rm=TRUE)
        t.low <- apply(betas.tissue[probes.leuko.low,],1,max, na.rm=TRUE)
    } else {
        t.high <- rep(0.0, length(probes.leuko.high))
        t.low <- rep(1.0, length(probes.leuko.low))
    }
    l.high <- as.numeric(as.matrix(ave.leuko[probes.leuko.high]))
    l.low <- as.numeric(as.matrix(ave.leuko[probes.leuko.low]))

    ## calculating Leukocyte Percentage
    leuko.estimate <- vapply(seq_len(ncol(betas.tissue)), function(i) {
        s.high <- betas.tissue[probes.leuko.high, i]
        s.low <- betas.tissue[probes.leuko.low, i]
        p <- c(
            (s.high - t.high) / (l.high - t.high),
            (s.low - t.low) / (l.low - t.low))

        if (sum(!is.na(p)) < 10) return(NA); # not enough data
        dd <- density(na.omit(p))
        dd$x[which.max(dd$y)]
    }, numeric(1))
    
    names(leuko.estimate) <- colnames(betas.tissue)
    leuko.estimate
}
