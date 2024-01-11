
#' filter data matrix by factor completeness
#' only works for discrete factors
#'
#' @param betas matrix data
#' @param fc factors, or characters
#' @return a boolean vector whether there is non-NA value for each tested
#' group for each probe
#' @examples
#' se0 <- sesameDataGet("MM285.10.SE.tissue")[1:100,]
#' se_ok <- checkLevels(SummarizedExperiment::assay(se0),
#'     SummarizedExperiment::colData(se0)$tissue)
#' sum(se_ok) # number of good probes
#' se1 <- se0[se_ok,]
#'
#' sesameDataGet_resetEnv()
#' @export
checkLevels <- function(betas, fc) {
    stopifnot(is(fc, "factor") || is(fc, "character"))
    apply(betas, 1, function(dt) {
        all(vapply(split(dt, fc), function(x) sum(!is.na(x))>0, logical(1)))
    })
}

model_contrasts <- function(mm, meta) {
    contrs <- names(attr(mm, "contrasts"))
    setNames(lapply(contrs, function(cont) {
        ## avoid X-prepended to levels that start with number
        x <- make.names(paste0("X",levels(factor(meta[[cont]]))))
        substr(x,2,nchar(x)) # remove the added X
    }), contrs)
}

#' Test differential methylation on each locus
#'
#' The function takes a beta value matrix with probes on the rows and
#' samples on the columns. It also takes a sample information data frame
#' (meta) and formula for testing. The function outputs a list of
#' coefficient tables for each factor tested.
#' @param betas beta values, matrix or SummarizedExperiment
#' rows are probes and columns are samples.
#' @param fm formula
#' @param meta data frame for sample information, column names
#' are predictor variables (e.g., sex, age, treatment, tumor/normal etc)
#' and are referenced in formula. Rows are samples.
#' When the betas argument is a SummarizedExperiment object, this
#' is ignored. colData(betas) will be used instead. The row order of the
#' data frame must match the column order of the beta value matrix.
#' @param BPPARAM number of cores for parallel processing, default to
#' SerialParam()
#' Use MulticoreParam(mc.cores) for parallel processing.
#' For Windows, try DoparParam or SnowParam.
#' @return a list of test summaries, summary.lm objects
#' @import stats
#' @import BiocParallel
#' @importFrom SummarizedExperiment assay
#' @importFrom SummarizedExperiment colData
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @examples
#' sesameDataCache() # in case not done yet
#' data <- sesameDataGet('HM450.76.TCGA.matched')
#' smry <- DML(data$betas[1:1000,], ~type, meta=data$sampleInfo)
#'
#' sesameDataGet_resetEnv()
#' @export
DML <- function(betas, fm, meta=NULL, BPPARAM=SerialParam()) {

    if(is(betas, "SummarizedExperiment")) {
        betas0 <- betas
        betas <- assay(betas0)
        meta <- colData(betas0)
    }
    stopifnot(nrow(meta) == ncol(betas))

    mm <- model.matrix(fm, meta)
    colnames(mm) <- make.names(colnames(mm))
    
    ## prepare holdout models
    contr2lvs <- model_contrasts(mm, meta)
    mm_holdout <- lapply(names(contr2lvs), function(cont) {
        mm[, !(colnames(mm) %in% paste0(cont, contr2lvs[[cont]]))] })
    names(mm_holdout) <- names(contr2lvs)

    ## fitting
    smry <- BiocParallel::bplapply(seq_len(nrow(betas)), function(i) {
        m0 <- lm(betas[i,]~.+0, data=as.data.frame(mm))
        sm <- summary(m0)
        sm$cov.unscaled <- NULL # reduce size of return
        sm$residuals <- NULL    # reduce size of return
        sm$terms <- NULL        # reduce size of return
        sm$Ftest <- do.call(cbind, lapply(mm_holdout, function(mm_) {
            m1 <- lm(betas[i,]~.+0, data=as.data.frame(mm_))
            anv <- anova(m1, m0)
            c(stat = anv[["F"]][2], pval = anv[["Pr(>F)"]][2])
        }))
        sm
    }, BPPARAM = BPPARAM)
    names(smry) <- rownames(betas)
    class(smry) <- "DMLSummary"
    attr(smry, "model.matrix") <- mm
    attr(smry, "fm") <- fm
    attr(smry, "contr2lvs") <- contr2lvs
    smry
}

#' Predict new data from DML
#'
#' This function is also important for investigating factor interactions.
#' 
#' @param betas beta values, matrix or SummarizedExperiment
#' rows are probes and columns are samples.
#' @param fm formula
#' @param pred new data for prediction, useful for studying effect size.
#' This argument is a data.frame to specify new data.
#' If the argument is NULL, all combinations of all contrasts will be used
#' as input. It might not work if there is a continuous variable input.
#' One may need to explicitly provide the input in a data frame.
#' @param meta data frame for sample information, column names
#' are predictor variables (e.g., sex, age, treatment, tumor/normal etc)
#' and are referenced in formula. Rows are samples.
#' When the betas argument is a SummarizedExperiment object, this
#' is ignored. colData(betas) will be used instead.
#' @param BPPARAM number of cores for parallel processing, default to
#' SerialParam()
#' Use MulticoreParam(mc.cores) for parallel processing.
#' For Windows, try DoparParam or SnowParam.
#' @return a SummarizedExperiment of predictions. The colData describes
#' the input of the prediction.
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @examples
#' data <- sesameDataGet('HM450.76.TCGA.matched')
#' 
#' ## use all contrasts as new input
#' res <- DMLpredict(data$betas[1:10,], ~type, meta=data$sampleInfo)
#'
#' ## specify new input
#' res <- DMLpredict(data$betas[1:10,], ~type, meta=data$sampleInfo,
#'   pred = data.frame(type=c("Normal","Tumour")))
#'
#' ## note that the prediction needs to be a factor of the same
#' ## level structure as the original training data.
#' pred = data.frame(type=factor(c("Normal"), levels=c("Normal","Tumour")))
#' res <- DMLpredict(data$betas[1:10,], ~type,
#'   meta=data$sampleInfo, pred = pred)
#'
#' @export
DMLpredict <- function(betas, fm, pred = NULL, meta = NULL,
    BPPARAM=SerialParam()) {

    if(is(betas, "SummarizedExperiment")) {
        betas0 <- betas
        betas <- assay(betas0)
        meta <- colData(betas0)
    }

    mm <- model.matrix(fm, meta)
    colnames(mm) <- make.names(colnames(mm))
    
    if (is.null(pred)) { # use all combinations of all contrasts
        contr2lvs <- model_contrasts(mm, meta)
        ## NOTE: conversion to factor is done by occurrence,
        ## not alphabetically. This is what we need. See ?expand.grid
        pred <- do.call(expand.grid, contr2lvs)
    }

    ## prepare prediction input
    mm_pred <- as.data.frame(model.matrix(fm, pred))
    colnames(mm_pred) <- make.names(colnames(mm_pred))
    stopifnot(all(colnames(mm_pred) %in% colnames(mm)))

    ## fitting and prediction
    res <- do.call(rbind, BiocParallel::bplapply(
        seq_len(nrow(betas)), function(i) {
            m0 <- lm(betas[i,]~.+0, data=as.data.frame(mm))
            predict(m0, mm_pred)
        }, BPPARAM = BPPARAM))
    rownames(res) <- rownames(betas)
    SummarizedExperiment(res, colData=pred)
}

#' Print DMLSummary object
#'
#' @param x a DMLSummary object
#' @param ... extra parameter for print
#' @return print DMLSummary result on screen
#' @examples
#' sesameDataCache() # in case not done yet
#' data <- sesameDataGet('HM450.76.TCGA.matched')
#' ## test the first 10
#' smry <- DML(data$betas[1:10,], ~type, meta=data$sampleInfo)
#' smry
#'
#' sesameDataGet_resetEnv()
#' @export
print.DMLSummary <- function(x, ...) {
    mm <- attr(x, "model.matrix")
    message(sprintf("DMLSummary Object with %d Loci, %d samples.\n",
        length(x), nrow(mm)))
    contrast_names <- paste0(names(attr(mm, "contrasts")), collapse=", ")
    message("Contrasts: ", contrast_names, "\n")
}

#' Extract slope information from DMLSummary
#' @param smry DMLSummary from DML command
#' @return a table of slope and p-value
#' @importFrom dplyr bind_rows
#' @importFrom dplyr bind_cols
#' @examples
#' sesameDataCache() # in case not done yet
#' data <- sesameDataGet('HM450.76.TCGA.matched')
#' smry <- DML(data$betas[1:10,], ~type, meta=data$sampleInfo)
#' slopes <- summaryExtractTest(smry)
#' 
#' sesameDataGet_resetEnv()
#' @export
summaryExtractTest <- function(smry) {

    contr2lvs <- attr(smry, "contr2lvs")
    ## exclude the sites with NAs, maybe we we can do better?
    smrylen <- vapply(smry, function(x) { nrow(x$coefficients); }, numeric(1))
    smry <- smry[smrylen == max(smrylen)]

    est <- do.call(bind_rows, lapply(smry, function(x) {
        x$coefficients[,'Estimate']; }))

    colnames(est) <- paste0("Est_", colnames(est))
    pvals <- do.call(bind_rows, lapply(smry, function(x) {
        x$coefficients[,"Pr(>|t|)"] }))
    colnames(pvals) <- paste0("Pval_", colnames(pvals))
    if(is.null(smry[[1]]$Ftest)) { # only continuous
        return(cbind(Probe_ID=names(smry),est,pvals))
    }
    f_pvals <- do.call(rbind, lapply(smry, function(x) {
        x$Ftest["pval",,drop=FALSE] }))
    colnames(f_pvals) <- paste0("FPval_", colnames(f_pvals))

    contr2lvs <- contr2lvs[vapply(
        contr2lvs, function(x) nchar(x[[1]]), numeric(1)) > 0]
    if (length(contr2lvs)>0) {
        ## this doesn't account for interaction terms
        effsize <- do.call(cbind, lapply(names(contr2lvs), function(cont) {
            lvs <- contr2lvs[[cont]]
            lvs <- lvs[2:length(lvs)]
            lvs <- lvs[paste0("Est_", cont, lvs) %in% colnames(est)]
            apply(est[, paste0("Est_", cont, lvs),drop=FALSE], 1, function(x) {
                max(x,0) - min(x,0) }) }))
        
        colnames(effsize) <- paste0("Eff_", names(contr2lvs))
    } else { effsize <- NULL; }
    bind_cols(Probe_ID=names(smry), est, pvals, f_pvals, effsize)
}

#' Compute effect size for different variables from prediction matrix
#'
#' The effect size is defined by the maximum variation of a variable with all
#' the other variables controled constant.
#'
#' @param pred predictions
#' @return a data.frame of effect sizes. Columns are different variables.
#' Rows are different probes.
#' @examples
#' data <- sesameDataGet('HM450.76.TCGA.matched')
#' res <- DMLpredict(data$betas[1:10,], ~type, meta=data$sampleInfo)
#' head(calcEffectSize(res))
#' @export
calcEffectSize <- function(pred) {
    vars <- colnames(colData(pred))
    if (length(vars) == 1) {
        eff <- data.frame(x = apply(
            assay(pred), 1, function(x) max(x) - min(x)))
        colnames(eff) <- vars[[1]]
        rownames(eff) <- rownames(pred)
        return(eff)
    }
    eff <- as.data.frame(do.call(cbind, lapply(vars, function(var) {
        other_vars <- vars[vars != var]
        col_indices <- seq_len(nrow(colData(pred)))
        Reduce(pmax, lapply(split(col_indices, colData(pred)[other_vars]),
            function(x) {
                apply(assay(pred)[,x], 1, function(x) max(x) - min(x))
            }))
    })))
    colnames(eff) <- vars
    rownames(eff) <- rownames(pred)
    eff
}

summaryExtractCf <- function(smry, contrast) {
    cf <- do.call(rbind, lapply(smry, function(x) {
        if (x$aliased[contrast]) { # missing fitting
            NA
        } else {
            x$coefficients[contrast,]
        }
    }))
    rownames(cf) <- names(smry)
    cf # probes x c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
}

dmr_merge_cpgs <- function(betas, probe.coords, dist.cutoff, seg.per.locus) {

    betas.noNA <- betas[!apply(betas, 1, function(x) all(is.na(x))),]
    cpg.ids <- intersect(rownames(betas.noNA), names(probe.coords))
    probe.coords <- GenomicRanges::sort(probe.coords[cpg.ids])
    betas.coord.srt <- betas.noNA[names(probe.coords),]
    
    cpg.ids <- rownames(betas.coord.srt)
    cpg.coords <- probe.coords[cpg.ids]
    cpg.chrm <- as.vector(GenomicRanges::seqnames(cpg.coords))
    cpg.start <- GenomicRanges::start(cpg.coords)
    cpg.end <- GenomicRanges::end(cpg.coords)

    n.cpg <- length(cpg.ids)

    ## euclidean distance suffcies
    beta.dist <- vapply(seq_len(n.cpg-1), function(i) sum(
        (betas.coord.srt[i,] - betas.coord.srt[i+1,])^2, na.rm=TRUE), 1)

    chrm.changed <- (cpg.chrm[-1] != cpg.chrm[-n.cpg])

    ## empirical cutoff based on quantiles
    if (is.null(dist.cutoff)) {
        dist.cutoff <- quantile(beta.dist, 1-seg.per.locus)
    }

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
    seg.ids <- c(seg.ids, seq.int(
        from=seg.ids[length(seg.ids)]+1, length.out=length(unmapped)))
    
    seg.chrm <- as.vector(tapply(cpg.chrm, seg.ids, function(x) x[1]))
    seg.start <- as.vector(tapply(cpg.start, seg.ids, function(x) x[1]))
    seg.end <- as.vector(tapply(cpg.end, seg.ids, function(x) x[length(x)]))

    list(id = seg.ids, chrm = seg.chrm,
        start = seg.start, end = seg.end, cpg.ids = cpg.ids)
}

dmr_combine_pval <- function(cf, segs) {
    ## mean of Estimate
    seg.est <- as.vector(tapply(
        cf[segs$cpg.ids, 'Estimate'], segs$id,
        function(x) mean(x, na.rm=TRUE) ))
    ## Stouffer's Z-score method
    seg.pval <- as.vector(tapply(
        cf[segs$cpg.ids, 'Pr(>|t|)'], segs$id,
        function(x) pnorm(sum(qnorm(x))/sqrt(length(x)))))
    seg.pval.adj <- p.adjust(seg.pval, method="BH")
    seg.ids.cf <- match(rownames(cf), segs$cpg.ids)
    s <- segs$id[seg.ids.cf]
    cf <- cbind(data.frame(
        Seg_ID = s,
        Seg_Chrm = segs$chrm[s],
        Seg_Start = segs$start[s],
        Seg_End = segs$end[s],
        Seg_Est = seg.est[s],
        Seg_Pval = seg.pval[s],
        Seg_Pval_adj = seg.pval.adj[s],
        Probe_ID = rownames(cf)
    ), as.data.frame(cf))
    message(sprintf(
        ' - %d significant segments.',
        sum(seg.pval<0.05, na.rm=TRUE)))
    message(sprintf(
        ' - %d significant segments (after BH).',
        sum(seg.pval.adj<0.05, na.rm=TRUE)))
    cf <- cf[order(cf$Seg_Est, cf$Seg_Chrm, cf$Seg_Start),]
    rownames(cf) <- NULL
    cf
}

DMGetProbeInfo <- function(platform, genome) {
    ## mft <- sesameDataGet(sprintf("%s.%s.manifest", platform, genome))
    mft <- sesameData_getManifestGRanges(platform, genome = genome)
    mft <- mft[GenomicRanges::seqnames(mft) != "*"]
    GenomicRanges::mcols(mft) <- NULL
    GenomicRanges::strand(mft) <- "*"
    mft <- sort(mft)
    mft
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
#' @param smry DML
#' @param contrast the pair-wise comparison or contrast
#' check colnames(attr(smry, "model.matrix")) if uncertain
#' @param dist.cutoff cutoff of beta value differences for two neighboring CGs
#' to be considered the same DMR (by default it's determined using the
#' quantile function on seg.per.locus)
#' @param seg.per.locus number of segments per locus
#' higher value leads to more segments
#' @param platform EPIC, HM450, MM285, ...
#' @param probe.coords GRanges object that defines CG coordinates
#' if NULL (default), then the default genome assembly is used.
#' Default genome is given by, e.g., sesameData_check_genome(NULL, "EPIC")
#' For additional mapping, download the GRanges object from
#' http://zwdzwd.github.io/InfiniumAnnotation
#' and provide the following argument
#' ..., probe.coords = sesameAnno_buildManifestGRanges("downloaded_file"),...
#' to this function.
#' @return coefficient table with segment ID and segment P-value
#' each row is a locus, multiple loci may share a segment ID if
#' they are merged to the same segment. Records are ordered by Seg_Est.
#' @importFrom SummarizedExperiment assay
#' @importFrom SummarizedExperiment colData
#' @examples
#'
#' sesameDataCache() # in case not done yet
#' 
#' data <- sesameDataGet('HM450.76.TCGA.matched')
#' smry <- DML(data$betas[1:1000,], ~type, meta=data$sampleInfo)
#' colnames(attr(smry, "model.matrix")) # pick a contrast from here
#' ## showing on a small set of 100 CGs
#' merged_segs <- DMR(data$betas[1:1000,], smry, "typeTumour", platform="HM450")
#'
#' sesameDataGet_resetEnv()
#' 
#' @export
DMR <- function(betas, smry, contrast,
    platform=NULL, probe.coords=NULL, dist.cutoff=NULL, seg.per.locus=0.5) {

    stopifnot(is(smry, "DMLSummary"))
    if(is(betas, "SummarizedExperiment")) {
        betas <- assay(betas)
    }

    ## sort by coordinates
    if (is.null(probe.coords)) {
        if (is.null(platform)) {
            platform <- inferPlatformFromProbeIDs(rownames(betas))
        }
        genome <- sesameData_check_genome(NULL, platform)
        probe.coords <- DMGetProbeInfo(platform, genome)
    }
    message("Merging correlated CpGs ... ", appendLF=FALSE)
    segs <- dmr_merge_cpgs(betas, probe.coords, dist.cutoff, seg.per.locus)
    message(sprintf('Generated %d segments.', segs$id[length(segs$id)]))
    message("Combine p-values ... ")
    cf <- summaryExtractCf(smry, contrast)
    ## make sure the beta values and coefficients match
    stopifnot(all(segs$cpg.ids %in% rownames(cf)))
    cf <- dmr_combine_pval(cf, segs)
    message("Done.")
    cf
}


#' List all contrasts of a DMLSummary
#' 
#' @param smry a DMLSummary object
#' @return a character vector of contrasts
#' @examples
#' data <- sesameDataGet('HM450.76.TCGA.matched')
#' smry <- DML(data$betas[1:10,], ~type, meta=data$sampleInfo)
#' dmContrasts(smry)
#'
#' sesameDataGet_resetEnv()
#' @export
dmContrasts <- function(smry) {
    stopifnot(is(smry, "DMLSummary"))
    colnames(attr(smry, "model.matrix"))
}
