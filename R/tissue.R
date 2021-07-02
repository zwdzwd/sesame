reference_plot_se = function(betas, se, color=c("blueYellow","fullJet")) {

    pkgTest("wheatmap")
    color = match.arg(color)
    if (color == "blueYellow") stop.points = c("blue","yellow")
    else stop.points = NULL

    cd = as_tibble(colData(se))
    rd = as_tibble(rowData(se))
    md = metadata(se)
    if (!is.null(betas) && is.null(dim(betas))) { # in case a vector
        betas = cbind(betas)
    }

    g = WHeatmap(assay(se), cmp=CMPar(stop.points=stop.points,
        dmin=0, dmax=1), name="b1") # reference
    if (!is.null(betas)) {          # query samples
        g = g + WHeatmap(betas[rownames(se),], RightOf("b1"),
            cmp=CMPar(stop.points=stop.points, dmin=0, dmax=1),
            name="b2", xticklabels=TRUE, xticklabels.n=ncol(betas))
        right = "b2"
    } else { # in case target is not given, plot just the reference
        right = "b1"
    }
    ## branch color bar (vertical)
    g = g + WColorBarV(rd$probebranches, RightOf(right, width=0.03),
        cmp=CMPar(label2color=md$branch_color), name="bh")
    ## tissue color bar (horizontal)
    g = g + WColorBarH(cd$Tissue_Corrected, TopOf("b1",height=0.03),
        cmp=CMPar(label2color=md$tissue_color), name="ti")
    ## legends
    g = g + WLegendV("ti", TopRightOf("bh", just=c('left','top'), h.pad=0.02),
        height=0.02)
    g = g + WLegendV('bh', Beneath(pad=0.06))
    g + WCustomize(mar.bottom=0.15, mar.right=0.06)
}

#' Compare mouse array data with mouse tissue references
#'
#' @param betas matrix of betas for the target sample
#' @param color either blueYellow or fullJet
#' @return grid object that contrast the target sample with
#' pre-built mouse tissue reference
#' @import wheatmap
#' @export
#' @examples
#' sesameDataCache("MM285") # if not done yet
#' b = sesameDataGet("MM285.10.tissue")$betas[,1:2]
#' compareMouseTissueReference(b)
#' @importFrom SummarizedExperiment assay
#' @importFrom SummarizedExperiment colData
#' @importFrom SummarizedExperiment rowData
compareMouseTissueReference = function(betas=NULL, color="blueYellow") {
    se = sesameDataGet("MM285.tissueSignature")
    reference_plot_se(betas, se, color=color)
}

#' Compare beta value against mouse blood reference
#'
#' @param betas matrix of betas for the target sample
#' @param color either blueYellow or fullJet
#' @return grid object that co-plots with a pre-built mouse blood reference
#' @export
#' @examples
#' sesameDataCache("MM285") # if not done yet
#' b = sesameDataGet("MM285.10.tissue")$betas[,10]
#' compareMouseBloodReference(b)
compareMouseBloodReference = function(betas=NULL, color="blueYellow") {
    se = sesameDataGet("MM285.bloodSignature")
    reference_plot_se(betas, se, color=color)
}

#' saveProbes saves probe information (delta_betas, auc, and coverage) as a 
#' tsv file for a given branch at in specified file directory (id_dir).
#'
#' @param branch String corresponding to the name of the branch for which the 
#' probe information will be saved.
#' @param probe_info List of probe information with named vectors delta_betas, 
#' auc, cand overage
#' @param id_dir String corresponding to the file directory to which the probe 
#' information is stored.
#'
#' @return None
#'
#' @import tibble
saveProbes = function(branch, probe_info, id_dir) {
    write.table(tibble(probeID=names(probe_info$delta_betas), delta_beta=probe_info$delta_betas, auc=probe_info$auc, coverage=probe_info$coverage, branchID=branch), file=sprintf("%s/%s.tsv", id_dir, branch), row.names=FALSE, quote=FALSE, sep="\t")
}

#' inferProbesHypo infers top hypomethylated probes according to a delta_max, 
#' auc_min, and coverage_min
#'
#' @param betas Matrix of beta value measurement for each sample (columns) 
#' across all probes of the array (rows)
#' @param branch_grouping Vector specifying the in-group (0), out-group (1), 
#' and ignored (2) samples as they appear in the column names of betas.
#' @param delta_max Numeric value corresponding to the maximum delta beta 
#' value for a probe to be considered.
#' @param auc_min Numeric value corresponding to the minimum AUC value 
#' required for a probe to be considered
#' @param coverage_min Numeric value corresponding to the minimum coverage 
#' requirement for a probe to be considered. Coverage is defined here as the 
#' proportion of samples without an NA value at a given probe.
#' @param scan_delta Local value indicating whether to vary the maximum delta 
#' beta when scanning for possible probes.
#' @param scan_auc Local value indicating whether to vary the minimum AUC 
#' threshold when scanning for possible probes.
#' @param verbose Logical value indicating whether to display the scanned  
#' results
#'
#' @return List containing delta_betas, auc, and coverage for a given 
#' branch and set of samples.
inferProbesHypo = function(betas, branch_grouping, delta_max=-0.6, auc_min=0.99, coverage_min=0.80, scan_delta=FALSE, scan_auc=FALSE, verbose=FALSE) {
    if (length(branch_grouping[branch_grouping==0]) == 1) {
        tmp = data.frame(betas[,branch_grouping==0])
        colnames(tmp) = colnames(betas)[branch_grouping==0]
        m0 = apply(tmp, 1,
               function(xx) mean(tail(sort(xx),n=5), na.rm=TRUE))
    } else {
        m0 = apply(betas[,branch_grouping==0],1,
               function(xx) mean(tail(sort(xx),n=5), na.rm=TRUE)) # .2
    }
    m1 = apply(betas[,branch_grouping==1],1,
               function(xx) mean(head(sort(xx),n=5), na.rm=TRUE)) # .8
    delta_betas = m0 - m1
    auc = apply(betas, 1, function(b1) {
        br = branch_grouping[branch_grouping %in% c(0,1)];
        b1 = b1[branch_grouping %in% c(0,1)];
        auc_wmw2(br, b1);})
    coverage = apply(betas[,branch_grouping==0 | branch_grouping==1],1,
               function(xx) sum(!is.na(xx))/length(xx))

    res = list(delta_betas=delta_betas, auc=auc, coverage=coverage)

    if (verbose) {
        viewInferedProbesHypo(res, delta_max, auc_min, scan_delta, scan_auc)
    }
    res
}

#' viewInferedProbesHypo view the results from the inferred top hypomethylated
#' probes according to a delta_max, auc_min, and coverage_min
#'
#' @param res List containing delta_betas, auc, and coverage for a given 
#' branch and set of samples.
#' @param delta_max Numeric value corresponding to the maximum delta beta 
#' value for a probe to be considered.
#' @param auc_min Numeric value corresponding to the minimum AUC value 
#' required for a probe to be considered
#' @param coverage_min Numeric value corresponding to the minimum coverage 
#' requirement for a probe to be considered. Coverage is defined here as the 
#' proportion of samples without an NA value at a given probe.
#' @param scan_delta Local value indicating whether to vary the maximum delta 
#' beta when scanning for possible probes.
#' @param scan_auc Local value indicating whether to vary the minimum AUC 
#' threshold when scanning for possible probes.
#'
#' @return None
viewInferedProbesHypo = function(res, delta_max=-0.6, auc_min=0.99, coverage_min=0.80, scan_delta=FALSE, scan_auc=FALSE) {
    if (scan_delta) {
        for (i in seq(-0.9,-0.1,by=0.05)) {
            probes = names(which(res$delta_betas <= i & res$auc >= auc_min & res$coverage >= coverage_min))
            message(sprintf("delta: %f Found %d probes", i, length(probes)))
        }
    }
    if (scan_auc) {
        for (i in seq(0.8,0.99,by=0.11)) {
            probes = names(which(res$delta_betas <= delta_max & res$auc >= i & res$coverage >= coverage_min))
            message(sprintf("auc: %f Found %d probes", i, length(probes)))
        }
    }
}

#' inferTopProbesHypo selects the top hypomethylated probes according to a 
#' delta_max, auc_min, and coverage_min with a defined sorting criteria 
# (sort_n or sort_delta_betas).
#'
#' @param res List containing delta_betas, auc, and coverage for a given 
#' branch and set of samples.
#' @param n Numeric value corresponding the number of probes to report (if 
#' sort_n is TRUE)
#' @param delta_max Numeric value corresponding to the maximum delta beta 
#' value for a probe to be considered.
#' @param auc_min Numeric value corresponding to the minimum AUC value 
#' required for a probe to be considered
#' @param coverage_min Numeric value corresponding to the minimum coverage 
#' requirement for a probe to be considered. Coverage is defined here as the 
#' proportion of samples without an NA value at a given probe.
#' @param sort_n Local value indicating whether to report the top n probes
#' @param sort_delta_beta Logical value indicating whether the report the top
#' probes satisfying delta_max.
#'
#' @return List containing delta_betas, auc, and coverage for the top selected
#' probes
inferTopProbesHypo = function(res, n=50, delta_max=-0.3, auc_min=0.99, coverage_min=0.80, sort_n=FALSE, sort_delta_beta=FALSE) {
    indx = which(res$auc >= auc_min & res$coverage >= coverage_min & res$delta_betas <= delta_max)
    if (sort_n) {
        # probes_n = names(head(sort(res$delta_betas[res$auc >= auc_min & res$coverage >= coverage_min]), n=n))
        
        indx_names = names(head(sort(res$delta_betas[indx]), n=n))
        probes_n = list(delta_betas=res$delta_betas[indx_names], auc=res$auc[indx_names], coverage=res$coverage[indx_names])
        if (length(probes_n$delta_betas) == 0){
            message(sprintf("Found 0 probes (max delta: 0)"))
        } else {
            message(sprintf("Found %d probes (max delta: %f)", length(probes_n$delta_betas),  max(probes_n$delta_betas)))
        }
    } 
    if (sort_delta_beta) {
        probes_delta_beta = list(delta_betas=res$delta_betas[indx], auc=res$auc[indx], coverage=res$coverage[indx])
        if (length(probes_delta_beta$delta_betas) == 0){
            message(sprintf("Found 0 probes (max delta: 0)"))
        } else {
            message(sprintf("Found %d probes (max delta: %f)", length(probes_delta_beta$delta_betas),  max(probes_delta_beta$delta_betas)))
        }
    }
    if (sort_n & !sort_delta_beta) {
        return(probes_n)
    }
    if (sort_delta_beta & !sort_n) {
        return(probes_delta_beta)
    }
    if (length(probes_n$delta_betas) > length(probes_delta_beta$delta_betas)) {
        return(probes_n)
    }
    return(probes_delta_beta)
}

#' buildReferenceHypo builds and saves a reference for the top hypomethylated
#' probes from using the data in betas and the meta data provided. The 
#' intermediate files are stored in the specified file directory (id_dir).
#'
#' @param betas Matrix of beta value measurement for each sample (columns) 
#' across all probes of the array (rows)
#' @param meta DataFrame corresponding the the meta of the analysis where each 
#' branch column name contains a vector specifying the in-group (0), out-group 
#' (1), and ignored (2) samples as they appear in the column names of betas.
#' @param id_dir String corresponding to the file directory to which the probe 
#' information is stored.
#' @param verbose Logical value indicating whether to display the scanned  
#' results
#' @param return_results Logical value indicating whether to return the 
#' results from the function
#'
#' @return List containing delta_betas, auc, and coverage for a given 
#' branch and set of samples (if return_results).
#'
#' @export
buildReferenceHypo = function(betas, meta, id_dir, verbose=FALSE, return_results=FALSE) {
    res = list()
    branch2probes_hypo = list()

    branches = names(meta)[5:ncol(meta)]

    for (branch in branches) {
        res[[branch]] = inferProbesHypo(betas, meta[[branch]], scan_delta=T, delta_max=-0.3, verbose=verbose)
        saveRDS(res[[branch]], sprintf("%s/fit_%s.rds", id_dir, branch))

        res[[branch]] = readRDS(sprintf("%s/fit_%s.rds", id_dir, branch))
        branch2probes_HM450_hypo[[branch]] = inferTopProbesHypo(res[[branch]],sort_delta_beta=T, delta_max=-0.3)
        saveProbes(branch, branch2probes_HM450_hypo[[branch]], id_dir=id_dir)
    }
    return(res)
}

#' inferProbesHyper infers top hypermethylated probes according to a 
#' delta_min, auc_min, and coverage_min
#'
#' @param betas Matrix of beta value measurement for each sample (columns) 
#' across all probes of the array (rows)
#' @param branch_grouping Vector specifying the in-group (0), out-group (1), 
#' and ignored (2) samples as they appear in the column names of betas.
#' @param delta_min Numeric value corresponding to the minimum delta beta 
#' value for a probe to be considered.
#' @param auc_min Numeric value corresponding to the minimum AUC value 
#' required for a probe to be considered
#' @param coverage_min Numeric value corresponding to the minimum coverage 
#' requirement for a probe to be considered. Coverage is defined here as the 
#' proportion of samples without an NA value at a given probe.
#' @param scan_delta Local value indicating whether to vary the minimum delta 
#' beta when scanning for possible probes.
#' @param scan_auc Local value indicating whether to vary the minimum AUC 
#' threshold when scanning for possible probes.
#' @param verbose Logical value indicating whether to display the scanned  
#' results
#'
#' @return List containing delta_betas, auc, and coverage for a given 
#' branch and set of samples.
inferProbesHyper = function(betas, branch_grouping, delta_min=0.6, auc_min = 0.99, coverage_min = 0.80, scan_delta = FALSE, scan_auc = FALSE) {
    if (length(branch_grouping[branch_grouping==0]) == 1) {
        tmp = data.frame(betas[,branch_grouping==0])
        colnames(tmp) = colnames(betas)[branch_grouping==0]
        m0 = apply(tmp, 1,
               function(xx) mean(head(sort(xx),n=5), na.rm=TRUE))
    } else {
        m0 = apply(betas[,branch_grouping==0],1,
               function(xx) mean(head(sort(xx),n=5), na.rm=TRUE)) # .2
    }
    m1 = apply(betas[,branch_grouping==1],1,
               function(xx) mean(tail(sort(xx),n=5), na.rm=TRUE)) # .8
    delta_betas = m0 - m1
    auc = apply(betas, 1, function(b1) {
        br = branch_grouping[branch_grouping %in% c(0,1)];
        b1 = b1[branch_grouping %in% c(0,1)];
        auc_wmw2(br, -b1);})
    coverage = apply(betas[,branch_grouping==0 | branch_grouping==1],1,
               function(xx) sum(!is.na(xx))/length(xx))

    res = list(delta_betas = delta_betas, auc = auc, coverage = coverage)

    viewInferedProbesHyper(res, delta_min, auc_min, scan_delta, scan_auc, coverage_min)

    res
}

#' viewInferedProbesHyper view the results from the inferred top 
#' hypermethylated probes according to a delta_min, auc_min, and coverage_min
#'
#' @param res List containing delta_betas, auc, and coverage for a given 
#' branch and set of samples.
#' @param delta_min Numeric value corresponding to the minimum delta beta 
#' value for a probe to be considered.
#' @param auc_min Numeric value corresponding to the minimum AUC value 
#' required for a probe to be considered
#' @param coverage_min Numeric value corresponding to the minimum coverage 
#' requirement for a probe to be considered. Coverage is defined here as the 
#' proportion of samples without an NA value at a given probe.
#' @param scan_delta Local value indicating whether to vary the minimum delta 
#' beta when scanning for possible probes.
#' @param scan_auc Local value indicating whether to vary the minimum AUC 
#' threshold when scanning for possible probes.
#'
#' @return None
viewInferedProbesHyper = function(res, delta_min=0.6, auc_min=0.99, scan_delta=FALSE, scan_auc=FALSE, coverage_min=0.80) {
    if (scan_delta) {
        for (i in seq(0.1, 0.9, by=0.05)) {
            probes = names(which(res$delta_betas >= i & res$auc >= auc_min & res$coverage >= coverage_min))
            message(sprintf("delta: %f Found %d probes", i, length(probes)))
        }
    }
    if (scan_auc) {
        for (i in seq(0.01, 0.20, by=0.01)) {
            probes = names(which(res$delta_betas >= delta_min & res$auc <= i & res$coverage >= coverage_min))
            message(sprintf("auc: %f Found %d probes", i, length(probes)))
        }
    }
}

#' inferTopProbesHyper selects the top hypermethylated probes according to a 
#' delta_min, auc_min, and coverage_min with a defined sorting criteria 
# (sort_n or sort_delta_betas).
#'
#' @param res List containing delta_betas, auc, and coverage for a given 
#' branch and set of samples.
#' @param n Numeric value corresponding the number of probes to report (if 
#' sort_n is TRUE)
#' @param delta_min Numeric value corresponding to the minimum delta beta 
#' value for a probe to be considered.
#' @param auc_min Numeric value corresponding to the minimum AUC value 
#' required for a probe to be considered
#' @param coverage_min Numeric value corresponding to the minimum coverage 
#' requirement for a probe to be considered. Coverage is defined here as the 
#' proportion of samples without an NA value at a given probe.
#' @param sort_n Local value indicating whether to report the top n probes
#' @param sort_delta_beta Logical value indicating whether the report the top
#' probes satisfying delta_min
#'
#' @return List containing delta_betas, auc, and coverage for the top selected
#' probes
inferTopProbesHyper = function(res, auc_min=0.99, n=50, coverage_min=0.80, delta_min=0.3, sort_n=FALSE, sort_delta_beta=FALSE) {
    indx = which(res$auc >= auc_min & res$coverage >= coverage_min & res$delta_betas >= delta_min)
    if (sort_n) {
        indx_names = names(tail(sort(res$delta_betas[indx]), n=n))
        probes_n = list(delta_betas=res$delta_betas[indx_names], auc=res$auc[indx_names], coverage=res$coverage[indx_names])
        if (length(probes_n$delta_betas) == 0){
            message(sprintf("Found 0 probes (min delta: 0)"))
        } else {
            message(sprintf("Found %d probes (min delta: %f)", length(probes_n$delta_betas),  min(probes_n$delta_betas)))
        }
    } 
    if (sort_delta_beta) {
        probes_delta_beta = list(delta_betas=res$delta_betas[indx], auc=res$auc[indx], coverage=res$coverage[indx])
        if (length(probes_delta_beta$delta_betas) == 0){
            message(sprintf("Found 0 probes (min delta: 0)"))
        } else {
            message(sprintf("Found %d probes (min delta: %f)", length(probes_delta_beta$delta_betas),  min(probes_delta_beta$delta_betas)))
        }

    }
    if (sort_n & !sort_delta_beta) {
        return(probes_n)
    }
    if (sort_delta_beta & !sort_n) {
        return(probes_delta_beta)
    }
    if (length(probes_n) > length(probes_delta_beta)) {
        return(probes_n)
    }
    return(probes_delta_beta)
}

#' buildReferenceHyper builds and saves a reference for the top hypermethylated
#' probes from using the data in betas and the meta data provided. The 
#' intermediate files are stored in the specified file directory (id_dir).
#'
#' @param betas Matrix of beta value measurement for each sample (columns) 
#' across all probes of the array (rows)
#' @param meta DataFrame corresponding the the meta of the analysis where each 
#' branch column name contains a vector specifying the in-group (0), out-group 
#' (1), and ignored (2) samples as they appear in the column names of betas.
#' @param id_dir String corresponding to the file directory to which the probe 
#' information is stored.
#' @param verbose Logical value indicating whether to display the scanned  
#' results
#' @param return_results Logical value indicating whether to return the 
#' results from the function
#'
#' @return List containing delta_betas, auc, and coverage for a given 
#' branch and set of samples (if return_results).
#'
#' @export
buildReferenceHyper = function(betas, meta, id_dir, verbose=FALSE, return_results=FALSE) {
    res = list()
    branch2probes_hyper = list()

    branches = names(meta)[5:ncol(meta)]

    for (branch in branches) {
        res[[branch]] = inferProbesHyper(betas, meta[[branch]], scan_delta=T, delta_min=0.3, verbose=verbose)
        saveRDS(res[[branch]], sprintf("%s/fit_%s.rds", id_dir, branch))

        res[[branch]] = readRDS(sprintf("%s/fit_%s.rds", id_dir, branch))
        branch2probes_hyper[[branch]] = inferTopProbesHyper(res[[branch]],sort_delta_beta=T, delta_min=0.3)
        saveProbes(branch, branch2probes_hyper[[branch]], id_dir=id_dir)
    }

    return(res)
}

#' readBranchProbesDeltaBeta reads all tsv files in the given file directory 
#' (id_dir) and returns a list of branches with named vectors of delta_beta 
#' for each probe 
#'
#' @param id_dir String corresponding to the file directory to which the probe 
#' information is stored.
#' results
#'
#' @return List of branches with named vectors of delta_beta for each probe 
readBranchProbesDeltaBeta = function(id_dir) {
    with(do.call(rbind, lapply(list.files(id_dir, ".tsv"), function(x) read_tsv(sprintf("%s/%s", id_dir, x), col_names=TRUE, , col_types = cols()))), split(setNames(delta_beta, probeID), branchID))
}

#' readBranchProbesAUC reads all tsv files in the given file directory 
#' (id_dir) and returns a list of branches with named vectors of auc 
#' for each probe 
#'
#' @param id_dir String corresponding to the file directory to which the probe 
#' information is stored.
#' results
#'
#' @return List of branches with named vectors of auc for each probe 
readBranchProbesAUC = function(id_dir) {
    with(do.call(rbind, lapply(list.files(id_dir, ".tsv"), function(x) read_tsv(sprintf("%s/%s", id_dir, x), col_names=TRUE, , col_types = cols()))), split(setNames(auc, probeID), branchID))
}

#' readBranchProbesCoverage reads all tsv files in the given file directory 
#' (id_dir) and returns a list of branches with named vectors of coverage 
#' for each probe 
#'
#' @param id_dir String corresponding to the file directory to which the probe 
#' information is stored.
#' results
#'
#' @return List of branches with named vectors of coverage for each probe 
readBranchProbesCoverage = function(id_dir) {
    with(do.call(rbind, lapply(list.files(id_dir, ".tsv"), function(x) read_tsv(sprintf("%s/%s", id_dir, x), col_names=TRUE, , col_types = cols()))), split(setNames(coverage, probeID), branchID))
}

#' buildTissueSEReference builds a Summarized Experiment formatted reference.
#'
#' @param betas Matrix of beta value measurement for each sample (columns) 
#' across all probes of the array (rows)
#' @param meta DataFrame corresponding the the meta of the analysis where each 
#' branch column name contains a vector specifying the in-group (0), out-group 
#' (1), and ignored (2) samples as they appear in the column names of betas.
#' @param id_dir String corresponding to the file directory to which the probe 
#' information is stored.
#' @param tissue_color Named vector with the tissue corresponding to color 
#' HEX code
#' @param branchID_color  Named vector with the branch corresponding to color
#' HEX code
#'
#' @return Summarized Experiment with probe selection (row data), sample 
#' selection (column data), meta data, and the betas (assay)
#'
#' @export
buildTissueSEReference = function(betas, meta, id_dir, tissue_color, branchID_color) {

    branch2probes_delta_beta = readBranchProbesDeltaBeta(id_dir=id_dir)
    branch2probes_delta_beta_ = unlist(setNames(branch2probes_delta_beta, NULL))
    branch2probes_auc = readBranchProbesAUC(id_dir=id_dir)
    branch2probes_auc_ = unlist(setNames(branch2probes_auc, NULL))
    branch2probe_coverage = readBranchProbesCoverage(id_dir=id_dir)
    branch2probe_coverage_ = unlist(setNames(branch2probe_coverage, NULL))

    bt = clusterWithRowGroupingW(betas, group2probes = lapply(branch2probes_delta_beta, names))
    
    bt = clusterWithColumnGroupingW(bt, grouping = meta$tissue, ordered_groups = names(tissue_color)) #[c(47, 20, 24, 31, 14, 30, 12, 11, 10, 29, 1:8, 22, 21, 27, 32, 33, 17, 18, 37, 38, 16, 25, 26, 39, 40:46, 36, 68, 69)]

    rd = tibble(probeID=rownames(bt), 
                branchID=rep(names(branch2probes_delta_beta), sapply(branch2probes_delta_beta, length)),
                delta_beta=as.numeric(branch2probes_delta_beta_[match(rownames(bt), names(branch2probes_delta_beta_))]),
                auc=as.numeric(branch2probes_auc_[match(rownames(bt), names(branch2probes_auc_))]),
                coverage=as.numeric(branch2probe_coverage_[match(rownames(bt), names(branch2probe_coverage_))]))
    cd = meta[match(colnames(bt), meta$sample),]
    se = SummarizedExperiment(assays=list(betas=bt), rowData=rd, colData=cd)
    metadata(se)$tissue_color = tissue_color
    metadata(se)$branchID_color = branchID_color
    se
}

#' mergeTissueSEReferences merges a hypomethylated Summarized Experiment 
#' reference with a hypermethylated Summarized Experiment reference and 
#' validates that some of their data matches.
#'
#' @param hypoSE Summarized Experiment with hypomethylated probe selection 
#' (row data), sample selection (column data), meta data, and the betas (assay)
#' @param hyperSE Summarized Experiment with hypermethylated probe selection 
#' (row data), sample selection (column data), meta data, and the betas (assay)
#'
#' @return Summarized Experiment with hypomethylated and hypermethylated probe
#' selection (row data), sample selection (column data), meta data, and the 
#' betas (assay)
#'
#' @export
mergeTissueSEReferences = function(hypoSE, hyperSE) {
    cd_hypo = colData(hypoSE)
    cd_hyper = colData(hyperSE)
    if (!all(cd_hypo == cd_hypo)) {
        cat("The column data between the two Summarized Experiments are not the same.")
        return(NULL)
    }

    meta_hypo = metadata(hypoSE)
    meta_hyper = metadata(hyperSE)
    if (length(setdiff(meta_hypo, meta_hyper)) != 0) {
        cat("The meta data between the two Summarized Experiments are not the same.")
        return(NULL)
    }

    # Infer hypo/hyper from sign of delta beta
    rd = rbind(rowData(hypoSE), rowData(hyperSE))
    data = rbind(assay(hypoSE), assay(hyperSE))

    se = SummarizedExperiment(assays=list(betas=data), rowData=rd, colData=cd_hypo)
    metadata(se) = metadata(hypoSE)

    return(se)
}

#' inferTissue infers the tissue of a single sample (as identified through the 
#' branchIDs in the row data of the reference) by reporting independent 
#' composition through cell type deconvolution.
#'
#' @param reference Summarized Experiment with either hypomethylated or 
#' hypermethylated probe selection (row data), sample selection (column data), 
#' meta data, and the betas (assay)
#' @param sample Named vector with probes and their corresponding beta value 
#' measurement
#' @param report_first Logical value indicating whether to only report the 
#' tissue with the highest composition
#' @param ignore_branches Vector of branches for which their relative 
#' composition will not be reported.
#' @param abs_delta_beta_min Numerical value indicating the absolute minimum 
#' required delta beta for the probe selection criteria
#'
#' @return Named vector of tissues and their relative composition. The sum of 
#' the result is not guaranteed to be equal to one.
#'
#' @export
inferTissue = function(reference, sample, report_first=TRUE, ignore_branches=c("ML-Hematopoiesis", "MLH-Lymphoid", "MLH-Myeloid", "MLHL-T"), abs_delta_beta_min=0.3) {
    rd = rowData(reference)
    branches = unique(rd$branchID)
    inferred.tissue = unlist(setNames(lapply(branches, 
        function(branch) {
            probes_hypo = rd$probeID[rd$branchID == branch & rd$delta_beta < 0 & abs(rd$delta_beta) >= abs_delta_beta_min]
            probes_hyper = rd$probeID[rd$branchID == branch & rd$delta_beta > 0 & abs(rd$delta_beta) >= abs_delta_beta_min]
            betas_hypo = c(0)
            betas_hyper = c(0)
            if(length(probes_hypo) > 0)
                betas_hypo = sample[probes_hypo]
            if(length(probes_hyper) > 0)
                betas_hyper = sample[probes_hyper]
            return((sum(betas_hyper, na.rm=T) + sum(1 - betas_hypo, na.rm=T)) / (length(betas_hypo) + length(betas_hyper)))
            }), branches))
    inferred.tissue = inferred.tissue[!(names(inferred.tissue) %in% ignore_branches)]
    if(report_first) return(inferred.tissue[1])
    return(inferred.tissue)
}   

#' inferTissue infers the tissue of a matrix of samples (as identified through
#' the branchIDs in the row data of the reference) by reporting independent 
#' composition through cell type deconvolution.
#'
#' @param reference Summarized Experiment with either hypomethylated or 
#' hypermethylated probe selection (row data), sample selection (column data), 
#' meta data, and the betas (assay)
#' @param betas Matrix of beta value measurement for each sample (columns) 
#' across all probes of the array (rows)
#' @param report_first Logical value indicating whether to only report the 
#' tissue with the highest composition
#' @param ignore_branches Vector of branches for which their relative 
#' composition will not be reported.
#' @param abs_delta_beta_min Numerical value indicating the absolute minimum 
#' required delta beta for the probe selection criteria
#'
#' @return Matrix of tissues inference (rows) across many samples (columns)
#' and their relative composition. The sum of the result for any one given 
#' sample is not guaranteed to be equal to one.
#'
#' @export
inferTissues = function(refrence, betas, report_first=TRUE, ignore_branches=c("ML-Hematopoiesis", "MLH-Lymphoid", "MLH-Myeloid", "MLHL-T"), abs_delta_beta_min=0.3) {
    res = apply(betas, 2, 
        function(sample) { inferTissue(reference=reference, sample=sample, report_first=F)
            })
    return(res)
}


#' inferTissue infers the tissue of a matrix of samples (as identified through
#' the branchIDs in the row data of the reference) by reporting independent 
#' composition through cell type deconvolution.
#'
#' @param results Matrix of tissues inference (rows) across many samples 
#' (columns) and their relative composition. The sum of the result for any one 
#' given sample is not guaranteed to be equal to one.
#' @param betas Matrix of beta value measurement for each sample (columns) 
#' across all probes of the array (rows)
#' @param reference Summarized Experiment with either hypomethylated or 
#' hypermethylated probe selection (row data), sample selection (column data), 
#' meta data, and the betas (assay)
#'
#' @return Summarized experiment with meta data of the inferred samples 
#' (column data), meta data, and results of tissue inference (assay).
#'
#' @export
buildResultsSE = function(results, meta, reference) {
    cd = meta[match(colnames(results), meta$sample),]
    se = SummarizedExperiment(assays=list(results=results), colData=cd)
    metadata(se)$tissue_color = metadata(reference)$tissue_color
    metadata(se)$branchID_color = metadata(reference)$branchID_color
    se
}


#' inferTissue infers the tissue of a matrix of samples (as identified through
#' the branchIDs in the row data of the reference) by reporting independent 
#' composition through cell type deconvolution.
#'
#' @param resultsSE Summarized experiment with meta data of the inferred 
#' samples (column data), meta data, and results of tissue inference (assay).
#'
#' @return ggplot2 object
#'
#' @export
plotInferTissueResults = function(resultsSE) {
    cd = as_tibble(colData(resultsSE))
    md = metadata(resultsSE)
    g = WHeatmap(assay(resultsSE), cmp=CMPar(dmax=1, dmin=0, colorspace.name = 'diverge_hcl'), name="main", 
               xticklabels=T, xticklabels.n=ncol(assay(resultsSE)),  
               xticklabel.fontsize=8, #xticklabel.rotat = 90, 
               yticklabels=T, yticklabels.n=nrow(assay(resultsSE)),  
               yticklabel.fontsize=10, yticklabel.pad = 0.04)
    g = g + WLegendV("main", BottomRightOf(h.pad=0.05))
    # tissue labels (horizontal)
    g = g + WColorBarH(cd$tissue, TopOf("main",height=0.03),
                     cmp=CMPar(label2color=md$tissue_color), name="tissue_label",
                     xticklabels=T, xticklabel.side='t', xticklabel.fontsize=10,
                     #label.space = 1,
                     #xticklabel.space=0.0001,
                     label.use.data=TRUE, 
                     label.pad=0.2
     )
    # branchID labels (vertical)
    g = g + WColorBarV(rownames(res), LeftOf("main"),
                     cmp=CMPar(label2color=md$branchID_color), name="branch_label",
                     yticklabels=T, yticklabel.side='t', yticklabel.fontsize=10,
                     #label.space = 1,
                     #xticklabel.space=0.0001,
                     #label.use.data=TRUE, 
                     label.pad=0.2
    )
    g = g + WCustomize(mar.top=0.23, mar.left = 0.12, mar.bottom = 0.12)
    g
}


