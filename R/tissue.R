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

#' inferTissue1 infers the tissue of a single sample (as identified through 
#' the branchIDs in the row data of the reference) by reporting independent 
#' composition through cell type deconvolution.
#'
#' @param sample Named vector with probes and their corresponding beta value 
#' measurement
#' @param reference Summarized Experiment with either hypomethylated or 
#' hypermethylated probe selection (row data), sample selection (column data), 
#' meta data, and the betas (assay)
#' @param platform String representing the array type of the betas and 
#' reference
#' @param ignore_branches Vector of branches for which their relative 
#' composition will not be reported.
#' @param abs_delta_beta_min Numerical value indicating the absolute minimum 
#' required delta beta for the probe selection criteria
#' @param auc_min Numeric value corresponding to the minimum AUC value 
#' required for a probe to be considered
#' @param coverage_min Numeric value corresponding to the minimum coverage 
#' requirement for a probe to be considered. Coverage is defined here as the 
#' proportion of samples without an NA value at a given probe.
#' @param n Numerical value indicating the number of probes to at most use for 
#' any one given branch
#'
#' @return Summarized experiment with meta data of the inferred samples 
#' (column data), meta data, and results of tissue inference (assay).
#'
#' @export
inferTissue = function(sample, reference=NA, platform=NA, ignore_branches=c("ML-Hematopoiesis", "MLH-Lymphoid", "MLH-Myeloid", "MLHL-T"), abs_delta_beta_min=0.3, auc_min=0.99, coverage_min=0.80, n=NA) {

    if (!is.numeric(sample)) {
        cat("Provided sample is not a named numeric vector.\n")
        return(NULL)
    }

    if (is.na(reference)) {
        if (is.na(platform)) {
            platform = inferPlatform(names(sample))
        }
        reference = readRDS(url(paste("http://zhouserver.research.chop.edu/InfiniumAnnotation/current/", platform, "/", platform, ".reference.signature.rds", sep="")))
        #sesameDataGet(paste(platform, ".reference.signature.rds" , sep=""))
    }
    
    rd = rowData(reference)
    branches = unique(rd$branchID)
    results = unlist(setNames(lapply(branches, 
        function(branch) {
            rd_branch = rd[rd$branchID == branch & abs(rd$delta_beta) >= abs_delta_beta_min, ]

            if (is.na(n)) {
                n = nrow(rd_branch)
            }

            probes = rownames(head(rd_branch[order(-abs(rd_branch$delta_beta)), ], n=n))

            rd_ = rd[probes, ]

            probes_hypo = rd_$probeID[rd$delta_beta < 0]
            probes_hyper = rd$probeID[rd$delta_beta > 0]

            betas_hypo = c(0)
            betas_hyper = c(0)
            if(length(probes_hypo) > 0)
                betas_hypo = sample[probes_hypo]
            if(length(probes_hyper) > 0)
                betas_hyper = sample[probes_hyper]
            return((sum(betas_hyper, na.rm=T) + sum(1 - betas_hypo, na.rm=T)) / (length(betas_hypo) + length(betas_hyper)))
            }), branches))
    results = results[!(names(results) %in% ignore_branches)]

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



