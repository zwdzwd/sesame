#' listDatabaseSets prints which database sets are available for a given release
#'
#' @return One list of vectors corresponding to aggregated database sets.
#'
#' @examples
#' listDatabaseSets()
#'
#' @export
listDatabaseSets = function() {
    meta = sesameData:::df_master
    meta = meta[grepl('KYCG', meta$Title), ]
    
    x = apply(meta, 1, function(row) {
        cat(sprintf("Accession: %s (n: %s)\n", 
                    format(row["Title"], width = 50, justify = "l"), 
                    row["N"]))})
}

#' getDatabaseSets retrieves database sets from a meta data sheet by querying 
#' the group, platform, reference columns. The data is returned as a list where the
#' names correspond to chosen database sets.
#'
#' @param titles vector containing the characters associated with the
#' selected database sets; only non-NA locations will be returned. Optional.
#' (Default: c("20210810_MM285_TFBS_ENCODE").
#' @param group string representing the group for which the database sets will
#' be returned. Optional. (Default: NA).
#' @param platform string representing the platform (EPIC, HM450, HM27, MM285)
#' for which database sets will be returned. Optional. (Default: NA).
#' @param reference string representing the reference (hg19, hg38, mm9, mm10) 
#' for which the database sets will be returned. Optional. (Default NA).
#' @param verbose Logical value indicating whether intermediate outputs will be
#' displayed to console. Optional. (Default: TRUE).
#'
#' @return One list of vectors corresponding to aggregated database sets.
#'
#' @examples
#' getDatabaseSets()
#'
#' @export
getDatabaseSets = function(titles=NA, group=NA, 
                           platform=NA, reference=NA, 
                           verbose=TRUE) {
    meta = sesameData:::df_master
    # meta = meta[meta$kyCG, ]
    meta = meta[grepl('KYCG', meta$Title), ]
    
    if (any(!is.na(titles))) {
        meta = meta[which(meta$Title %in% titles), ]
    }
    
    if (!is.na(group)) {
        meta = meta[grepl(group, meta$Title, ignore.case=TRUE), ]
    }
    
    platform = sprintf('Platform%s', platform)
    if (platform %in% colnames(meta))
        meta = meta[as.logical(meta[[platform]]), ]
    
    if (!is.na(reference)) {
        meta = meta[which(meta$Reference %in% reference), ]
    }
    
    if (verbose) {
        print(sprintf("Retrieving %d databaseSets...", sum(meta$N)))
    }
    
    databaseSets = flattenlist(lapply(unlist(na.omit(meta$Title)), 
                                      function(title) {
                                          sesameDataGet(title, verbose=verbose)
                                      })
    )
    
    return(databaseSets)
}


#' flattenlist flattens a multidimensional list into a single dimensional list.
#'
#' @param x Multidimensional list.
#'
#' @return A single dimensional list.
#'
#' @import methods
#'
#' @examples
#' flattenlist(list(a=list(1,2,3), b=list(4,5,6)))
flattenlist = function(x) {
    morelists = vapply(x, function(x_) is(x_, 'list'), TRUE)
    out = c(x[!morelists], unlist(x[morelists], recursive=FALSE))
    if(sum(morelists)){
        Recall(out)
    } else{
        return(out)
    }
}


#' compareDatbaseSetOverlap tests for the pariwise overlap between given
#' list of database sets using a distance metric.
#'
#' @param databaseSets List of vectors corresponding to the database sets of
#' interest with associated meta data as an attribute to each element. Optional.
#' (Default: NA)
#' @param metric String representing the similarity metric to use. Optional.
#' (Default: "Jaccard").
#' @param verbose Logical value indicating whether to display intermediate
#' text output about the type of test. Optional. (Default: FALSE)
#'
#' @return An upper triangular matrix containing a metric (Jaccard) comparing
#' the pairwise distances between database sets.
#'
#' @examples
#' databaseSets = list(a=c("a", "b"), b=c("a", "e", "f"), c=c("q", "a"))
#' compareDatbaseSetOverlap(databaseSets)
#'
#' @export
compareDatbaseSetOverlap = function(databaseSets=NA,
                                    metric="Jaccard",
                                    verbose=FALSE) {
    ndatabaseSets = length(databaseSets)
    names = names(databaseSets)
    m = matrix(0, nrow=ndatabaseSets, ncol=ndatabaseSets)
    colnames(m) = names
    rownames(m) = names
    for (i in seq(ndatabaseSets - 1)) { 
        for (j in seq(i + 1, ndatabaseSets)) {
            cat(i, j, '\n')
            m[i, j] = length(intersect(databaseSets[[i]], databaseSets[[j]])) /
                length(union(databaseSets[[i]], databaseSets[[j]]))
        }
    }
    return(m)
}


#' getDatabaseSetOverlap tests for the overlap of set of probes (querySet) in a
#' single given feature (database set)
#'
#' @param querySet Vector of probes corresponding to a single database set
#' of interest.
#' @param databaseSets List of vectors corresponding to the database sets of
#' interest with associated meta data as an attribute to each element.
#' @param platform String corresponding to the type of platform to use. Either
#' MM285, EPIC, HM450, or HM27. If it is not provided, it will be inferred
#' from the query set probeIDs (Default: NA).
#' @param verbose Logical value indicating whether to display intermediate
#' text output about the type of test. Optional. (Default: FALSE)
#'
#' @return A sparse data.frame containing all of the meta data from all database
#' sets.
#'
#' @examples
#' querySet=c("cg29176188_TC21", "cg29176794_TC21")
#' databaseSet=c("cg29176188_TC21", "cg29176794_TC21")
#' getDatabaseSetOverlap(querySet, databaseSet)
#'
#' @export
getDatabaseSetOverlap = function(querySet,
                                 databaseSets,
                                 platform=NA,
                                 verbose=TRUE) {
    if (all(is.na(databaseSets))) {
        if (verbose) {
            cat("The databaseSets were not defined.", 
                "Loading in databaseSets based on platform.")
        }
        if (is.na(platform)) {
            if (verbose) {
                cat("The platform was not defined.",
                    "Inferring platform from probeIDs.")
            }
            platform = inferPlatformFromProbeIDs(querySet)
        }
        databaseSets = getDatabaseSets(platform=platform, verbose=verbose)
    }
    
    metadata = as.data.frame(
        do.call(rbind,
                lapply(databaseSets,
                       function(databaseSet) {
                           rowmeta = attr(databaseSet, "meta")
                           if (!is.null(rowmeta)) 
                               rowmeta = c(meta=TRUE, rowmeta)
                           else
                               rowmeta = c(meta=FALSE, rowmeta)
                           nQ = length(querySet)
                           nD = length(databaseSet)
                           overlap = length(intersect(querySet, databaseSet))
                           rowmeta = c(rowmeta, nQ=nQ, nD=nD, overlap=overlap)
                           return(rowmeta)
                       })
        ))
    
    metadata$meta = as.logical(metadata$meta)
    
    return(metadata)
}


#' testEnrichment1 tests for the enrichment of set of probes (query set) in a
#' single given feature (database set)
#'
#' @param querySet Vector of probes of interest (e.g., significant probes)
#' @param databaseSet Vector corresponding to the database sets of
#' interest with associated meta data as an attribute to each element.
#' @param universeSet Vector of probes in the universe set containing all of
#' the probes to be considered in the test.
#' @param estimate.type String indicating the estimate to report. (Default:
#' "ES")
#' @param p.value.adj Logical value indicating whether to report the adjusted
#' p-value. (Default: FALSE)
#' @param verbose Logical value indicating whether to display intermediate
#' text output about the type of test. Optional. (Default: FALSE)
#'
#' @import utils
#'
#' @return One list containing features corresponding the test estimate,
#' p-value, and type of test.
testEnrichment1 = function(querySet, databaseSet, universeSet,
                           estimate.type="ES", p.value.adj=FALSE,
                           verbose=FALSE) {
    if (is.numeric(querySet)) { # a named vector of continuous value
        if(is.numeric(databaseSet)) { # numeric db
            if (verbose) {
                cat("Query set: Continuous\t",
                    "Database set: Continuous\t[Spearman test]\n")
            }
            results = testEnrichmentSpearman(
                querySet=querySet,
                databaseSet=databaseSet)
        } else {
            if (verbose) {
                cat("Query set: Continuous\t",
                    "Database set: Discrete\t\t[FGSEA test]\n")
            }
            results = testEnrichmentFGSEA(
                querySet=querySet,
                databaseSet=databaseSet,
                p.value.adj=p.value.adj,
                estimate.type=estimate.type)
        }
    } else { # categorical query
        if(is.numeric(databaseSet)) { # numeric db
            ## do fgsea(switched arguments)
            if (verbose) {
                cat("Query set: Discrete\t", 
                    "Database set: Continuous\t[FGSEA test]\n")
            }
            results = testEnrichmentFGSEA(
                querySet=databaseSet,
                databaseSet=querySet,
                p.value.adj=p.value.adj,
                estimate.type=estimate.type)
        } else { # categorical db
            if (verbose) {
                cat("Query set: Discrete\t",
                    "Database set: Discrete\t\t[Fisher exact test]\n")
            }
            results = testEnrichmentFisher(
                querySet=querySet,
                databaseSet=databaseSet,
                universeSet=universeSet)
        }
    }
    return(results)
}


#' testEnrichmentAll tests for the enrichment of set of probes (query set) in
#' a number of features (database sets).
#'
#' @param querySet Vector of probes of interest (e.g., significant probes)
#' @param databaseSets List of vectors corresponding to the database sets of
#' interest with associated meta data as an attribute to each element. Optional.
#' (Default: NA)
#' @param universeSet Vector of probes in the universe set containing all of
#' the probes to be considered in the test. If it is not provided, it will be
#' inferred from the provided platform. (Default: NA).
#' @param platform String corresponding to the type of platform to use. Either
#' MM285, EPIC, HM450, or HM27. If it is not provided, it will be inferred
#' from the query set probeIDs (Default: NA).
#' @param estimate.type String indicating the estimate to report. (Default:
#' "ES")
#' @param p.value.adj Logical value indicating whether to report the adjusted
#' p-value. (Default: FALSE).
#' @param n.fdr Integer corresponding to the number of comparisons made. 
#' Optional. (Default: NA).
#' @param verbose Logical value indicating whether to display intermediate
#' text output about the type of test. Optional. (Default: FALSE).
#'
#' @return One list containing features corresponding the test estimate,
#' p-value, and type of test.
#'
#' @examples
#' testEnrichmentAll(c("cg0000029"))
#'
#' @export
testEnrichmentAll = function(querySet, databaseSets=NA, universeSet=NA,
                             platform=NA, estimate.type="ES", p.value.adj=FALSE,
                             n.fdr=NA, verbose=FALSE) {
    if (all(is.na(universeSet))) {
        if (verbose) {
            cat("The universeSet was not defined.", 
                "Loading in universeSet based on platform.")
        }
        if (is.na(platform)) {
            if (verbose) {
                cat("The platform was not defined.",
                    "Inferring platform from probeIDs.")
            }
            platform = inferPlatformFromProbeIDs(querySet)
        }
        
        manifests = c("MM285.mm10.manifest",
                      "EPIC.hg19.manifest",
                      "HM450.hg19.manifest",
                      "HM27.hg19.manifest")
        
        manifest = manifests[grepl(platform, manifests)]
        
        universeSet = tryCatch({
            names(sesameDataGet(manifest))
        },
        error = function (condition) {
            print("ERROR:")
            print(paste("  Message:",conditionMessage(condition)))
            print(paste("  Call: ",conditionCall(condition)))
            return(NULL)
        },
        finally= function() {
            print(sprintf("Invalid platform [%s]", platform))
            return(NULL)
        })
    }
    
    if (all(is.na(databaseSets))) {
        if (verbose) {
            cat("The databaseSets were not defined.", 
                "Loading in databaseSets based on platform.")
        }
        if (is.na(platform)) {
            if (verbose) {
                cat("The platform was not defined.",
                    "Inferring platform from probeIDs.")
            }
            platform = inferPlatformFromProbeIDs(querySet)
        }
        databaseSets = getDatabaseSets(platform=platform, verbose=verbose)
    }
    
    results = data.frame(
        do.call(rbind,
                lapply(databaseSets,
                       function(databaseSet) testEnrichment1(
                           querySet=querySet,
                           databaseSet=databaseSet,
                           universeSet=universeSet,
                           p.value.adj=p.value.adj,
                           estimate.type=estimate.type,
                           verbose=verbose)
                )
        ))
    
    if (is.na(n.fdr))
        results$p.adjust.fdr = p.adjust(results$p.value, method='fdr')
    else
        results$p.adjust.fdr = p.adjust(results$p.value, method='fdr', n=n.fdr)
    
    metadata = data.frame(
        do.call(rbind,
                lapply(databaseSets,
                       function(databaseSet) {
                           output = attr(databaseSet, "meta")
                           if (!is.null(output)) 
                               return(append(c(meta=TRUE), output))
                           return(c(meta=FALSE))
                       })
        )
    )
    
    rank = list()
    rank$estimate.rank = rank(-results$estimate, ties.method='first')
    rank$p.value.rank = rank(results$p.value, ties.method='first')
    rank$overlap.rank = rank(results$overlap, ties.method='first')
    rank$max.rank = apply(data.frame(rank), 1, max)
    rank$mean.rank = apply(data.frame(rank), 1, mean)
    rank = data.frame(rank, row.names=row.names(results))
    output = cbind(results, rank, metadata)
    
    output = output[order(output$p.value, decreasing=FALSE), ]
    return(output)
}


#' testEnrichmentGene tests for the enrichment of set of probes
#' (querySet) in gene regions.
#'
#' @param querySet Vector of probes of interest (e.g., probes belonging to a
#' given platform)
#' @param platform String corresponding to the type of platform to use. Either
#' MM285, EPIC, HM450, or HM27. If it is not provided, it will be inferred
#' from the query set querySet (Default: NA)
#' @param verbose Logical value indicating whether to display intermediate
#' text output about the type of test. Optional. (Default: FALSE)
#'
#' @return One list containing features corresponding the test estimate,
#' p-value, and type of test.
#'
#' @examples
#' testEnrichmentGene(c("cg0000029"), platform="EPIC")
#'
#' @export
testEnrichmentGene = function(querySet, platform=NA, verbose=FALSE) {
    if (is.na(platform)) {
        if (verbose) {
            print("The platform was not defined. Inferring platform from probeIDs.")
        }
        platform = inferPlatformFromProbeIDs(querySet)
    }
    
    probeID2gene = tryCatch({
        sesameDataGet(
            sprintf('%s.probeID2gene.20210913', platform))
    },
    error = function (condition) {
        print("ERROR:")
        print(paste("  Message:",conditionMessage(condition)))
        print(paste("  Call: ",conditionCall(condition)))
        return(NULL)
    },
    finally= function() {
        print(sprintf("Invalid platform [%s]", platform))
        return(NULL)
    })
    
    databaseSetNames = probeID2gene$genesUniq[match(querySet, 
                                                    probeID2gene$probeID)]
    
    databaseSetNames = na.omit(unique(
        unlist(lapply(databaseSetNames,
                      function(databaseSetName) {
                          strsplit(databaseSetName, ";")
                      }))))
    
    if (length(databaseSetNames) == 0) return(NULL)
    
    databaseSets = getDatabaseSets(group="gene", platform=platform)
    
    n = length(databaseSets)
    
    databaseSets = databaseSets[names(databaseSets) %in% databaseSetNames]
    
    return(testEnrichmentAll(querySet, databaseSets, platform=platform, n.fdr=n))
}


#' testEnrichmentFisher uses Fisher's exact test to estimate the association
#' between two categorical variables.
#'
#' @param querySet Vector of probes of interest (e.g., significant probes)
#' @param databaseSet Vectors corresponding to the database set of
#' interest with associated meta data as an attribute to each element.
#' @param universeSet Vector of probes in the universe set containing all of
#' the probes to be considered in the test. (Default: NULL)
#'
#' @import stats
#'
#' @return A DataFrame with the estimate/statistic, p-value, and name of test
#' for the given results.
testEnrichmentFisher = function(querySet, databaseSet, universeSet) {
    test = "fisher"
    if (length(intersect(querySet, databaseSet)) == 0) {
        return(data.frame(estimate=0,
                          p.value=1,
                          test=test,
                          nQ=length(querySet),
                          nD = length(databaseSet),
                          overlap=0
        ))
    }
    
    mtx = matrix(c(
        length(intersect(querySet, databaseSet)),
        length(setdiff(databaseSet, querySet)),
        length(setdiff(querySet, databaseSet)),
        length(setdiff(universeSet, union(databaseSet, querySet)))),
        nrow = 2,
        dimnames = list(
            querySet = c("Q_in","Q_out"),
            databaseSet = c("D_in","D_out")))
    
    res = fisher.test(mtx)
    
    result = data.frame(
        estimate = calcFoldChange(mtx),
        p.value = res$p.value,
        test = test,
        nQ = length(querySet),
        nD = length(databaseSet),
        overlap = length(intersect(querySet, databaseSet))
    )
    return(result)
}

#' calcFoldChange calculates fold change given a 2x2 matrix of counts.
#'
#' @param mtx 2x2 matrix of values corresponding to overlapping counts between
#' two sets of a categorical variable.
#'
#' @return A numerical value corresponding to the fold change enrichment,
calcFoldChange = function(mtx){
    num = mtx[1, 1] / (mtx[1, 1] + mtx[1, 2])
    den = (mtx[1, 1] + mtx[2, 1]) / sum(mtx)
    num / den
} 


#' testEnrichmentFGSEA uses the FGSEA test to estimate the association of a
#' categorical variable against a continuous variable.
#'
#' @param querySet Vector of probes of interest (e.g., significant probes)
#' @param databaseSet Vector of probes corresponding to a single database set
#' of interest.
#' @param p.value.adj Logical value indicating whether to report the adjusted
#' p-value. (Default: FALSE).
#' @param estimate.type String indicating the estimate to report. Optional.
#' (Default: "ES").
#' 
#' 
#' @import fgsea
#'
#' @return A DataFrame with the estimate/statistic, p-value, and name of test
#' for the given results.
testEnrichmentFGSEA = function(querySet, databaseSet, p.value.adj=FALSE,
                               estimate.type="ES") {
    test="fgsea"
    overlap = length((intersect(names(querySet), databaseSet)))
    if (overlap == 0) {
        return(data.frame(estimate=0,
                          p.value=1,
                          test=test,
                          nQ=length(databaseSet),
                          nD=length(querySet),
                          overlap=overlap
        ))
    }
    res = fgsea(pathways=list(pathway=databaseSet), 
                stats=querySet)
    
    if (p.value.adj) {
        p.value = res$padj
    } else {
        p.value = res$pval
    }
    
    if (estimate.type == "log2err") {
        estimate = res$log2err
    } else if (estimate.type == "NES") {
        estimate = res$NES
    } else if (estimate.type == "leadingEdge") {
        estimate = res$leadingEdge
    } else if (estimate.type == "ES") {
        estimate = res$ES
    } else {
        print(sprintf("Incorrect estimate.type: [%s].", estimate.type))
        return(NULL)
    }
    
    result = data.frame(
        estimate = estimate,
        p.value = p.value,
        test = test,
        nQ=length(databaseSet),
        nD=length(querySet),
        overlap=overlap
    )
    return(result)
}


#' testEnrichmentSpearman uses the Spearman statistical test to estimate the 
#' association between two continuous variables.
#'
#' @param querySet Vector of probes of interest (e.g., significant probes)
#' @param databaseSet List of vectors corresponding to the database set of
#' interest with associated meta data as an attribute to each element.
#'
#' @import stats
#'
#' @return A DataFrame with the estimate/statistic, p-value, and name of test
#' for the given results.
testEnrichmentSpearman = function(querySet, databaseSet) {
    test = "spearman"
    if (length(intersect(names(querySet), names(databaseSet))) == 0) {
        return(data.frame(estimate=0,
                          p.value=1,
                          test=test,
                          nQ=length(querySet),
                          nD=length(databaseSet),
                          overlap=0
        ))
    }
    
    databaseSet = databaseSet[match(names(querySet), names(databaseSet))]
    
    res = cor.test(
        querySet,
        databaseSet,
        method = test
    )
    result = data.frame(
        estimate = res$estimate[[1]],
        p.value = res$p.value,
        test = test
    )
    return(result)
}


#' calcDatabaseSetStatistics1 calculates features of x
#'
#' @param x Vector of numeric values
#'
#' @return Vector with ~20 different engineered features
#'
#' @import stats
calcDatabaseSetStatistics1 = function(x) {
    a = data.frame(mean=apply(x, 2, mean, na.rm=TRUE),
                   median=apply(x, 2, median, na.rm=TRUE),
                   var=apply(x, 2, var, na.rm=TRUE),
                   sd=apply(x, 2, sd, na.rm=TRUE),
                   skew=apply(x, 2, var, na.rm=TRUE),
                   iqr=apply(x, 2, IQR, na.rm=TRUE),
                   range=apply(x, 2, max, na.rm=TRUE) - apply(x, 2, min, na.rm=TRUE),
                   min=apply(x, 2, min, na.rm=TRUE),
                   max=apply(x, 2, max, na.rm=TRUE))
    b = apply(x, 2, quantile, na.rm=TRUE, probs=seq(0, 1, 0.1))
    return(cbind(a, t(b)))
}


#' calcDatabaseSetStatisticsAll builds dataset for a given betas matrix 
#' composed of engineered features from the given database sets
#'
#' @param betas matrix of beta values where probes are on the rows and samples
#' are on the columns
#' @param databaseSets List of vectors corresponding to probe locations for
#' which the features will be extracted
#' 
#' @examples 
#' betas = getBetas("MM285", dev=TRUE, verbose=TRUE)
#' databaseSetNames = c('20210630_MM285_mm10_CpGDensity',
#' '20210630_MM285_mm10_CGI', 20210816_MM285_mm10_distToTSS',
#' '20210210_MM285_design', '20210630_MM285_mm10_probe_type')
#' databaseSets = getDatabaseSets(databaseSetNames, dev=TRUE)
#' calcDatabaseSetStatisticsAll(betas, databaseSets)
#' 
#' @return Vector for a given sample columns are features across different
#' databaseSets
#' 
#' @export
calcDatabaseSetStatisticsAll = function(betas, databaseSets) {
    a = do.call(cbind, 
                lapply(names(databaseSets),
                       function(databaseSetName) {
                           databaseSet = databaseSets[[databaseSetName]]
                           if (length(databaseSet) >= nrow(betas)) return(FALSE)
                           if (is.numeric(databaseSet)) {
                               probes = names(databaseSet)
                           } else {
                               probes = databaseSet
                           }
                           
                           statistics = suppressWarnings(
                               calcDatabaseSetStatistics1(
                                   betas[na.omit(match(probes, rownames(betas))), ]))
                           names(statistics) = unlist(lapply(names(statistics), function(colname) {
                               paste(databaseSetName, colname, sep="-")
                           }))
                           return(statistics)
                       }))
    b = a[, !grepl("FALSE", colnames(a))]
    c = b[, !apply(b, 2, function(x) {any(is.na(x) | is.infinite(x))})]
    return(c)
}


#' skew determines the skew of a distribution x, taken from the Moments package
#'
#' @param x Vector of numeric values
#' @param na.rm Logical value corresponding to whether NA will be ignored
#'
#' @return Numeric value quantifying the skew of the distribution x
skew = function (x, na.rm = FALSE) {
    if (is.matrix(x))
        apply(x, 2, skew, na.rm = na.rm)
    else if (is.vector(x)) {
        if (na.rm)
            x <- x[!is.na(x)]
        n <- length(x)
        (sum((x - mean(x))^3)/n)/(sum((x - mean(x))^2)/n)^(3/2)
    }
    else if (is.data.frame(x))
        vapply(x, skew, na.rm = na.rm)
    else skew(as.vector(x), na.rm = na.rm)
}


#' plotVolcano creates a volcano plot of -log2(p.value) and log(estimate)
#' given data with fields estimate and p.value.
#'
#' @param data DataFrame where each field is a database name with two fields
#' for the estimate and p.value.
#' @param title String representing the title label. Optional. (Default: NA)
#' @param subtitle String representing the subtitle label. Optional. (Default:
#' NA)
#'
#' @return ggplot volcano plot
#'
#' @import ggplot2
#' @import ggrepel
#'
#' @examples
#' data=data.frame(estimate=c(runif(10)), p.value=c(runif(10)))
#' plotVolcano(data)
#'
#' @export
plotVolcano = function(data, title=NA, subtitle=NA, n.fdr=FALSE) {
    options(ggrepel.max.overlaps = 10)
    
    if ("Target" %in% colnames(data))
        data["label"] = unlist(data[["Target"]])
    else
        data["label"] = rownames(data)
    
    if (is.na(title)) {
        title = "Volcano plot"
    }
    title = gsub('(.{1,80})(\\s|$)', '\\1\n', title)
    
    if (is.na(subtitle)) {
        subtitle = ''
    }
    
    if (n.fdr) {
        data$p.value = data$p.adjust.fdr
    }
    
    subtitle = gsub('(.{1,80})(\\s|$)', '\\1\n', subtitle)
    
    # TODO: repalce with column specifying sig vs non sig
    
    if (any(data$p.value <= 0.05)) {
        g = ggplot(data=data, aes(x=log2(estimate), y=-log10(p.value),
                                  color = cut(p.value, c(-Inf, 0.05))))
    } else {
        g = ggplot(data=data, aes(x=estimate, y=p.value))
    }
    g = g + geom_point() + 
        xlab("log2 Fold Change")
    
    if(is.na(n.fdr)) {
        g = g + 
            ylab("-log10 p-value") +
            scale_colour_discrete(
                name = "Significance (p < 0.05)",
                labels=c("Significant", "Not Significant")
            )
    } else {
        g = g + 
            ylab("-log10 q-value") +
            scale_colour_discrete(
                name = "Significance (q < 0.05)",
                labels=c("Significant", "Not Significant")
            )
    }
    g = g + labs(
        title = title,
        subtitle = subtitle,
        fill = "pvalue"
    ) +
        theme(
            plot.title = element_text(size=16, face = "bold"),
            axis.text = element_text(size=12),
            axis.title = element_text(size=12),
            legend.title = element_text(size=12),
            legend.text = element_text(size=12)
        ) +
        geom_text_repel(
            data = subset(data, p.value < 0.05),
            aes(label = label),
            size = 5,
            box.padding = unit(0.35, "lines"),
            point.padding = unit(0.3, "lines")
        )
    g
}


#' plotLollipop creates a lollipop plot of log(estimate) given data with fields
#' estimate.
#'
#' @param data DataFrame where each field is a database name with its estimate.
#' @param n Integer representing the number of top enrichments to report.
#' Optional. (Default: 10)
#' @param title String representing the title label. Optional. (Default: NA)
#' @param subtitle String representing the subtitle label. Optional. (Default:
#' NA)
#'
#' @return ggplot lollipop plot
#'
#' @import ggplot2
#'
#' @examples
#' data=data.frame(estimate=c(runif(10, 0, 10)))
#' plotLollipop(data)
#'
#' @export
plotLollipop = function(data, n=10, title=NA, subtitle=NA) {
    # data = data[which(as.logical(data$meta)), ]
    
    if ("Target" %in% colnames(data))
        data["label"] = unlist(data[["Target"]])
    else
        data["label"] = rownames(data)
    
    data = head(data[order(data$estimate, decreasing=TRUE), ], n=n)
    
    if (is.na(title)) {
        title = 'Lollipop Plot'
    }
    
    if (is.na(subtitle)) {
        subtitle = ''
    }
    
    ggplot(data, aes(x=label, 
                     y=log2(estimate), 
                     label=sprintf('%.2f',log2(estimate)))) +
        geom_hline(yintercept=0) +
        geom_segment(aes(y=0, 
                         x=reorder(label, -estimate), 
                         yend=log2(estimate), xend=label), color='black') +
        geom_point(aes(fill=pmax(-1.5,log2(estimate))), 
                   stat='identity', 
                   size=10, 
                   alpha=0.95, 
                   shape=21) +
        scale_fill_gradientn(name='Fold Change',
                             colours=c('#2166ac','#333333','#b2182b'),
                             # limits=c(min(log2(data$estimate + 1)),
                             #         max(log2(data$estimate + 1)) )
        ) +
        geom_text(color='white', size=3) +
        labs(title=title, subtitle=subtitle) +
        geom_label(aes(x=label,
                       y=ifelse(estimate>1,
                                log2(estimate) + 0.8,
                                log2(estimate) - 0.5),
                       label=label),
                   alpha=0.8) +
        # new_scale("fill") +
        # scale_fill_manual(values=hmm_colors) +
        ylab("Log2 Enrichment") +
        theme(axis.title.x = element_blank(),
              axis.text.x = element_blank(),
              axis.ticks.x = element_blank())
}

#' createGeneNetwork creates databaseSet network using the given similarity
#' metric.
#'
#' @param databaseSets Vector of probes corresponding to a single database set
#' of interest.
#' @param metric String representing the similarity score to use. Optional.
#' (Default: "Jaccard").
#'
#' @return ggplot lollipop plot
#'
#' @import RCy3
#' @import reshape2
#'
#' @examples
#' databaseSets = list(a=c("a", "b"), b=c("a", "e", "f"), c=c("q", "a"))
#' createDatabaseSetNetwork(databaseSets)
#'
#' @export
createDatabaseSetNetwork = function(databaseSets, 
                                    title="Database Interaction Network", 
                                    collection="DatabaseSets") {
    m = getDatabaseSetPairwiseDistance(databaseSets, metric="jaccard")
    saveRDS(m, "/Users/ethanmoyer/Dropbox/Ongoing_knowYourCpG/data/databaseSetNetwork.rds")
    
    m_ = m
    m = m_[seq(50), seq(50)]
    
    m_melted = melt(m); colnames(m_melted) = c("gene1", "gene2", "metric")
    m_melted = m_melted[m_melted$metric != 0, ]
    
    # Used for additional attributes like color, size, name. This is for GSM
    nodes <- data.frame(id=colnames(m),
                        # group=c("A","A","B","B"), # categorical strings
                        # score=as.integer(c(20,10,15,5)), # integers
                        stringsAsFactors=FALSE)
    # This is for Target
    edges <- data.frame(source=m_melted$gene1,
                        target=m_melted$gene2,
                        # interaction=NULL, Maybe for positive/negative assocation
                        weight=m_melted$metric, # numeric
                        stringsAsFactors=FALSE)
    
    # Return nodes and edges
    return(list(nodes=nodes, edges=edges))
    
}

