
#' getDatabaseSetOverlap tests for the overlap of set of probes (query) in a
#' single given feature (database set)
#'
#' @param query Vector of probes corresponding to a single database set
#' of interest.
#' @param databases List of vectors corresponding to the database sets of
#' interest with associated meta data as an attribute to each element.
#' @param platform String corresponding to the type of platform to use. Either
#' MM285, EPIC, HM450, or HM27. If it is not provided, it will be inferred
#' from the query set probeIDs (Default: NA).
#' @param verbose Logical value indicating whether to display intermediate
#' text output about the type of test. Optional. (Default: FALSE)
#'
#' @return A sparse data.frame containing all of the meta data from all database
#' sets.
#' @export
#' @examples
#'
#' library(SummarizedExperiment)
#' MM285.tissueSignature <- sesameDataGet('MM285.tissueSignature')
#' df <- rowData(MM285.tissueSignature)
#' query <- df$Probe_ID[df$branch == "B_cell"]
#' databaseNames <- c('KYCG.MM285.seqContextN.20210630', 
#' 'KYCG.MM285.designGroup.20210210')
#' databases <- do.call(c, lapply(databaseNames, sesameDataGet))
#' getDatabaseSetOverlap(query, databases)
#' sesameDataClearCache()
#' 
getDatabaseSetOverlap <- function(
    query, databases, platform=NA, verbose=TRUE) {
    
    if (all(is.na(databases))) {
        if (verbose) {
            message("The databases were not defined. ", 
                "Loading in databases based on platform.")
        }
        if (is.na(platform)) {
            if (verbose) {
                message("The platform was not defined. ",
                    "Inferring platform from probeIDs.")
            }
            if (is.numeric(query)) {
                platform <- inferPlatformFromProbeIDs(names(query))
            } else {
                platform <- inferPlatformFromProbeIDs(query)
            }
        }
        databaseNames <- sesameData:::df_master$Title[
            grepl(paste("KYCG.", platform, sep=''), 
                sesameData:::df_master$Title)]
        databases <- do.call(c, lapply(databaseNames, sesameDataGet))
    }
    
    metadata <- as.data.frame(
        do.call(rbind,
            lapply(databases,
                function(database) {
                    rowmeta <- attr(database, "meta")
                    if (!is.null(rowmeta)) {
                        rowmeta <- c(meta=TRUE, rowmeta)
                    } else {
                        rowmeta <- c(meta=FALSE, rowmeta)
                    }
                    nQ <- length(query)
                    nD <- length(database)
                    overlap <- length(intersect(query, database))
                    rowmeta <- c(rowmeta, nQ=nQ, nD=nD, overlap=overlap)
                    return(rowmeta)
                })
        ))
    
    metadata$meta <- as.logical(metadata$meta)
    
    metadata
}


#' testEnrichment1 tests for the enrichment of set of probes (query set)
#' in a single given feature (database set)
#'
#' @param query Vector of probes of interest (e.g., significant probes)
#' @param database Vector corresponding to the database sets of
#' interest with associated meta data as an attribute to each element.
#' @param universe Vector of probes in the universe set containing all of
#' the probes to be considered in the test.
#' @param estimate.type String indicating the estimate to report. (Default:
#' "ES")
#' @import utils
#' @return One list containing features corresponding the test estimate,
#' p-value, and type of test.
testEnrichment1 <- function(
    query, database, universe, estimate.type="ES") {
    
    if (is.numeric(query)) { # a named vector of continuous value
        if(is.numeric(database)) { # numeric db
            results <- testEnrichmentSpearman(
                query=query,
                database=database)
        } else {
            results <- testEnrichmentFGSEA(
                query = query,
                database = database,
                estimate.type=estimate.type)
        }
    } else { # categorical query
        if(is.numeric(database)) { # numeric db
            results <- testEnrichmentFGSEA(
                query = database,
                database = query,
                estimate.type=estimate.type)
        } else { # categorical db
            results <- testEnrichmentFisher(
                query = query,
                database = database,
                universe = universe)
        }
    }
    results
}

inferPlatformFromQuery <- function(query) {
    if (is.numeric(query)) {
        inferPlatformFromProbeIDs(names(query))
    } else {
        inferPlatformFromProbeIDs(query)
    }
}

#' testEnrichment tests for the enrichment of set of probes (query set) in
#' a number of features (database sets).
#'
#' @param query Vector of probes of interest (e.g., significant probes)
#' @param databases List of vectors corresponding to the database sets of
#' interest with associated meta data as an attribute to each element. Optional.
#' (Default: NA)
#' @param universe Vector of probes in the universe set containing all of
#' the probes to be considered in the test. If it is not provided, it will be
#' inferred from the provided platform. (Default: NA).
#' @param platform String corresponding to the type of platform to use. Either
#' MM285, EPIC, HM450, or HM27. If it is not provided, it will be inferred
#' from the query set probeIDs (Default: NA).
#' @param estimate.type String indicating the estimate to report. (Default:
#' "ES")
#' @param return.meta Logical value indicating whether to return meta data 
#' columns for those database sets containing sparse meta data information.
#' @param verbose Logical value indicating whether to display intermediate
#' text output about the type of test. Optional. (Default: FALSE).
#'
#' @return One list containing features corresponding the test estimate,
#' p-value, and type of test.
#' @importFrom dplyr bind_rows
#' @examples
#' 
#' library(SummarizedExperiment)
#' df <- rowData(sesameDataGet('MM285.tissueSignature'))
#' query <- df$Probe_ID[df$branch == "B_cell"]
#' testEnrichment(query, 'KYCG.MM285.designGroup.20210210')
#' sesameDataClearCache()
#'
#' @export
testEnrichment <- function(
    query, databases = NA, universe=NA,
    platform=NA, estimate.type="ES", return.meta=FALSE, verbose=FALSE) {

    if (all(is.na(universe))) { # infer uset from platform if not given
        if (is.na(platform)) {     # infer platform from probe ID
            platform <- inferPlatformFromQuery(query)
        }
        
        manifests <- c("MM285.mm10.manifest", "EPIC.hg19.manifest",
            "HM450.hg19.manifest", "HM27.hg19.manifest")
        manifest <- manifests[grepl(platform, manifests)]
        stopifnot(length(manifest) == 1 && all(manifest %in% manifests))
        universe <- names(sesameDataGet(manifest))
    }

    if (all(is.na(databases))) { # db not give, load a default set
        if (is.na(platform)) { platform <- inferPlatformFromQuery(query); }
        query_names <- grep("(chromHMM)|(designGroup|probeType)",
            sesameData:::df_master$Title[
                grepl(paste("KYCG.", platform, sep=''), 
                    sesameData:::df_master$Title)], value=TRUE)
    } else if (is.character(databases)) { # db given in names
        query_names <- guess_dbnames(databases)
    } else {
        query_names <- NULL
    }

    if (!is.null(query_names)) { # needs sesameDataGet
        dblist <- lapply(query_names, sesameDataGet)
        gpnames <- do.call(c, lapply(seq_along(dblist), function(ii) {
            rep(names(dblist)[ii], length(dblist[[ii]]))}))
        dbnames <- do.call(c, lapply(dblist, names))
        databases <- do.call(c, dblist)
    } else { # given explicitly
        gpnames <- rep("", length(databases))
        dbnames <- names(databases)
    }
    
    res <- data.frame(do.call(rbind, lapply(
        databases, function(db) {
            testEnrichment1(
                query = query,
                database = db,
                universe = universe,
                estimate.type = estimate.type)
        })))

    res$db <- dbnames
    res$group <- gpnames
    res$fdr <- p.adjust(res$p.value, method='fdr')
    rownames(res) <- NULL

    meta <- do.call(bind_rows, lapply(databases, function(db) {
        m1 <- attr(db, "meta")
        if (is.null(m1)) {
            data.frame(hasMeta = FALSE)
        } else {
            c(m1, hasMeta = TRUE)
        }
    }))
    res <- cbind(res, meta)
    res[order(-res$estimate), ]
}

#' testEnrichmentGene tests for the enrichment of set of probes
#' (query) in gene regions.
#'
#' @param query Vector of probes of interest (e.g., probes belonging to a
#' given platform)
#' @param databases List of vectors corresponding to the database sets of
#' interest with associated meta data as an attribute to each element. Optional.
#' (Default: NA)
#' @param platform String corresponding to the type of platform to use. Either
#' MM285, EPIC, HM450, or HM27. If it is not provided, it will be inferred
#' from the query set query (Default: NA)
#' @param verbose Logical value indicating whether to display intermediate
#' text output about the type of test. Optional. (Default: FALSE)
#'
#' @return One list containing features corresponding the test estimate,
#' p-value, and type of test.
#'
#' @examples
#' 
#' library(SummarizedExperiment)
#' MM285.tissueSignature <- sesameDataGet('MM285.tissueSignature')
#' df <- rowData(MM285.tissueSignature)
#' query <- df$Probe_ID[df$branch == "B_cell"]
#' testEnrichmentGene(query, platform="MM285", verbose=FALSE)
#'
#' @export
testEnrichmentGene <- function(
    query, databases=NA, platform=NA, verbose=FALSE) {

    if (is.na(platform)) {
        if (verbose)
            message("The platform was not defined.',
                'Inferring platform from probeIDs.")
        if (is.numeric(query)) {
            platform <- inferPlatformFromProbeIDs(names(query))
        } else {
            platform <- inferPlatformFromProbeIDs(query)
        }
    }
    if (is.na(databases)) {
        databases <- sesameDataGet(
            sprintf('KYCG.%s.gene.20210923', platform))
    }
    
    probeID2gene <- attr(databases, 'probeID2gene')
    databaseNames <- probeID2gene$genesUniq[match(query, 
        probeID2gene$probeID)]
    
    databaseNames <- na.omit(unique(
        unlist(lapply(databaseNames,
            function(databaseName) {
                strsplit(databaseName, ";")
            }))))
    
    if (length(databaseNames) == 0) return(NULL)
    
    databases <- databases[names(databases) %in% databaseNames]
    testEnrichment(query, databases, platform=platform)
}


#' testEnrichmentFisher uses Fisher's exact test to estimate the association
#' between two categorical variables.
#'
#' @param query Vector of probes of interest (e.g., significant probes)
#' @param database Vectors corresponding to the database set of
#' interest with associated meta data as an attribute to each element.
#' @param universe Vector of probes in the universe set containing all of
#' the probes to be considered in the test. (Default: NULL)
#' 
#' @import stats
#' 
#' @return A DataFrame with the estimate/statistic, p-value, and name of test
#' for the given results.
testEnrichmentFisher <- function(query, database, universe) {
    q_and_d <- length(intersect(query, database))
    ## if (q_and_d == 0) {
    ##     return(data.frame(
    ##         estimate = NA,
    ##         p.value = 1,
    ##         test = "Fisher",
    ##         nQ = length(query),
    ##         nD = length(database),
    ##         overlap = 0
    ##     ))
    ## }
    l_d <- length(database)
    l_q <- length(query)
    d_min_q <- l_d - q_and_d
    q_min_d <- l_q - q_and_d
    min_q_d <- length(universe) - l_q - l_d + q_and_d
    
    res <- fisher.test(matrix(c(
        q_and_d, d_min_q, q_min_d, min_q_d), nrow = 2))
    
    data.frame(
        estimate = (
            q_and_d / (q_and_d + q_min_d) / (q_and_d + d_min_q) *
                (q_and_d + q_min_d + d_min_q + min_q_d)),
        p.value = res$p.value,
        test = "Fisher",
        nQ = length(query),
        nD = length(database),
        overlap = q_and_d)
}

#' testEnrichmentFGSEA uses the FGSEA test to estimate the association of a
#' categorical variable against a continuous variable.
#'
#' @param query Vector of probes of interest (e.g., significant probes)
#' @param database Vector of probes corresponding to a single database set
#' of interest.
#' @param p.value.adj Logical value indicating whether to report the adjusted
#' p-value. (Default: FALSE).
#' @param estimate.type String indicating the estimate to report. Optional.
#' (Default: "ES").
#' 
#' @import fgsea
#'
#' @return A DataFrame with the estimate/statistic, p-value, and name of test
#' for the given results.
testEnrichmentFGSEA <- function(query, database, p.value.adj=FALSE,
    estimate.type="ES") {
    test <- "fgsea"
    overlap <- length((intersect(names(query), database)))
    if (overlap == 0) {
        return(data.frame(estimate=0,
            p.value=1,
            test=test,
            nQ=length(database),
            nD=length(query),
            overlap=overlap
        ))
    }
    res <- fgsea(pathways=list(pathway=database), 
        stats=na.omit(query))
    
    if (p.value.adj) {
        p.value <- res$padj
    } else {
        p.value <- res$pval
    }
    
    if (estimate.type == "log2err") {
        estimate <- res$log2err
    } else if (estimate.type == "NES") {
        estimate <- res$NES
    } else if (estimate.type == "leadingEdge") {
        estimate <- res$leadingEdge
    } else if (estimate.type == "ES") {
        estimate <- res$ES
    } else {
        message(sprintf("Incorrect estimate.type: [%s].", estimate.type))
        return(NULL)
    }
    
    data.frame(
        estimate = estimate,
        p.value = p.value,
        test = test,
        nQ = length(database),
        nD = length(query),
        overlap = overlap
    )
}


#' testEnrichmentSpearman uses the Spearman statistical test to estimate
#' the association between two continuous variables.
#'
#' @param query Vector of probes of interest (e.g., significant probes)
#' @param database List of vectors corresponding to the database set of
#' interest with associated meta data as an attribute to each element.
#'
#' @import stats
#'
#' @return A DataFrame with the estimate/statistic, p-value, and name
#' of test for the given results.
testEnrichmentSpearman <- function(query, database) {
    test <- "spearman"
    if (length(intersect(names(query), names(database))) == 0) {
        return(data.frame(
            estimate = 0,
            p.value = 1,
            test = test,
            nQ = length(query),
            nD = length(database),
            overlap = 0
        ))
    }
    
    database <- database[match(names(query), names(database))]
    
    res <- cor.test(query, database, method = "spearman")
    data.frame(
        estimate = res$estimate[[1]],
        p.value = res$p.value,
        test = test,
        overlap = length(query)
    )
}

guess_dbnames <- function(nms) {
    df <- sesameDataList()
    df <- df[grep("KYCG", df$Title),]
    vapply(nms, function(nm) {
        if (nm %in% df$Title) {
            return(nm)
        } else if (length(grep(nm, df$Title)) == 1) {
            return(grep(nm, df$Title, value=TRUE))
        } else if (length(grep(nm, df$Title)) == 0) {
            res <- df$Title[apply(do.call(cbind, lapply(
                strsplit(nm, "\\.")[[1]], function(q1) grepl(q1, df$Title))),
                1, all)]
            if (length(res) == 1) {
                return(res[1])
            }
        }
        return(nm)
    }, character(1))
}

#' dbStats builds dataset for a given betas matrix 
#' composed of engineered features from the given database sets
#'
#' @param betas matrix of beta values where probes are on the rows and
#' samples are on the columns
#' @param dbs List of vectors corresponding to probe locations for
#' which the features will be extracted
#' @param fun aggregation function, default to mean
#' @param na.rm whether to remove NA
#' @return matrix with samples on the rows and database set on the columns
#' @export
#' @examples
#' 
#' library(SummarizedExperiment)
#' se <- sesameDataGet('MM285.20Kx467.SE')
#' stats <- dbStats(assay(se), c('KYCG.MM285.probeType.20210630'))
#' head(stats)
#' sesameDataClearCache()
#' 
dbStats <- function(betas, dbs, fun = mean, na.rm = TRUE) {
    if (is(betas, "numeric")) { betas <- cbind(sample = betas); }
    if (is.character(dbs)) {
        nms <- guess_dbnames(dbs)
        dbs <- do.call(c, lapply(nms, sesameDataGet))
    }
    stats <- do.call(cbind, lapply(names(dbs), function(db_nm) {
        db <- dbs[[db_nm]]
        betas1 <- betas[db[db %in% rownames(betas)],,drop=FALSE]
        if (nrow(betas1) == 0) { return(rep(NA, ncol(betas))); }
        apply(betas1, 2, fun, na.rm = na.rm)
    }))
    colnames(stats) <- names(dbs)
    rownames(stats) <- colnames(betas)
    stats
}

#' createGeneNetwork creates database network using the Jaccard index.
#'
#' @param databases Vector of probes corresponding to a single database set
#' of interest.
#' @return ggplot lollipop plot
#' @import reshape2
#'
#' @export
#' @examples
#'
#' databaseNames <- c('KYCG.MM285.seqContextN.20210630')
#' databases <- do.call(c, lapply(databaseNames, sesameDataGet))
#' createDatabaseSetNetwork(databases)
#' sesameDataClearCache()
#'
createDatabaseSetNetwork <- function(databases) {
    m <- compareDatbaseSetOverlap(databases, metric="jaccard")
    
    m_melted <- melt(m)
    colnames(m_melted) <- c("gene1", "gene2", "metric")
    m_melted <- m_melted[m_melted$metric != 0, ]
    
    nodes <- data.frame(id=colnames(m),
        stringsAsFactors=FALSE)
    edges <- data.frame(source=m_melted$gene1,
        target=m_melted$gene2,
        weight=m_melted$metric, # numeric
        stringsAsFactors=FALSE)
    
    list(nodes=nodes, edges=edges)
}

#' calculates the pariwise overlap between given list
#' of database sets using a distance metric.
#'
#' @param databases List of vectors corresponding to the database sets of
#' interest with associated meta data as an attribute to each element. Optional.
#' (Default: NA)
#' @param metric String representing the similarity metric to use. Optional.
#' (Default: "Jaccard").
#' @param verbose Logical value indicating whether to display intermediate
#' text output about the type of test. Optional. (Default: FALSE)
#' @return An upper triangular matrix containing a metric (Jaccard) comparing
#' the pairwise distances between database sets.
#' @export
#' @examples
#' 
#' databaseNames <- c('KYCG.MM285.seqContextN.20210630')
#' databases <- do.call(c, lapply(databaseNames, sesameDataGet))
#' compareDatbaseSetOverlap(databases)
#' sesameDataClearCache()
#'
compareDatbaseSetOverlap <- function(
    databases=NA, metric = "Jaccard", verbose = FALSE) {
    
    ndatabases <- length(databases)
    names <- names(databases)
    m <- matrix(0, nrow=ndatabases, ncol=ndatabases)
    colnames(m) <- names
    rownames(m) <- names
    for (i in seq(ndatabases - 1)) { 
        for (j in seq(i + 1, ndatabases)) {
            message(i, " ", j, '\n')
            m[i, j] <- length(intersect(databases[[i]], databases[[j]])) /
                length(union(databases[[i]], databases[[j]]))
        }
    }
    m
}
