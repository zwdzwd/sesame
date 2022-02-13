
databases_getMeta <- function(dbs) {
    meta <- do.call(bind_rows, lapply(dbs, function(db) {
        m1 <- attributes(db)
        m1 <- m1[names(m1) != "names"] # a special attribute for continuous db
        if ("meta" %in% names(m1)) { # backward compatibility, to delete
            m1 <- c(m1[!(names(m1) %in% c("meta"))], m1$meta)
        } else {
            m1 <- m1[!(names(m1) %in% c("meta"))]
        }
        
        if (is.null(m1)) {
            data.frame(hasMeta = FALSE)
        } else {
            c(m1, hasMeta = TRUE)
        }
    }))
    meta[,colnames(meta)!="hasMeta"]
}

## testEnrichment1 <- function(query, database, universe) {
    
##     if (is.numeric(query)) { # a named vector of continuous value
##         if(is.numeric(database)) { # numeric db
##             res <- testEnrichmentSpearman(query=query, database=database)
##         } else {
##             res <- testEnrichmentGSEA(query = query, database = database)
##         }
##     } else if (is.character(query)) { # categorical query
##         if(is.numeric(database)) { # numeric db
##             res <- testEnrichmentGSEA(query = database, database = query)
##         } else { # categorical db
##             res <- testEnrichmentFisher(query = query, database = database,
##                 universe = universe)
##         }
##     } else {
##         stop("Query is neither numerical or categorical.")
##     }
##     res
## }

queryCheckPlatform <- function(platform, query = NULL) {
    if (is.null(platform)) {
        stopifnot(!is.null(query))
        if (is.numeric(query)) {
            platform <- inferPlatformFromProbeIDs(names(query))
        } else {
            platform <- inferPlatformFromProbeIDs(query)
        }
    }
    platform
}

inferUniverse <- function(platform) {
    mfts <- c(
        "MM285.address", "EPIC.address",
        "Mammal40.address", "HM450.address", "HM27.address")
    mft <- mfts[grepl(platform, mfts)]
    stopifnot(length(mft) == 1 && all(mft %in% mfts))
    sesameDataGet(mft)$ordering$Probe_ID
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
#' @param alternative "two.sided", "greater", or "less"
#' @param platform String corresponding to the type of platform to use. Either
#' MM285, EPIC, HM450, or HM27. If it is not provided, it will be inferred
#' from the query set probeIDs (Default: NA).
#' @param silent output message? (Default: FALSE)
#' @return One list containing features corresponding the test estimate,
#' p-value, and type of test.
#' @importFrom dplyr bind_rows
#' @examples
#' 
#' library(SummarizedExperiment)
#' df <- rowData(sesameDataGet('MM285.tissueSignature'))
#' query <- df$Probe_ID[df$branch == "B_cell"]
#' res <- testEnrichment(query, "chromHMM")
#' sesameDataGet_resetEnv()
#'
#' @export
testEnrichment <- function(
    query, databases = NULL, universe = NULL, alternative = "greater",
    platform = NULL, silent = FALSE) {

    platform <- queryCheckPlatform(platform, query)
    if (is.null(universe)) {
        universe <- inferUniverse(platform) }
    
    if (is.null(databases)) {
        dbs <- c(KYCG_getDBs(KYCG_listDBGroups( # by default, all dbs + gene
            platform, type="categorical")$Title, silent = silent),
            KYCG_buildGeneDBs(query, platform, silent = silent))
    } else if (is.character(databases)) {
        dbs <- KYCG_getDBs(databases, platform = platform, silent = silent)
    } else {
        dbs <- databases
    }
    ## there shouldn't be empty databases, but just in case
    dbs <- dbs[vapply(dbs, length, integer(1)) > 0]
    if (!silent) {
        message(sprintf("Testing against %d database(s)...", length(dbs)))
    }
    
    res <- do.call(bind_rows, lapply(dbs, function(db) {
        testEnrichmentFisher(query = query, database = db,
            universe = universe, alternative = alternative)}))

    ## adjust p.value after merging
    res$FDR <- p.adjust(res$p.value, method='fdr')
    rownames(res) <- NULL

    ## bind meta data
    res <- cbind(res, databases_getMeta(dbs))
    res[order(res$log10.p.value, -abs(res$estimate)), ]
}

#' Aggregate test enrichment results
#'
#' @param result_list a list of results from testEnrichment
#' @param column the column name to aggregate (Default: estimate)
#' @param return_df whether to return a merged data frame
#' @return a matrix for all results
#' @importFrom reshape2 melt
#' @examples
#'
#' ## pick some big TFBS-overlapping CpG groups
#' cg_lists <- KYCG_getDBs("MM285.TFBS")
#' queries <- cg_lists[(sapply(cg_lists, length) > 40000)]
#' result_list <- lapply(queries, testEnrichment, "MM285.chromHMM")
#' mtx <- aggregateTestEnrichments(result_list)
#' 
#' @export
aggregateTestEnrichments <- function(
    result_list, column = "estimate", return_df = FALSE) {
    mtx <- do.call(cbind, lapply(result_list[[1]]$dbname, function(db) {
        vapply(result_list,
            function(x) x$estimate[x$dbname == db], numeric(1))}))
    colnames(mtx) <- result_list[[1]]$dbname
    if (return_df) {
        melt(mtx, value.name = column, varnames = c("query", "db"))
    } else {
        mtx
    }
}

#' build gene-probe association database
#'
#' @param query the query probe list. If NULL, use all the probes
#' on the platform
#' @param platform HM450, EPIC, MM285, Mammal40, will infer from
#' query if not given
#' @param max_distance probe-gene distance for association
#' @param silent suppress messages
#' @return gene databases
#' @importFrom GenomicRanges findOverlaps
#' @importFrom S4Vectors subjectHits
#' @importFrom S4Vectors queryHits
#' @examples
#' query <- c("cg04707299", "cg13380562", "cg00480749")
#' dbs <- KYCG_buildGeneDBs(query)
#' testEnrichment(query, dbs)
#' @export
KYCG_buildGeneDBs <- function(
    query = NULL, platform = NULL, max_distance = 10000, silent = FALSE) {
    
    platform <- queryCheckPlatform(platform, query)
    genes <- sesameData_txnToGeneGRanges(
        sesameData_getTxnGRanges(
            sesameData_check_genome(NULL, platform)))
    all_probes <- sesameData_getManifestGRanges(platform)
    if (!is.null(query)) {
        probes <- all_probes[names(all_probes) %in% query] }

    ## skip non-overlapping genes, strand always ignored
    genes <- subsetByOverlaps(
        genes, probes + max_distance, ignore.strand = TRUE)
    hits <- findOverlaps(
        genes, all_probes + max_distance, ignore.strand = TRUE)
    dbs <- split(names(all_probes)[subjectHits(hits)],
        names(genes)[queryHits(hits)])
    gene_names <- genes[names(dbs)]$gene_name
    res <- lapply(seq_along(dbs), function(i) {
        d1 <- dbs[[i]];
        attr(d1, "group") <- sprintf("KYCG.%s.gene.00000000", platform);
        attr(d1, "dbname") <- names(dbs)[i];
        attr(d1, "gene_name") <- gene_names[i];
        d1;})
    names(res) <- names(dbs)
    message(sprintf("Building %d gene DBs for %s...", length(res), platform))
    res
}

#' testEnrichmentFisher uses Fisher's exact test to estimate the association
#' between two categorical variables.
#'
#' Estimates log2 Odds ratio
#'
#' @param query Vector of probes of interest (e.g., significant probes)
#' @param database Vectors corresponding to the database set of
#' interest with associated meta data as an attribute to each element.
#' @param universe Vector of probes in the universe set containing all of
#' @param alternative greater or two.sided (default: greater)
#' the probes to be considered in the test. (Default: NULL)
#' 
#' @import stats
#' 
#' @return A DataFrame with the estimate/statistic, p-value, and name of test
#' for the given results.
testEnrichmentFisher <- function(query, database, universe,
    alternative = "greater") {

    l_d <- length(database)
    l_q <- length(query)
    qd <- length(intersect(query, database))
    dnq <- l_d - qd
    qnd <- l_q - qd
    nqd <- length(universe) - l_q - l_d + qd

    if (alternative == "two.sided") {
        pval_g <- phyper(qd-1, qd + qnd, nqd + dnq, dnq + qd,
            lower.tail = FALSE, log.p = TRUE) / log(10)
        pval_l <- phyper(qd, qd + qnd, nqd + dnq, dnq + qd,
            lower.tail = TRUE, log.p = TRUE) / log(10)
        log10.p.value <- pmin(pmin(pval_g, pval_l) + log(2), 0) / log(10)
        ## log10.p.value <- log10(fisher.test(matrix(c(
        ##     qd, dnq, qnd, nqd), nrow = 2))$p.value)
    } else if (alternative == "greater") {
        log10.p.value <- phyper(qd-1, qd + qnd, nqd + dnq, dnq + qd,
            lower.tail = FALSE, log.p = TRUE) / log(10)
    } else if (alternative == "less") {
        log10.p.value <- phyper(qd, qd + qnd, nqd + dnq, dnq + qd,
            lower.tail = TRUE, log.p = TRUE) / log(10)
    } else { stop("alternative must be either greater, less or two-sided.") }
    
    odds_ratio <- qd / qnd / dnq * nqd # can be NaN if 0
    data.frame(
        estimate = log2(odds_ratio), p.value = 10**(log10.p.value),
        log10.p.value = log10.p.value, test = "Log2(OR)",
        nQ = length(query), nD = length(database), overlap = qd)
}

calcES <- function(dCont, dDisc) {
    presence <- names(dCont) %in% dDisc
    s <- ifelse(presence, 1/sum(presence), -1/sum(!presence))
    cs <- cumsum(s)
    list(es_max = max(cs), es_min = min(cs))
}

calcES_Significance <- function(dCont, dDisc, permut=100, precise=FALSE) {

    dCont <- sort(dCont)
    dContName <- names(dCont)
    dDiscN <- length(dDisc)
    dContN <- length(dCont)
    s <- rep(-1/(dContN-dDiscN), dContN)

    ess <- do.call(rbind, lapply(seq_len(permut), function(i) {
        s[sample.int(dContN, dDiscN)] <- 1/dDiscN
        cs <- cumsum(s)
        data.frame(es_max = max(cs), es_min = min(cs))
    }))

    es <- calcES(dCont, dDisc)
    res <- list(
        es_small = es$es_max,
        es_large = -es$es_min,
        pv_small = 1-ecdf(ess$es_max)(es$es_max),
        pv_large = ecdf(ess$es_min)(es$es_min))
    
    if (res$pv_small < 0.01 || res$pv_large < 0.01) {
        if (permut < 1000 && precise) {
            res <- calcES_Significance(dCont, dDisc, permut = 1000)
        } else {
            ## high precisions are approximated by Gaussian
            ## TODO: should also report log.p=TRUE in the future
            if (res$pv_small == 0) {
                res$pv_small <- pnorm(
                    es$es_max, mean=mean(ess$es_max),
                    sd=sd(ess$es_max), lower.tail=FALSE)
            }
            if (res$pv_large == 0) {
                res$pv_large <- pnorm(
                    es$es_max, mean=mean(ess$es_max),
                    sd=sd(ess$es_max), lower.tail=TRUE)
            }
        }
    }
    res
}

testEnrichmentGSEA1 <- function(query, database, precise=FALSE) {
    test <- "GSEA Enrichment Score"
    overlap <- intersect(names(query), database)

    if (length(overlap) != length(database)) {
        warning("Not every data in database has query.")
        warning(sprintf("Using %d in %d data for testing.",
            length(overlap), length(database)))
        warning("Be careful interpreting results.")
    }
    
    if (length(overlap) == 0 || length(overlap) == length(query)) {
        return(data.frame(
            estimate = 0, p.value = 1,
            log10.p.value = 0, test = test,
            nQ = length(database), nD = length(query),
            overlap = length(overlap))) }

    res <- calcES_Significance(query, overlap, precise=precise)

    if (res$es_large > res$es_small) {
        ## negative sign represent enrichment for the large end
        data.frame(
            estimate = -res$es_large, p.value = res$pv_large,
            log10.p.value = log10(res$pv_large), test = test,
            nQ = length(database), nD = length(query),
            overlap = length(overlap))
    } else {
        data.frame(
            estimate = res$es_small, p.value = res$pv_small,
            log10.p.value = log10(res$pv_small), test = test,
            nQ = length(database), nD = length(query),
            overlap = length(overlap))
    }
}
    
#' testEnrichmentGSEA uses the GSEA test to estimate the association of a
#' categorical variable against a continuous variable.
#'
#' estimate represent enrichment score and negative estimate indicate a
#' test for depletion
#'
#' @param query query, if numerical, expect categorical database, if
#' categorical expect numerical database
#' @param databases database, numerical or categorical, but needs to be
#' different from query
#' @param platform EPIC, MM285, ..., infer if not given
#' @param silent suppress message (default: FALSE)
#' @param precise whether to compute precise p-value (up to numerical limit)
#' of interest.
#' @return A DataFrame with the estimate/statistic, p-value, and name of test
#' for the given results.
#' @examples
#' query <- KYCG_getDBs("KYCG.MM285.designGroup")[["TSS"]]
#' res <- testEnrichmentGSEA(query, "MM285.seqContextN")
#' @export
testEnrichmentGSEA <- function(query, databases = NULL,
    platform = NULL, silent = FALSE, precise = FALSE) {

    platform <- queryCheckPlatform(platform, query)
    if (is.null(databases)) {
        if (is.character(query)) {
            dbs <- KYCG_getDBs(KYCG_listDBGroups( # by default, all dbs
                platform, type="numerical")$Title, silent = silent)
        } else if (is.numeric(query)) {
            dbs <- KYCG_getDBs(KYCG_listDBGroups( # all dbs, could be slow
                platform, type="categorical")$Title, silent = silent)
        }
    } else if (is.character(databases)) {
        dbs <- KYCG_getDBs(databases, platform = platform, silent = silent)
    } else {
        dbs <- databases
    }
    ## there shouldn't be empty databases, but just in case
    dbs <- dbs[vapply(dbs, length, integer(1)) > 0]
    if (!silent) {
        message(sprintf("Testing against %d database(s)...", length(dbs))) }

    if (is.character(query) && all(vapply(dbs, is.numeric, logical(1)))) {
        res <- do.call(bind_rows, lapply(dbs, function(db) {
            testEnrichmentGSEA1(
                query = db, database = query, precise = precise)}))
    } else if (
        is.numeric(query) && all(vapply(dbs, is.character, logical(1)))) {
        res <- do.call(bind_rows, lapply(dbs, function(db) {
            testEnrichmentGSEA1(
                query = query, database = db, precise = precise)}))
    } else { stopifnot("query must be numerical or categorical"); }
    
    ## adjust p.value after merging
    res$FDR <- p.adjust(res$p.value, method='fdr')
    rownames(res) <- NULL

    ## bind meta data
    res <- cbind(res, databases_getMeta(dbs))
    res[order(res$p.value, -abs(res$estimate)), ]
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
    test <- "Spearman's rho"
    if (length(intersect(names(query), names(database))) == 0) {
        return(data.frame(
            estimate = 0, p.value = 1,
            log10.p.value = 0, test = test,
            nQ = length(query), nD = length(database),
            overlap = 0))
    }
    
    database <- database[match(names(query), names(database))]
    
    res <- cor.test(query, database, method = "spearman")
    data.frame(
        estimate = res$estimate[[1]], p.value = res$p.value,
        log10.p.value = log10(res$estimate[[1]]), test = test,
        nQ = length(query), nD = length(database),
        overlap = length(query))
}

guess_dbnames <- function(nms, allow_multi=FALSE) {
    df <- sesameDataList()
    df <- df[grep("KYCG", df$Title),]
    do.call(c, lapply(nms, function(nm) {
        if (nm %in% df$Title) {
            return(nm)
        } else if (length(grep(nm, df$Title)) >= 1) {
            ret <- grep(nm, df$Title, value=TRUE)
            if (!allow_multi) { ret <- ret[1]; }
            return(ret)
        } else if (length(grep(nm, df$Title)) == 0) {
            res <- df$Title[apply(do.call(cbind, lapply(
                strsplit(nm, "\\.")[[1]], function(q1) grepl(q1, df$Title))),
                1, all)]
            if (length(res) == 1) {
                return(res[1])
            }
        }
        return(nm)
    }))
}

#' List database group names
#'
#' @param filter keywords for filtering
#' @param type categorical, numerical (default: all)
#' @return a list of db group names
#' @examples
#' head(KYCG_listDBGroups("chromHMM"))
#' @export
KYCG_listDBGroups <- function(
    filter = NULL, type = NULL) {
    
    gps <- sesameDataList("KYCG", full=TRUE)[,c("Title","Description")]
    gps$type <- vapply(strsplit(
        gps$Description, " "), function(x) x[2], character(1))
    gps$Description <- str_replace(
        gps$Description, "KYCG categorical database holding ", "")
    if (!is.null(filter)) {
        gps <- gps[grepl(filter, gps$Title),]
    }
    if (!is.null(type)) {
        gps <- gps[gps$type %in% type,]
    }
    gps
}

#' Get databases by full or partial names of the database group(s)
#'
#' @param group_nms database group names
#' @param db_names name of the database, fetech only the given databases
#' @param platform EPIC, HM450, MM285, ... If given, will restrict to
#' that platform.
#' @param summary return a summary of database instead of db itself
#' @param allow_multi allow multiple groups to be returned for
#' @param silent no messages
#' each query.
#' @return a list of databases
#' @examples
#' dbs <- KYCG_getDBs("MM285.chromHMM")
#' dbs <- KYCG_getDBs(c("MM285.chromHMM", "MM285.probeType"))
#' @export
KYCG_getDBs <- function(group_nms, db_names = NULL, platform = NULL,
    summary = FALSE, allow_multi = FALSE, silent = FALSE) {
    
    if (!is.character(group_nms)) {
        return(group_nms)
    }
    group_nms <- guess_dbnames(group_nms, allow_multi = TRUE)
    if (!is.null(platform)) {
        group_nms <- grep(platform, group_nms, value = TRUE) }
    group_nms_ <- paste(group_nms, sep="\n")
    if (!silent) {
        message("Selected the following database groups:")
        invisible(lapply(seq_along(group_nms_), function(i) {
            message(sprintf("%d. %s", i, group_nms_[i]))
        }))}
    res <- do.call(c, lapply(unname(group_nms), function(nm) {
        dbs <- sesameDataGet(nm)
        setNames(lapply(seq_along(dbs), function(ii) {
            db <- dbs[[ii]]
            attr(db, "group") <- nm
            attr(db, "dbname") <- names(dbs)[ii]
            db
        }), names(dbs))}))

    if (summary) {
        do.call(bind_rows, lapply(res, attributes))
    } else if (is.null(db_names)) {
        res
    } else {
        stopifnot(all(db_names %in% names(res)))
        res[db_names]
    }
}

#' dbStats builds dataset for a given betas matrix 
#' composed of engineered features from the given database sets
#'
#' @param betas matrix of beta values where probes are on the rows and
#' samples are on the columns
#' @param databases List of vectors corresponding to probe locations for
#' which the features will be extracted
#' @param fun aggregation function, default to mean
#' @param na.rm whether to remove NA
#' @param f_min min fraction of non-NA for aggregation function to apply
#' @param n_min min number of non-NA for aggregation function to apply,
#' overrides f_min
#' @return matrix with samples on the rows and database set on the columns
#' @examples
#' 
#' library(SummarizedExperiment)
#' se <- sesameDataGet('MM285.20Kx467.SE')
#' head(dbStats(assay(se), "MM285.probeType")[,1:3])
#' sesameDataGet_resetEnv()
#' 
#' @export
dbStats <- function(
    betas, databases, fun = mean, na.rm = TRUE, n_min = NULL, f_min = 0.1) {

    if (is(betas, "numeric")) { betas <- cbind(sample = betas); }
    if (is.character(databases)) {
        dbs <- KYCG_getDBs(databases)
    } else {
        dbs <- databases
    }
    stats <- do.call(cbind, lapply(names(dbs), function(db_nm) {
        db <- dbs[[db_nm]]
        betas1 <- betas[db[db %in% rownames(betas)],,drop=FALSE]
        n_probes <- nrow(betas1)
        if (n_probes == 0) { return(rep(NA, ncol(betas))); }
        nacnt <- colSums(!is.na(betas1), na.rm = TRUE)
        stat1 <- apply(betas1, 2, fun, na.rm = na.rm)
        if(is.null(n_min)) {
            n_min1 <- n_probes * f_min
        } else {
            n_min1 <- n_min
        }
        stat1[nacnt < n_min1] <- NA
        stat1
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
createDBNetwork <- function(databases) {
    m <- compareDatbaseSetOverlap(databases, metric="jaccard")

    ## databaseNames <- c('KYCG.MM285.seqContextN.20210630')
    ## databases <- do.call(c, lapply(databaseNames, sesameDataGet))
    ## createDatabaseSetNetwork(databases)
    ## sesameDataGet_resetEnv()
    
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
#' @return An upper triangular matrix containing a metric (Jaccard) comparing
#' the pairwise distances between database sets.
compareDatbaseSetOverlap <- function(
    databases=NA, metric = "Jaccard") {

    ## databaseNames <- c('KYCG.MM285.seqContextN.20210630')
    ## databases <- do.call(c, lapply(databaseNames, sesameDataGet))
    ## compareDatbaseSetOverlap(databases)
    ## sesameDataGet_resetEnv()

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
