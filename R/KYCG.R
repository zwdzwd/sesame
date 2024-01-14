
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

queryCheckPlatform <- function(platform, query = NULL, silent = FALSE) {
    if (is.null(platform)) {
        stopifnot(!is.null(query))
        if (is.numeric(query)) {
            platform <- inferPlatformFromProbeIDs(names(query), silent = silent)
        } else {
            platform <- inferPlatformFromProbeIDs(query, silent = silent)
        }
    }
    platform
}

inferUniverse <- function(platform) {
    mfts <- c(
        "MM285.address", "EPIC.address", "EPICv2.address",
        "Mammal40.address", "HM450.address", "HM27.address")
    mft <- paste0(platform, ".address")
    if (!(mft %in% mfts)) {
        stop("Platform unsupported. Please provide universe set explicitly.")
    }
    stopifnot()
    sesameDataGet(mft)$ordering$Probe_ID
}

subsetDBs <- function(dbs, universe) {
    dbs <- lapply(dbs, function(db) {
        db1 <- intersect(db, universe)
        attributes(db1) <- attributes(db)
        db1
    })
    dbs <- dbs[length(dbs) > 0]
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
#' @param include_genes include gene link enrichment testing
#' @param platform String corresponding to the type of platform to use. Either
#' MM285, EPIC, HM450, or HM27. If it is not provided, it will be inferred
#' from the query set probeIDs (Default: NA).
#' @param silent output message? (Default: FALSE)
#' @return A data frame containing features corresponding to the test estimate,
#' p-value, and type of test.
#' @importFrom dplyr bind_rows
#' @examples
#' 
#' library(SummarizedExperiment)
#' df <- rowData(sesameDataGet('MM285.tissueSignature'))
#' query <- df$Probe_ID[df$branch == "B_cell"]
#' res <- testEnrichment(query, "chromHMM", platform="MM285")
#' sesameDataGet_resetEnv()
#'
#' @export
testEnrichment <- function(
    query, databases = NULL, universe = NULL, alternative = "greater",
    include_genes = FALSE, platform = NULL, silent = FALSE) {

    platform <- queryCheckPlatform(platform, query, silent = silent)
    
    if (is.null(databases)) {
        dbs <- c(KYCG_getDBs(KYCG_listDBGroups( # by default, all dbs + gene
            platform, type="categorical")$Title, silent = silent))
    } else if (is.character(databases)) {
        dbs <- KYCG_getDBs(databases, platform = platform, silent = silent)
    } else {
        dbs <- databases
    }

    if (include_genes) {
        dbs <- c(dbs, KYCG_buildGeneDBs(query, platform, silent = silent))
    }
    
    ## there shouldn't be empty databases, but just in case
    dbs <- dbs[vapply(dbs, length, integer(1)) > 0]
    if (!silent) {
        message(sprintf("Testing against %d database(s)...", length(dbs)))
    }

    if (is.null(universe)) {
        universe <- inferUniverse(platform)
    } else { # subset the dbs by universe
        dbs <- subsetDBs(dbs, universe) }
    
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

#' Convenient function for testing enrichment of gene linkage
#'
#' @param query probe set of interest
#' @param platform string corresponding to the type of platform to use. Either
#' MM285, EPIC, HM450, or HM27. If it is not provided, it will be inferred
#' from the query set probe IDs.
#' @param silent whether to output message
#' @param ... addition argument provided to testEnrichment
#' @return A data frame containing features corresponding to the test estimate,
#' p-value, and type of test etc.
#' @examples
#' query <- c("cg04707299", "cg13380562", "cg00480749")
#' testEnrichment(query, platform = "EPIC")
#' @export
testEnrichmentGene <- function(query, platform = NULL, silent = FALSE, ...) {
    platform <- queryCheckPlatform(platform, query, silent = silent)
    dbs_gene <- KYCG_buildGeneDBs(query, platform)
    testEnrichment(query, databases = dbs_gene, platform = platform, ...)
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
#' @param genome hg38, mm10, ..., will infer if not given.
#' For additional mapping, download the GRanges object from
#' http://zwdzwd.github.io/InfiniumAnnotation
#' and provide the following argument
#' ..., genome = sesameAnno_buildManifestGRanges("downloaded_file"),...
#' to this function.
#' @param max_distance probe-gene distance for association
#' @param silent suppress messages
#' @return gene databases
#' @importFrom GenomicRanges findOverlaps
#' @importFrom S4Vectors subjectHits
#' @importFrom S4Vectors queryHits
#' @examples
#' query <- c("cg04707299", "cg13380562", "cg00480749")
#' dbs <- KYCG_buildGeneDBs(query, platform = "EPIC")
#' testEnrichment(query, dbs, platform = "EPIC")
#' @export
KYCG_buildGeneDBs <- function(
    query = NULL, platform = NULL,
    genome = NULL, max_distance = 10000, silent = FALSE) {
    
    platform <- queryCheckPlatform(platform, query, silent = silent)
    genes <- sesameData_getTxnGRanges(
        sesameData_check_genome(NULL, platform), merge2gene=TRUE)
    all_probes <- sesameData_getManifestGRanges(platform, genome = genome)
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

    nD <- length(database)
    nQ <- length(query)
    nDQ <- length(intersect(query, database))
    nU <- length(universe)

    testEnrichmentFisherN(nD, nQ, nDQ, nU, alternative = alternative)
}

testEnrichmentFisherN <- function(
    nD, nQ, nDQ, nU, alternative = "greater") {
    
    nDmQ <- nD - nDQ
    nQmD <- nQ - nDQ
    nUmDQ <- nU - nQ - nD + nDQ

    if (alternative == "two.sided") {
        pvg <- phyper(
            nDQ-1, nDQ + nQmD, nUmDQ + nDmQ, nDmQ + nDQ,
            lower.tail = FALSE, log.p = TRUE) / log(10)
        pvl <- phyper(
            nDQ, nDQ + nQmD, nUmDQ + nDmQ, nDmQ + nDQ,
            lower.tail = TRUE, log.p = TRUE) / log(10)
        log10.p.value <- pmin(pmin(pvg, pvl) + log(2), 0) / log(10)
        ## log10.p.value <- log10(fisher.test(matrix(c(
        ##     nDQ, nDmQ, nQmD, nUmDQ), nrow = 2))$p.value)
    } else if (alternative == "greater") {
        log10.p.value <- phyper(
            nDQ-1, nDQ + nQmD, nUmDQ + nDmQ, nDmQ + nDQ,
            lower.tail = FALSE, log.p = TRUE) / log(10)
    } else if (alternative == "less") {
        log10.p.value <- phyper(
            nDQ, nDQ + nQmD, nUmDQ + nDmQ, nDmQ + nDQ,
            lower.tail = TRUE, log.p = TRUE) / log(10)
    } else {
        stop("alternative must be either greater, less or two-sided.")
    }
    
    odds_ratio <- nDQ / nQmD / nDmQ * nUmDQ # can be NaN if 0
    odds_ratio[odds_ratio == Inf] <- .Machine$double.xmax
    odds_ratio[odds_ratio == 0] <- .Machine$double.xmin
    data.frame(
        estimate = log2(odds_ratio),
        p.value = 10**(log10.p.value),
        log10.p.value = log10.p.value,
        test = "Log2(OR)",
        nQ = nQ, nD = nD, overlap = nDQ,
        cf_Jaccard = nDQ / (nD + nQmD),
        cf_overlap = nDQ / pmin(nD, nQ), # Szymkiewicz–Simpson
        cf_NPMI = (log2(nD)+log2(nQ)-2*log2(nU))/(log2(nDQ)-log2(nU))-1,
        cf_SorensenDice = 2 * nDQ/(nD + nQ))
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

    ## es <- calcES(dCont, dDisc)
    presence <- names(dCont) %in% dDisc
    s <- ifelse(presence, 1/sum(presence), -1/sum(!presence))
    cs <- cumsum(s)
    es_max <- max(cs)
    es_min <- min(cs)
    res <- list(es_small = es_max,
                es_large = -es_min,
                pv_small = 1-ecdf(ess$es_max)(es_max),
                pv_large = ecdf(ess$es_min)(es_min))

    ## if significant, try to be more precise
    if (res$pv_small < 0.01 || res$pv_large < 0.01) {
        if (permut < 1000 && precise) {
            res <- calcES_Significance(dCont, dDisc, permut = 1000)
        } else { # approximated by Gaussian (TODO: also report log.p=TRUE)
            if (res$pv_small == 0) {
                res$pv_small <- pnorm(
                    es_max, mean=mean(ess$es_max),
                    sd=sd(ess$es_max), lower.tail=FALSE) }
            if (res$pv_large == 0) {
                res$pv_large <- pnorm(
                    es_max, mean=mean(ess$es_max),
                    sd=sd(ess$es_max), lower.tail=TRUE)
            }}}
    
    res
}

testEnrichmentSEA1 <- function(query, database, precise=FALSE, full=FALSE) {
    
    test <- "Set Enrichment Score"
    overlap <- intersect(names(query), database)
    if (length(overlap) != length(database)) {
        warning("Not every data in database has query.")
        warning(sprintf("Using %d in %d data for testing.",
            length(overlap), length(database)))
    }
    
    if (length(overlap) == 0 || length(overlap) == length(query)) {
        return(data.frame(
            estimate = 0, p.value = 1,
            log10.p.value = 0, test = test,
            nQ = length(database), nD = length(query),
            overlap = length(overlap))) }

    res <- calcES_Significance(query, overlap, precise=precise)

    if (res$es_large > res$es_small) {
        df <- data.frame( ## negative sign represent enriching for large values
            estimate = -res$es_large, p.value = res$pv_large,
            log10.p.value = log10(res$pv_large), test = test,
            nQ = length(database), nD = length(query),
            overlap = length(overlap))
    } else {
        df <- data.frame(
            estimate = res$es_small, p.value = res$pv_small,
            log10.p.value = log10(res$pv_small), test = test,
            nQ = length(database), nD = length(query),
            overlap = length(overlap))
    }
    
    if (full) {
        list(res = df, dCont = query, dDisc = overlap)
    } else {
        df
    }
}
    
#' uses the GSEA-like test to estimate the association of a
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
#' @param prepPlot return the raw enrichment scores and presence vectors
#' for plotting
#' @return A DataFrame with the estimate/statistic, p-value, and name of test
#' for the given results.
#' @examples
#' query <- KYCG_getDBs("KYCG.MM285.designGroup")[["TSS"]]
#' res <- testEnrichmentSEA(query, "MM285.seqContextN")
#' @export
testEnrichmentSEA <- function(query, databases,
    platform = NULL, silent = FALSE, precise = FALSE, prepPlot = FALSE) {

    platform <- queryCheckPlatform(platform, query, silent = silent)
    stopifnot(!is.null(databases))
    if (is.character(databases)) { # infer database from string
        dbs <- KYCG_getDBs(databases, platform = platform, silent = silent)
    } else {
        dbs <- databases
    }
    dbs <- dbs[vapply(dbs, length, integer(1)) > 0] # no empty db
    if (!silent) {
        message(sprintf("Testing against %d database(s)...", length(dbs))) }

    if (is.character(query) && all(vapply(dbs, is.numeric, logical(1)))) {
        res <- lapply(dbs, function(db) {
            testEnrichmentSEA1(query = db, database = query,
                                precise = precise, full=prepPlot)})
    } else if (
        is.numeric(query) && all(vapply(dbs, is.character, logical(1)))) {
        res <- lapply(dbs, function(db) {
            testEnrichmentSEA1(query = query, database = db,
                                precise = precise, full=prepPlot)})
    } else { stop("query and db must be one numerical and one categorical"); }

    if (prepPlot) { return(res)
    } else { res <- do.call(bind_rows, res) }
    
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

guess_dbnames <- function(nms, platform = NULL,
    allow_multi = FALSE, type = NULL, silent = FALSE,
    ignore.case = FALSE) {
    
    gps <- KYCG_listDBGroups(type = type)
    nms <- do.call(c, lapply(nms, function(nm) {
        if (nm %in% gps$Title) {
            return(nm)
        } else if (length(grep(nm, gps$Title, ignore.case=ignore.case)) >= 1) {
            ret <- grep(nm, gps$Title, value=TRUE, ignore.case=ignore.case)
            if (!allow_multi) { ret <- ret[1]; }
            return(ret)
        } else if (length(grep(nm, gps$Title, ignore.case=ignore.case)) == 0) {
            res <- gps$Title[apply(do.call(cbind, lapply(
                strsplit(nm, "\\.")[[1]],
                function(q1) grepl(q1, gps$Title, ignore.case=ignore.case))),
                1, all)]
            if (length(res) == 1) {
                return(res[1])
            }
        }
        return(nm)
    }))
    if (!is.null(platform)) {
        nms <- grep(platform, nms, value = TRUE)
    }
    if (!silent) {
        message("Selected the following database groups:")
        invisible(lapply(seq_along(nms), function(i) {
            message(sprintf("%d. %s", i, nms[i]))
        }))
    }
    nms
}

#' List database group names
#'
#' @param filter keywords for filtering
#' @param path file path to downloaded knowledgebase sets
#' @param type categorical, numerical (default: all)
#' @return a list of db group names
#' @examples
#' head(KYCG_listDBGroups("chromHMM"))
#' ## or KYCG_listDBGroups(path = "~/Downloads")
#' @export
KYCG_listDBGroups <- function(filter = NULL, path = NULL, type = NULL) {

    if (is.null(path)) {
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
    } else {
        gps <- basename(list.files(path, recursive = TRUE))
    }
    gps
}

## #' A convenience function for downloading knowledgebase sets
## #'
## #' @param platform EPICv2, EPIC, HM450 etc.
## #' @param fdr directory to which feature files will be downloaded
## #' @return untarred feature folders
## #' @export
## KYCG_downloadDBs <- function(platform, fdr) {
##     fdr <- path.expand(fdr)
##     dir.create(fdr)
##     URL <- paste0("https://zhouserver.research.chop.edu/",
##         sprintf("InfiniumAnnotation/%s/annotations.tar.gz", platform))
##     download.file(URL, paste0(fdr,"/tmp.tar.gz"))
##     untar(paste0(fdr,"/tmp.tar.gz"), exdir=paste0(fdr,"/"))
##     unlink(paste0(fdr,"/tmp.tar.gz"))
##     db_groups <- list.files(fdr, recursive=TRUE)
##     message(sprintf("Downloaded %d knowledgebase groups to %s",
##         length(db_groups), fdr))
## }

#' Load database groups
#'
#' @param in_paths folder that contains all databases
#' @param group_use_filename whether to use file name for groups
#' @return a list of db group names
#' @examples
#'
#' ## download regulatory annotations from
#' ## http://zwdzwd.github.io/InfiniumAnnotation
#' ## unzip the file
#' if (FALSE) {
#' dbs <- KYCG_loadDBs(path_to_unzipped_folder)
#' }
#' @export
KYCG_loadDBs <- function(in_paths, group_use_filename=FALSE) {
    if (length(in_paths)==1 && dir.exists(in_paths)) {
        groupnms <- grep(".gz$",
            list.files(in_paths, recursive=TRUE), value=TRUE)
        in_paths <- file.path(in_paths, groupnms)
    } else {
        groupnms <- basename(in_paths)
    }
    do.call(c, lapply(seq_along(groupnms), function(i) {
        tbl <- read.table(in_paths[i],sep="\t",header=TRUE)
        dbs <- split(tbl$Probe_ID, tbl$Knowledgebase)
        lapply(names(dbs), function(gp_dbname) {
            gp_dbname_lst <- str_split(gp_dbname, ";", n=2)[[1]]
            db1 <- dbs[[gp_dbname]];
            ## group names come from file instead of file name
            if (group_use_filename) {
                attr(db1, "group") <- sub(".gz$","",groupnms[i])
            } else {
                attr(db1, "group") <- gp_dbname_lst[[1]]
            }
            if (length(gp_dbname_lst) > 1) {
                attr(db1, "dbname") <- gp_dbname_lst[[2]]
            } else {
                attr(db1, "dbname") <- attr(db1, "group")
            }
            db1;})
    }))
}

#' Get databases by full or partial names of the database group(s)
#'
#' @param group_nms database group names
#' @param db_names name of the database, fetech only the given databases
#' @param platform EPIC, HM450, MM285, ... If given, will restrict to
#' that platform.
#' @param summary return a summary of database instead of db itself
#' @param allow_multi allow multiple groups to be returned for
#' @param ignore.case ignore case or not
#' @param type numerical, categorical, default: all
#' @param silent no messages
#' each query.
#' @return a list of databases, return NULL if no database is found
#' @examples
#' dbs <- KYCG_getDBs("MM285.chromHMM")
#' dbs <- KYCG_getDBs(c("MM285.chromHMM", "MM285.probeType"))
#' @export
KYCG_getDBs <- function(group_nms, db_names = NULL, platform = NULL,
    summary = FALSE, allow_multi = FALSE,
    ignore.case = FALSE, type = NULL, silent = FALSE) {
    
    if (!is.character(group_nms)) {
        return(group_nms)
    }

    group_nms <- guess_dbnames(group_nms, platform = platform,
        allow_multi = TRUE, type = type, silent = silent,
        ignore.case = ignore.case)
    ## group_nms <- group_nms[sesameDataHas(group_nms)]
    group_nms <- group_nms[group_nms %in% sesameDataList()$Title]
    if (length(group_nms) == 0) {
        return(NULL)
    }
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

#' Annotate Probe IDs using KYCG databases
#'
#' see sesameData_annoProbes if you'd like to annotate by genomic coordinates
#' (in GRanges)
#' @param query probe IDs in a character vector
#' @param databases character or actual database (i.e. list of probe IDs)
#' @param db_names specific database (default to all databases)
#' @param platform EPIC, MM285 etc. will infer from probe IDs if not given
#' @param indicator return the indicator matrix instead of a concatenated
#' annotation (in the case of have multiple annotations)
#' @param sep delimiter used in paste
#' @param silent suppress message
#' @return named annotation vector, or indicator matrix
#' @examples
#' query <- names(sesameData_getManifestGRanges("MM285"))
#' anno <- KYCG_annoProbes(query, "designGroup", silent = TRUE)
#' @export
KYCG_annoProbes <- function(query, databases, db_names = NULL,
    platform = NULL, sep = ",", indicator = FALSE, silent = FALSE) {

    platform <- queryCheckPlatform(platform, query, silent = silent)
    if (is.character(databases)) {
        dbs <- KYCG_getDBs(databases, db_names = db_names,
            platform = platform, silent = silent, type = "categorical")
    } else {
        dbs <- databases
        names(dbs) <- vapply(dbs, function(db) {
            paste0(attr(db, "group"), "-", attr(db, "dbname"))}, character(1))
    }

    ind <- do.call(cbind, lapply(dbs, function(db) {
        query %in% db
    }))
    if (indicator) {
        rownames(ind) <- query
        colnames(ind) <- names(dbs)
        return(ind)
    } else {
        anno <- apply(ind, 1, function(x) paste(names(dbs)[x], collapse=sep))
        anno <- ifelse(anno == "", NA, anno)
        names(anno) <- query
        return(anno)
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
#' @param long produce long-form result
#' @return matrix with samples on the rows and database set on the columns
#' @examples
#' 
#' library(SummarizedExperiment)
#' se <- sesameDataGet('MM285.467.SE.tissue20Kprobes')
#' head(dbStats(assay(se), "MM285.chromHMM")[,1:3])
#' sesameDataGet_resetEnv()
#' 
#' @importFrom reshape2 melt
#' @export
dbStats <- function(
    betas, databases, fun = mean, na.rm = TRUE, n_min = NULL,
    f_min = 0.1, long = FALSE) {

    if (is(betas, "numeric")) { betas <- cbind(sample = betas); }
    if (is.character(databases)) {
        dbs <- KYCG_getDBs(databases)
    } else {
        dbs <- databases
    }
    stats <- do.call(cbind, lapply(dbs, function(db) {
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
    if (!is.null(names(dbs))) {
        colnames(stats) <- names(dbs)
    } else { # use attributes instead of names
        colnames(stats) <- vapply(dbs,
            function(x) attr(x, "dbname"), character(1))
    }
    rownames(stats) <- colnames(betas)
    if (long) {
        stats <- melt(stats, varnames = c("query", "db"), value.name = "value")
    }
    stats
}

#' createGeneNetwork creates database network using the Jaccard index.
#'
#' @param databases Vector of probes corresponding to a single database set
#' of interest.
#' @return ggplot lollipop plot
#' @importFrom reshape2 melt
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

