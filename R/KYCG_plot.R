#' plot enrichment test result
#'
#' @param df test enrichment result data frame
#' @param fdr_max maximum fdr for capping
#' @param n_label number of database to label
#' @param min_estimate minimum estimate
#' @param short_label use short label
#' @return grid object
#' @importFrom stringr str_replace
#' @importFrom tibble rownames_to_column
#' @import ggplot2
#' @examples
#' query <- KYCG_getDBs("MM285.designGroup")[["PGCMeth"]]
#' res <- testEnrichment(query, platform="MM285")
#' KYCG_plotEnrichAll(res)
#' 
#' @export
KYCG_plotEnrichAll <- function(
    df, fdr_max = 25, n_label = 15, min_estimate = 0, short_label = TRUE) {

    gp_size <- sort(table(df$group))
    gp_width <- log(2+gp_size)
    e1 <- df[order(factor(df$group, levels=names(gp_size)), df$dbname),]
    e1$inc <- (gp_width / gp_size)[e1$group]
    e1$inc1 <- c(0,ifelse(e1$group[-1] != e1$group[-nrow(e1)], 1, 0))
    e1$inc2 <- cumsum(e1$inc + e1$inc1)

    if (length(grep("^KYCG", e1$group))>0) {
        e1$group <- str_replace(e1$group,"KYCG.","")
        e1$group <- vapply(strsplit(e1$group, "\\."),
            function(x) paste0(x[2:(length(x)-1)], collapse="."), character(1))
    }
    if ("gene_name" %in% colnames(e1)) {
        e1$dbname[e1$group == "gene"] <- e1$gene_name[e1$group == "gene"] }

    e2 <- e1[e1$estimate > min_estimate & e1$FDR < 0.01 ,]
    e2$FDR[e2$FDR < 10**-fdr_max] <- 10**-(fdr_max*1.1)

    e3 <- rownames_to_column(as.data.frame(do.call(rbind, lapply(
        split(e1$inc2, e1$group), function(x)
            c(beg=min(x), middle=mean(x), end=max(x))))), "group")

    inc2 <- FDR <- estimate <- group <- dbname <- beg <- middle <- NULL
    if (short_label) {
        e2$dbname <- vapply(
            strsplit(e2$dbname, ";"), function(x) {
                if(length(x)>1) { x[[2]]; } else { x[[1]]; }}, character(1)) }
    requireNamespace("ggrepel")
    ggplot(e2, aes(inc2, -log10(FDR))) +
        geom_point(aes(size=estimate, color=group), alpha=0.5) +
        ggrepel::geom_text_repel(data = e2[head(order(e2$FDR), n = n_label),],
            aes(label=dbname, color=group), size = 3,
            ## box.padding = unit(0.35, "lines"),
            ## point.padding = unit(0.3, "lines"),
            direction="y", nudge_y=0.2, max.overlaps=100) +
        annotate("text", -1, fdr_max*0.96,
            label="Values above this line are capped.",
            hjust=0, vjust=1, color="grey60") +
        geom_hline(yintercept = fdr_max, linetype="dotted", color="grey60") +
        geom_segment(aes(x = beg, y = 0, xend = end, yend = 0, color=group),
            size=3, data=e3) +
        geom_text(data=e3,aes(middle, -1, label=group, color=group),
            vjust=1, hjust=1, angle=30) + scale_color_discrete(guide="none") +
        ylim(-6, fdr_max*1.2) + xlab("") +
        scale_size_continuous(guide=guide_legend(title="log2(OR)")) +
        coord_cartesian(clip="off") + theme_minimal() +
        theme(axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            panel.grid.minor.x = element_blank())
}

#' @importFrom dplyr slice_min
#' @importFrom dplyr ungroup
preparePlotDF <- function(
    df, n, order_by, short_label = FALSE, label_by = "dbname") {
    ## suppress R CMD CHECK no visible binding warning
    db1 <- FDR <- NULL
    
    stopifnot("estimate" %in% colnames(df) && "FDR" %in% colnames(df))
    df1 <- df[df$nD >0,]
    df1$FDR[df1$FDR==0] <- .Machine$double.xmin # cap FDR

    if ("group" %in% colnames(df1) && !short_label) {
        gp <- sprintf("%s~", vapply(str_split(
            df1$group, "\\."), function(x) {
                if(length(x)>3) {x[3]} else {x[1]}}, character(1)))
    } else {
        gp <- ""
    }
    ## TODO: make everything use dbname
    if (label_by %in% colnames(df1)) {
        df1$db1 <- paste0(gp, df1[[label_by]])
    } else if ("feat" %in% colnames(df1)) { # genome-wide data
        df1$db1 <- paste0(gp, df1$feat)
    }
    if (length(unique(df1$db1)) != nrow(df1)) {
        df1 <- df1 %>% group_by(db1) %>% slice_min(
            order_by=FDR, n=1, with_ties=FALSE) %>% ungroup()
    }

    ord <- df1[[order_by]]
    if (order_by == "estimate") { ord <- -ord; }
    df1 <- df1[order(ord, -df1$estimate),]
    df1 <- head(df1, n=n)

    df1$db1 <- factor(df1$db1, levels=rev(df1$db1))
    df1
}

#' Bar plot to show most enriched CG groups from testEnrichment
#'
#' The input data frame should have an "estimate" and
#' a "FDR" columns.
#' 
#' Top CG groups are determined by estimate (descending order).
#'
#' @param df KYCG result data frame
#' @param y the column to be plotted on y-axis
#' @param n number of CG groups to plot
#' @param order_by the column by which CG groups are ordered
#' @param label whether to label significant bars
#' @return grid plot object
#'
#' @import ggplot2
#' @examples
#' KYCG_plotBar(data.frame(
#'   estimate=runif(10,0,10), FDR=runif(10,0,1), nD=10,
#'   overlap=as.integer(runif(10,0,30)), group="g", dbname=seq_len(10)))
#' @export
KYCG_plotBar <- function(df, y = "-log10(FDR)",
    n = 20, order_by = "FDR", label = FALSE) {

    stopifnot("estimate" %in% colnames(df) && "FDR" %in% colnames(df))
    df1 <- preparePlotDF(df, n, order_by)
    if (y == "-log10(FDR)") {
        df1[["-log10(FDR)"]] <- -log10(df1$FDR)
    }

    p <- ggplot(df1, aes_string("db1", y)) +
        geom_bar(stat="identity") +
        coord_flip() + ylab(y) + xlab("CpG Group")
    
    if (label) {
        ## only significant ones are labeled
        df1_label <- df1[df1$FDR < 0.05,]
        df1_label$pos_label <- df1_label[[y]]/2
        df1_label$label <- sprintf("N=%d", df1_label$overlap)
        p <- p + geom_label(aes_string(
            x="db1", y="pos_label", label="label"),
            data = df1_label, alpha=0.6, hjust=0.5)
    }
    p
    ## p2 <- ggplot(df1, aes_string("db1", "estimate")) +
    ##     geom_bar(stat="identity") +
    ##     coord_flip() + ylab("Log2(OR)") + xlab("") +
    ##     theme(axis.text.y = element_blank())
    ## WGG(p1) + WGG(p2, RightOf(width=0.5, pad=0))
}



#' Dot plot to show most enriched CG groups from testEnrichment
#'
#' The input data frame should have an "estimate" and
#' a "FDR" columns.
#' 
#' Top CG groups are determined by estimate (descending order).
#'
#' @param df KYCG result data frame
#' @param y the column to be plotted on y-axis
#' @param n number of CG groups to plot
#' @param order_by the column by which CG groups are ordered
#' @param size_by the column by which CG group size plot
#' @param color_by the column by which CG groups are colored
#' @param label_by the column for label
#' @param short_label omit group in label
#' @param title plot title
#' @return grid plot object (by ggplot)
#'
#' @import ggplot2
#' @examples
#' KYCG_plotDot(data.frame(
#'   estimate=runif(10,0,10), FDR=runif(10,0,1), nD=runif(10,10,20),
#'   overlap=as.integer(runif(10,0,30)), group="g", dbname=seq_len(10)))
#' @export
KYCG_plotDot <- function(df, y = "-log10(FDR)",
    n = 20, order_by = "FDR", title = "Enriched Databases",
    label_by = "dbname", size_by = "overlap", color_by = "estimate",
    short_label = FALSE) {

    df1 <- preparePlotDF(
        df, n, order_by, short_label = short_label, label_by = label_by)

    if (y == "-log10(FDR)") {
        df1[["-log10(FDR)"]] <- -log10(df1$FDR)
    }
    ggplot(df1) +
        geom_point(aes(.data[["db1"]], .data[[y]],
            size = .data[[size_by]], color = .data[[color_by]])) +
        coord_flip() + ggtitle(title) +
        scale_color_gradient(low="blue",high="red") +
        ylab(y) + xlab("")
}

#' creates a volcano plot of -log2(p.value) and log(estimate)
#' given data with fields estimate and p.value.
#'
#' @param df DataFrame where each field is a database name with two fields
#' for the estimate and p.value.
#' @param label_by column in df to be used as the label (default: dbname)
#' @param alpha Float representing the cut-off alpha value for the plot. 
#' Optional. (Default: 0.05)
#' @return ggplot volcano plot
#' @import ggplot2
#' @examples
#' 
#' KYCG_plotVolcano(data.frame(
#'   estimate=runif(10,0,10), FDR=runif(10,0,1), nD=runif(10,10,20),
#'   overlap=as.integer(runif(10,0,30)), group="g", dbname=seq_len(10)))
#'
#' @export
KYCG_plotVolcano <- function(df, label_by="dbname", alpha=0.05) {
    ## suppress R CMD CHECK no visible binding warning
    estimate <- FDR <- label <- NULL

    ## volcano plot cannot plot extreme effect size
    df <- df[abs(df$estimate) < 1000,]
    df[["-log10(FDR)"]] <- -log10(df$FDR)
    df$Significance <- ifelse(
        df$FDR < alpha, "Significant", "Not significant")
    ## TODO: replace with column specifying sig vs non sig
    g <- ggplot(data = df,
        aes_string(x = "estimate", y = "-log10(FDR)", color = "Significance"))
    g <- g + geom_point() + xlab("log2(OR)")
    g <- g + ylab("-log10 FDR") +
        scale_colour_manual(
            name = sprintf("Significance (q < %s)", alpha),
            values = c("Significant" = "red", "Not significant" = "black"))
    requireNamespace("ggrepel")
    g <- g + ggrepel::geom_text_repel(
        data = df[df$FDR < alpha & df$estimate > 0,],
        aes_string(label = label_by), size = 5,
        box.padding = unit(0.35, "lines"),
        point.padding = unit(0.3, "lines"),
        show.legend = FALSE)
    g
}


#' creates a lollipop plot of log(estimate) given data with
#' fields estimate.
#'
#' @param df DataFrame where each row is a database name with its
#' estimate.
#' @param label_column column in df to be used as the label (default: dbname)
#' @param n Integer representing the number of top enrichments to report.
#' Optional. (Default: 10)
#' @return ggplot lollipop plot
#' @import ggplot2
#' @examples
#' 
#' KYCG_plotLollipop(data.frame(
#'   estimate=runif(10,0,10), FDR=runif(10,0,1), nD=runif(10,10,20),
#'   overlap=as.integer(runif(10,0,30)), group="g",
#'   dbname=as.character(seq_len(10))))
#' 
#' @export
KYCG_plotLollipop <- function(df, label_column="dbname", n=20) {
    ## suppress R CMD CHECK no visible binding warning
    estimate <- label <- NULL

    df$label <- df[[label_column]]
    df <- head(df[order(df$estimate, decreasing=TRUE), ], n=n)

    allest <- df$estimate[!is.infinite(df$estimate)]
    cap <- max(allest) * 1.4
    cap_line <- max(allest) * 1.2
    df$estimate[df$estimate == Inf] <- cap
    
    ggplot(df, aes_string(x = "label", y = "estimate", label = "label")) +
        geom_hline(yintercept = 0) +
        geom_segment(aes(
            x=reorder(label, -estimate), y=0,
            yend=estimate, xend=label), color='black') +
        geom_point(fill="black", stat='identity', size=15,
            alpha=0.95, shape=21) +
        scale_fill_gradientn(name='Log2(OR)',
            colours=c('#2166ac','#333333','#b2182b')) +
        geom_text(color='white', size=3) +
        ylab("Log2(OR)") +
        geom_hline(yintercept = cap_line, linetype = "dashed") +
        ylim(min(min(allest)*1.3,0), max(max(allest)*1.5,0)) +
        theme_minimal() +
        theme(axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank())
}

#' create a waterfall plot of log(estimate) given test enrichment
#'
#' @param df data frame where each row is a database with test
#' enrichment result
#' @param order_by the column by which CG groups are ordered
#' @param size_by the column by which CG group size plot
#' @param n_label number of datapoints to label
#' @param label_by column in df to be used as the label (default: dbname)
#' @return grid
#' @import ggplot2
#' @examples
#'
#' library(SummarizedExperiment)
#' df <- rowData(sesameDataGet('MM285.tissueSignature'))
#' query <- df$Probe_ID[df$branch == "fetal_brain" & df$type == "Hypo"]
#' results <- testEnrichment(query, "TFBS", platform="MM285")
#' KYCG_plotWaterfall(results)
#' 
#' @export
KYCG_plotWaterfall <- function(df,
    order_by="Log2(OR)", size_by="-log10(FDR)",
    label_by="dbname", n_label=10) {

    df$label <- df[[label_by]]
    if (size_by == "-log10(FDR)" ||
        order_by == "-log10(FDR)" ||
        label_by == "-log10(FDR)") {
        df[["-log10(FDR)"]] <- -log10(df$FDR)
    }
    if (df$test[[1]] == "Log2(OR)" && (
        size_by == "Log2(OR)" || order_by == "Log2(OR)" ||
        label_by == "Log2(OR)")) {
        df[["Log2(OR)"]] <- df$estimate
        message(sprintf("%d extremes are capped.",
            sum(abs(df[["Log2(OR)"]]) > 1000)))
        ## cap extremes
        df[["Log2(OR)"]][df[["Log2(OR)"]] > 1000] <- 1000
        df[["Log2(OR)"]][df[["Log2(OR)"]] < -1000] <- -1000
        ## df <- df[abs(df$estimate) < 1000,] # skip extremes
    }

    df <- df[order(df[[order_by]]),]
    df$index <- seq_len(nrow(df))
    
    requireNamespace("ggrepel")
    ggplot(df, aes(.data[["index"]], .data[[order_by]])) +
        geom_point(aes(size=.data[[size_by]]), alpha=0.6) +
        geom_hline(yintercept=0, linetype="dashed", color="grey60") +
        theme_minimal() + ylab(order_by) + xlab("Databases") +
        ggrepel::geom_text_repel(
            data = df[head(order(df$log10.p.value),
                n = min(n_label, nrow(df)*0.5)),],
            aes_string(label="label"), nudge_x=-nrow(df)/10,
            max.overlaps=999)
}

#' Plot meta gene or other meta genomic features
#'
#' @param result_list one or a list of testEnrichment
#' @return a grid plot object
#' @examples
#' cg_lists <- KYCG_getDBs("MM285.TFBS")
#' queries <- cg_lists[(sapply(cg_lists, length) > 40000)]
#' result_list <- lapply(queries, testEnrichment,
#'     "MM285.metagene", silent=TRUE, platform="MM285")
#' 
#' KYCG_plotMetaEnrichment(result_list)
#' @export
KYCG_plotMetaEnrichment <- function(result_list) {

    if (is.data.frame(result_list)) { # 1 testEnrichment result
        result_list <- list(result_list)
    }

    stopifnot(all(c("dbname", "label") %in% colnames(result_list[[1]])))
    df <- aggregateTestEnrichments(result_list, return_df = TRUE)
    
    ggplot(df) +
        annotate("rect", xmin = -1, xmax = 10, ymin = -Inf,
            ymax = Inf, fill = "grey80", alpha = .5, color = NA) +
        geom_line(aes_string("db", "estimate", color="query")) +
        scale_x_continuous(breaks=as.integer(result_list[[1]]$db),
            labels=result_list[[1]]$label) +
        annotate("text", x=min(as.integer(result_list[[1]]$db)),
            y=0.05, label="Enrichment", hjust = 0, vjust = 0) +
        annotate("text", x=min(as.integer(result_list[[1]]$db)),
            y=-0.05, label="Depletion", hjust = 0, vjust = 1) +
        geom_hline(yintercept = 0.0, linetype="dashed") +
        ylab("Log2 Fold Enrichment") + xlab("") +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
}

#' Plot meta gene or other meta genomic features
#'
#' @param betas a named numeric vector or a matrix
#' (row: probes; column: samples)
#' @param platform if not given and x is a SigDF, will be inferred
#' the meta features
#' @importFrom reshape2 melt
#' @return a grid plot object
#' @examples
#' sdf <- sesameDataGet("EPIC.1.SigDF")
#' KYCG_plotMeta(getBetas(sdf))
#' @export
KYCG_plotMeta <- function(betas, platform = NULL) {

    if (!is.matrix(betas)) {
        betas <- cbind(sample=betas)
    }
    if (is.null(platform)) {
        platform <- inferPlatformFromProbeIDs(rownames(betas))
    }
    stopifnot(!is.null(platform))

    dbs <- KYCG_getDBs(sprintf("%s.metagene", platform))
    df <- dbStats(betas, dbs, long=TRUE)
    dflabel <- data.frame(
        ord = as.integer(names(dbs)),
        reg = vapply(dbs, function(x) attr(x, "label"), character(1)))
    
    ggplot(df) +
        annotate("rect", xmin = -1, xmax = 10, ymin = -Inf,
            ymax = Inf, fill = "grey80", alpha = .5, color = NA) +
        geom_line(aes_string("db", "value", group="query")) +
        scale_x_continuous(breaks=dflabel$ord, labels=dflabel$reg) +
        ylab("Mean DNA Methylation Level") + xlab("") +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
}


#' Plot point range for a list of enrichment testing results
#' against the same set of databases
#'
#' @param result_list a list of testEnrichment resultsx
#' @return grid plot object
#' @importFrom reshape2 melt
#' @importFrom dplyr summarize
#' @examples
#'
#' ## pick some big TFBS-overlapping CpG groups
#' cg_lists <- KYCG_getDBs("MM285.TFBS")
#' queries <- cg_lists[(sapply(cg_lists, length) > 40000)]

#' result_list <- lapply(queries, testEnrichment,
#'     "MM285.chromHMM", platform="MM285")
#' KYCG_plotPointRange(result_list)
#' 
#' @export
KYCG_plotPointRange <- function(result_list) {

    ord <- mean_betas <- state <- est <- NULL
    
    mtx <- aggregateTestEnrichments(result_list)
    df <- melt(mtx, varnames = c("sample","state"), value.name = "est")
    df <- summarize(group_by(df, state), 
        ave = mean(pmax(-4, est),na.rm=TRUE),
        sd = sd(pmax(-10,est),na.rm=TRUE))
    df$ymin <- df$ave - df$sd
    df$ymax <- df$ave + df$sd
    
    df$state <- factor(df$state, levels = df$state[order(df$ave)])
    
    ggplot(df) +
        geom_pointrange(aes_string("state", "ave", ymin="ymin", ymax="ymax")) +
        geom_hline(yintercept=0, linetype='dashed') +
        ylab("Log2 Fold Enrichment") + xlab("") +
        scale_y_continuous(position="right") +
        annotate("text", x=0.5, y=0.5,
            label="Enrichment", angle=-90, hjust=1) +
        annotate("text", x=0.5, y=-0.5,
            label="Depletion", angle=90, hjust=0) +
        coord_flip()
}

#' KYCG_plotManhattan makes a manhattan plot to summarize EWAS results
#'
#' @param vals named vector of values (P,Q etc), vector name is Probe ID.
#' @param platform String corresponding to the type of platform to use for
#' retrieving GRanges coordinates 
#' of probes. Either MM285, EPIC, HM450, or HM27. If it is not provided, it
#' will be inferred from the query set probeIDs (Default: NA).
#' @param genome hg38, mm10, ..., will infer if not given.
#' For additional mapping, download the GRanges object from
#' http://zwdzwd.github.io/InfiniumAnnotation
#' and provide the following argument
#' ..., genome = sesameAnno_buildManifestGRanges("downloaded_file"),...
#' to this function.
#' @param title title for plot
#' @param label_min Threshold above which data points will be labelled with
#' Probe ID
#' @param col color
#' @param ylabel y-axis label
#' @return a ggplot object
#' @examples
#' 
#' ## see vignette for examples
#' sesameDataGet_resetEnv()
#' 
#' @export
KYCG_plotManhattan <- function(
    vals, platform = NULL, genome = NULL, title = NULL,
    label_min = 100, col = c("wheat1", "sienna3"), ylabel="Value") {

    stopifnot(is(vals, "numeric"))
    if (is.null(platform)) { platform <- queryCheckPlatform(
        platform, query=vals, silent = FALSE) }
    genome <- sesameData_check_genome(genome, platform)
    gr <- sesameData_getManifestGRanges(platform, genome=genome)
    seqLength <- sesameData_getGenomeInfo(genome)$seqLength
    
    v <- vals[names(gr)]
    gr <- gr[!is.na(v)]
    SummarizedExperiment::mcols(gr)$val <- v[!is.na(v)]
    cumLength <- setNames(c(0,cumsum(
        as.numeric(seqLength))[-length(seqLength)]), names(seqLength))
    midLength <- cumLength + seqLength/2
    SummarizedExperiment::mcols(gr)$pos <- cumLength[
        as.character(seqnames(gr))] + end(gr)
    SummarizedExperiment::mcols(gr)$Probe_ID <- names(gr)
    
    df <- as_tibble(gr)
    df$seqnames <- factor(df$seqnames, levels=names(seqLength))
    requireNamespace("ggrepel")
    ggplot(df, aes_string(x="pos", y="val")) + 
        geom_point(aes_string(color="seqnames"), alpha=0.8, size=1.3) + 
        ggrepel::geom_text_repel(data=df[df$val > label_min,],
            aes_string(label="Probe_ID")) +
        scale_color_manual(values = rep(col, length(seqLength))) +
        scale_x_continuous(labels = names(midLength), breaks= midLength) +
        scale_y_continuous(expand = c(0, 0)) +  
        theme_bw() +
        theme( 
            legend.position="none",
            panel.border = element_blank(),
            panel.grid.major.x = element_blank(),
            panel.grid.minor.x = element_blank()
        ) + labs(title=title) + xlab("Chromosome") + ylab(ylabel)
}

#' Plot Set Enrichment
#'
#' @param result result object as returned from an element of the list of
#' testEnrichmentSEA(..., prepPlot=TRUE)
#' @param n_sample number of CpGs to sample
#' @param n_presence number of overlap to sample for the plot
#' @return grid object for plot
#' @examples
#' query <- KYCG_getDBs("KYCG.MM285.designGroup")[["VMR"]]
#' db <- KYCG_getDBs("MM285.seqContextN", "distToTSS")
#' res <- testEnrichmentSEA(query, db, prepPlot = TRUE)
#' KYCG_plotSetEnrichment(res[[1]])
#' 
#' @export
KYCG_plotSetEnrichment <- function(
    result, n_sample = 1000, n_presence = 200) {

    stopifnot("dDisc" %in% names(result))
    dCont <- sort(result$dCont)
    dDisc <- result$dDisc
    presence <- names(dCont) %in% dDisc
    cs <- cumsum(ifelse(presence, 1/sum(presence), -1/sum(!presence)))
    index <- as.integer(seq(1, length(cs), length.out=n_sample))

    pos <- which(presence)
    if(length(pos) > n_presence) { pos <- sample(pos, n_presence) }

    WGG(ggplot(data.frame(index=index, cs=cs[index])) +
        geom_segment(data=data.frame(pos=pos),
            aes_string(x = "pos", xend = "pos", y = -0.02, yend = 0.02),
            color="grey50") +
        geom_line(aes_string(x="index", y="cs"), color="darkred") +
        xlab("") + ylab("ES(S)")) +
    WGG(ggplot(data.frame(index=index, var=dCont[index]),
        aes_string(x="index", y="var")) +
        geom_area() +
        xlab("CpGs") + ylab("Phenotype Var"), Beneath(height=0.5))
}
