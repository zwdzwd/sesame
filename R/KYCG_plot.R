preparePlotDF <- function(df, n_min, n_max, max_fdr) {
    df <- df[df$nD >0,]
    df$FDR[df$FDR==0] <- .Machine$double.xmin
    df <- df[order(df$FDR),]
    if (sum(df$FDR < max_fdr) < n_min) {
        df1 <- head(df, n=n_min)
    } else {
        df <- df[df$estimate > 0,] # enrichment only, exclude depletion
        df1 <- df[df$FDR < max_fdr,]
        df1 <- head(df1, n=n_max)
    }
    
    gp <- vapply(str_split(df1$group, "\\."), function(x) x[3], character(1))
    if ("Target" %in% colnames(df1)) {
        df1$db1 <- sprintf("%s: %s (%s)", gp, df1$Target, df1$dbname)
    } else {
        df1$db1 <- sprintf("%s: %s", gp, df1$dbname)
    }
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
#' @param n_min minimum number of databases to report
#' @param n_max maximum number of databases to report
#' @param max_fdr maximum FDR
#' @return grid plot object
#'
#' @import ggplot2
#' @examples
#' KYCG_plotBar(data.frame(
#'   estimate=runif(10,0,10), FDR=runif(10,0,1), nD=10,
#'   overlap=as.integer(runif(10,0,30)), group="g", dbname=seq_len(10)))
#' @export
KYCG_plotBar <- function(df, n_min = 10, n_max = 30, max_fdr = 0.05) {

    db1 <- FDR <- overlap <- estimate <- NULL
    stopifnot("estimate" %in% colnames(df) && "FDR" %in% colnames(df))

    df1 <- preparePlotDF(df, n_min, n_max, max_fdr)
    p1 <- ggplot(df1, aes(db1, -log10(FDR))) + geom_bar(stat="identity") +
        coord_flip() + ylab("-log10(P-value)") + xlab("CpG Group") +
        geom_label(aes(x=db1, y=-log10(FDR)/2,
            label=sprintf("%d CGs", overlap)),
            data = df1[df1$FDR < 0.05,], alpha=0.6, hjust=0.5)
    p2 <- ggplot(df1, aes(db1, estimate)) + geom_bar(stat="identity") +
        coord_flip() + ylab("Enrichment Score") + xlab("") +
        theme(axis.text.y = element_blank())
    WGG(p1) + WGG(p2, RightOf(width=0.5, pad=0))
}



#' Dot plot to show most enriched CG groups from testEnrichment
#'
#' The input data frame should have an "estimate" and
#' a "FDR" columns.
#' 
#' Top CG groups are determined by estimate (descending order).
#'
#' @param df KYCG result data frame
#' @param n_min minimum number of databases to report
#' @param n_max maximum number of databases to report
#' @param max_fdr maximum FDR
#' @return grid plot object
#'
#' @import ggplot2
#' @examples
#' KYCG_plotDot(data.frame(
#'   estimate=runif(10,0,10), FDR=runif(10,0,1), nD=runif(10,10,20),
#'   overlap=as.integer(runif(10,0,30)), group="g", dbname=seq_len(10)))
#' @export
KYCG_plotDot <- function(df, n_min = 10, n_max = 30, max_fdr = 0.05) {

    db1 <- FDR <- overlap <- estimate <- NULL
    stopifnot("estimate" %in% colnames(df) && "FDR" %in% colnames(df))

    df1 <- preparePlotDF(df, n_min, n_max, max_fdr)
    ggplot(df1) +
        geom_point(aes(db1, -log10(FDR), size=overlap, color=estimate)) +
        coord_flip() + ggtitle("Enriched Databases") +
        scale_color_gradient(low="blue",high="red") +
        ylab("-log10(FDR)") + xlab("")
}

#' creates a volcano plot of -log2(p.value) and log(estimate)
#' given data with fields estimate and p.value.
#'
#' @param data DataFrame where each field is a database name with two fields
#' for the estimate and p.value.
#' @param alpha Float representing the cut-off alpha value for the plot. 
#' Optional. (Default: 0.05)
#' @return ggplot volcano plot
#' @import ggplot2
#' @import ggrepel
#' @examples
#' 
#' KYCG_plotVolcano(data.frame(
#'   estimate=runif(10,0,10), FDR=runif(10,0,1), nD=runif(10,10,20),
#'   overlap=as.integer(runif(10,0,30)), group="g", dbname=seq_len(10)))
#'
#' @export
KYCG_plotVolcano <- function(data, alpha=0.05) {
    ## suppress R CMD CHECK no visible binding warning
    estimate <- FDR <- label <- NULL
    
    if ("Target" %in% colnames(data)) {
        data$label <- data$Target
    } else {
        data$label <- data$dbname
    }

    data <- data[data$estimate > -Inf,]
    ## TODO: replace with column specifying sig vs non sig
    if (any(data$FDR <= alpha)) {
        g <- ggplot(data=data, aes(x=estimate, y=-log10(FDR),
            color = cut(FDR, c(-Inf, alpha))))
    } else {
        g <- ggplot(data=data, aes(x=estimate, y=FDR))
    }
    g <- g + geom_point() + xlab("log2 Fold Change")
    
    g <- g + 
        ylab("-log10 q-value") +
        scale_colour_discrete(
            name = sprintf("Significance (q < %s)", alpha),
            labels=c("Significant", "Not Significant"))

    options(ggrepel.max.overlaps = 10)
    g + geom_text_repel(
        data = subset(data, FDR < 0.0005 & estimate > 0),
        aes(label = label), size = 5,
        box.padding = unit(0.35, "lines"),
        point.padding = unit(0.3, "lines"),
        show.legend = FALSE)
}


#' creates a lollipop plot of log(estimate) given data with
#' fields estimate.
#'
#' @param df DataFrame where each field is a database name with its
#' estimate.
#' @param n Integer representing the number of top enrichments to report.
#' Optional. (Default: 10)
#' @return ggplot lollipop plot
#'
#' @import ggplot2
#'
#' @examples
#' 
#' KYCG_plotLollipop(data.frame(
#'   estimate=runif(10,0,10), FDR=runif(10,0,1), nD=runif(10,10,20),
#'   overlap=as.integer(runif(10,0,30)), group="g",
#'   dbname=as.character(seq_len(10))))
#' 
#' @export
KYCG_plotLollipop <- function(df, n=10) {
    ## suppress R CMD CHECK no visible binding warning
    estimate <- label <- NULL
    
    if ("Target" %in% colnames(df))
        df$label <- df$Target
    else
        df$label <- df$dbname
    
    df <- head(df[order(df$estimate, decreasing=TRUE), ], n=n)
    
    ggplot(df, aes(x = label, 
        y = estimate, label = sprintf('%.2f', estimate))
    ) + geom_hline(yintercept = 0
    ) + geom_segment(aes(y = 0, 
        x = reorder(label, -estimate), 
        yend = estimate, xend=label), color='black'
    ) + geom_point(aes(fill=pmax(-2, 2)),
        stat='identity', size=10, alpha=0.95, shape=21
    ) + scale_fill_gradientn(name='Fold Change',
        colours=c('#2166ac','#333333','#b2182b')
    ) + geom_text(color='white', size=3
    ) + geom_label(aes(x = label,
        y = ifelse(estimate > 0, estimate + 0.8, estimate - 0.5),
        label = label), alpha=0.8
    ) + ylab("Log2 Enrichment"
    ) + theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
}

#' Plot meta gene or other meta genomic features
#'
#' @param result_list one or a list of testEnrichment
#' @return a grid plot object
#' @examples
#' cg_lists <- KYCG_getDBs("MM285.TFBS")
#' queries <- cg_lists[(sapply(cg_lists, length) > 40000)]
#' result_list <- lapply(queries, testEnrichment,
#'     "MM285.metagene", silent=TRUE)
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

    meta <- KYCG_getDBs(sprintf("%s.metagene", platform))
    mb <- do.call(cbind, lapply(
        meta, function(m1) apply(betas[m1,,drop=FALSE],2, mean, na.rm=TRUE)))
    df <- melt(mb, varnames=c("query","db"), value.name="mean_betas")
    dflabel <- data.frame(
        ord = as.integer(names(meta)),
        reg = vapply(meta, function(x) attr(x, "label"), character(1)))
    
    ggplot(df) +
        annotate("rect", xmin = -1, xmax = 10, ymin = -Inf,
            ymax = Inf, fill = "grey80", alpha = .5, color = NA) +
        geom_line(aes_string("db", "mean_betas")) +
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

#' result_list <- lapply(queries, testEnrichment, "MM285.chromHMM")
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
    
    df$state <- factor(df$state, levels = df$state[order(df$ave)])
    
    ggplot(df) +
        geom_pointrange(aes(state, ave, ymin=ave-sd, ymax=ave+sd)) +
        geom_hline(yintercept=0, linetype='dashed') +
        ylab("Log2 Fold Enrichment") + xlab("") +
        scale_y_continuous(position="right") +
        annotate("text", x=0.5, y=0.5,
            label="Enrichment", angle=-90, hjust=1) +
        annotate("text", x=0.5, y=-0.5,
            label="Depletion", angle=90, hjust=0) +
        coord_flip()
}

#' plot enrichment test result
#'
#' @param df test enrichment result data frame
#' @return grid object
#' @importFrom stringr str_replace
#' @importFrom tibble rownames_to_column
#' @import ggplot2
#' 
#' @export
KYCG_plotEnrichAll <- function(df, fdr_max = 25) {

    gp_size <- sort(table(df$group))
    gp_width <- log(gp_size)
    e1 <- df[order(factor(df$group, levels=names(gp_size)), df$dbname),]
    e1$inc <- (gp_width / gp_size)[e1$group]
    e1$inc1 <- c(0,ifelse(e1$group[-1] != e1$group[-nrow(e1)], 1, 0))
    e1$inc2 <- cumsum(e1$inc + e1$inc1)

    e1$group <- str_replace(e1$group,"KYCG.","")
    e1$group <- vapply(strsplit(e1$group, "\\."),
        function(x) paste0(x[2:(length(x)-1)], collapse="."), character(1))

    e2 <- e1[e1$estimate > 1 & e1$FDR < 0.01 ,]
    e2$FDR[e2$FDR < 10**-fdr_max] <- 10**-(fdr_max*1.1)

    e3 <- rownames_to_column(as.data.frame(do.call(rbind, lapply(
        split(e1$inc2, e1$group), function(x)
            c(beg=min(x), mid=mean(x), end=max(x))))), "group")

    ggplot(e2, aes(inc2, -log10(FDR))) +
        geom_point(aes(size=estimate, color=group), alpha=0.5) +
        geom_text_repel(data = e2[order(e2$FDR)[1:15],],
            aes(label=dbname, color=group), size = 3,
            ## box.padding = unit(0.35, "lines"),
            ## point.padding = unit(0.3, "lines"),
            direction="y", nudge_y=0.2, max.overlaps=100) +
        annotate("text", 0, fdr_max*0.96,
            label="Values above this line are capped.",
            hjust=0, vjust=1, color="grey60") +
        geom_hline(yintercept = fdr_max, linetype="dotted", color="grey60") +
        geom_segment(aes(x = beg, y = 0, xend = end, yend = 0, color=group),
            size=3, data=e3) +
        geom_text(data=e3,aes(mid, -1, label=group, color=group),
            vjust=1, hjust=1, angle=30) + scale_color_discrete(guide="none") +
        ylim(-6, max(-log10(e2$FDR))*1.1) + xlab("") +
        scale_size_continuous(guide=guide_legend(title="log2(FC)")) +
        coord_cartesian(clip="off") + theme_minimal() +
        theme(axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            panel.grid.minor.x = element_blank())
}
