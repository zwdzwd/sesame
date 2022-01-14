preparePlotDF <- function(df, n_min, n_max, max_fdr) {
    df <- df[df$nD >0,]
    df$FDR[df$FDR==0] <- .Machine$double.xmin
    df <- df[order(df$FDR),]
    if (sum(df$FDR < max_fdr) < n_min) {
        df1 <- head(df, n=n_min)
    } else {
        df <- df[df$estimate > 1,] # enrichment only, exclude depletion
        df1 <- df[df$FDR < max_fdr,]
        df1 <- head(df1, n=n_max)
    }
    
    gp <- vapply(str_split(df1$group, "\\."), function(x) x[3], character(1))
    if ("Target" %in% colnames(df1)) {
        df1$db1 <- sprintf("%s: %s (%s)", gp, df1$Target, df1$db)
    } else {
        df1$db1 <- sprintf("%s: %s", gp, df1$db)
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
#'   overlap=as.integer(runif(10,0,30)), group="g", db=seq_len(10)))
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
#'   overlap=as.integer(runif(10,0,30)), group="g", db=seq_len(10)))
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
#' data=data.frame(estimate=c(runif(10)), FDR=c(runif(10)))
#' KYCG_plotVolcano(data)
#'
#' @export
KYCG_plotVolcano <- function(data, alpha=0.05) {
    ## suppress R CMD CHECK no visible binding warning
    estimate <- FDR <- label <- NULL
    
    if ("Target" %in% colnames(data)) {
        data["label"] <- unlist(data[["Target"]])
    } else {
        data["label"] <- rownames(data)
    }
    
    ## TODO: replace with column specifying sig vs non sig
    if (any(data$FDR <= alpha)) {
        g <- ggplot(data=data, aes(x=log2(estimate), y=-log10(FDR),
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
        data = subset(data, FDR < 0.0005),
        aes(label = label),
        size = 5,
        box.padding = unit(0.35, "lines"),
        point.padding = unit(0.3, "lines"))
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
#' KYCG_plotLollipop(data.frame(estimate=c(runif(10, 0, 10))))
#'
#' @export
KYCG_plotLollipop <- function(df, n=10) {
    ## suppress R CMD CHECK no visible binding warning
    estimate <- label <- NULL
    
    if ("Target" %in% colnames(df))
        df["label"] <- unlist(df[["Target"]])
    else
        df["label"] <- rownames(df)
    
    df <- head(df[order(df$estimate, decreasing=TRUE), ], n=n)
    
    ggplot(df, aes(x=label, 
        y=log2(estimate), 
        label=sprintf('%.2f',log2(estimate)))
    ) + geom_hline(yintercept=0
    ) + geom_segment(aes(y=0, 
        x=reorder(label, -estimate), 
        yend=log2(estimate), xend=label), color='black'
    ) + geom_point(aes(fill=pmax(-1.5,log2(estimate))), 
        stat='identity', 
        size=10, 
        alpha=0.95, 
        shape=21
    ) + scale_fill_gradientn(name='Fold Change',
        colours=c('#2166ac','#333333','#b2182b')
    ) + geom_text(color='white', size=3
    ) + geom_label(aes(x=label,
        y=ifelse(estimate>1,
            log2(estimate) + 0.8,
            log2(estimate) - 0.5),
        label=label),
        alpha=0.8
    ) + ylab("Log2 Enrichment"
    ) + theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
}

