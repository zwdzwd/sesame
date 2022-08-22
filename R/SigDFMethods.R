#' SigDF validation from a plain data frame
#'
#' @param df a \code{data.frame} with Probe_ID, MG, MR, UG, UR, col and mask
#' @param platform a string to specify the array platform
#' @param ctl optional control probe data frame
#' @return a \code{SigDF} object
#' @examples
#' sesameDataCache() # if not done yet
#' sdf <- sesameDataGet('EPIC.1.SigDF')
#' @export
SigDF <- function(df, platform = "EPIC", ctl=NULL) {

    df <- df[,c("Probe_ID", "MG","MR","UG","UR","col","mask")]

    ## in case following the manifest
    if (is.factor(df$col) && length(levels(df$col)) == 2) {
        df$col <- as.character(df$col)
        df$col[is.na(df$col)] <- "2"
        df$col <- factor(df$col, levels=c("G","R","2"))
    }
    
    sdf <- structure(df, class=c("SigDF", "data.frame"))
    attr(sdf, "platform") <- platform
    attr(sdf, "controls") <- ctl
    rownames(sdf) <- NULL
    sdf
}

#' Convenience function to output platform attribute of SigDF
#'
#' @param sdf a SigDF object
#' @param verbose print more messages
#' @return the platform string for the SigDF object
#' @examples
#' sesameDataCache()
#' sdf <- sesameDataGet('EPIC.1.SigDF')
#' sdfPlatform(sdf)
#' 
#' @export
sdfPlatform <- function(sdf, verbose = FALSE) {
    if ("platform" %in% attributes(sdf)) {
        attr(sdf, "platform")
    } else {
        inferPlatformFromProbeIDs(sdf$Probe_ID, silent = !verbose)
    }
}

sdfMsg <- function(sdf, verbose, msg, ...) {
    msg <- sprintf(msg, ...)
    msg <- sprintf("[%s] %s", Sys.time(), msg)
    attr(sdf, "msg") <- c(attr(sdf, "msg"), msg)
    if (verbose) {
        message(msg)
    }
    sdf
}

#' remove masked probes from SigDF
#'
#' @param sdf input SigDF object
#' @return a SigDF object without masked probes
#' @export
#' @examples
#' sesameDataCache()
#' sdf <- sesameDataGet("EPIC.1.SigDF")
#' sdf <- pOOBAH(sdf)
#'
#' sdf_noMasked <- noMasked(sdf)
#'
#' @export
noMasked <- function(sdf) { # filter masked probes
    sdf[!sdf$mask,,drop=FALSE]
}

InfIR <- function(sdf) {
    sdf[sdf$col == "R",,drop=FALSE]
}

InfIG <- function(sdf) {
    sdf[sdf$col == "G",,drop=FALSE]
}

InfI <- function(sdf) {
    sdf[sdf$col != "2",,drop=FALSE]
}

InfII <- function(sdf) {
    sdf[sdf$col == "2",,drop=FALSE]
}

oobG <- function(sdf) {
    dR <- InfIR(sdf)
    c(dR$MG, dR$UG)
}

oobR <- function(sdf) {
    dG <- InfIG(sdf)
    c(dG$MR, dG$UR)
}

#' get the controls attributes
#'
#' @param sdf a \code{SigDF}
#' @param verbose print more messages
#' @return the controls data frame
#' @examples
#' sesameDataCache() # if not done yet
#' sdf <- sesameDataGet('EPIC.1.SigDF')
#' head(controls(sdf))
#' @export
controls <- function(sdf, verbose = FALSE) {
    stopifnot(is(sdf, "SigDF"))
    if (sesameDataHas(sprintf(
        "%s.address", sdfPlatform(sdf, verbose = verbose)))) {
        df <- sesameDataGet(sprintf(
            "%s.address", sdfPlatform(sdf, verbose = verbose)))$controls
        if (is.null(df)) {
            return(sdf[grepl("^ctl", sdf$Probe_ID),])
        } else {
            cbind(df, sdf[
                match(paste0("ctl_",df$Address), sdf$Probe_ID),
                c("MG","MR","UG","UR")])
        }
    } else { # custom array
        return(sdf[grepl("^ctl", sdf$Probe_ID),])
    }
}

#' write SigDF to table file
#'
#' @param sdf the \code{SigDF} to output
#' @param ... additional argument to write.table
#' @return write SigDF to table file
#' @examples
#' sesameDataCache() # if not done yet
#' sdf <- sesameDataGet('EPIC.1.SigDF')
#' sdf_write_table(sdf, file=sprintf("%s/sigdf.txt", tempdir()))
#' @export
sdf_write_table <- function(sdf, ...) {
    write.table(sdf, row.names=FALSE, ...)
}

#' read a table file to SigDF
#'
#' @param fname file name
#' @param platform array platform (will infer if not given)
#' @param verbose print more information
#' @param ... additional argument to read.table
#' @return read table file to SigDF
#' @examples
#' sesameDataCache() # if not done yet
#' sdf <- sesameDataGet('EPIC.1.SigDF')
#' fname <- sprintf("%s/sigdf.txt", tempdir())
#' sdf_write_table(sdf, file=fname)
#' sdf2 <- sdf_read_table(fname)
#' @export
sdf_read_table <- function(fname, platform = NULL, verbose=FALSE, ...) {
    df <- read.table(fname, header=TRUE, ...)
    sdf <- structure(df, class=c("SigDF", "data.frame"))
    sdf$col <- factor(sdf$col, levels=c("G","R","2"))
    sdf$mask <- as.logical(sdf$mask)
    if (is.null(platform)) {
        attr(sdf, "platform") <- inferPlatformFromProbeIDs(
            sdf$Probe_ID, silent = !verbose)
    }
    sdf
}
