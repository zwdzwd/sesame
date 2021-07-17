#' SigDF constructor from a plain data frame
#'
#' @param df a \code{data.frame}
#' @param platform a string to specify the array platform
#' @param ctl optional control probe data frame
#' @return a \code{SigDF} object
#' @examples
#' sesameDataCache("EPIC") # if not done yet
#' df <- as.data.frame(sesameDataGet('EPIC.1.SigDF'))
#' @export
SigDF = function(df, platform = "EPIC", ctl=NULL) {
    sdf = structure(df, class=c("SigDF", "data.frame"))
    sdf$col = factor(sdf$col, levels=c("G","R","2"))
    attr(sdf, "platform") = platform
    attr(sdf, "controls") = ctl
    rownames(sdf) = NULL
    sdf
}

#' Print SigDF object
#'
#' @param x a SigDF object
#' @param ... extra parameter for print
#' @return print SigDF result on screen
#' @examples
#' sesameDataCache("EPIC") # if not done yet
#' sdf = sesameDataGet('EPIC.1.SigDF')
#' sdf
#' @export
print.SigDF = function(x, ...) {
    stopifnot(is(x, "SigDF"))
    cat(paste0(
        sprintf("SigDF - %s\n", platform(x)),
        sprintf(" - %d Infinium-I Probes\n", nrow(InfI(x))),
        sprintf(" - %d Infinium-II Probes\n", nrow(InfII(x))),
        sprintf(" - %d Control Probes\n",
            ifelse(is.null(controls(x)), 0, nrow(controls(x)))),
        sprintf(" - %d Number of Masked Probes\n", sum(x$mask))))
    print(cbind("-"="-",
        Row=c(1,2,nrow(x)-1,nrow(x)),
        as.data.frame(x[c(1,2,nrow(x)-1,nrow(x)),])), row.names=FALSE)
}

#' Convenience function to output platform attribute of SigDF
#'
#' @param sdf a SigDF object
#' @examples
#' sesameDataCache("EPIC")
#' sdf = = sesameDataGet('EPIC.1.SigDF')
#' platform(sdf)
#' @export
platform = function(sdf) {
    stopifnot(is(sdf, "SigDF"))
    attr(sdf, "platform")
}

noMasked = function(sdf) { # filter masked probes
    sdf[!sdf$mask,,drop=FALSE]
}

InfIR = function(sdf) {
    sdf[sdf$col == "R",,drop=FALSE]
}

InfIG = function(sdf) {
    sdf[sdf$col == "G",,drop=FALSE]
}

InfI = function(sdf) {
    sdf[sdf$col != "2",,drop=FALSE]
}

InfII = function(sdf) {
    sdf[sdf$col == "2",,drop=FALSE]
}

oobG = function(sdf) {
    dR = InfIR(sdf)
    c(dR$MG, dR$UG)
}

oobR = function(sdf) {
    dG = InfIG(sdf)
    c(dG$MR, dG$UR)
}

#' get the controls attributes
#'
#' @param sdf a \code{SigDF}
#' @return the controls data frame
#' @examples
#' sesameDataCache("EPIC") # if not done yet
#' sdf = sesameDataGet('EPIC.1.SigDF')
#' head(controls(sdf))
#' @export
controls = function(sdf) {
    stopifnot(is(sdf, "SigDF"))
    attr(sdf, "controls")
}

#' write SigDF to table file
#'
#' @param sdf the \code{SigDF} to output
#' @param ... additional argument to write.table
#' @return write SigDF to table file
#' @examples
#' sesameDataCache("EPIC") # if not done yet
#' sdf = sesameDataGet('EPIC.1.SigDF')
#' sdf_write_table(sdf, file=sprintf("%s/sigdf.txt", tempdir()))
#' @export
sdf_write_table = function(sdf, ...) {
    write.table(sdf, row.names=FALSE, ...)
}

#' read a table file to SigDF
#'
#' @param fname file name
#' @param platform array platform (will infer if not given)
#' @param ... additional argument to read.table
#' @return read table file to SigDF
#' @examples
#' sesameDataCache("EPIC") # if not done yet
#' sdf = sesameDataGet('EPIC.1.SigDF')
#' fname = sprintf("%s/sigdf.txt", tempdir())
#' sdf_write_table(sdf, file=fname)
#' sdf2 = sdf_read_table(fname)
#' @export
sdf_read_table = function(fname, platform = NULL, ...) {
    df = read.table(fname, header=TRUE, ...)
    sdf = structure(df, class=c("SigDF", "data.frame"))
    sdf$col = factor(sdf$col, levels=c("G","R","2"))
    sdf$mask = as.logical(sdf$mask)
    if (is.null(platform)) {
        attr(sdf, "platform") = inferPlatformFromProbeIDs(sdf$Probe_ID)
    }
    sdf
}
