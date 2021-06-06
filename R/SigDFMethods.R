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
    cat("SigDF -", platform(x),
        sprintf("\n - %d Infinium-I Probes\n", nrow(InfI(x))),
        sprintf(" - %d Infinium-II Probes\n", nrow(InfII(x))),
        sprintf(" - %d Control Probes\n", nrow(controls(x))),
        sprintf(" - %d Number of Masked Probes\n", sum(x$mask)))
    print(cbind("- "=" - ",
        Row=c(1,2,nrow(x)-1,nrow(x)),
        as.data.frame(x[c(1,2,nrow(x)-1,nrow(x)),])), row.names=FALSE)
}

platform = function(sdf) {
    stopifnot(is(sdf, "SigDF"))
    attr(sdf, "platform")
}

controls = function(sdf) {
    stopifnot(is(sdf, "SigDF"))
    attr(sdf, "controls")
}

noMask = function(sdf) {
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
