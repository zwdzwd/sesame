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
