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

IR = function(sdf) {
    sdf[sdf$col == "R",,drop=FALSE]
}

IG = function(sdf) {
    sdf[sdf$col == "G",,drop=FALSE]
}

II = function(sdf) {
    sdf[sdf$col == "2",,drop=FALSE]
}

oobG = function(sdf) {
    with(IR(sdf), c(MG, UG))
}

oobR = function(sdf) {
    with(IG(sdf), c(MR, UR))
}
