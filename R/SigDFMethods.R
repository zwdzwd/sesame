platform = function(sdf) {
    stopifnot(is(sdf, "SigDF"))
    attr(sdf, "platform")
}

controls = function(sdf) {
    stopifnot(is(sdf, "SigDF"))
    attr(sdf, "platform")
}

IR = function(sdf) {
    sdf[sdf$col == "R",]
}

IG = function(sdf) {
    sdf[sdf$col == "G",]
}

II = function(sdf) {
    sdf[sdf$col == "2",]
}
