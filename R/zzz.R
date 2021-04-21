#' @import utils
.onAttach <- function(libname, pkgname) {
    packageStartupMessage('
Loading SeSAMe package.
Please cache the annotation data for your array platform (e.g. EPIC) by calling
sesameDataCache("EPIC")
or
sesameDataCacheAll()
This needs to be done only once per SeSAMe installation.
')
}
