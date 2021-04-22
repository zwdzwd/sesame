#' @import utils
.onAttach <- function(libname, pkgname) {
    packageStartupMessage('
----------------------------------------------------------
| SEnsible Step-wise Analysis of DNA MEthylation (SeSAMe)
| --------------------------------------------------------
| Please cache the annotation data for your array platform
| (e.g. EPIC) by calling "sesameDataCache("EPIC")"
| or "sesameDataCacheAll()". This needs to be done only
| once per SeSAMe installation.
----------------------------------------------------------
')
}
