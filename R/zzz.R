.onAttach <- function(libname, pkgname) {
    packageStartupMessage('
----------------------------------------------------------
| SEnsible Step-wise Analysis of DNA MEthylation (SeSAMe)
| --------------------------------------------------------
| Please cache auxiliary data by "sesameDataCache()".
| This needs to be done only once per SeSAMe installation.
----------------------------------------------------------
')
}
