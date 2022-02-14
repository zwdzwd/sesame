
#' Convert signal M and U to SigDF
#'
#' This overcomes the issue of missing IDAT files. However,
#' out-of-band signals will be missing or faked.
#' 
#' @param sigM methylated signal, a numeric vector
#' @param sigU unmethylated signal, a numirc vector
#' @param Probe_IDs probe ID vector
#' @param oob assumed value for out-of-band signals
#' @return SigDF
#' @examples
#' sigM <- c(11436, 6068, 2864)
#' sigU <- c(1476, 804, 393)
#' probes <- c("cg07881041", "cg23229610", "cg03513874")
#' sdf <- parseGEOsignalMU(sigM, sigU, probes)
#' @export
parseGEOsignalMU <- function(sigM, sigU, Probe_IDs, oob=500) {
    platform <- inferPlatformFromProbeIDs(Probe_IDs)
    addr <- sesameDataGet(paste0(platform, ".address"))$ordering
    M <- sigM[match(addr$Probe_ID, Probe_IDs)]
    U <- sigU[match(addr$Probe_ID, Probe_IDs)]
    col <- ifelse(is.na(addr$col), "2", as.character(addr$col))
    MG <- ifelse(col == "2", NA, ifelse(col == "G", M, oob))
    MR <- ifelse(col == "2", NA, ifelse(col == "R", M, oob))
    UG <- ifelse(col == "2", M, ifelse(col == "G", U, oob))
    UR <- ifelse(col == "2", U, ifelse(col == "R", U, oob))
    sdf <- data.frame(Probe_ID = addr$Probe_ID,
        MG = MG, MR = MR, UG = UG, UR = UR,
        col = factor(col, levels=c("G","R","2")), mask = addr$mask)
    class(sdf) <- c("SigDF", class(sdf))
    sdf
}

