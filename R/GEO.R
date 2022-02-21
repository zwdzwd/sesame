
#' Convert signal M and U to SigDF
#'
#' This overcomes the issue of missing IDAT files. However,
#' out-of-band signals will be missing or faked (sampled from a
#' normal distribution).
#' 
#' @param sigM methylated signal, a numeric vector
#' @param sigU unmethylated signal, a numirc vector
#' @param Probe_IDs probe ID vector
#' @param oob.mean assumed mean for out-of-band signals
#' @param oob.sd assumed standard deviation for out-of-band signals
#' @param platform platform code, will infer if not given
#' @return SigDF
#' @examples
#' sigM <- c(11436, 6068, 2864)
#' sigU <- c(1476, 804, 393)
#' probes <- c("cg07881041", "cg23229610", "cg03513874")
#' sdf <- parseGEOsignalMU(sigM, sigU, probes, platform = "EPIC")
#' @export
parseGEOsignalMU <- function(
    sigM, sigU, Probe_IDs, oob.mean = 500, oob.sd = 300, platform = NULL) {

    if (is.null(platform)) {
        platform <- inferPlatformFromProbeIDs(Probe_IDs) }
    addr <- sesameDataGet(paste0(platform, ".address"))$ordering
    M <- sigM[match(addr$Probe_ID, Probe_IDs)]
    U <- sigU[match(addr$Probe_ID, Probe_IDs)]
    col <- ifelse(is.na(addr$col), "2", as.character(addr$col))
    oobs <- pmax(50,rnorm(length(col), mean = oob.mean, sd = oob.sd))
    MG <- ifelse(col == "2", NA, ifelse(col == "G", M, oobs))
    MR <- ifelse(col == "2", NA, ifelse(col == "R", M, oobs))
    UG <- ifelse(col == "2", M, ifelse(col == "G", U, oobs))
    UR <- ifelse(col == "2", U, ifelse(col == "R", U, oobs))
    sdf <- data.frame(Probe_ID = addr$Probe_ID,
        MG = MG, MR = MR, UG = UG, UR = UR,
        col = factor(col, levels=c("G","R","2")), mask = addr$mask)
    class(sdf) <- c("SigDF", class(sdf))
    sdf
}

