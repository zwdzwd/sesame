
normalizeSetM <- function(input, ref, U) {
    bn <- normalize.quantiles.use.target(matrix(input), ref)
    U * bn / (1-bn)
}

calcMode <- function(x) {
    dd <- density(na.omit(x))
    dd$x[which.max(dd$y)]
}

match1To2_1state <- function(sdf) {
    dR <- InfIR(sdf)
    dR$betas <- getBetas(dR, mask=F)
    dG <- InfIG(sdf)
    dG$betas <- getBetas(dG, mask=F)
    d2 <- InfII(sdf)
    d2$betas <- getBetas(d2, mask=F)

    b2 <- normalize.quantiles.use.target(matrix(dG$betas), d2$betas)
    dG$MG <- dG$UG * b2 / (1 - b2)
    b2 <- normalize.quantiles.use.target(matrix(dR$betas), d2$betas)
    dR$MR <- dR$UR * b2 / (1 - b2)
    sdf2 <- rbind(dR, dG, d2)
    sdf2[order(sdf2$Probe_ID),]
}

#' normalize Infinium I probe betas to Infinium II
#'
#' This is designed to counter tail inflation in Infinium I probes.
#'
#' @param sdf SigDF
#' @param min_dbeta the default algorithm perform 2-state
#' quantile-normalization of the unmethylated and methylated modes
#' separately. However, when the two modes are too close, we fall back
#' to a one-mode normalization. The threshold defines the maximum
#' inter-mode distance.
#' @return SigDF
#' @examples
#'
#' library(RPMM)
#' sdf <- sesameDataGet("MM285.1.SigDF")
#' sesameQC_plotBetaByDesign(sdf)
#' sesameQC_plotBetaByDesign(matchDesign(sdf))
#'
#' @export
matchDesign <- function(sdf, min_dbeta = 0.3) {
    dR <- InfIR(sdf)
    dG <- InfIG(sdf)
    d2 <- InfII(sdf)

    b2 <- getBetas(d2, mask = FALSE)
    m2 <- as.integer(betaMix2States(b2))

    if (sum(m2==1, na.rm=TRUE) < 10 || 
            sum(m2==2, na.rm=TRUE) < 10 || 
            abs(calcMode(b2[m2 == 1]) - calcMode(b2[m2 == 2])) < min_dbeta) {
        return(match1To2_1state(sdf)) }

    bR <- getBetas(dR, mask = FALSE)
    mR <- as.integer(betaMix2States(bR))
    bG <- getBetas(dG, mask = FALSE)
    mG <- as.integer(betaMix2States(bG))

    dR$MR[mR==1] <- normalizeSetM(bR[mR==1], b2[m2==1], dR$UR[mR==1])
    dR$MR[mR==2] <- normalizeSetM(bR[mR==2], b2[m2==2], dR$UR[mR==2])
    dG$MG[mG==1] <- normalizeSetM(bG[mG==1], b2[m2==1], dG$UG[mG==1])
    dG$MG[mG==2] <- normalizeSetM(bG[mG==2], b2[m2==2], dG$UG[mG==2])
    sdf2 <- rbind(dR, dG, d2)
    sdf2[order(sdf2$Probe_ID),]
}

betaMix2States <- function(x, n_samples = 10000, th_init = 0.5) {
    pkgTest("RPMM")
    if (sum(!is.na(x)) > n_samples) {
        x1 <- sample(na.omit(x), n_samples)
    } else {
        x1 <- na.omit(x)
    }
    m <- matrix(0, nrow = length(x1), ncol = 2) # membership matrix
    m[x1 <= th_init, 1] <- 1
    m[x1 > th_init, 2] <- 1
    
    fitres <- blc(matrix(x1), m, maxiter = 5, tol = 0.001, verbose = FALSE)
    m1 <- apply(fitres$w, 1, which.max)
    th = mean(max(x1[m1 == 1]), min(x1[m1 == 2]))
    m2 <- cut(x, breaks=c(0, th, 1), include.lowest = TRUE)
    names(m2) <- names(x)
    m2
}
