
normalizeSetM <- function(input, ref, U) {
    bn <- normalize.quantiles.use.target(matrix(input), ref)
    U * bn / (1-bn)
}

calcMode <- function(x) {
    dd <- density(na.omit(x))
    dd$x[which.max(dd$y)]
}

valleyDescent <- function(x1, x2) {

    m1 <- calcMode(x1)
    m2 <- calcMode(x2)
    dd <- density(na.omit(c(x1, x2)))
    dfunc <- approxfun(dd$x, dd$y)
    lo <- min(m1, m2)
    hi <- max(m1, m2)
    va <- min(dfunc(c(x1[x1 >= lo & x1 <= hi], x2[x2 >= lo & x2 <= hi])),
        na.rm=TRUE)
    va / min(dfunc(c(lo, hi)), na.rm=TRUE)
}

match1To2_1state <- function(sdf) {
    dR <- noMasked(InfIR(sdf))
    bR <- getBetas(dR)
    dG <- noMasked(InfIG(sdf))
    bG <- getBetas(dG)
    d2 <- noMasked(InfII(sdf))
    b2 <- getBetas(d2)

    dG$MG <- normalizeSetM(bG, b2, dG$UG)
    dR$MR <- normalizeSetM(bR, b2, dR$UR)
    sdf2 <- rbind(dR, dG, d2)
    sdf2 <- rbind(sdf2, sdf[!(sdf$Probe_ID %in% sdf2$Probe_ID),])
    sdf2[order(sdf2$Probe_ID),]
}

match1To2_3states <- function(sdf) {
    dR <- noMasked(InfIR(sdf))
    bR <- getBetas(dR)
    dG <- noMasked(InfIG(sdf))
    bG <- getBetas(dG)
    d2 <- noMasked(InfII(sdf))
    b2 <- getBetas(d2)

    mR <- as.integer(betaMix3States(bR))
    mG <- as.integer(betaMix3States(bG))
    m2 <- as.integer(betaMix3States(b2))
    
    dR$MR[mR==1] <- normalizeSetM(bR[mR==1], b2[m2==1], dR$UR[mR==1])
    dR$MR[mR==2] <- normalizeSetM(bR[mR==2], b2[m2==2], dR$UR[mR==2])
    dR$MR[mR==3] <- normalizeSetM(bR[mR==3], b2[m2==3], dR$UR[mR==3])
    dG$MG[mG==1] <- normalizeSetM(bG[mG==1], b2[m2==1], dG$UG[mG==1])
    dG$MG[mG==2] <- normalizeSetM(bG[mG==2], b2[m2==2], dG$UG[mG==2])
    dG$MG[mG==3] <- normalizeSetM(bG[mG==3], b2[m2==3], dG$UG[mG==3])
    sdf2 <- rbind(dR, dG, d2)
    sdf2 <- rbind(sdf2, sdf[!(sdf$Probe_ID %in% sdf2$Probe_ID),])
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
    dR <- noMasked(InfIR(sdf))
    dG <- noMasked(InfIG(sdf))
    d2 <- noMasked(InfII(sdf))

    b2 <- getBetas(d2)
    m2 <- as.integer(betaMix2States(b2))

    ## message(calcMode(b2[m2 == 1]), " ", calcMode(b2[m2 == 2]))
    ## message(valleyDescent(b2[m2 == 1], b2[m2 == 2]))
    if (sum(m2==1, na.rm=TRUE) > 100 &&
            sum(m2==2, na.rm=TRUE) > 100 &&
            abs(calcMode(b2[m2 == 1]) - calcMode(b2[m2 == 2])) > 0.7) {
        return(match1To2_3states(sdf)) }
    
    if (sum(m2==1, na.rm=TRUE) < 10 || 
            sum(m2==2, na.rm=TRUE) < 10 ||
            valleyDescent(b2[m2==1], b2[m2==2]) >= 0.8 ||
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
    sdf2 <- rbind(sdf2, sdf[!(sdf$Probe_ID %in% sdf2$Probe_ID),])
    sdf2[order(sdf2$Probe_ID),]
}

betaMix2States <- function(x, n_samples = 10000, th_init = 0.5) {
    if (sum(!is.na(x)) > n_samples) {
        x1 <- sample(na.omit(x), n_samples)
    } else {
        x1 <- na.omit(x)
    }
    m <- matrix(0, nrow = length(x1), ncol = 2) # membership matrix
    m[x1 <= th_init, 1] <- 1
    m[x1 > th_init, 2] <- 1
    
    fitres <- RPMM::blc(
        matrix(x1), m, maxiter = 5, tol = 0.001, verbose = FALSE)
    m1 <- apply(fitres$w, 1, which.max)
    th <- mean(max(x1[m1 == 1]), min(x1[m1 == 2]))
    m2 <- cut(x, breaks=c(0, th, 1), include.lowest = TRUE)
    names(m2) <- names(x)
    m2
}

betaMix3States <- function(
    x, n_samples = 10000, th_init1 = 0.2, th_init2 = 0.7) {
    
    if (sum(!is.na(x)) > n_samples) {
        x1 <- sample(na.omit(x), n_samples)
    } else {
        x1 <- na.omit(x)
    }
    m <- matrix(0, nrow = length(x1), ncol = 3) # membership matrix
    m[x1 <= th_init1, 1] <- 1
    m[x1 > th_init1 & x1 <= th_init2, 2] <- 1
    m[x1 > th_init2, 3] <- 1

    fitres <- RPMM::blc(
        matrix(x1), m, maxiter = 5, tol = 0.001, verbose = FALSE)
    m1 <- apply(fitres$w, 1, which.max)
    th1 <- mean(max(x1[m1 == 1]), min(x1[m1 == 2]))
    th2 <- mean(max(x1[m1 == 2]), min(x1[m1 == 3]))
    m2 <- cut(x, breaks=c(0, th1, th2, 1), include.lowest = TRUE)
    names(m2) <- names(x)
    m2
}
