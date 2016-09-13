
#' Dye bias correction by matching green to red
#'
#' @param sset a \code{SignalSet}
#' @param in.place to modify sset in place
#' @return a modified \code{SignalSet}
#' @importFrom preprocessCore normalize.quantiles.use.target
#' @importFrom stats approx
#' @export
dyeBiasCorrTypeINorm <- function(sset, in.place=FALSE) {

  library(preprocessCore)
  if (!in.place)
    sset <- sset$clone()

  maxIG <- max(sset$IG)
  minIG <- min(sset$IG)
  maxIR <- max(sset$IR)
  minIR <- min(sset$IR)

  IR1 <- sort(as.numeric(sset$IR))
  IR2 <- sort(as.vector(normalize.quantiles.use.target(matrix(IR1), as.vector(sset$IG))))
  IRmid <- (IR1 + IR2) / 2.0
  maxIRmid <- max(IRmid)
  minIRmid <- min(IRmid)
  fitfunRed <- function(data) {
    insupport    <- data <= maxIR & data >= minIR & (!is.na(data))
    oversupport  <- data > maxIR & (!is.na(data))
    undersupport <- data < minIR & (!is.na(data))
    data[insupport]    <- approx(x=IR1, y=IRmid, xout=data[insupport])$y
    data[oversupport]  <- data[oversupport] - maxIR + maxIRmid
    data[undersupport] <- minIRmid/ minIR * data[undersupport]
    data
  }

  IG1 <- sort(as.numeric(sset$IG))
  IG2 <- sort(as.vector(normalize.quantiles.use.target(matrix(IG1), as.vector(sset$IR))))
  IGmid <- (IG1 + IG2) / 2.0
  maxIGmid <- max(IGmid)
  minIGmid <- min(IGmid)
  fitfunGrn <- function(data) {
    insupport    <- data <= maxIG & data >= minIG & (!is.na(data))
    oversupport  <- data > maxIG & (!is.na(data))
    undersupport <- data < minIG & (!is.na(data))
    data[insupport]    <- approx(x=IG1, y=IGmid, xout=data[insupport])$y
    data[oversupport]  <- data[oversupport] - maxIG + maxIGmid
    data[undersupport] <- minIGmid/ minIG * data[undersupport]
    data
  }

  ## fit type II
  sset$II[,'U'] <- fitfunRed(sset$II[,'U'])
  sset$II[,'M'] <- fitfunGrn(sset$II[,'M'])

  ## IR
  IR <- fitfunRed(sset$IR)
  ## dim(IRmid) <- dim(sset$IR)
  ## dimnames(IRmid) <- dimnames(sset$IR)
  sset$IR <- IR
  ## IG
  IG <- fitfunGrn(sset$IG)
  ## dim(IGmid) <- dim(sset$IG)
  ## dimnames(IGmid) <- dimnames(sset$IG)
  sset$IG <- IG

  ## fit control
  sset$ctl[,'R'] <- fitfunRed(sset$ctl[,'R'])
  sset$ctl[,'G'] <- fitfunGrn(sset$ctl[,'G'])

  ## fit oob
  sset$oobR <- fitfunRed(sset$oobR)
  sset$oobG <- fitfunGrn(sset$oobG)

  sset
}

## M, U, green matched to red distribution
dyeBiasCorrTypeINormG2R <- function(sset, in.place=FALSE) {

  library(preprocessCore)
  if (!in.place)
    sset <- sset$clone()

  maxIG <- max(sset$IG)
  minIG <- min(sset$IG)
  maxIR <- max(sset$IR)
  minIR <- min(sset$IR)

  IG1 <- sort(as.numeric(sset$IG))
  IG2 <- sort(as.vector(normalize.quantiles.use.target(matrix(IG1), as.vector(sset$IR))))
  fitfun <- function(xx) approx(x=IG1, y=IG2, xout=xx)$y

  ## fit type II
  insupport <- sset$II[,'M'] <= maxIG & sset$II[,'M'] >= minIG
  oversupport <- sset$II[,'M'] > maxIG
  undersupport <- sset$II[,'M'] < minIG
  sset$II[insupport,'M'] <- fitfun(sset$II[insupport,'M'])
  sset$II[oversupport,'M'] <- maxIR + (sset$II[oversupport,'M'] - maxIG) * (maxIR-minIR) / (maxIG-minIG)
  sset$II[undersupport,'M'] <- minIR / minIG * sset$II[undersupport,'M']

  ## fit IG
  IG.fit <- fitfun(sset$IG)
  dim(IG.fit) <- dim(sset$IG)
  dimnames(IG.fit) <- dimnames(sset$IG)
  sset$IG <- IG.fit

  ## fit control
  insupport <- sset$ctl[,'G'] <= maxIG & sset$ctl[,'G'] >= minIG & (!is.na(sset$ctl[,'G']))
  oversupport <- sset$ctl[,'G'] > maxIG & (!is.na(sset$ctl[,'G']))
  undersupport <- sset$ctl[,'G'] < minIG & (!is.na(sset$ctl[,'G']))
  sset$ctl[insupport,'G'] <- fitfun(sset$ctl[insupport,'G'])
  sset$ctl[oversupport,'G'] <- maxIR + (sset$ctl[oversupport,'G'] - maxIG) * (maxIR-minIR) / (maxIG-minIG)
  sset$ctl[undersupport,'G'] <- minIR / minIG * sset$ctl[undersupport,'G']

  ## fit oob
  insupport <- sset$oobG <= maxIG & sset$oobG >= minIG & (!is.na(sset$oobG))
  oversupport <- sset$oobG > maxIG & (!is.na(sset$oobG))
  undersupport <- sset$oobG < minIG & (!is.na(sset$oobG))
  sset$oobG[insupport] <- fitfun(sset$oobG[insupport])
  sset$oobG[oversupport] <- maxIR + (sset$oobG[oversupport] - maxIG) * (maxIR-minIR) / (maxIG-minIG)
  sset$oobG[undersupport] <- minIR / minIG * sset$oobG[undersupport]

  sset
}


## red matched to green distribution
dyeBiasCorrTypeINormR2G <- function(sset, in.place=FALSE) {

  library(preprocessCore)
  if (!in.place)
    sset <- sset$clone()

  maxIG <- max(sset$IG)
  minIG <- min(sset$IG)
  maxIR <- max(sset$IR)
  minIR <- min(sset$IR)

  IR1 <- sort(as.numeric(sset$IR))
  IR2 <- sort(as.vector(normalize.quantiles.use.target(matrix(IR1), as.vector(sset$IG))))
  fitfun <- function(xx) approx(x=IR1, y=IR2, xout=xx)$y

  ## fit type II
  insupport <- sset$II[,'U'] <= maxIR & sset$II[,'U'] >= minIR
  oversupport <- sset$II[,'U'] > maxIR
  undersupport <- sset$II[,'U'] < minIR
  sset$II[insupport,'U'] <- fitfun(sset$II[insupport,'U'])
  sset$II[oversupport,'U'] <- maxIG + (sset$II[oversupport,'U'] - maxIR) * (maxIG-minIG) / (maxIR-minIR)
  sset$II[undersupport,'U'] <- minIG / minIR * sset$II[undersupport,'U']

  ## fit IR
  IR.fit <- fitfun(sset$IR)
  dim(IR.fit) <- dim(sset$IR)
  dimnames(IR.fit) <- dimnames(sset$IR)
  sset$IR <- IR.fit

  ## fit control
  insupport <- sset$ctl[,'R'] <= maxIR & sset$ctl[,'R'] >= minIR & (!is.na(sset$ctl[,'R']))
  oversupport <- sset$ctl[,'R'] > maxIR & (!is.na(sset$ctl[,'R']))
  undersupport <- sset$ctl[,'R'] < minIR & (!is.na(sset$ctl[,'R']))
  sset$ctl[insupport,'R'] <- fitfun(sset$ctl[insupport,'R'])
  sset$ctl[oversupport,'R'] <- maxIG + (sset$ctl[oversupport,'R'] - maxIR) * (maxIG-minIG) / (maxIR-minIR)
  sset$ctl[undersupport,'R'] <- minIG / minIR * sset$ctl[undersupport,'R']

  ## fit oob
  insupport <- sset$oobR <= maxIR & sset$oobR >= minIR & (!is.na(sset$oobR))
  oversupport <- sset$oobR > maxIR & (!is.na(sset$oobR))
  undersupport <- sset$oobR < minIR & (!is.na(sset$oobR))
  sset$oobR[insupport] <- fitfun(sset$oobR[insupport])
  sset$oobR[oversupport] <- maxIG + (sset$oobR[oversupport] - maxIR) * (maxIG-minIG) / (maxIR-minIR)
  sset$oobR[undersupport] <- minIG / minIR * sset$oobR[undersupport]

  sset
}

## M + U green matched to red distribution
dyeBiasCorrTypeINormMpU <- function(sset, in.place=FALSE) {
  if (!in.place)
    sset <- sset$clone()
  maxIG <- max(rowSums(sset$IG))
  minIG <- min(rowSums(sset$IG))
  maxIR <- max(rowSums(sset$IR))
  minIR <- min(rowSums(sset$IR))

  IG1 <- sort(as.numeric(rowSums(sset$IG)))
  IG2 <- sort(as.vector(normalize.quantiles.use.target(matrix(IG1), as.vector(rowSums(sset$IR)))))
  fitfun <- function(xx) approx(x=IG1, y=IG2, xout=xx)$y

  ## fit type II
  insupport <- sset$II[,'M'] <= maxIG & sset$II[,'M'] >= minIG
  oversupport <- sset$II[,'M'] > maxIG
  undersupport <- sset$II[,'M'] < minIG
  sset$II[insupport,'M'] <- fitfun(sset$II[insupport,'M'])
  sset$II[oversupport,'M'] <- maxIR + (sset$II[oversupport,'M'] - maxIG) * (maxIR-minIR) / (maxIG-minIG)
  sset$II[undersupport,'M'] <- minIR / minIG * sset$II[undersupport,'M']

  ## fit IG
  IG.fit <- fitfun(sset$IG)
  dim(IG.fit) <- dim(sset$IG)
  dimnames(IG.fit) <- dimnames(sset$IG)
  sset$IG <- IG.fit

  ## fit control
  insupport <- sset$ctl[,'G'] <= maxIG & sset$ctl[,'G'] >= minIG & (!is.na(sset$ctl[,'G']))
  oversupport <- sset$ctl[,'G'] > maxIG & (!is.na(sset$ctl[,'G']))
  undersupport <- sset$ctl[,'G'] < minIG & (!is.na(sset$ctl[,'G']))
  sset$ctl[insupport,'G'] <- fitfun(sset$ctl[insupport,'G'])
  sset$ctl[oversupport,'G'] <- maxIR + (sset$ctl[oversupport,'G'] - maxIG) * (maxIR-minIR) / (maxIG-minIG)
  sset$ctl[undersupport,'G'] <- minIR / minIG * sset$ctl[undersupport,'G']

  ## fit oob
  insupport <- sset$oobG <= maxIG & sset$oobG >= minIG & (!is.na(sset$oobG))
  oversupport <- sset$oobG > maxIG & (!is.na(sset$oobG))
  undersupport <- sset$oobG < minIG & (!is.na(sset$oobG))
  sset$oobG[insupport] <- fitfun(sset$oobG[insupport])
  sset$oobG[oversupport] <- maxIR + (sset$oobG[oversupport] - maxIG) * (maxIR-minIR) / (maxIG-minIG)
  sset$oobG[undersupport] <- minIR / minIG * sset$oobG[undersupport]

  sset
}
