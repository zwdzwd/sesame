
#' build control matrix for funnorm
#' 
#' build control matrix for funnorm
#'
#' @param dmp an object of class SignalSet
#' @return a vector with control summaries
#' @export
BuildControlMatrix450k <- function(dmp) {
  
  ctls <- split(dmp$ctl, dmp$ctl$type)

  cm <- NULL

  ## bisulfite conversion type II
  cm <- c(cm, bisulfite2=mean(ctls[['BISULFITE CONVERSION II']]$R,na.rm=T))

  ## bisulfite conversion type I
  cm <- c(cm, bisulfite1=mean(
                ctls[['BISULFITE CONVERSION I']][sprintf('BS.Conversion.I.C%s', 1:3),'G'] +
                  ctls[['BISULFITE CONVERSION I']][sprintf('BS.Conversion.I.C%s', 4:6),'R']))

  ## staining
  cm <- c(cm, stain.red=ctls[['STAINING']]['DNP..High.', 'R'],
          stain.green=ctls[['STAINING']]['Biotin..High.','G'])

  ## extension
  cm <- c(cm,
          extRed1=ctls[['EXTENSION']]['Extension..A.','R'],
          extRed2=ctls[['EXTENSION']]['Extension..T.','R'],
          extGrn1=ctls[['EXTENSION']]['Extension..C.','G'],
          extGrn2=ctls[['EXTENSION']]['Extension..G.','G'])

  ## hybridization
  d <- ctls[['HYBRIDIZATION']]$G
  cm <- c(cm, setNames(d, paste0('hybe',1:length(d))))

  ## target removal
  d <- ctls[['TARGET REMOVAL']]$G
  cm <- c(cm, setNames(d, paste0('targetrem',1:length(d))))

  ## non-polymorphic
  cm <- c(cm,
          nonpolyRed1=ctls[['NON-POLYMORPHIC']]['NP..A.','R'],
          nonpolyRed2=ctls[['NON-POLYMORPHIC']]['NP..T.','R'],
          nonpolyGrn1=ctls[['NON-POLYMORPHIC']]['NP..C.','G'],
          nonpolyGrn2=ctls[['NON-POLYMORPHIC']]['NP..G.','G'])

  ## specificity type II
  d <- ctls[['SPECIFICITY II']]
  cm <- c(cm,
          structure(d$G, names=paste0('spec2Grn', 1:dim(d)[1])),
          structure(d$R, names=paste0('spec2Red', 1:dim(d)[1])))
  cm <- c(cm, spec2.ratio = mean(d$G,na.rm=TRUE) / mean(d$R,na.rm=TRUE))

  ## specificity type I green
  d <- ctls[['SPECIFICITY I']][sprintf('GT.Mismatch.%s..PM.',1:3),]
  cm <- c(cm, structure(d$G, names=paste0('spec1Grn',1:dim(d)[1])))
  cm <- c(cm, spec1.ratio1 = mean(d$R, na.rm=TRUE)/mean(d$G, na.rm=TRUE))

  ## specificity type I red
  d <- ctls[['SPECIFICITY I']][sprintf('GT.Mismatch.%s..PM.',4:6),]
  cm <- c(cm, structure(d$R, names=paste0('spec1Red',1:dim(d)[1])))
  cm <- c(cm, spec1.ratio2 = mean(d$G, na.rm=TRUE)/mean(d$R, na.rm=TRUE))

  ## average specificity ratio
  cm <- c(cm, spec1.ratio = unname((cm['spec1.ratio1']+cm['spec1.ratio2'])/2.0))

  ## normalization
  cm <- c(cm, c(normA=mean(ctls[['NORM_A']]$R, na.rm=TRUE), 
                normT=mean(ctls[['NORM_T']]$R, na.rm=TRUE), 
                normC=mean(ctls[['NORM_C']]$G, na.rm=TRUE), 
                normG=mean(ctls[['NORM_G']]$G, na.rm=TRUE)))

  cm <- c(cm, dyebias=(cm['normC']+cm['normG']) / (cm['normA']+cm['normT']))

  ## out-of-band probe quantiles
  cm <- c(cm, oob.ratio=median(dmp$oobG) / median(dmp$oobR))
  cm <- c(cm, structure(quantile(dmp$oobG, c(0.01,0.5,0.99)), names=paste0('oob', c(1,50,99))))

  cm

}

#' build quantiles for funnorm
#' 
#' build quantiles for funnorm
#'
#' quantiles are stratified by IR autosomal, IG autosomal, II autosomal, X all design and Y all design and by methylated and unmethylated channel.
#' @param dmp a SignalSet
#' @param n number of grid points in quantiles.
#' @return a list of quantiles in each category
#' @export
BuildQuantiles450k <- function(dmp, n=500) {

  data(hm450.hg19.probe2chr)
  
  probs <- seq(from=0, to=1, length.out=n)

  auto.IR <- dmp$IR[!(hm450.hg19.probe2chr[rownames(dmp$IR)] %in% c('chrX','chrY')),]
  quantile.auto.IR.M <- quantile(auto.IR[,'M'], probs)
  quantile.auto.IR.U <- quantile(auto.IR[,'U'], probs)

  auto.IG <- dmp$IG[!(hm450.hg19.probe2chr[rownames(dmp$IG)] %in% c('chrX','chrY')),]
  quantile.auto.IG.M <- quantile(auto.IG[,'M'], probs)
  quantile.auto.IG.U <- quantile(auto.IG[,'U'], probs)

  auto.II <- dmp$II[!(hm450.hg19.probe2chr[rownames(dmp$II)] %in% c('chrX','chrY')),]
  quantile.auto.II.M <- quantile(auto.II[,'M'], probs)
  quantile.auto.II.U <- quantile(auto.II[,'U'], probs)

  all <- rbind(dmp$IR, dmp$IG, dmp$II)
  X.all <- all[hm450.hg19.probe2chr[rownames(all)] == 'chrX',]
  quantile.X.all.M <- quantile(X.all[,'M'], probs)
  quantile.X.all.U <- quantile(X.all[,'U'], probs)
  Y.all <- all[hm450.hg19.probe2chr[rownames(all)] == 'chrY',]
  quantile.Y.all.M <- quantile(Y.all[,'M'], probs)
  quantile.Y.all.U <- quantile(Y.all[,'U'], probs)

  list(auto.IR.M=quantile.auto.IR.M,
       auto.IR.U=quantile.auto.IR.U,
       auto.IG.M=quantile.auto.IG.M,
       auto.IG.U=quantile.auto.IG.U,
       auto.II.M=quantile.auto.II.M,
       auto.II.U=quantile.auto.II.U,
       X.all.M=quantile.X.all.M,
       X.all.U=quantile.X.all.U,
       Y.all.M=quantile.Y.all.M,
       Y.all.U=quantile.Y.all.U)

}

#' funnorm regress
#'
#' regress out control effect from quantile distribution
#'
#' @param cms control matrices from all samples
#' @param qntiles quantile distributions from all samples
#' @param k number of principal components from control matrices for regression
#' @import preprocessCore
#' @return new quantile distriubtions from all samples
#' @export
FunnormRegress <- function(cms, qntiles, genders=NULL, k=2) {

  mm <- do.call(rbind, cms)
  
  ## Impute missing value using the mean across samples
  mm <- apply(mm, 2, function(col) {
    col[is.na(col)] <- mean(col, na.rm=TRUE)
    col
  })
  
  ## Scale and rescale
  mm <- scale(mm)
  mm[mm > 3] <- 3
  mm[mm < (-3)] <- -3
  mm <- scale(mm)

  ## Normalize quantiles - autosomes
  mprint("Normalizing autosomes.")
  FunnormRegress1 <- function(qmat, mm, k=2) {
    stopifnot(identical(colnames(qmat), rownames(mm)))
    qmat[1,] <- 0
    qmat[nrow(qmat),] <- qmat[nrow(qmat)-1,] + 1000
    global.mean.quantiles <- rowMeans(qmat)
    res <- qmat - global.mean.quantiles
    mm.pcs <- prcomp(mm)$x[,1:k,drop=FALSE]
    design <- model.matrix(~mm.pcs)
    fits <- lm.fit(x = design, y = t(res))
    qntiles.n <- global.mean.quantiles + t(fits$residuals)
  }

  qntiles.IR.M <- FunnormRegress1(sapply(qntiles, function(x) x$auto.IR.M), mm, k=k)
  qntiles.IG.M <- FunnormRegress1(sapply(qntiles, function(x) x$auto.IG.M), mm, k=k)
  qntiles.II.M <- FunnormRegress1(sapply(qntiles, function(x) x$auto.II.M), mm, k=k)
  qntiles.IR.U <- FunnormRegress1(sapply(qntiles, function(x) x$auto.IR.U), mm, k=k)
  qntiles.IG.U <- FunnormRegress1(sapply(qntiles, function(x) x$auto.IG.U), mm, k=k)
  qntiles.II.U <- FunnormRegress1(sapply(qntiles, function(x) x$auto.II.U), mm, k=k)

  ## X chromosome
  mprint("Normalizing X-chromosome.")
  if ((!is.null(genders)) &&
      sum(genders == 0)>10 && sum(genders == 1)>10) {
    qntiles.X.all.M <- cbind(
      FunnormRegress1(sapply(
        qntiles[genders[names(qntiles)]==0],
        function(x) x$X.all.M), mm, k=k),
      FunnormRegress1(sapply(
        qntiles[genders[names(qntiles)]==1],
        function(x) x$X.all.M), mm, k=k))
    qntiles.X.all.U <- cbind(
      FunnormRegress1(sapply(
        qntiles[genders[names(qntiles)]==0],
        function(x) x$X.all.U), mm, k=k),
      FunnormRegress1(sapply(
        qntiles[genders[names(qntiles)]==1],
        function(x) x$X.all.U), mm, k=k))
  } else {
    qntiles.X.all.M <- FunnormRegress1(sapply(
      qntiles, function(x) x$X.all.M), mm, k=k)
    qntiles.X.all.U <- FunnormRegress1(sapply(
      qntiles, function(x) x$X.all.U), mm, k=k)
  }

  ## Y chromosome using conventional quantile normalization
  mprint("Normalizing Y-chromosome.")
  if ((!is.null(genders)) &&
      any(genders == 0) && any(genders == 1)) {
    norm.Y.all.M <- cbind(
      preprocessCore::normalize.quantiles(
        sapply(qntiles[genders[names(qntiles)]==0],
               function(x) x$Y.all.M)),
      preprocessCore::normalize.quantiles(
        sapply(qntiles[genders[names(qntiles)]==1],
               function(x) x$Y.all.M)))
    norm.Y.all.U <- cbind(
      preprocessCore::normalize.quantiles(
        sapply(qntiles[genders[names(qntiles)]==0],
               function(x) x$Y.all.U)),
      preprocessCore::normalize.quantiles(
        sapply(qntiles[genders[names(qntiles)]==1],
               function(x) x$Y.all.U)))
  } else {
    norm.Y.all.M <- preprocessCore::normalize.quantiles(
      sapply(qntiles, function(x) x$Y.all.M))
    norm.Y.all.U <- preprocessCore::normalize.quantiles(
      sapply(qntiles, function(x) x$Y.all.U))
  }
}

#' Infer gender from signals
#'
#' Infer gender from signals
#' @param dmps a list of \code{SignalSet}s
#' @return named vector of 0 (for female) and 1 (for male)
#' @export
InferGenders <- function(dmps) {
  XYmedian <- t(sapply(dmps, function(dmp) {
    all.signals <- c(apply(dmp$IR,1,max), apply(dmp$IG,1,max), apply(dmp$II,1,max))
    c(median(all.signals[hm450.hg19.probe2chr[names(all.signals)] == 'chrX']),
      median(all.signals[hm450.hg19.probe2chr[names(all.signals)] == 'chrY']))
  }))
  colnames(XYmedian) <- c('x','y')
  genders <- ifelse(XYmedian[,'y'] < 500, 0,
                    ifelse(XYmedian[,'y'] < 500 & XYmedian[,'x'] > 5000, 0, 1))
  genders
}
