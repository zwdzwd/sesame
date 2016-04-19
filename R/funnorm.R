#' QuantileSet
#'
#' Construct a QuantileSet object
#' 
#' This class holds the information for funnorm quantile normalization
#' 
#' @param auto.IG.M quantiles for type I, green channel, autosomal probes "methylated" allele
#' @param auto.IG.U quantiles for type I, green channel, autosomal probes "unmethylated" allele
#' @param auto.IR.M quantiles for type I red channel, autosomal probes "methylated" allele
#' @param auto.IR.U quantiles for type I, red channel, autosomal probes "unmethylated" allele
#' @param auto.II.M quantiles for type II, autosomal, probes "methylated" allele
#' @param auto.II.U quantiles for type II, autosomal, probes "unmethylated" allele
#' @param X.all.M quantiles for X chromosome probes, both type I and II, "methylated" allele
#' @param X.all.U quantiles for X chromosome probes, both type I and II, "unmethylated" allele
#' @return an object of class \code{QuantileSet}
#' @export
QuantileSet <- function(
  auto.IG.M=NULL, auto.IG.U=NULL,
  auto.IR.M=NULL, auto.IR.U=NULL,
  auto.II.M=NULL, auto.II.U=NULL,
  X.all.M=NULL, X.all.U=NULL) {
  
  structure(list(
    auto.IG.M=auto.IG.M, auto.IG.U=auto.IG.U,
    auto.IR.M=auto.IR.M, auto.IR.U=auto.IR.U,
    auto.II.M=auto.II.M, auto.II.U=auto.II.U,
    X.all.M=X.all.M, X.all.U=X.all.U), class="QuantileSet")
}


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
  cm <- c(cm, bisulfite2=mean(ctls[['BISULFITE CONVERSION II']]$R, na.rm=TRUE))

  ## bisulfite conversion type I
  cm <- c(cm, bisulfite1=mean(
                ctls[['BISULFITE CONVERSION I']][sprintf('BS.Conversion.I.C%s', 1:3),'G'] +
                  ctls[['BISULFITE CONVERSION I']][sprintf('BS.Conversion.I.C%s', 4:6),'R'],
                na.rm=TRUE))

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
#' @return an object of \code{QuantileSet}
#' @export
BuildQuantiles450k <- function(dmp, n=500) {

  data(hm450.hg19.probe2chr)
  
  probs <- seq(from=0, to=1, length.out=n)

  auto.IG <- dmp$IG[!(hm450.hg19.probe2chr[rownames(dmp$IG)] %in% c('chrX','chrY')),]
  auto.IG.M <- quantile(auto.IG[,'M'], probs)
  auto.IG.U <- quantile(auto.IG[,'U'], probs)

  auto.IR <- dmp$IR[!(hm450.hg19.probe2chr[rownames(dmp$IR)] %in% c('chrX','chrY')),]
  auto.IR.M <- quantile(auto.IR[,'M'], probs)
  auto.IR.U <- quantile(auto.IR[,'U'], probs)

  auto.II <- dmp$II[!(hm450.hg19.probe2chr[rownames(dmp$II)] %in% c('chrX','chrY')),]
  auto.II.M <- quantile(auto.II[,'M'], probs)
  auto.II.U <- quantile(auto.II[,'U'], probs)

  all <- rbind(dmp$IG, dmp$IR, dmp$II)
  X.all <- all[hm450.hg19.probe2chr[rownames(all)] == 'chrX',]
  X.all.M <- quantile(X.all[,'M'], probs)
  X.all.U <- quantile(X.all[,'U'], probs)
  ## Y.all <- all[hm450.hg19.probe2chr[rownames(all)] == 'chrY',]
  ## quantile.Y.all.M <- quantile(Y.all[,'M'], probs)
  ## quantile.Y.all.U <- quantile(Y.all[,'U'], probs)

  QuantileSet(auto.IG.M, auto.IG.U,
              auto.IR.M, auto.IR.U,
              auto.II.M, auto.II.U,
              X.all.M, X.all.U)
}

#' Interpolate signal
#'
#' Interpolate signal using quantile distribution
#'
#' @param dmp an object of class \code{SignalSet}
#' @param qs an object of class \code{QuantileSet}
InterpolateSignalUseQntiles <- function(dmp, qntiles) {

  data(hm450.hg19.probe2chr)
  auto.IG <- dmp$IG[!(hm450.hg19.probe2chr[rownames(dmp$IG)] %in% c('chrX','chrY')),]
  auto.IR <- dmp$IR[!(hm450.hg19.probe2chr[rownames(dmp$IR)] %in% c('chrX','chrY')),]
  auto.II <- dmp$II[!(hm450.hg19.probe2chr[rownames(dmp$II)] %in% c('chrX','chrY')),]

  all <- rbind(dmp$IG, dmp$IR, dmp$II)
  category <- do.call(c, lapply(c('IG','IR','II'), function(nm) rep(nm, length(dmp[[nm]]))))
  X.all <- all[hm450.hg19.probe2chr[rownames(all)] == 'chrX',]

  lapply(c('auto.IG','auto.IR','auto.II','X.all'), function(nm) {
    lapply(c('M','U'), function(a) {
      assign(paste0(nm,'.',a), get(nm)[,a])
    })
  })

  ## actual interpolation for each category
  separate.n <- sapply(names(qntiles), function(category)
    .InterpolateSignalUseQntiles1(get(category), qntiles[[category]]), USE.NAMES=TRUE, simplify=FALSE)

  ## build back a SignalSet
  s <- do.call(SignalSet, sapply(c('IR','IG','II'), function(nm) {
    d <- cbind(get(paste0('auto.', nm,'.M')), get(paste0('auto.', nm, '.U')))
    x <- cbind(separate.n[['X.all.M']], separate.n[['X.all.U']])
    d <- rbind(d, x[category[names(x)] == nm])
    d
  }, USE.NAMES=TRUE, simplify=FALSE))
}

.InterpolateSignalUseQntiles1 <- function(signal, qntiles) {
  n <- length(qntiles)
  qntiles.densified <- do.call(c, lapply(1:(n-1), function(i) {
    seq(qntiles[i], qntiles[i+1], (qntiles[i+1]-qntiles[i])/n)[-(n+1)]
  }))
  preprocessCore.normalize.quantiles.use.target(signal, qntiles.densified)
}

#' funnorm regress
#'
#' regress out control effect from quantile distribution
#'
#' @param cms control matrices from all samples
#' @param qntiles quantile distributions from all samples
#' @param k number of principal components from control matrices for regression
#' @import preprocessCore
#' @return a list of new quantile distriubtions from all samples
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
  stopifnot(identical(names(qntiles),rownames(mm)))
  
  ## Quantile normalize - autosomes
  mprint("Normalizing autosomes.")
  FunnormRegress1 <- function(qmat, mm, k=2) {
    stopifnot(identical(colnames(qmat), rownames(mm)))
    qmat[1,] <- 0
    qmat[nrow(qmat),] <- qmat[nrow(qmat)-1,] + 1000
    global.mean.quantiles <- rowMeans(qmat, na.rm=TRUE)
    res <- qmat - global.mean.quantiles
    mm.pcs <- prcomp(mm)$x[,1:k,drop=FALSE]
    design <- model.matrix(~mm.pcs)
    fits <- lm.fit(x = design, y = t(res))
    qmat.n <- global.mean.quantiles + t(fits$residuals)
  }

  auto.names <- sprintf('auto.I%s.%s', rep(c('R','G','I'),2), rep(c('M','U'),each=3))
  qmat.n <- lapply(
    auto.names,
    function(target.nm) {
      FunnormRegress1(
        vapply(qntiles, function(x) x[[target.nm]],
               numeric(length(qntiles[[1]][[1]]))), mm, k=k)
    })
  names(qmat.n) <- auto.names

  ## Quantile normalize X chromosome
  fmle <- genders[names(qntiles)]==0
  male <- genders[names(qntiles)]==1
  mprint("Normalizing X-chromosome.")
  if ((!is.null(genders)) &&
      sum(genders == 0)>10 && sum(genders == 1)>10) {
    qmat.n[['X.all.M']] <- cbind(
      FunnormRegress1(sapply(
        qntiles[fmle], function(x) x$X.all.M), mm[fmle,], k=k),
      FunnormRegress1(sapply(
        qntiles[male], function(x) x$X.all.M), mm[male,], k=k))
    qmat.n[['X.all.U']] <- cbind(
      FunnormRegress1(sapply(
        qntiles[fmle], function(x) x$X.all.U), mm[fmle,], k=k),
      FunnormRegress1(sapply(
        qntiles[male], function(x) x$X.all.U), mm[male,], k=k))
  } else {
    qmat.n[['X.all.M']] <- FunnormRegress1(sapply(
      qntiles, function(x) x$X.all.M), mm, k=k)
    qmat.n[['X.all.U']] <- FunnormRegress1(sapply(
      qntiles, function(x) x$X.all.U), mm, k=k)
  }

  ## Wrap up and return
  qntiles.n <- sapply(
    names(qntiles), function(sample) {
      do.call(QuantileSet, lapply(qmat.n, function(m) m[,sample]))},
    USE.NAMES=TRUE, simplify=FALSE)
}

#' Quantile normalize
#'
#' Quantile normalize a list of \code{SignalSet}s
#'
#' This is a wrapper of the preprocessCore's quantile normalization
#' 
#' @param dmps a list of objects of class \code{SignalSet}
#' @param genders specify whether to normalize separately by gender
#' @return dmps.n a list of objects of class \code{SignalSet} after normalization
#' @import preprocessCore
#' @export
QuantileNormalize <- function(dmps, genders=NULL) {

  M <- lapply(dmps, function(dmp) {
    c(dmp$IR[,'M'], dmp$IG[,'M'], dmp$II[,'M'])
  })

  U <- lapply(dmps, function(dmp) {
    c(dmp$IR[,'U'], dmp$IG[,'U'], dmp$II[,'U'])
  })

  if ((!is.null(genders)) &&
      any(genders == 0) && any(genders == 1)) {
    M.n <- cbind(
      preprocessCore::normalize.quantiles(M[genders[colnames(M)]==0]),
      preprocessCore::normalize.quantiles(M[genders[colnames(M)]==0]))
    U.n <- cbind(
      preprocessCore::normalize.quantiles(U[genders[colnames(U)]==0]),
      preprocessCore::normalize.quantiles(U[genders[colnames(U)]==0]))
  } else {
    M.n <- preprocessCore::normalize.quantiles(M)
    U.n <- preprocessCore::normalize.quantiles(U)
  }
  
  all.n <- cbind(M=M.n, U=U.n, rownames=rownames(M))
  nr <- vapply(c('IR','IG','II'), function(x) nrow(dmp[[x]]), numeric(1))
  SignalSet(
    IR = all.n[1:nr[1],], 
    IG = all.n[(nr[1]+1):(nr[1]+nr[2]),], 
    II = all.n[(nr[1]+nr[2]+1):(nr[1]+nr[2]+nr[3]),])
}

#' Infer gender from signals
#'
#' Infer gender from signals
#' @param dmps a list of \code{SignalSet}s
#' @return named vector of 0 (for female) and 1 (for male)
#' @export
InferGenders <- function(dmps) {
  XYmedian <- t(vapply(dmps, function(dmp) {
    all.signals <- c(apply(dmp$IR,1,max), apply(dmp$IG,1,max), apply(dmp$II,1,max))
    c(median(all.signals[hm450.hg19.probe2chr[names(all.signals)] == 'chrX']),
      median(all.signals[hm450.hg19.probe2chr[names(all.signals)] == 'chrY']))
  }, numeric(2)))
  colnames(XYmedian) <- c('x','y')
  genders <- ifelse(XYmedian[,'y'] < 500, 0,
                    ifelse(XYmedian[,'y'] < 500 & XYmedian[,'x'] > 5000, 0, 1))
  genders
}
