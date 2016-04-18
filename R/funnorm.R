#' build control matrix for funnorm
#' 
#' build control matrix for funnorm
#'
#' @param dmp an object of class SignalSet
#' @return a vector with control summaries
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
cm <- c(cm, stain.red=ctls[['STAINING']]['DNP..High.', 'R'], stain.green=ctls[['STAINING']]['Biotin..High.','G'])

## extension
cm <- c(cm,
extRed1=ctls[['EXTENSION']]['Extension..A.','R'],
extRed2=ctls[['EXTENSION']]['Extension..T.','R'],
extGrn1=ctls[['EXTENSION']]['Extension..C.','G'],
extGrn2=ctls[['EXTENSION']]['Extension..G.','G'])

## hybridization
d <- ctls[['HYBRIDIZATION']]$G
cm <- c(cm, structure(d, names=paste0('hybe',1:length(d))))

## target removal
d <- ctls[['TARGET REMOVAL']]$G
cm <- c(cm, structure(d, names = paste0('targetrem',1:length(d))))

## non-polymorphic
cm <- c(cm,
  nonpolyRed1=ctls[['NON-POLYMORPHIC']]['NP..A.','R'],
  nonpolyRed2=ctls[['NON-POLYMORPHIC']]['NP..T.','R'],
  nonpolyGrn1=ctls[['NON-POLYMORPHIC']]['NP..C.','G'],
  nonpolyGrn2=ctls[['NON-POLYMORPHIC']]['NP..G.','G'])

## specificity type II
d <- ctls[['SPECIFICITY II']]
cm <- c(cm, structure(d$G, names=paste0('spec2Grn', 1:dim(d)[1])), structure(d$R, names=paste0('spec2Red', 1:dim(d)[1])))
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
#' quantiles are stratified by IR autosomal, IG autosomal, II autosomal, X all design and Y all design.
#' @param dmp a SignalSet
#' @param n number of grid points in quantiles.
#' @return a list of quantiles in each category
BuildQuantiles450k <- function(dmp, n=500) {

library()
hm450.hg19 <- get450k()

library(matrixStats)
probs <- seq(from=0, to=1, length.out=n)

auto.IR <- dmp$IR[as.vector(!(seqnames(hm450.hg19[rownames(dmp$IR),]) %in% c('chrX','chrY'))),]
quantile.auto.IR <- quantile(auto.IR, probs)

auto.IG <- dmp$IG[as.vector(!(seqnames(hm450.hg19[rownames(dmp$IG),]) %in% c('chrX','chrY'))),]
quantile.auto.IG <- quantile(auto.IG, probs)

auto.II <- dmp$II[as.vector(!(seqnames(hm450.hg19[rownames(dmp$II),]) %in% c('chrX','chrY'))),]
quantile.auto.II <- quantile(auto.II, probs)

all <- rbind(dmp$IR, dmp$IG, dmp$II)
X.all <- all[as.vector(seqnames(hm450.hg19[rownames(all),]) == 'chrX'),]
quantile.X.all <- quantile(X.all, probs)

Y.all <- all[as.vector(seqnames(hm450.hg19[rownames(all),]) == 'chrY'),]
quantile.Y.all <- quantile(Y.all, probs)

list(auto.IR=quantile.auto.IR, auto.IG=quantile.auto.IG, auto.II=quantile.auto.II, X.all=quantile.X.all, Y.all=quantile.Y.all)

}
