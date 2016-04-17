#' build control matrix for funnorm
#' 
#' build control matrix for funnorm
#'
#' \param 
BuildControlMatrix450k <- function(dmp) {

ctls <- split(dmp$ctl, dmp$ctl$type)

ctl.mat <- NULL

## bisulfite conversion type II
ctl.mat <- c(ctl.mat, bisulfite2=mean(red.ctls[['BISULFITE CONVERSION II']]$R,na.rm=T))

## bisulfite conversion type I
ctl.mat <- c(ctl.mat, bisulfite1=mean(
 ctls[['BISULFITE CONVERSION I']][sprintf('BS.Conversion.I.C%s', 1:3),'G'] +
 ctls[['BISULFITE CONVERSION I']][sprintf('BS.Conversion.I.C%s', 4:6),'R']))

## staining
ctl.mat <- c(ctl.mat, stain.red=ctls[['STAINING']]['DNP..High.', 'R'], stain.green=ctls[['STAINING']]['Biotin..High.','G'])

## extension
ctl.mat <- c(ctl.mat, 
extRed1=ctls[['EXTENSION']]['Extension..A.','R'],
extRed2=ctls[['EXTENSION']]['Extension..T.','R'],
extGrn1=ctls[['EXTENSION']]['Extension..C.','G'],
extGrn2=ctls[['EXTENSION']]['Extension..G.','G'])

## hybridization
d <- ctls[['HYBRIDIZATION']]$G
ctl.mat <- c(ctl.mat, structure(d, names=paste0('hybe',1:length(d))))

## target removal
d <- ctls[['TARGET REMOVAL']]$G
ctl.mat <- c(ctl.mat, structure(d, names = paste0('targetrem',1:length(d))))

## non-polymorphic
ctl.mat <- c(ctl.mat,
nonpolyRed1=ctls[['NON-POLYMORPHIC']]['NP..A.','R'],
nonpolyRed2=ctls[['NON-POLYMORPHIC']]['NP..T.','R'],
nonpolyGrn1=ctls[['NON-POLYMORPHIC']]['NP..C.','G'],
nonpolyGrn2=ctls[['NON-POLYMORPHIC']]['NP..G.','G'])

## specificity type II
d <- ctls[['SPECIFICITY II']]
cm <- c(cm, structure(d$G, names=paste0('spec2Grn', 1:dim(d)[1])), structure(d$R, names=paste0('spec2Red', 1:dim(d)[1])))
cm <- c(cm, spec2.ratio = mean(d,na.rm=TRUE) / mean(d,na.rm=TRUE))

## specificty type I green
d <- ctls[['SPECIFICITY I']][sprintf('GT.Mismatch.%s..PM.',1:3),]
cm <- c(cm, structure(d$G, names=paste0('spec1Grn',1:dim(d)[1])))
cm <- c(cm, spec1.ratio.grn = mean(d$R, na.rm=TRUE)/mean(d$G, na.rm=TRUE))

## specifity type I red
d <- ctls[['SPECIFICITY I']][sprintf('GT.Mismatch.%s..PM.',4:6),]
cm <- c(cm, structure(d$R, names=paste0('spec1Red',1:dim(d)[1])))
cm <- c(cm, spec1.ratio.red = mean(d$G, na.rm=TRUE)/mean(d$R, na.rm=TRUE))

## normalization
cm <- c(cm, c(normA=mean(ctls[['NORM_A']]$R, na.rm=TRUE), 
  normT=mean(ctls[['NORM_T']]$R, na.rm=TRUE), 
  normC=mean(ctls[['NORM_C']]$G, na.rm=TRUE), 
  normG=mean(ctls[['NORM_G']]$G, na.rm=TRUE)))

dyebias <- structure((normC+normG) / (normA+normT)

## out of band probes
cm <- c(cm, oob.ratio=median(dmp$oobG) / median(dmp$oobR))
cm <- c(cm, structure(quantile(dmp$oobG, c(0.01,0.5,0.99)), names=paste0('oob', c(1,50,99))))

}
