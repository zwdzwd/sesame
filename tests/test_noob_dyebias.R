create.validation <- function() {
  library(methylumi)
  d <- getwd()
  setwd('data/tcga.random6/')
  map <- read.csv('samples.csv', stringsAsFactors=F)
  MSET <- methylumIDAT(map, parallel=T)
  MSET2 <- methylumi.bgcorr(MSET)
  MSET3 <- stripMethyLumiSet(MSET2)
  mset <- normalizeMethyLumiSet(MSET3)
  setwd(d)

  mset <- processMethylumi()
  betas.save <- betas(mset)
  save(betas.save, file='data/tcga.random6/test_noob_dyebias.Rout.rda')
}

# create.validation

library(devtools)
load_all("..",export_all=FALSE)

dms <- ReadIDATsFromSampleSheet(
  "data/tcga.random6/samples.csv", base.dir='data/tcga.random6')
dmps <- lapply(dms, ChipAddressToProbe)
pvals <- lapply(dmps, DetectPValue)
dmps.noob <- lapply(dmps, BackgroundCorrectionNoob)
dmps.noob.dye <- DyeBiasCorrectionMostBalanced(dmps.noob)
betas <- mapply(SignalToBeta, dmps.noob.dye, pvals)

load('data/tcga.random6/test_noob_dyebias.Rout.rda')
MPrint('Are testing the same as validation: ', all.equal(betas, betas.save))

MPrint('Head old:')
print(betas.save[1:5,])

MPrint('Head new:')
print(betas[1:5,])

MPrint('Tail old:')
print(betas.save[(nrow(betas.save)-5):nrow(betas.save),])

MPrint('Tail new:')
print(betas[(nrow(betas)-5):nrow(betas),])
