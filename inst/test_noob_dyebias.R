#!/usr/bin/env Rscript

create.validation <- function() {
  library(methylumi)
  d <- getwd()
  setwd('../../test/data/tcga.random6/')
  map <- read.csv('samples.csv', stringsAsFactors=F)
  MSET <- methylumIDAT(map, parallel=T)
  MSET2 <- methylumi.bgcorr(MSET)
  MSET3 <- stripMethyLumiSet(MSET2)
  mset <- normalizeMethyLumiSet(MSET3)
  setwd(d)

  mset <- processMethylumi()
  betas.save <- betas(mset)
  save(betas.save, file='../../test/data/tcga.random6/test_noob_dyebias.Rout.rda')
}

# create.validation

library(sesame)
library(parallel)

options(mc.cores=8)
ssets <- readIDATsFromSheet("../../test/data/tcga.random6/samples.csv", base.dir='../../test/data/tcga.random6', mc=T)
ssets <- mclapply(ssets, noob)
ssets <- mclapply(ssets, dyeBiasCorr)
betas <- signalToBeta(ssets, na.mask=F, mc=T)

load('../../test/data/tcga.random6/test_noob_dyebias.Rout.rda')
message('Are testing the same as validation: ', all.equal(betas, betas.save))

message('Head old:')
print(betas.save[1:5,])

message('Head new:')
print(betas[1:5,])

message('Tail old:')
print(betas.save[(nrow(betas.save)-5):nrow(betas.save),])

message('Tail new:')
print(betas[(nrow(betas)-5):nrow(betas),])

message('is.na old:', sum(is.na(betas.save)))
message('is.na new:', sum(is.na(betas)))
