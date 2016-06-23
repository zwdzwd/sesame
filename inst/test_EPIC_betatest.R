#!/usr/bin/env Rscript

library(sesame)
library(parallel)
options(mc.cores=4)

ssets <- readIDATsFromDir('../../test/data/EPIC.beta.test/', mc=T)
ssets <- mclapply(ssets, noob)
ssets <- mclapply(ssets, dyeBiasCorr)

sset <- ssets[[1]]
cat(sset$inferSex())
