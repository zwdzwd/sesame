library(devtools)
load_all('..', export_all=FALSE)

dms <- ReadIDATsFromDir('data/EPIC.beta.test/')
ssets <- lapply(dms, ChipAddressToSignal)
pvals <- lapply(ssets, DetectPValue)
ssets.noob <- lapply(ssets, BackgroundCorrectionNoob)
ssets.noob.dye <- DyeBiasCorrectionMostBalanced(ssets.noob)
betas <- mapply(SignalToBeta, ssets.noob.dye, pvals)

