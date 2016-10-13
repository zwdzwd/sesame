
## map probe to chromosome
#' hm27 probe2chromosome
'hm27.hg19.probe2chr'
#' hm450 probe2chromosome
'hm450.hg19.probe2chr'
#' EPIC probe2chromosome
'EPIC.hg19.probe2chr'

#' EPIC female clean chrY probes
'EPIC.female.clean.chrY.probes'
#' hm450 female clean chrY probes
'hm450.female.clean.chrY.probes'

#' sex inference model
'sex.inference.model'

#' EPIC female X linked probes
'EPIC.female.xlinked.chrX.probes'
#' hm450 female X linked probes
'hm450.female.xlinked.chrX.probes'

## the design for each probes
#' hm27 data ordering
'hm27.ordering'
#' hm450 data ordering
'hm450.ordering'
#' EPIC data ordering
'EPIC.ordering'

## the control probes
#' hm27 controls
'hm27.controls'
#' hm450 controls
'hm450.controls'
#' EPIC controls
'EPIC.controls'

## mask for snp and repeats
#' hm27 masks
'hm27.mask'
#' hm450 masks
'hm450.mask'
#' EPIC masks
'EPIC.mask'

#' EPIC mapped probes hg19
'EPIC.mapped.probes.hg19'
#' EPIC mapped probes hg38
'EPIC.mapped.probes.hg38'
#' hm450 mapped probes hg19
'hm450.mapped.probes.hg19'
#' hm450 mapped probes hg38
'hm450.mapped.probes.hg38'

#' hg19 chromosome info
'hg19.chrominfo'
#' hg38 chromosome info
'hg38.chrominfo'

#' ccs probes for ethnicity inference
'ethnicity.ccs.probes'
#' rs probes for ethnicity inference
'ethnicity.rs.probes'
#' model for ethnicity inference
'ethnicity.model'

## #' EPIC masks
## 'EPIC.mask'

cacheEnv <- new.env()

getBuiltInData <- function(nm, platform) {
  datanm <- paste0(platform,'.',nm)
  if (!exists(datanm, envir=cacheEnv)) {
    return(get(datanm, envir=cacheEnv))
  }

  base.dir <- 'http://zwdzwd.io/sesame/20161013/'

  x <- readRDS(gzcon(url(paste(base.dir, datanm, sep='/'))))
  assign(datanm, x, envir=cacheEnv)

  x
}

