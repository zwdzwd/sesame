
## chromosome information (deprecated)
## hm27.hg19.probe2chr
## hm450.hg19.probe2chr
## EPIC.hg19.probe2chr

## for sex inference
## sex.inference.model
## EPIC.female.clean.chrY.probes
## hm450.female.clean.chrY.probes
## EPIC.female.xlinked.chrX.probes
## hm450.female.xlinked.chrX.probes

## production probes
## hm27.ordering
## hm450.ordering
## EPIC.ordering

## the control probes
## hm27.controls
## hm450.controls
## EPIC.controls

## mask for snp and repeats
## hm27.mask
## hm450.mask
## EPIC.mask

## for cnv segmentation
## EPIC.mapped.probes.hg19
## EPIC.mapped.probes.hg38
## hm450.mapped.probes.hg19
## hm450.mapped.probes.hg38
## hg19.chrominfo
## hg38.chrominfo

### ccs probes for ethnicity inference
## ethnicity.ccs.probes
## rs probes for ethnicity inference
## ethnicity.rs.probes
## model for ethnicity inference
## ethnicity.model

## EPIC masks
## EPIC.mask

cacheEnv <- new.env()

getBuiltInData <- function(nm, platform='') {
  if (platform != '') {
    datanm <- paste0(platform,'.',nm)
  } else {
    datanm <- nm
  }
  if (exists(datanm, envir=cacheEnv)) {
    return(get(datanm, envir=cacheEnv))
  }

  base.dir <- 'http://zwdzwd.io/sesame/20161013/'

  message('Downloading ', datanm, ' from ', base.dir, '.')
  x <- readRDS(gzcon(url(paste0(base.dir, datanm, '.rds'))))
  assign(datanm, x, envir=cacheEnv)

  x
}

