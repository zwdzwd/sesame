
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

getBuiltInData <- function(nm, platform='', subdir='') {
  if (platform != '') {
    datanm <- paste0(platform,'.',nm)
  } else {
    datanm <- nm
  }
  if (exists(datanm, envir=cacheEnv)) {
    return(get(datanm, envir=cacheEnv))
  }

  base.dir <- 'http://zwdzwd.io/sesame/20161013/'
  if (subdir != '') {
    base.dir <- paste0(base.dir, subdir, '/');
  }

  seshome <- Sys.getenv('SESAMEHOME')
  dir.create(seshome, showWarnings = FALSE) # make sure directory exists
  localpath <- paste0(seshome, '/', datanm, '.rds')
  if (seshome != "") {
    if (file.exists(localpath)) {
      x <- readRDS(localpath)
    } else {
      x <- NULL
    }
  }

  if (is.null(x)) {
    message('Caching ', datanm, '... ', appendLF=FALSE)
    x <- readRDS(gzcon(url(paste0(base.dir, datanm, '.rds'))))
    if (seshome != "") {
      saveRDS(x, file=localpath)
    }
  }
  assign(datanm, x, envir=cacheEnv)
  message('Done.')
  x
}

#' Cache all the built in data from remote
#'
#' need environment SESAMEHOME be set
#'
#' @importFrom XML htmlTreeParse
#' 
#' @export
cacheBuiltInData <- function() {
  seshome <- Sys.getenv('SESAMEHOME')
  if (seshome == '') {
    stop("SESAMEHOME is not set. Abort caching.")
  }
  dir.create(seshome, showWarnings = FALSE)
  base.dir <- 'http://zwdzwd.io/sesame/20161013/'
  dwfiles <- unname(unlist(htmlTreeParse(base.dir, useInternalNodes=T)["//a/@href"]))
  for (dwfile in dwfiles[grep('*.rds$', dwfiles)]) {
    download.file(paste0(base.dir,'/',dwfile), paste0(seshome,'/',dwfile))
  }

  subname <- 'cellref'
  subdir <- paste0(base.dir, '/', subname, '/')
  dwfiles <- unname(unlist(htmlTreeParse(subdir, useInternalNodes=T)["//a/@href"]))
  dir.create(paste0(seshome, '/', subname), showWarnings = FALSE)
  for (dwfile in dwfiles[grep('*.rds$', dwfiles)]) {
    download.file(paste0(base.dir,'/',subname,'/',dwfile), paste0(seshome,'/',subname,'/',dwfile))
  }

  subname <- 'examples'
  subdir <- paste0(base.dir, '/', subname, '/')
  dwfiles <- unname(unlist(htmlTreeParse(subdir, useInternalNodes=T)["//a/@href"]))
  dir.create(paste0(seshome, '/', subname), showWarnings = FALSE)
  for (dwfile in dwfiles[grep('*.rds$', dwfiles)]) {
    download.file(paste0(base.dir,'/',subname,'/',dwfile), paste0(seshome,'/',subname,'/',dwfile))
  }
}

#' Retrieve SeSAMe examples
#'
#' @param nm name
#' @param platform optional, hm450, EPIC or hm27
#' @return example object
#' @export
sesameGetExample <- function(nm, platform='') {
  getBuiltInData(nm, platform, 'examples');
}
