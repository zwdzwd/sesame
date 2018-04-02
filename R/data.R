## content in this file is obsolete

## chromosome information (deprecated)
## HM27.hg19.probe2chr
## HM450.hg19.probe2chr
## EPIC.hg19.probe2chr

## for sex inference
## sex.inference.model
## EPIC.female.clean.chrY.probes
## HM450.female.clean.chrY.probes
## EPIC.female.xlinked.chrX.probes
## HM450.female.xlinked.chrX.probes

## production probes
## HM27.ordering
## HM450.ordering
## EPIC.ordering

## the control probes
## HM27.controls
## HM450.controls
## EPIC.controls

## mask for snp and repeats
## HM27.mask
## HM450.mask
## EPIC.mask

## for cnv segmentation
## EPIC.mapped.probes.hg19
## EPIC.mapped.probes.hg38
## HM450.mapped.probes.hg19
## HM450.mapped.probes.hg38
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
        platform <- toupper(platform) # HM450, HM27 or EPIC
        datanm <- paste0(platform,'.',nm)
    } else {
        datanm <- nm
    }
    if (exists(datanm, envir=cacheEnv)) {
        return(get(datanm, envir=cacheEnv))
    }

    x <- NULL
    base.dir <- 'http://zwdzwd.io/sesame/current/'
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
        message('Done.')
    }
    assign(datanm, x, envir=cacheEnv)
    x
    
}

## Cache all the built in data from remote
##
## need environment SESAMEHOME be set
##
## @importFrom XML htmlTreeParse
## @importFrom utils download.file
## @return the object NULL
## @examples
## \dontrun{
## cacheBuiltInData() # download annotation data to $SESAMEHOME
## }
##
## cat("Data will be deposited to", Sys.getenv('SESAMEHOME'), "\n")
cacheBuiltInData <- function() {

    pkgTest('XML')
    pkgTest('utils')
    seshome <- Sys.getenv('SESAMEHOME')
    if (seshome == '') {
        stop("SESAMEHOME is not set. Abort caching.")
    }
    dir.create(seshome, showWarnings = FALSE)
    base.dir <- 'http://zwdzwd.io/sesame/current/'

    dwfiles <- unname(
        unlist(XML::htmlTreeParse(base.dir, useInternalNodes=TRUE)["//a/@href"]))

    for (dwfile in dwfiles[grep('*.rds$', dwfiles)]) {
        download.file(
            paste0(base.dir,'/',dwfile),
            paste0(seshome,'/',dwfile))
    }

    subname <- 'cellref'
    subdir <- paste0(base.dir, '/', subname, '/')
    dwfiles <- unname(
        unlist(XML::htmlTreeParse(subdir, useInternalNodes=TRUE)["//a/@href"]))
    
    dir.create(paste0(seshome, '/', subname), showWarnings = FALSE)
    for (dwfile in dwfiles[grep('*.rds$', dwfiles)]) {
        download.file(
            paste0(base.dir,'/',subname,'/',dwfile),
            paste0(seshome,'/',subname,'/',dwfile))
    }

    subname <- 'examples'
    subdir <- paste0(base.dir, '/', subname, '/')
    dwfiles <- unname(
        unlist(XML::htmlTreeParse(subdir, useInternalNodes=TRUE)["//a/@href"]))
    
    dir.create(paste0(seshome, '/', subname), showWarnings = FALSE)
    for (dwfile in dwfiles[grep('*.rds$', dwfiles)]) {
        download.file(
            paste0(base.dir,'/',subname,'/',dwfile),
            paste0(seshome,'/',subname,'/',dwfile))
    }
    NULL
}

## betas <- SeSAMeGetExample('HM450.betas.TCGA-2L-AAQA-01A-21D-A38H-05')
SeSAMeGetExample <- function(nm, platform='') {
    getBuiltInData(nm, platform, 'examples');
}
