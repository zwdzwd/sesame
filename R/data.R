
## map probe to chromosome
'hm450.hg19.probe2chr'
'EPIC.hg19.probe2chr'

## the design for each probes
'hm450.ordering'
'EPIC.ordering'

## the control probes
'hm450.controls'
'EPIC.controls'

## mask for snp and repeats
'hm450.mask'

GetBuiltInData <- function(nm, platform) {
  datanm <- paste0(platform,'.',nm)
  data(list=datanm)
  get(datanm)
}

