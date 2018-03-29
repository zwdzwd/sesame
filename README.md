# SeSAMe - SEnsible Step-wise Analysis of Methylation data [![Travis-CI Build Status](https://travis-ci.org/zwdzwd/sesame.svg?branch=master)](https://travis-ci.org/zwdzwd/sesame)
                  
## Topics

   0. Introduction
   1. Install
   2. Usage
   3. Bugs
   4. About

## Introduction

SeSAMe is an R package for processing Infinium DNA methylation data. SeSAMe currently supports EPIC, HM450 and HM27 platforms.

## Install
    
You need devtools from install.packages('devtools'). Then,

```R
library(devtools)
install_github("zwdzwd/sesame")
```

optional: download data for offline processing. Otherwise, sesame needs to access internet for data needed

```R
Sys.setenv(SESAMEHOME='<my.download.path>')
cacheBuiltInData()
```
put `export SESAMEHOME=<my.download.path>` to .bashrc to avoid future setting of SESAMEHOME

## Usage

### Data IO
    
To read IDATs from directory
```R
library(sesame)
ssets <- readIDATsFromDir('IDATs/')
```
        
A simple list of "SignalSet"s are returned.
        
Other options of data import include

+ readIDATs - from path
+ readIDATsFromSheet - from data.frame "barcode" column
        
### Background subtraction
    
noob:
```R
noob(sset)
```
      
### Dye bias correction

linear scaling
```R
dyeBiasCorr(sset)
```
      
quantile interpolation with Type-I probes
```R
dyeBiasCorrTypeINorm(sset)
```
      
### Get betas

```R
getBetas(sset)
```
use option quality.mask = TRUE/FALSE to mask probes with mapping issues, SNPs 	and non-uniqueness, described in Zhou 2016 NAR.
use option nondetection.mask = TRUE/FALSE to mask nondetection calls.

Both masks are recommended to ensure data quality and defaulted to TRUE.
        
### Sample/experiment QC

#### sex
```R
inferSex(sset)
```

#### ethnicity
```R
inferEthnicity(sset)
```

#### age
```R
betas <- SeSAMeGetExample('HM450.betas.TCGA-2L-AAQA-01A-21D-A38H-05')
predictAgeHorvath353(betas)
```

#### mean intensity
```R
meanIntensity(sset)
```

#### bisulfite conversion control using [GCT score](https://academic.oup.com/nar/article/45/4/e22/2290930)
```R
bisConversionControl(sset)
```

### visualize probes

```R
betas <- sesameGetExample('hm450.betas.TCGAnormalPAAD')
```

visualize probes from a gene
```R
visualizeGene('DNMT1', betas, platform='hm450')
```

visualize probes from arbitrary region
```R
visualizeRegion('chr19',10260000,10380000, betas, platform='hm450')
```

visualize by probe names
```R
visualizeProbes(c("cg02382400", "cg03738669"), betas, platform='hm450')
```

### Differential methylation

test differential methylation on each locus
```R
betas <- sesameGetExample('HM450.betas.76matchedTCGAchr20')
sample.info <- sesameGetExample('HM450.sampleinfo.76matchedTCGAchr20')
cf <- DMR(betas, sample.info, ~patient+type, platform='hm450')
```

top differentially methylated probes on factor "type"
```R
cf1 <- cf$typeTumour
head(topLoci(cf1))
```

top segments
```R
head(topSegments(cf1))
```
      
visualize top DMR
```R
visualizeProbes(rownames(cf1[cf1[,'Seg.ID']==topSegments(cf1)$Seg.ID[1],]), betas, upstream=5000, dwstream=5000, platform='hm450',heat.height=3.5)
```


### CNV
    
uses LNCaP EPIC data from GSE86833
```R
example.sset <- SeSAMeGetExample('EPIC.sset.LNCaP.Rep1')
segs <- cnSegmentation(example.sset)
```

To visualize segments,
```R
visualizeSegments(segs)
```

### cell composition deconvolution

Use blood set as example,

```R
g <- diffRefSet(getRefSet(c('CD4T','CD19B','CD14Monocytes','CD56NK', 'granulocytes'), platform='hm450'))
betas <- SeSAMeGetExample('HM450.betas.TCGA-2L-AAQA-01A-21D-A38H-05')
estimateCellComposition(g, betas[rownames(g)])$frac
```

## Bugs
    
Bug reports are appreciated. Register issues at the sesame [issue tracker](http://github.com/zwdzwd/sesame/issues)
    
    
## About

W Zhou, PW Laird, H Shen, in submission
    
