## for hm450

library('FDb.InfiniumMethylation.hg19')

hm450.hg19 <- get450k()

hm450.hg19.df <- as.data.frame(hm450.hg19)
hm450.hg19.probe2chr <- setNames(hm450.hg19.df$seqnames, rownames(hm450.hg19.df))

save(hm450.hg19.probe2chr, file='../data/hm450.hg19.probe2chr.rda')

save(hm450.ordering, file='../data/hm450.ordering.rda')
save(hm450.controls, file='../data/hm450.controls.rda')
save(hm450.mask, file='../data/hm450.mask.rda')

## for EPIC

## sensible and correct sourceseq (same as the probe sequence direction), coordinates for rs probes.

manifest <- read.csv('tests/data//EPIC.official.manifest/MethylationEPIC_v-1-0_B2.csv',header=T,stringsAsFactors=F, row.names=1)

ct <- which(rownames(manifest) == '[Controls]')
manifest.ordering <- manifest[1:(ct-1),]
manifest.ordering.cg <- manifest.ordering[grep('^cg', manifest.ordering$Name),]
manifest.ordering.ch <- manifest.ordering[grep('^ch', manifest.ordering$Name),]
manifest.ordering.rs <- manifest.ordering[grep('^rs', manifest.ordering$Name),]
context35 <- sapply(gregexpr('CG', paste(substr(manifest.ordering.cg$Forward_Sequence, 26, 60),
                                         substr(manifest.ordering.cg$Forward_Sequence, 65, 99))), length)

getup <- function(x.order) {
  upstream <- DNAStringSet(paste0(substr(x.order$Forward_Sequence, 1, 60),substr(x.order$Forward_Sequence, 62, 63)))
  dwstream <- DNAStringSet(paste0(substr(x.order$Forward_Sequence, 62, 63),substr(x.order$Forward_Sequence, 65, 124)))
  upstream.r <- reverseComplement(upstream)
  dwstream.r <- reverseComplement(dwstream)
  up <- substr(upstream, 13,62)
  up.r <- substr(upstream.r, 1,50)
  dw <- substr(dwstream,1,50)
  dw.r <- substr(dwstream.r,13,62)
  sourceseq <- DNAStringSet(x.order$SourceSeq)
  isup <- (sourceseq==up | sourceseq==up.r)
  isdw <- (sourceseq==dw | sourceseq==dw.r)
  sourceseq <- DNAStringSet(ifelse(isup, substr(upstream,13,62), substr(dwstream.r,13,62)))
  list(up=isup, dw=isdw, sourceseq=sourceseq)
}

isup <- getup(manifest.ordering.cg)
unique(isup$up+isup$dw)
gr <- GRanges(seqnames=Rle(paste0('chr',manifest.ordering.cg$CHR)),
              ranges=IRanges(manifest.ordering.cg$MAPINFO, end=manifest.ordering.cg$MAPINFO+1, names=manifest.ordering.cg$Name),
              strand='*',
              addressA=as.numeric(manifest.ordering.cg$AddressA_ID),
              addressB=as.numeric(manifest.ordering.cg$AddressB_ID),
              channel=sapply(manifest.ordering.cg$Color_Channel, function(x) if (x=='') 'Both' else x),
              designtype=manifest.ordering.cg$Infinium_Design_Type,
              probetype='cg',
              orientation=Rle(as.factor(ifelse(isup$up,'up','down'))),
              sourceseq=isup$sourceseq,
              probecpgcnt=sapply(gregexpr('CG',isup$sourceseq), length),
              context35=context35,
              probeStart=ifelse(isup$up, manifest.ordering.cg$MAPINFO-48, manifest.ordering.cg$MAPINFO),
              probeEnd=ifelse(isup$up, manifest.ordering.cg$MAPINFO+1, manifest.ordering.cg$MAPINFO+49)
              )

context35 <- sapply(gregexpr('CG', paste(substr(manifest.ordering.ch$Forward_Sequence, 26, 60),
                                         substr(manifest.ordering.ch$Forward_Sequence, 65, 99))), length)

grepup <- function(x.order) {
  upstream <- DNAStringSet(paste0(substr(x.order$Forward_Sequence, 1, 60),substr(x.order$Forward_Sequence, 62, 63)))
  dwstream <- DNAStringSet(paste0(substr(x.order$Forward_Sequence, 62, 63),substr(x.order$Forward_Sequence, 65, 124)))
  upstream.r <- reverseComplement(upstream)
  dwstream.r <- reverseComplement(dwstream)
  up <- substr(upstream, 13,62)
  up.r <- substr(upstream.r, 1,50)
  dw <- substr(dwstream,1,50)
  dw.r <- substr(dwstream.r,13,62)
  sourceseq <- substr(x.order$SourceSeq,1,48)
  isup <- sapply(seq_along(sourceseq), function(i) length(grep(sourceseq[i], upstream[i]))+length(grep(sourceseq[i], upstream.r[i])))
  isdw <- sapply(seq_along(sourceseq), function(i) length(grep(sourceseq[i], dwstream[i]))+length(grep(sourceseq[i], dwstream.r[i])))
  sourceseq <- DNAStringSet(ifelse(isup, substr(upstream.r, 13,62), substr(dwstream.r,13,62)))
  list(up=isup, dw=isdw, sourceseq=sourceseq)
}

isup <- grepup(manifest.ordering.ch)
unique(isup$up+isup$dw)
gr.ch <- GRanges(seqnames=Rle(paste0('chr',manifest.ordering.ch$CHR)),
                 ranges=IRanges(manifest.ordering.ch$MAPINFO, end=manifest.ordering.ch$MAPINFO, names=manifest.ordering.ch$Name),
                 strand='*',
                 addressA=as.numeric(manifest.ordering.ch$AddressA_ID),
                 addressB=as.numeric(manifest.ordering.ch$AddressB_ID),
                 channel=sapply(manifest.ordering.ch$Color_Channel, function(x) if (x=='') 'Both' else x),
                 designtype=manifest.ordering.ch$Infinium_Design_Type,
                 probetype='ch',
                 orientation=Rle(as.factor(ifelse(isup$up,'up','down'))),
                 sourceseq=isup$sourceseq,
                 probecpgcnt=sapply(gregexpr('CG',isup$sourceseq), length),
                 context35=context35,
                 probeStart=ifelse(isup$up, manifest.ordering.ch$MAPINFO-48, manifest.ordering.ch$MAPINFO),
                 probeEnd=ifelse(isup$up, manifest.ordering.ch$MAPINFO+1, manifest.ordering.ch$MAPINFO+49)
                 )


library(SNPlocs.Hsapiens.dbSNP.20120608)
snps <- SNPlocs.Hsapiens.dbSNP.20120608
library(SNPlocs.Hsapiens.dbSNP.20110815)
library(SNPlocs.Hsapiens.dbSNP144.GRCh37)
snps <- SNPlocs.Hsapiens.dbSNP144.GRCh37
rsnames <- rownames(manifest.ordering.rs)
names(rsnames) <- rsnames
## replace rs13369115 by rs10155413
rsnames[rsnames == 'rs13369115'] <- 'rs10155413'
snplocs <- snpsById(snps, rsnames, ifnotfound='warning')

gr.snp <- GRanges(seqnames=gsub('ch','chr', as.character(seqnames(snplocs))),
                  ranges=ranges(snplocs),
                  strand='*',
                  addressA=as.numeric(manifest.ordering.rs$AddressA_ID),
                  addressB=as.numeric(manifest.ordering.rs$AddressB_ID),
                  channel=sapply(manifest.ordering.rs$Color_Channel, function(x) if (x=='') 'Both' else x),
                  designtype=manifest.ordering.rs$Infinium_Design_Type,
                  probetype='rs',
                  orientation=NA,
                  sourceseq=DNAStringSet(manifest.ordering.rs$AlleleA_ProbeSeq),
                  probecpgcnt=NA,
                  context35=NA,
                  probeStart=NA,
                  probeEnd=NA
                  )
names(gr.snp) <- names(rsnames)

EPIC.manifest <- c(gr,gr.ch,gr.snp)

EPIC.manifest <- sortSeqlevels(EPIC.manifest)
seqlevels(EPIC.manifest)
EPIC.manifest <- sort(EPIC.manifest)
save(EPIC.manifest, file='data/EPIC.manifest.rda')

setwd('/primary/home/wanding.zhou/tools/sesame/sesame')

EPIC.ordering <- with(manifest.ordering, data.frame(Probe_ID=Name, M=as.numeric(AddressB_ID), U=as.numeric(AddressA_ID), DESIGN=Infinium_Design_Type, COLOR_CHANNEL=Color_Channel, col=ifelse(Color_Channel=="Red","R", ifelse(Color_Channel=='Grn', 'G', NA)), row.names=rownames(manifest.ordering), stringsAsFactors=FALSE))
EPIC.ordering$COLOR_CHANNEL[EPIC.ordering$COLOR_CHANNEL==""] <- 'Both'
EPIC.ordering$COLOR_CHANNEL <- as.factor(EPIC.ordering$COLOR_CHANNEL)
EPIC.ordering$col <- as.factor(EPIC.ordering$col)
save(EPIC.ordering, file='data/EPIC.ordering.rda')

manifest.controls <- manifest[(ct+1):nrow(manifest),]
EPIC.controls <- with(manifest.controls, data.frame(Address=rownames(manifest.controls), Type=Name, Color_Channel=AddressA_ID, Name=AlleleA_ProbeSeq))

save(EPIC.controls, file='data/EPIC.controls.rda')

## zcat hg19_dbsnp.vcf.gz | awk 'NR==FNR{if(match($1,/^\"(rs[^\"]*)\"/,c)){a[c[1]]=1}}NR!=FNR{if($3 in a){print $3,"chr"$1,$2}}' IlluminaEpicManifest.csv - >rs2chr.csv
## zcat hg19

rs2chr <- read.table('tests/data/EPIC.official.manifest/rs2chr.csv',stringsAsFactors=FALSE)
rschr <- rep(NA, nrow(manifest.ordering))
rschr[match(rs2chr$V1, rownames(manifest.ordering))] <- rs2chr$V2
EPIC.hg19.probe2chr <- setNames(with(manifest.ordering, as.factor(ifelse(CHR=="", as.character(rschr), as.character(paste0('chr',CHR))))),rownames(manifest.ordering))
more.missing <- which(is.na(EPIC.hg19.probe2chr))
load('data/hm450.hg19.probe2chr.rda')
EPIC.hg19.probe2chr[more.missing] <- hm450.hg19.probe2chr[names(more.missing)]
any(is.na(EPIC.hg19.probe2chr)) # no NA
save(EPIC.hg19.probe2chr, file='data/EPIC.hg19.probe2chr.rda')

## d <- data.frame(a=as.character(EPIC.hg19.probe2chr[match(names(hm450.hg19.probe2chr),names(EPIC.hg19.probe2chr))]), b=as.character(hm450.hg19.probe2chr),row.names=names(hm450.hg19.probe2chr))

## for 27k
library(methylumi)

hm27.hg19 <- get27k()
hm27.hg19.df <- as.data.frame(hm27.hg19)
hm27.hg19.probe2chr <- setNames(hm27.hg19.df$seqnames, rownames(hm27.hg19.df))

save(hm27.hg19.probe2chr, file='../data/hm27.hg19.probe2chr.rda')
data(hm27.ordering)
save(hm27.ordering, file='~/tools/sesame/sesame/data/hm27.ordering.rda')
data(hm27.controls)
save(hm27.controls, file='~/tools/sesame/sesame/data/hm27.controls.rda')

load('tests//data/tcga.27k.random6/probesToMask.27k.rda')
hm27.mask <- rownames(as.data.frame(toMask))
save(hm27.mask, file='../data/hm27.mask.rda')
