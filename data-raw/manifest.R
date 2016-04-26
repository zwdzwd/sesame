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

manifest <- read.csv('~/projects/shen-secondary/2016_04_20_EPIC_beta_test/IlluminaEpicManifest.csv',header=T,stringsAsFactors=F, row.names=1)

ct <- which(rownames(manifest) == '[Controls]')
manifest.ordering <- manifest[1:(ct-1),]
EPIC.ordering <- with(manifest.ordering, data.frame(Probe_ID=Name, M=AddressB_ID, U=AddressA_ID, DESIGN=Infinium_Design_Type, COLOR_CHANNEL=Color_Channel, col=ifelse(Color_Channel=="Red","R", ifelse(Color_Channel=='Grn', 'G', NA)), row.names=rownames(manifest.ordering), stringsAsFactors=FALSE))
EPIC.ordering$COLOR_CHANNEL[EPIC.ordering$COLOR_CHANNEL==""] <- 'Both'
EPIC.ordering$COLOR_CHANNEL <- as.factor(EPIC.ordering$COLOR_CHANNEL)
EPIC.ordering$col <- as.factor(EPIC.ordering$col)
save(EPIC.ordering, file='~/tools/biscuitr/biscuitr/data/EPIC.ordering.rda')

manifest.controls <- manifest[(ct+1):nrow(manifest),]
EPIC.controls <- with(manifest.controls, data.frame(Address=rownames(manifest.controls), Type=Name, Color_Channel=AddressA_ID, Name=AlleleA_ProbeSeq))
save(EPIC.controls, file='~/tools/biscuitr/biscuitr/data/EPIC.controls.rda')

## zcat hg19_dbsnp.vcf.gz | awk 'NR==FNR{if(match($1,/^\"(rs[^\"]*)\"/,c)){a[c[1]]=1}}NR!=FNR{if($3 in a){print $3,"chr"$1,$2}}' IlluminaEpicManifest.csv - >rs2chr.csv

rs2chr <- read.table('~/projects/shen-secondary/2016_04_20_EPIC_beta_test/rs2chr.csv',stringsAsFactors=FALSE)
rschr <- rep(NA, nrow(manifest.ordering))
rschr[match(rs2chr$V1, rownames(manifest.ordering))] <- rs2chr$V2
EPIC.hg19.probe2chr <- setNames(with(manifest.ordering, as.factor(ifelse(CHR=="", as.character(rschr), as.character(paste0('chr',CHR))))),rownames(manifest.ordering))
more.missing <- which(is.na(EPIC.hg19.probe2chr))
EPIC.hg19.probe2chr[more.missing] <- hm450.hg19.probe2chr[names(more.missing)]
any(is.na(EPIC.hg19.probe2chr))
save(EPIC.hg19.probe2chr, file='~/tools/biscuitr/biscuitr/data/EPIC.hg19.probe2chr.rda')

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
