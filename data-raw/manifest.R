library('FDb.InfiniumMethylation.hg19')

hm450.hg19 <- get450k()

hm450.hg19.df <- as.data.frame(hm450.hg19)
hm450.hg19.probe2chr <- setNames(hm450.hg19.df$seqnames, rownames(hm450.hg19.df))

save(hm450.hg19.probe2chr, file='../data/hm450.hg19.probe2chr.rda')

save(hm450.ordering, file='../data/hm450.ordering.rda')
save(hm450.controls, file='../data/hm450.controls.rda')
