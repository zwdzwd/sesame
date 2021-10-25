## TODO need to export this function
parseGEOSignalABFile <- function(
    path, platform='HM450', drop=TRUE, parallel=TRUE,
    id_column=TRUE, format="signalMU") {

    df = read.table(gzfile(path), header=TRUE, check.names = FALSE)
    samples = colnames(df)
    if (id_column) {
        samples = samples[2:length(samples)]
    }
    if (format == "signalAB") {
        s1 = ".Signal_A"
        s2 = ".Signal_B"
    } else if (format == "signalMU") {
        s1 = "meth"
        s2 = "nometh"
        samples = sub("_pval","",samples)
    }
    samples = sub(s1, "", samples)
    samples = sub(s2, "", samples)
    samples = unique(samples)

    # make sure each sample has both Signal_A and Signal_B
    stopifnot(all(vapply(samples, function(s) exists(
        paste0(s,'.Signal_A'), df), logical(1))))
    stopifnot(all(vapply(samples, function(s) exists(
        paste0(s,'.Signal_B'), df), logical(1))))
    
    ord0 = sesameDataGet(paste0(platform, '.address'))$ordering
    sdfs = lapply(samples, function(sample) {
        ## sample = "C_318"
        auxM = df[[paste0(sample,"_",s1)]][match(ord0$Probe_ID, df[[1]])]
        auxU = df[[paste0(sample,"_",s2)]][match(ord0$Probe_ID, df[[1]])]
        ord0$MG = NA
        ord0$MR = NA
        ord0$UG = NA
        ord0$UR = NA
        idx = ord0$DESIGN == "I" & ord0$col == "G"
        ord0$MG[idx] = auxM[idx]
        ord0$UG[idx] = auxU[idx]

        idx = ord0$DESIGN == "I" & ord0$col == "R"
        ord0$MR[idx] = auxM[idx]
        ord0$MU[idx] = auxU[idx]

        idx = ord0$DESIGN == "II"
        ord0$UG[idx] = auxM[idx]
        ord0$UR[idx] = auxU[idx]

        sdf = SigDF(ord0, platform)
        auxP = df[[paste0(sample,"_pval")]][match(ord0$Probe_ID, df[[1]])]
        sdf$mask[auxP >= 0.05] = TRUE
        sdf
    })
    names(sdfs) = samples
    sdfs
}

## path = "~/Downloads/GSE57362_signal_intensities.txt.gz"
## saveRDS(sdfs, file="~/Dropbox/tmp/GSE57362_SigDFs.rds")

