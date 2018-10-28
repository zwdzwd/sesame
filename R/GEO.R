
#' Parse GEO signal-A/B File into a SigSet List
#' 
#' This function is meant to be a convenience function for 
#' parsing data from Signal_A and Signal_B file provided by
#' GEO. In many cases, this function generates a "partial"
#' SigSet due to lack of out-of-band signal and control 
#' probe measurement in those Signal_A/B files.
#' The detection p-value is based on a fixed normal distribution
#' rather than from negative control or OOB probes.
#' 
#' @param path path to Signal-A/B file downlaoded from GEO. 
#' The file can remain gzipped.
#' @param platform HM450, EPIC or HM27
#' @param drop whether to reduce to SigSet when there is only
#' one sample.
#' @param parallel whether to use multiple cores.
#' @return a SigSetList or a SigSet
#' @examples
#' path = '/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/hm450/GSE36369/GSE36369_NonEBV_SignalA_SignalB_3samples_1k.txt.gz'
#' @export
parseGEOSignalABFile <- function(
    path, platform='HM450', drop=TRUE, parallel=TRUE) {
    
    df <- read.table(gzfile(path), header=TRUE, check.names = FALSE)
    samples <- colnames(df)
    samples <- unique(sub('.Signal_[AB]','',samples[samples != 'TargetID']))
    
    # make sure each sample has both Signal_A and Signal_B
    stopifnot(all(vapply(samples, function(s) exists(
        paste0(s,'.Signal_A'), df), logical(1))))
    stopifnot(all(vapply(samples, function(s) exists(
        paste0(s,'.Signal_B'), df), logical(1))))
    
    ord0 <- sesameDataGet(paste0(platform, '.address'))$ordering
    load.sset <- function(sample) {
        sset <- SigSet(platform)
        ## by GEO convention, A is unmethylated and B is methylated
        ## IG
        ord <- ord0[((ord0$DESIGN=='I')&(ord0$col=='G')),]
        aux <- !is.na(match(df$TargetID, ord$Probe_ID))
        IG(sset) <- as.matrix(data.frame(
            U=df[aux, paste0(sample, '.Signal_B')],
            M=df[aux, paste0(sample, '.Signal_A')]), 
            row.names=df$TargetID[aux])
        ## IR
        ord <- ord0[((ord0$DESIGN=='I')&(ord0$col=='R')),]
        aux <- !is.na(match(df$TargetID, ord$Probe_ID))
        IR(sset) <- as.matrix(data.frame(
            U=df[aux, paste0(sample, '.Signal_B')],
            M=df[aux, paste0(sample, '.Signal_A')]),
            row.names=df$TargetID[aux])
        ## II
        ord <- ord0[ord0$DESIGN=='II',]
        aux <- !is.na(match(df$TargetID, ord$Probe_ID))
        II(sset) <- as.matrix(data.frame(
            U=df[aux, paste0(sample, '.Signal_B')],
            M=df[aux, paste0(sample, '.Signal_A')]),
            row.names=df$TargetID[aux])
        sset <- detectionPfixedNorm(sset)
        sset
    }
    
    if (drop && length(samples)==1) {
        load.sset(samples[1])
    } else if (parallel) {
        SigSetList(mclapply(samples, load.sset, verbose=TRUE))
    } else {
        SigSetList(lapply(samples, load.sset))
    }
}
