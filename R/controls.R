#' Build control summary matrix
#'
#' The function takes a \code{SigSet} as input and outputs the control matrix
#' summary vector. This vector summarizes one single QC metric for the array
#' control. This includes bisulfite control, stain signal extension efficiency
#' and more.
#' 
#' @param sset an object of class SigSet
#' @return a vector with control summaries
#' @examples
#' sset <- makeExampleSeSAMeDataSet()
#' control.summary <- buildControlMatrix450k(sset)
#' 
#' @export
buildControlMatrix450k <- function(sset) {
    
    ctls <- split(sset@ctl, sset@ctl$type)

    cm <- NULL

    ## bisulfite conversion type II
    cm <- c(
        cm, bisulfite2=mean(ctls[['BISULFITE CONVERSION II']]$R, na.rm=TRUE))

    ## bisulfite conversion type I
    cm <- c(
        cm,
        bisulfite1=mean(
        ctls[['BISULFITE CONVERSION I']][sprintf(
            'BS.Conversion.I.C%s', seq_len(3)),'G'] +
                ctls[['BISULFITE CONVERSION I']][sprintf(
                    'BS.Conversion.I.C%s', 4:6),'R'],
        na.rm=TRUE))

    ## staining
    cm <- c(
        cm, stain.red=ctls[['STAINING']]['DNP..High.', 'R'],
        stain.green=ctls[['STAINING']]['Biotin..High.','G'])

    ## extension
    cm <- c(
        cm,
        extRed1=ctls[['EXTENSION']]['Extension..A.','R'],
        extRed2=ctls[['EXTENSION']]['Extension..T.','R'],
        extGrn1=ctls[['EXTENSION']]['Extension..C.','G'],
        extGrn2=ctls[['EXTENSION']]['Extension..G.','G'])

    ## hybridization
    d <- ctls[['HYBRIDIZATION']]$G
    cm <- c(cm, setNames(d, paste0('hybe',seq_along(d))))

    ## target removal
    d <- ctls[['TARGET REMOVAL']]$G
    cm <- c(cm, setNames(d, paste0('targetrem',seq_along(d))))

    ## non-polymorphic
    cm <- c(
        cm,
        nonpolyRed1=ctls[['NON-POLYMORPHIC']]['NP..A.','R'],
        nonpolyRed2=ctls[['NON-POLYMORPHIC']]['NP..T.','R'],
        nonpolyGrn1=ctls[['NON-POLYMORPHIC']]['NP..C.','G'],
        nonpolyGrn2=ctls[['NON-POLYMORPHIC']]['NP..G.','G'])

    ## specificity type II
    d <- ctls[['SPECIFICITY II']]
    cm <- c(
        cm,
        structure(d$G, names=paste0('spec2Grn', seq_len(dim(d)[1]))),
        structure(d$R, names=paste0('spec2Red', seq_len(dim(d)[1]))))
    cm <- c(cm, spec2.ratio = mean(d$G,na.rm=TRUE) / mean(d$R,na.rm=TRUE))

    ## specificity type I green
    d <- ctls[['SPECIFICITY I']][sprintf('GT.Mismatch.%s..PM.',seq_len(3)),]
    cm <- c(cm, structure(d$G, names=paste0('spec1Grn',seq_len(dim(d)[1]))))
    cm <- c(cm, spec1.ratio1 = mean(d$R, na.rm=TRUE)/mean(d$G, na.rm=TRUE))

    ## specificity type I red
    d <- ctls[['SPECIFICITY I']][sprintf('GT.Mismatch.%s..PM.',4:6),]
    cm <- c(cm, structure(d$R, names=paste0('spec1Red',seq_len(dim(d)[1]))))
    cm <- c(cm, spec1.ratio2 = mean(d$G, na.rm=TRUE)/mean(d$R, na.rm=TRUE))

    ## average specificity ratio
    cm <- c(cm, spec1.ratio = unname(
                (cm['spec1.ratio1']+cm['spec1.ratio2'])/2.0))

    ## normalization
    cm <- c(cm, c(
        normA=mean(ctls[['NORM_A']]$R, na.rm=TRUE), 
        normT=mean(ctls[['NORM_T']]$R, na.rm=TRUE), 
        normC=mean(ctls[['NORM_C']]$G, na.rm=TRUE), 
        normG=mean(ctls[['NORM_G']]$G, na.rm=TRUE)))

    cm <- c(cm, dyebias=(cm['normC']+cm['normG']) / (cm['normA']+cm['normT']))

    ## out-of-band probe quantiles
    if (is.null(sset@oobG)) {
        cm <- c(
            cm,
            oob.ratio=NA,
            structure(rep(NA,3), names=paste0('oob',c(1,50,99))))
    } else {
        cm <- c(cm, oob.ratio=median(sset@oobG) / median(sset@oobR))
        cm <- c(cm, structure(
            quantile(sset@oobG, c(0.01,0.5,0.99)),
            names=paste0('oob', c(1,50,99))))
    }

    cm
}
