#' Infer Species
#'
#' @param sset a \code{SigSet}
#' @param df_as a data.frame of alignment score for each probe.
#' @param topN Top n positive and negative probes used to infer species.
#' @param threshold.pos pvalue < threshold.pos are considered to be positive probes (default is 0.01).
#' @param threshold.neg pvalue > threshold.neg are considered to be negative probes (default is 0.2).
#' @param ret.max whether to return the species with maximal AUC.
#' @param balance whether to balance the postive and negative probes size (default is TRUE).
#' @param threshold.sucess.rate threshold of success rate to determine mouse species.
#' @return a list of auc, pvalue, species (NCBI official species names) and taxid.
#' We infer species based on probes pvalues and alignment score.
#' AUC was calculated for each specie, y_true is 1 or 0 
#' for pval < threshold.pos or pval > threshold.neg, respeceively,
#'
#' Our function works on a single sample.
#' @examples
#' df_as <- sesameDataGet('Mammal40.alignmentScore')
#' manifests=sesameDataGet("Mammal40.address")
#' mft_mm=cbind(manifests$ordering[,1:3],manifests$species$mus_musculus)
#' idat <- searchIDATprefixes(system.file("extdata/", package = "sesameData"))['GSM4411953']
#' sset <- readIDATpair(idat,manifest=mft_mm, platform='Mammal40')
#' inferred.species <- inferSpecies(sset,df_as)
#' @export

inferSpecies <- function(sset,df_as=NULL,topN=3000,
			threshold.pos=0.01,threshold.neg=0.1,ret.max=T,
			balance=T,threshold.sucess.rate = 0.8,platform='MM285') {
    if (is.null(df_as)) {
	# Load alignment score (df_as) from .rda according to the specific species.
    	df_as=sesameDataGet(paste(platform,'alignmentScore',sep='.'))
	}
    pvalue=pval(sset)
    # Only keep the overlapped probes between matrix of pvalue and alignment score.
    pvalue <- pvalue[intersect(names(pvalue),rownames(df_as))]
    # Get positive probes (pvalue <= 0.01) and sort from the smallest to the largest
    pos_probes=sort(pvalue[pvalue <= threshold.pos],decreasing=F)
    # Get negative probes (pvalue >= 0.1) and sort from the largest to the smallest.
    neg_probes=sort(pvalue[pvalue >= threshold.neg],decreasing=T)
    success.rate = length(pvalue[pvalue<=0.05]) / length(pvalue)
	
    # If success.rate larger than 0.8, then directly assign species 'mouse' to this sample
    if (success.rate >= threshold.sucess.rate && sset@platform=='MM285') {
	return(list(auc=1,taxid='10090',species='Mus musculus'))
	}
	
    # balance means keep the same number of positive and negative probes.
    if (balance) {
	topN <- min(length(neg_probes),length(pos_probes))
    }
    if (length(pos_probes) > topN){
	pos_probes <- pos_probes[1:topN]
    }
    if (length(neg_probes) > topN){
	neg_probes <- neg_probes[1:topN]
    }
	
    # for positive probes (pvalue <= 0.01), y_true <- 1, for negative probes (pvalue > 0.1), y_true <- 0
    y_true=structure(c(rep(1,length(pos_probes)),rep(0,length(neg_probes))),
		     names=c(names(pos_probes),names(neg_probes)))
	
    # y_pred is the alignment score.
    df_as=df_as[c(names(pos_probes),names(neg_probes)),]
	
    if (length(y_true) == 0){
	    if (ret.max==T){
		    return(list(auc=NA,species=NA,taxid=NA))
	    }
	    else {
		    return(as.data.frame(t(sapply(colnames(df_as),function() {
			    list(auc=NA,species=NA,taxid=NA)
		    }))))
	    }
    }
        
    # Calculate AUC based on y_true and y_pred (df_as)
    auc = sapply(colnames(df_as),function(s) {
	    labels <- as.logical(y_true)
	    n1 <- sum(labels)
	    n2 <- sum(!labels)
	    R1 <- sum(rank(df_as[,s])[labels])
	    U1 <- R1 - n1 * (n1 + 1)/2
	    auc <- U1/(n1 * n2)
	    return(auc)
    })
        
    # if ret.max, then only the species with the maximal AUC would be returned, else, a vector of species and AUC would be returned.
    if (ret.max == T){
	    l=as.list(auc[which.max(auc)])
	    l['species']=unlist(strsplit(names(l)[1],"\\|"))[1]
	    l['taxid']=unlist(strsplit(names(l)[1],"\\|"))[2]
		names(l)[1] <- 'auc'
	    return(l)
    }
    else {
	    return(auc)
    }
}
