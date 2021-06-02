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
#' sset <- sesameDataGet('EPIC.1.LNCaP')$sset
#' df_as <- sesameDataGet('Mammal40.alignmentScore')
#' res=inferSpecies(sset,df_as)
#' res=inferSpecies(sset,df_as,ret.max=F)
#' result=do.call(cbind, mclapply(idats, function(x) {
#' 		sset=readIDATpair(x,manifest=mft, platform='mammal')
#'      	res=inferSpecies(sset,df_as,ret.max=T)
#'    },mc.cores=8)
#'    )
#' @export

inferSpecies <- function(sset,df_as=NULL,topN=3000,
			threshold.pos=0.01,threshold.neg=0.1,ret.max=T,
			balance=T,threshold.sucess.rate = 0.8) {
    if (is.null(df_as)) {
    	df_as=sesameDataGet(paste(sset@platform,'alignmentScore',sep='.'))
	}
    pvalue=pval(sset)
    pvalue <- pvalue[intersect(names(pvalue),rownames(df_as))]
    pos_probes=sort(pvalue[pvalue <= threshold.pos],decreasing=F)
    neg_probes=sort(pvalue[pvalue >= threshold.neg],decreasing=T)
	success.rate = length(pvalue[pvalue<=0.05]) / length(pvalue)
	if (success.rate >= threshold.sucess.rate && sset@platform=='Mouse') {
		return(list(auc=1,taxid='10090',species='Mus musculus'))
	}
    if (balance) {
	    topN <- min(length(neg_probes),length(pos_probes))
    }
    if (length(pos_probes) > topN){
	    pos_probes <- pos_probes[1:topN]
    }
    if (length(neg_probes) > topN){
	    neg_probes <- neg_probes[1:topN]
    }
    y_true=structure(c(rep(1,length(pos_probes)),rep(0,length(neg_probes))),
		     names=c(names(pos_probes),names(neg_probes)))
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
        
    
    auc = sapply(colnames(df_as),function(s) {
		labels <- as.logical(y_true)
	    n1 <- sum(labels)
	    n2 <- sum(!labels)
	    R1 <- sum(rank(df_as[,s])[labels])
	    U1 <- R1 - n1 * (n1 + 1)/2
	    auc <- U1/(n1 * n2)
	    #z <- (R1 - n1*(n1+n2+1)/2) / sqrt(n1*n2*(n1+n2+1))
	    #p1 <- pnorm(z, lower.tail=FALSE)
	    #p <- wilcox.test(df_as[names(pos_probes),s],df_as[names(neg_probes),s],alternative='greater')
	    return(auc) #p.value=p$p.value
    })
        
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
