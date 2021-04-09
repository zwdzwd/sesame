#' Infer Species
#'
#' @param sset a \code{SigSet}
#' @param df_as a data.frame of alignment score for each probe.
#' @param infer.human whether to use the fraction of negative probes to infer human.
#' @param topN Top n positive and negative probes used to infer species.
#' @return a list of auc, pvalue, species (NCBI official species names) and taxid.
#' We infer species based on probes pvalues and alignment score.
#' AUC was calculated for each specie, y_true is 1 or 0 
#' for pval < 0.01 or pval > 0.2, respeceively,
#' y_pred is normalized alignment score (AS / 50).
#'
#' Our function works on a single sample.
#' @importFrom PRROC roc.curve
#' @import sesameData
#' @examples
#' sset <- sesameDataGet('EPIC.1.LNCaP')$sset
#' df_as <- sesameDataGet() #TODO
#' inferSpecies(sset)
#' 
#' @examples
#' library(tidyverse)
#' library(sesame)
#' library(parallel)
#' library(readxl)
#' indir="/mnt/isilon/zhou_lab/projects/20191212_GEO_datasets/IDATs/Mammal40/"
#' AS_dir="/mnt/isilon/zhoulab/labprojects/20200228_Mouse_Array_Project/20210324_SeSAme_Wubin/data/AS_A_probe_species_matrix.txt.gz"
#' samplesheet="/mnt/isilon/zhoulab/labsamplesheets/2021/20210325_Mammal40_Merged_GEO_Steve_Wubin_samplesheets.xlsx"
#' 
#' sample_info=read_excel(samplesheet)
#' df_as=read.table(AS_dir,header=T,sep='\t')
#' rownames(df_as) <- do.call(paste,c(df_as[c("Species","taxid")],sep="|"))
#' df_as <- t(df_as[,!names(df_as) %in% c("Species","taxid")])
#' df_as[1:3,1:3]
#' 
#' sample_info$IDs <- do.call(paste,c(sample_info[c("NCBI_scientific_name","taxid")],sep="|"))
#' overlapped_species=intersect(unique(sample_info$IDs),colnames(df_as))
#' sample_info=sample_info[sample_info$IDs %in% overlapped_species,]
#' df_as=df_as[,overlapped_species]
#' df_as=df_as / 50
#' idats=searchIDATprefixes(indir)
#' idats=idats[sample_info$Basename]
#' idats=idats[! is.na(idats)]
#' 
#' mft <- readRDS('/mnt/isilon/zhoulab/labprojects/20201124_Horvath_mammal_array/mammal_array_manifest_original.rds')
#' sset=readIDATpair(idats[1],manifest=mft, platform='mammal')
#' # res=inferSpecies(sset,df_as)
#' res=inferSpecies(sset,df_as,ret.max=F)
#' result=do.call(cbind, mclapply(idats, function(x) {
#' 		sset=readIDATpair(x,manifest=mft, platform='mammal')
#'      res=inferSpecies(sset,df_as,ret.max=T)
#'    },mc.cores=8)
#'    )

#' result=do.call(cbind, mclapply(idats, function(x) {
#'       sset=readIDATpair(x,manifest=mft, platform='mammal')
#'       res=inferSpecies(sset,df_as,ret.max=F)
#'       return(res[,'auc'])
#'    },mc.cores=8)
#'    )
#'
#' write.table(as.data.frame(result),file="AUC_for_species.txt",sep='\t',quote=FALSE,col.names=T)
#' @export

inferSpecies <- function(sset,df_as=NULL,infer.human=T,topN=3000,ret.max=T) {
    pvalue=pval(sset)
    pvalue <- pvalue[intersect(names(pvalue),rownames(df_as))]
    pos_probes=sort(pvalue[pvalue< 0.01],decreasing=F)
    neg_probes=sort(pvalue[pvalue > 0.2],decreasing=T)
    n_neg=length(neg_probes)
    if (infer.human & (n_neg / length(pvalue) < 0.1)){
        return('Homo sapiens|9606')
        }
    if (length(pos_probes) > topN){
        pos_probes <- pos_probes[1:topN]
        }
    if (n_neg > topN){
        neg_probes <- neg_probes[1:topN]
        }
    #probes <- c(names(pos_probes),names(neg_probes))
    #y_true <- sapply(pvalue[probes],function(x) {if (x<0.01) {return(1)} else {return(0)}})
    y_true=structure(c(rep(1,length(pos_probes)),rep(0,length(neg_probes))),
						names=c(names(pos_probes),names(neg_probes)))
	df_as=df_as[c(names(pos_probes),names(neg_probes)),]
    if (length(y_true) == 0){
        if (ret.max==T){
            return(list(auc=NA,p.value=NA,species=NA,taxid=NA))
            }
        else {
			return(as.data.frame(t(sapply(colnames(df_as),function(species) {
					return(list(auc=NA,p.value=NA,species=NA,taxid=NA))
				}))))
            }
        }
        
    
    res = as.data.frame(t(sapply(colnames(df_as),function(s) {
			# x=df_as[,s]
			# if(length(x) == 0){
				# return(list(auc=NA,p.value=1,species=NA,taxid=NA))
				# }
			PRROC_obj <- roc.curve(scores.class0 = df_as[,s], weights.class0=y_true,curve=F)
			p <- wilcox.test(df_as[names(pos_probes),s],df_as[names(neg_probes),s],alternative='greater')
			return(list(auc=PRROC_obj$auc,p.value=p$p.value))
		})))
        

    if (ret.max == T){
		r=t(res[which.max(res$auc),])
		l=r[,1]
		l['species']=unlist(strsplit(colnames(r)[1],"\\|"))[1]
		l['taxid']=unlist(strsplit(colnames(r)[1],"\\|"))[2]
		return(l)
        }
    else {
		return(res)
        }
}
