
#' Infer Species
#'
#' @param sset a \code{SigSet}
#' @param df_as a data.frame of alignment score for each probe.
#' @param hm.pos.por whether to use the postive probes portion to infer human.
#' @param topN Top n positive and negative probes used to infer species.
#' @return species (NCBI official species names)
#' We infer species based on probes pvalues and alignment score.
#' We calculate AUC for each species, y_true is whether pval < 0.01
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
#' overlapped_species=intersect(unique(IDs),colnames(df_as))
#' sample_info=sample_info[sample_info$IDs %in% overlapped_species,]
#' df_as=df_as[,overlapped_species]
#' df_as=df_as / 50
#' idats=searchIDATprefixes(indir)
#' idats=idats[sample_info$Basename]
#' idats=idats[! is.na(idats)]
#' 
#' mft <- readRDS('/mnt/isilon/zhoulab/labprojects/20201124_Horvath_mammal_array/mammal_array_manifest_original.rds')
#' sset=readIDATpair(idats[1],manifest=mft, platform='mammal')
#' 
#' @export
inferSpecies <- function(sset,df_as=NULL,hm.pos.por=T,topN=3000,ret.max=T) {
    pvalue=pval(sset)
    pvalue <- pvalue[intersect(names(pvalue),rownames(df_as))]
    N=length(pvalue)
    pos_probes=sort(pvalue[pvalue< 0.01],decreasing=F)
    neg_probes=sort(pvalue[pvalue > 0.2],decreasing=T)
    n_pos=length(pos_probes)
    n_neg=length(neg_probes)
    if (hm.pos.por & (n_neg / N < 0.1)){
        return('Homo sapiens|9606')
        }
    if (n_pos > topN){
        pos_probes=pos_probes[1:topN]
        }
    if (n_neg > topN){
        neg_probes=neg_probes[1:topN]
        }
    probes=c(names(pos_probes),names(neg_probes))
    y_true=sapply(pvalue[probes],function(x) {if (x<0.01) {return(1)} else {return(0)}})
    #y_true=c(rep(1,length(pos_probes)),rep(0,length(neg_probes)))
    if (length(y_true) == 0){
        if (ret.max==T){
            return(NA)
            }
        else {
            return(rep(NA,dim(df_as)[2]))
            }
        }
        
    
    res = as.data.frame(t(sapply(colnames(df_as),function(species) {
            x=df_as[probes,species]
            if(length(x) == 0){
                return(list(auc=NA,p.value=1))
                }
            PRROC_obj <- roc.curve(scores.class0 = x, weights.class0=y_true,curve=F)
            p <- wilcox.test(x[names(pos_probes)],x[names(neg_probes)],alternative='greater')
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
