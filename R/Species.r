
#' Infer Species
#'
#' @param sset a \code{SigSet}
#' @param df_as a data.frame of alignment score for each probe.
#' @param hm.pos.por whether to use the postive probes portion to infer human.
#' @param topN Top n positive and negative probes used to infer species.
#' @param plot Whether to plot ROC curve.
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
#' df_as=t(read.table(AS_dir,header=T,row.names=1,sep='\t'))
#' overlapped_species=intersect(unique(sample_info$NCBI_scientific_name),colnames(df_as))
#' sample_info=sample_info[sample_info$NCBI_scientific_name %in% overlapped_species,]
#' df_as=df_as[,overlapped_species]
#' df_as=df_as / 50
#' idats=searchIDATprefixes(indir)
#' idats=idats[sample_info$Basename]
#' idats=idats[! is.na(idats)]
#' 
#' mft <- readRDS('/mnt/isilon/zhoulab/labprojects/20201124_Horvath_mammal_array/mammal_array_manifest_original.rds')
#' sset=readIDATpair(idats[1],manifest=mft, platform='mammal')
#' 
#' result=do.call(cbind, mclapply(idats, function(x) {
#'     sset=readIDATpair(x,manifest=mft, platform='mammal')
#'     AUC=inferSpecies(sset,df_as,hm.pos.por=F,topN=2000,plot=F,ret.max=F)
#' },mc.cores=8)
#' )
#' write.table(as.data.frame(result),file="AUC_for_species.txt",sep='\t',quote=FALSE,col.names=T)

#' @export
inferSpecies <- function(sset,df_as=NULL,hm.pos.por=F,topN=2000,plot=F,ret.max=T) {
    pvalue=pval(sset)
    N=length(pvalue)
    pos_probes=sort(pvalue[pvalue< 0.01],decreasing=F)
    n_pos=length(pos_probes)
    if( hm.pos.por & (n_pos / N > 0.8)){
        return("Homo sapiens")
        }
    neg_probes=sort(pvalue[pvalue > 0.2],decreasing=T)
    if(n_pos > topN){
        pos_probes=pos_probes[1:topN]
        }
    if(length(neg_probes) > topN){
        neg_probes=neg_probes[1:topN]
        }
    probes=intersect(row.names(df_as),c(names(pos_probes),names(neg_probes)))
    y_true=sapply(pvalue[probes],function(x) {if (x<0.01) {return(1)} else {return(0)}})
    #y_true=c(rep(1,length(pos_probes)),rep(0,length(neg_probes)))
    if(length(y_true) == 0){
        if (ret.max==T){
            return(NA)
            }
        else {
            return(rep(NA,dim(df_as)[2]))
            }
        }
        
    AUC=sapply(as.data.frame(df_as[probes,]),function(x) {
            if(length(x) == 0){
                return(NA)
                }
            PRROC_obj <- roc.curve(scores.class0 = x, weights.class0=y_true,curve=plot)
            if (plot){
                plot(PRROC_obj)
                }
            return(PRROC_obj$auc)
        })
    if (ret.max == T){
        return(names(AUC)[which.max(AUC)])
        #What if there are two maximal AUC?
        }
    else {
        return(AUC)
        }
}

