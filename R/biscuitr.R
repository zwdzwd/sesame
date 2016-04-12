#' Analyze DNA methylation data
#' 
#' biscuitR allows analysing DNA methylation data in R.
#' 
#' This package complements the biscuit software.
#' @aliases biscuitr
#' @author Wanding Zhou
#' @references To appear
#' @seealso To appear
#' @examples To appear
"_PACKAGE"
#> [1] "_PACKAGE"


#samples <- read.csv("samples1k.csv", stringsAsFactors=FALSE)
#samples <- read.csv("samples.csv", stringsAsFactors=FALSE)
## options("mc.cores"=20)
## options("mc.preschedule"=FALSE)

#' Import all IDATs
#' 
#' Import all IDATs from a directory
#' 
#' @param dir directory name.
#' @author Wanding Zhou
#' @return a data frame with signal intensity at each address
#' 
#' @examples
#' batchReadIDATs(dir.name)
batchReadIDATs <- function() {
  library(illuminaio)
  
  for (j in 0:10) {
    i <- 1
    dms <- do.call('cbind', lapply(samples$barcode[(j*1000+1):min((j+1)*1000,length(samples$barcode))], function(x) {
      i <<- i+1;
      if (i%%100==0) {
        cat(i,"\n");
        gc();
      }
      ida.grn <- readIDAT(paste0(x,"_Grn.idat"))
      ida.red <- readIDAT(paste0(x,"_Red.idat"))
      dm1 <- cbind(cy3=ida.grn$Quants[,"Mean"], cy5=ida.red$Quants[,"Mean"])
      colnames(dm1) <- c(paste0(x,'.cy3'), paste0(x, '.cy5'))
      dm1
    }))
    save(dms, file=paste0("datam.",j,".rda"))
  }
}




