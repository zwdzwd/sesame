#' Make a simulated SeSAMe data set
#'
#' Constructs a simulated dataset.
#'
#' @param platform optional, HM450, EPIC or HM27
#' @param n optional, 
#' @return Object of class \code{SignalSet}
#' @examples
#' sset <- makeExampleSeSAMeDataSet()
#'
#' @export
makeExampleSeSAMeDataSet <- function(n=1000, platform='HM450') {

    dm.ordering <- get(paste0(platform, '.ordering'))
    sset <- SignalSet$new(platform)
    probes <- rownames(
        dm.ordering[dm.ordering$DESIGN=='I' &
                        dm.ordering$COLOR_CHANNEL=='Grn',])
    mt <- matrix(pmax(rnorm(length(probes)*2, 4000, 200),0), ncol=2)
    rownames(mt) <- probes
    colnames(mt) <- c('M','U')
    sset$IG <- mt
    mt <- matrix(pmax(rnorm(length(probes)*2, 400, 200),0), ncol=2)
    rownames(mt) <- probes
    colnames(mt) <- c('M','U')
    sset$oobR <- mt

    probes <- rownames(dm.ordering[(
        dm.ordering$DESIGN=='I' & dm.ordering$COLOR_CHANNEL=='Red'),])
    mt <- matrix(pmax(rnorm(length(probes)*2, 4000, 200),0), ncol=2)
    rownames(mt) <- probes
    colnames(mt) <- c('M','U')
    sset$IR <- mt
    mt <- matrix(pmax(rnorm(length(probes)*2, 400, 200),0), ncol=2)
    rownames(mt) <- probes
    colnames(mt) <- c('M','U')
    sset$oobG <- mt

    probes <- rownames(dm.ordering[dm.ordering$DESIGN=='II',])
    mt <- matrix(pmax(rnorm(length(probes)*2, 4000, 200),0), ncol=2)
    rownames(mt) <- probes
    colnames(mt) <- c('M','U')
    sset$II <- mt

    dm.controls <- get(paste0(platform, '.controls'))
    ctl <- as.data.frame(matrix(pmax(rnorm(
        2*nrow(dm.controls), 400, 300),0), ncol=2))
    rownames(ctl) <- make.names(dm.controls$Name,unique=TRUE)
    ctl <- cbind(ctl, dm.controls[, c("Color_Channel","Type")])
    colnames(ctl) <- c('G','R','col','type')
    sset$ctl <- ctl

    sset
}


#' Make a tiny toy simulated EPIC data set
#'
#' @return Object of class \code{SignalSet}
#' @examples
#' sset <- makeExampleTinyEPICDataSet()
#'
#' @export
makeExampleTinyEPICDataSet <- function() {
    
    platform <- 'EPIC'
    sset <- SignalSet$new(platform)
    probes <- c(
        "cg18478105", "cg01763666", "cg25813447",
        "cg07779434", "cg13417420", "cg24133276")
    mt <- matrix(as.integer(
        pmax(rnorm(length(probes)*2, 4000, 200),0)), ncol=2)
    rownames(mt) <- probes
    colnames(mt) <- c('M','U')
    sset$IG <- mt
    mt <- matrix(as.integer(
        pmax(rnorm(length(probes)*2, 400, 200),0)), ncol=2)
    rownames(mt) <- probes
    colnames(mt) <- c('M','U')
    sset$oobR <- mt
    
    probes <- c(
        "cg16619049", "cg01782097", "cg12712429",
        "cg24373735", "cg18865112", "cg22226438")
    mt <- matrix(as.integer(
        pmax(rnorm(length(probes)*2, 4000, 200),0)), ncol=2)
    rownames(mt) <- probes
    colnames(mt) <- c('M','U')
    sset$IR <- mt
    mt <- matrix(as.integer(
        pmax(rnorm(length(probes)*2, 400, 200),0)), ncol=2)
    rownames(mt) <- probes
    colnames(mt) <- c('M','U')
    sset$oobG <- mt
    
    probes <- c(
        "cg07881041", "cg23229610", "cg03513874",
        "cg05451842", "cg14797042", "cg09838562")
    mt <- matrix(as.integer(
        pmax(rnorm(length(probes)*2, 4000, 200),0)), ncol=2)
    rownames(mt) <- probes
    colnames(mt) <- c('M','U')
    sset$II <- mt
    
    dm.controls <- get(paste0(platform, '.controls'))
    ctl <- as.data.frame(matrix(as.integer(
        pmax(rnorm(2*nrow(dm.controls), 400, 300),0)), ncol=2))
    rownames(ctl) <- make.names(dm.controls$Name,unique=TRUE)
    ctl <- cbind(ctl, dm.controls[, c("Color_Channel","Type")])
    colnames(ctl) <- c('G','R','col','type')
    sset$ctl <- ctl
    
    sset
}
