#' Make a simulated SeSAMe data set
#'
#' Constructs a simulated dataset.
#'
#' @param platform optional, HM450, EPIC or HM27
#' @return Object of class \code{SignalSet}
#' @examples
#' sset <- makeExampleSeSAMeDataSet()
#'
#' @export
makeExampleSeSAMeDataSet <- function(platform='HM450') {

    dm.ordering <- getBuiltInData('ordering', platform)
    sset <- SignalSet$new(platform)
    probes <- rownames(dm.ordering[dm.ordering$DESIGN=='I' & dm.ordering$COLOR_CHANNEL=='Grn',])
    mt <- matrix(pmax(rnorm(length(probes)*2, 4000, 200),0), ncol=2)
    rownames(mt) <- probes
    colnames(mt) <- c('M','U')
    sset$IG <- mt
    mt <- matrix(pmax(rnorm(length(probes)*2, 400, 200),0), ncol=2)
    rownames(mt) <- probes
    colnames(mt) <- c('M','U')
    sset$oobR <- mt

    probes <- rownames(dm.ordering[dm.ordering$DESIGN=='I' & dm.ordering$COLOR_CHANNEL=='Red',])
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

    dm.controls <- getBuiltInData('controls', platform)
    ctl <- as.data.frame(matrix(pmax(rnorm(2*nrow(dm.controls), 400, 300),0), ncol=2))
    rownames(ctl) <- make.names(dm.controls$Name,unique=TRUE)
    ctl <- cbind(ctl, dm.controls[, c("Color_Channel","Type")])
    colnames(ctl) <- c('G','R','col','type')
    sset$ctl <- ctl

    sset
}
