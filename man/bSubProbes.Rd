% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{bSubProbes}
\alias{bSubProbes}
\title{subset beta value matrix by probes}
\usage{
bSubProbes(betas, probes)
}
\arguments{
\item{betas}{beta value matrix}

\item{probes}{probe set}
}
\value{
subsetted beta value matrix
}
\description{
subset beta value matrix by probes
}
\examples{
probes <- getAutosomeProbes('HM450')
betas <- sesameDataGet('HM450.1.TCGA.PAAD')$betas
betas <- bSubProbes(betas, probes)
}
