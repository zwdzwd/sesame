% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SigDFMethods.R
\name{controls}
\alias{controls}
\title{get the controls attributes}
\usage{
controls(sdf)
}
\arguments{
\item{sdf}{a \code{SigDF}}
}
\value{
the controls data frame
}
\description{
get the controls attributes
}
\examples{
sesameDataCache("EPIC") # if not done yet
sdf = sesameDataGet('EPIC.1.SigDF')
head(controls(sdf))
}
