% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{probeID_designType}
\alias{probeID_designType}
\title{Extract the probe type field from probe ID
This only works with the new probe ID system.
See https://github.com/zhou-lab/InfiniumAnnotation for illustration}
\usage{
probeID_designType(Probe_ID)
}
\arguments{
\item{Probe_ID}{Probe ID}
}
\value{
a vector of '1' and '2' suggesting Infinium-I and Infinium-II
}
\description{
Extract the probe type field from probe ID
This only works with the new probe ID system.
See https://github.com/zhou-lab/InfiniumAnnotation for illustration
}
\examples{
probeID_designType("cg36609548_TC21")
}
