% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sesame.R
\name{bisConversionControl}
\alias{bisConversionControl}
\title{Compute internal bisulfite conversion control}
\usage{
bisConversionControl(sdf)
}
\arguments{
\item{sdf}{a SigDF}
}
\value{
GCT score (the higher, the more incomplete conversion)
}
\description{
Compute GCT score for internal bisulfite conversion control. The function
takes a \code{SigSet} as input. The higher the GCT score, the more likely
the incomplete conversion.
}
\examples{
sesameDataCache("EPIC") # if not done yet
sdf <- sesameDataGet('EPIC.1.SigDF')
bisConversionControl(sdf)

}
