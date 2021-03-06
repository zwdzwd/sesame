% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/visualize.R
\name{visualizeRegion}
\alias{visualizeRegion}
\title{Visualize Region}
\usage{
visualizeRegion(
  chrm,
  plt.beg,
  plt.end,
  betas,
  platform = c("EPIC", "HM450", "MM285"),
  refversion = c("hg38", "hg19", "mm10"),
  sample.name.fontsize = 10,
  heat.height = NULL,
  draw = TRUE,
  show.sampleNames = TRUE,
  show.samples.n = NULL,
  show.probeNames = TRUE,
  cluster.samples = FALSE,
  nprobes.max = 1000,
  na.rm = FALSE,
  dmin = 0,
  dmax = 1
)
}
\arguments{
\item{chrm}{chromosome}

\item{plt.beg}{begin of the region}

\item{plt.end}{end of the region}

\item{betas}{beta value matrix (row: probes, column: samples)}

\item{platform}{EPIC, HM450, or MM285}

\item{refversion}{hg38, hg19, or mm10}

\item{sample.name.fontsize}{sample name font size}

\item{heat.height}{heatmap height (auto inferred based on rows)}

\item{draw}{draw figure or return betas}

\item{show.sampleNames}{whether to show sample names}

\item{show.samples.n}{number of samples to show (default: all)}

\item{show.probeNames}{whether to show probe names}

\item{cluster.samples}{whether to cluster samples}

\item{nprobes.max}{maximum number of probes to plot}

\item{na.rm}{remove probes with all NA.}

\item{dmin}{data min}

\item{dmax}{data max}
}
\value{
graphics or a matrix containing the captured beta values
}
\description{
The function takes a genomic coordinate (chromosome, start and end) and a
beta value matrix (probes on the row and samples on the column). It plots
the beta values as a heatmap for all probes falling into the genomic region.
If `draw=TRUE` the function returns the plotted grid graphics object.
Otherwise, the selected beta value matrix is returned.
`cluster.samples=TRUE/FALSE` controls whether hierarchical clustering is
applied to the subset beta value matrix.
}
\examples{
betas <- sesameDataGet('HM450.76.TCGA.matched')$betas
visualizeRegion('chr20', 44648623, 44652152, betas, 'HM450')
}
