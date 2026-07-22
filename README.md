# SeSAMe - SEnsible Step-wise Analysis of Methylation data 

[![last commit](https://img.shields.io/github/last-commit/zwdzwd/sesame.svg?style=flat-square)](https://github.com/zwdzwd/sesame/commits/master)
                  
SeSAMe is an R package for processing Infinium DNA methylation data. SeSAMe currently supports EPICv2, EPIC, HM450, HM27, MSA, MM285 and Mammal40 platforms and dynamically generated manifest.

To install from Github,
```R
BiocManager::install("zwdzwd/sesame")
```

See the package [Home Page on Bioconductor](https://bioconductor.org/packages/release/bioc/html/sesame.html) and the [Developmental Branch](https://bioconductor.org/packages/devel/bioc/html/sesame.html).

It also has a depended [data package](https://github.com/zwdzwd/sesameData) for annotation and example data.

## Two ways to run SeSAMe

SeSAMe comes as two implementations that share the same methods and are developed
together:

- **sesame** (this package) — the **reference implementation** of the SeSAMe
  methods and the interactive R toolkit for exploratory analysis and KYCG
  enrichment. Actively maintained on Bioconductor, and the oracle the C version
  is validated against.

- **[sesame-cli](https://zwdzwd.github.io/sesame-cli/)** — a standalone C program
  (no R, no Bioconductor) covering the full workflow: IDAT → betas → QC →
  differential methylation → copy number, with visualization through
  [cinderplot](https://github.com/zhou-lab/cinderplot). Built for large cohorts,
  pipeline deployment, and non-R environments. A good fit when you process many
  samples or run in a pipeline.

  ```sh
  conda install -c zhou-lab -c conda-forge sesame-cli
  ```

The two agree closely but are not bit-exact: sesame-cli reimplements the
preprocessing independently, and where the R code has a numerically fragile
construction it is reformulated. On probes with identical raw input, betas agree
to a median of ~6e-6 and ~1e-3 worst case. Every intentional divergence is
justified and recorded in the sesame-cli
[NUMERICS.md](https://github.com/zwdzwd/sesame-cli/blob/main/NUMERICS.md).

## Bugs
    
Bug reports are appreciated. Register issues at the SeSAMe [issue tracker](http://github.com/zwdzwd/sesame/issues).
    
    
## About

Please cite and reference [SeSAMe: reducing artifactual detection of DNA methylation by Infinium BeadChips in genomic deletions](https://doi.org/10.1093/nar/gky691) for more details.
