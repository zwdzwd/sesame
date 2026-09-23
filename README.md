# SeSAMe - SEnsible Step-wise Analysis of Methylation data 

[![last commit](https://img.shields.io/github/last-commit/zwdzwd/sesame.svg?style=flat-square)](https://github.com/zwdzwd/sesame/commits/master)
                  
SeSAMe is an R package for processing Infinium DNA methylation data. SeSAMe currently supports EPICv2, EPIC, HM450, HM27, MSA, MM285 and Mammal40 platforms and dynamically generated manifest.

To install from Github,
```R
BiocManager::install("zwdzwd/sesame")
```

See the package [Home Page on Bioconductor](https://bioconductor.org/packages/release/bioc/html/sesame.html) and the [Developmental Branch](https://bioconductor.org/packages/devel/bioc/html/sesame.html).

It also has a depended [data package](https://github.com/zwdzwd/sesameData) for annotation and example data.

## Two implementations, one method

SeSAMe comes as two implementations that share the same methods and are
developed together. They are **parallel, not sequential** — neither is an
upgrade of the other, and the version lines say which is which:

- **sesame** (this package) — **SeSAMe v1**, on the **1.x** series. The
  reference implementation of the SeSAMe methods and the interactive R toolkit
  for exploratory analysis and KYCG enrichment. Actively maintained on
  Bioconductor, and the oracle the C version is validated against. This package
  stays on 1.x.

  ```R
  BiocManager::install("sesame")
  ```

- **[SeSAMe2](https://zwdzwd.github.io/sesame-cli/)** — the **second
  implementation**, on the **2.x** series: a standalone C program (no R, no
  Bioconductor, no network) covering the full workflow — IDAT → betas → QC →
  differential methylation → copy number → SNP genotyping, with visualization
  through [cinderplot](https://github.com/zhou-lab/cinderplot). Built for large
  cohorts, pipeline deployment, and non-R environments.

  ```sh
  conda install -c zhou-lab -c conda-forge sesame yame
  ```

Both installables are called `sesame`, because they are one project and one
method. They never meet in a resolver: `BiocManager` reads Bioconductor's
index, conda reads the zhou-lab channel. So `sesame 2.0.0` from conda is not a
newer version of `sesame 1.31.x` from Bioconductor — it is the other
implementation.

R is SeSAMe2's permanent oracle, and the agreement is gated on every release.
Raw betas (`prep=""`) are bit-identical; across the full default pipeline,
probes with identical raw input agree to a median of ~6e-6 and ~1e-3 worst
case. Every intentional numerical divergence is recorded in SeSAMe2's
[NUMERICS.md](https://github.com/zwdzwd/sesame-cli/blob/main/NUMERICS.md), and
the few places where the two answer genuinely different questions — all of them
annotation-lineage differences, none a tolerance — in
[DIVERGENCES.md](https://github.com/zwdzwd/sesame-cli/blob/main/DIVERGENCES.md).

## Bugs
    
Bug reports are appreciated. Register issues at the SeSAMe [issue tracker](http://github.com/zwdzwd/sesame/issues).
    
    
## About

Please cite and reference [SeSAMe: reducing artifactual detection of DNA methylation by Infinium BeadChips in genomic deletions](https://doi.org/10.1093/nar/gky691) for more details.
