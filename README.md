<p align="center">
  <img src="https://raw.githubusercontent.com/zjwinn/HaploCatcher/main/man/figures/logo.png" title="HaploCatcher: A Predictive Haplotyping Package" width="200" alt="HaploCatcher hex logo" />
</p>

# HaploCatcher

*A Predictive Haplotyping Package*

## Introduction
This is a public repository for the R package 'HaploCatcher'. This package utilizes genome-wide molecular data, paired with historical locus/QTL/gene information, to train machine learning algorithms to produce predictive haplotype calls for lines which are genotyped with high density molecular platforms, yet not genotyped for the specific locus/QTL/gene in question. Any questions related to development, maintenance, and potential errors can be directed here.

## Installation
To install HaploCatcher in R, users can use the following function taken from the [devtools package](https://www.rdocumentation.org/packages/devtools/versions/2.4.5):
```r
# Install HaploCatcher
devtools::install_github("zjwinn/HaploCatcher")
```
It is recommended to install using the above code to get the most up-to-date version of HaploCatcher. To install from [CRAN](https://CRAN.R-project.org/package=HaploCatcher) instead, use the following code in R:
```r
# Install HaploCatcher from CRAN
install.packages("HaploCatcher")
```

## Getting Started
To understand all necessary inputs for HaploCatcher, users can follow the tutorial laid out in the following R vignette titled, ["An Intro to HaploCatcher."](https://cran.r-project.org/web//packages//HaploCatcher/vignettes/An_Intro_to_HaploCatcher.html) To call on the vignette directly from the R terminal, use the following R code after installing HaploCatcher:
```r
# Call vignette
vignette("An_Intro_to_HaploCatcher")
```

## News
See the [release notes](https://github.com/zjwinn/HaploCatcher/releases) for a full history of changes. The current release, **HaploCatcher 2.0.1**, is a major restructure for speed and maintainability that also adds user-definable case labels (`het_label` / `neg_label`) so you are no longer locked into the `gene` / `het_gene` / `non_gene` naming convention.
