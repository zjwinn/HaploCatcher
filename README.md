<img src="https://raw.githubusercontent.com/zjwinn/HaploCatcher/main/HaploCatcher_Logo.png" title="HaploCatcher: A Predictive Haplotyping Package" width="500" />

# Introduction
This is a public repository for the R package 'HaploCatcher'. This package utilizes genome-wide molecular data, paired with historical locus/QTL/gene information, to train machine learning algorithms to produce predictive haplotype calls for lines which are genotyped with high density molecular platforms, yet not genotyped for the specific locus/QTL/gene in question. Any questions related to development, maintenance, and potential errors can be directed here.

# Installation
To install HaploCatcher in R, users can use the following fucntion taken from the (devtools package)[https://www.rdocumentation.org/packages/devtools/versions/2.4.5] 
```r
# Install HaploCatcher
devtools::install_github("zjwinn/HaploCatcher")
```
It is recommended to install using the above code to get the most up-to-date version of HaploCatcher. To install from (CRAN)[https://cran.rstudio.com/] instead, use the following code in R:
```r
# Install HaploCatcher from CRAN
install.packages("HaploCatcher")
```

# Getting Started
To understand all necesary inputs for HaploCatcher, users can follow the tutorial laid out in the following R vignette titled, "(An Intro to HaploCatcher)[https://cran.r-project.org/web//packages//HaploCatcher/vignettes/An_Intro_to_HaploCatcher.html]." To call on the vignette directly from the R terminal, use the following R code after installing HaploCatcher:
```r
# Call vignette
vignette("An_Intro_to_HaploCatcher")
```
