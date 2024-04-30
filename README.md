<img src="https://raw.githubusercontent.com/zjwinn/HaploCatcher/main/HaploCatcher_Logo.png" title="HaploCatcher: A Predictive Haplotyping Package" width="500" />

# Introduction
This is a public repository for the R package 'HaploCatcher'. This package utilizes genome-wide molecular data, paired with historical locus/QTL/gene information, to train machine learning algorithms to produce predictive haplotype calls for lines which are genotyped with high density molecular platforms, yet not genotyped for the specific locus/QTL/gene in question. Any questions related to development, maintenance, and potential errors can be directed to: Dr. Zachary J. Winn [zwinn@outlook.com].

# Installation
To install HaploCatcher in R, users can use the following fucntion taken from the (devtools package)[https://www.rdocumentation.org/packages/devtools/versions/2.4.5] 

```r
#install HaploCatcher
devtools::install_github("zjwinn/HaploCatcher")
```
