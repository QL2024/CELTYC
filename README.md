---
title: "CELTYC"
author:
- name: "Qi Luo and Andrew E. Teschendorff"
  affiliation: 
  - CAS Key Lab of Computational Biology, SINH
date: "2024-04-02"
package: CELTYC
output:
  BiocStyle::html_document:
    toc_float: true
---

# Summary

CELTYC is an R package that implements cell type specific clustering pipeline given a DNA methylation data matrix. It includes the process of estimating cell type fractions within a tissue type, identifying disease-associated DMCs altered in specific cell types, and unsupervised clustering with the identified DMCTs for different cell types.

# Installation

To install:

```r
library(devtools)
devtools::install_github("QL2024/CELTYC")
```
