SpecClustPack
=============

A bunch of R functions related to spectral clustering.

Installation
===

The *SpecClustPack* package can be installed in R directly from GitHub by 
using devtools. 

```r
library(devtools)
install_github("SpecClustPack", "norbertbin")
```

Function to simulate from Stochastic Blockmodel
===
```r
blockPMat = matrix(c(.6,.3,.3,.6), nrow=2)
nMembers = c(5,5)
simSBM(blockpMat, nMembers)
```