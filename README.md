SpecClustPack
=============

A bunch of R functions related to spectral clustering.

Installation
===

The **SpecClustPack** package can be installed in R directly from GitHub by 
using devtools. 

```r
library(devtools)
install_github("SpecClustPack", "norbertbin")
```

Simulate from Stochastic Blockmodel
===
```r
blockPMat = matrix(c(.6,.2,.2,.6), nrow=2)
nMembers = c(5,5)

simSBM(blockPMat, nMembers)
```
##  10 x 10 sparse Matrix of class "dsCMatrix"
##                         
##  [1,] . . 1 1 1 . . . . 1
##  [2,] . . 1 1 1 . . 1 . .
##  [3,] 1 1 . . . . . . . .
##  [4,] 1 1 . . 1 . . . . .
##  [5,] 1 1 . 1 . . . . . .
##  [6,] . . . . . . 1 1 . 1
##  [7,] . . . . . 1 . 1 . 1
##  [8,] . 1 . . . 1 1 . . 1
##  [9,] . . . . . . . . . .
## [10,] 1 . . . . 1 1 1 . .
```