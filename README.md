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

adjMat = simSBM(blockPMat, nMembers)
adjMat
```

```
##  10 x 10 sparse Matrix of class "dsCMatrix"
##                         
##  [1,] . . 1 1 1 . . . . 1
##  [2,] . . 1 1 1 . . 1 . .
##  [3,] 1 1 . . . . . . . .
##  [4,] 1 1 . . 1 . . . 1 .
##  [5,] 1 1 . 1 . . . . . .
##  [6,] . . . . . . 1 1 . 1
##  [7,] . . . . . 1 . 1 . 1
##  [8,] . 1 . . . 1 1 . . 1
##  [9,] . . . 1 . . . . . .
## [10,] 1 . . . . 1 1 1 . .
```

Plot the SBM Probability Matrix
===
```r
plotSBM(blockPMat, nMembers)
```

Plot the Simulated Adjacency Matrix
===
```r
plotAdj(adjMat)
```

Run Spectral Clustering
===
By default, the **specClust** function uses regularized spectral 
clustering with row normalization, but can be adjusted by changing 
the *method* and *rowNorm* parameters. 
```r
clusters = specClust(adjMat, 2)
clusters
```

```
## [1] 1 1 1 1 1 2 2 2 1 2
```

Compute the Mis-clustering Rate
===
The function **misClustRate** computes the proportion of mis-clustered nodes
(up to identifiability) given the cluster sizes.
```r
misClustRate(clusters, nMembers)
```

```
## [1] 0.1
```

Estimate SBM Probabilities
===
The function **estSBM** estimates the block probability matrix given the
adjacency matrix and the cluster assignments.
```r
estSBM(adjMat, clusters)
```

```
##            [,1]       [,2]
## [1,] 1.00000000 0.08333333
## [2,] 0.08333333 0.53333333
```