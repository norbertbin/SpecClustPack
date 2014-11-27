#' Performs an approximate spectral clustering on an adjacency matrix by using
#' the Nystrom extension.
#'
#' @param adjMat A graph adjacency matrix.
#' @param nBlocks The number of clusters in the graph.
#' @param rowNorm If true, eigenvector rows should be normalized before
#' running kmeans (default is true).
#' @param nIter Number of kmeans iterations (default is 10).
#' @param verbose If true, return a list with cluster assignments, within
#' cluster sum of squares, and eigendecomposition (default is false). 
#'
#' @export
#' @return A vector of node cluster assignments. Or, if verbose is set to true
#' a list with cluster assignments and additional information.
#'
#' @keywords spectral clustering
partSpecClust <- function(adjMat, nBlocks, subSampleSize, rowNorm = T,
                          nIter = 10, verbose = F) {

    nNodes = dim(adjMat)[1]
    
    # get sorted node degrees with indices 
    deg = rowSums(adjMat)
    sortDeg = sort.int(deg, decreasing = T, index.return = T)
    topDegInd = sortDeg$ix[1:subSampleSize]
    
    # get inverse index mapping to re-order final results
    fullMapping = c(topDegInd, (1:nNodes)[-topDegInd])
    invMapping = match(1:nNodes, fullMapping)

    # get subsample adjacency matrix
    subAdjMat = adjMat[topDegInd, topDegInd]

    # compute approx. eigenvectors
    # for numerical stability get extra eigenvector from eigs
    subEigen = eigs(subAdjMat, nBlocks + 1)
    partBMat = adjMat[-topDegInd, topDegInd]
    tempU = partBMat %*% subEigen$vectors[,1:nBlocks] %*%
        solve(Diagonal(subEigen$values[1:nBlocks]))
    approxEV = rbind2(eigen$vectors[,1:nBlocks], tempU)

    if(rowNorm == T) {
        # project eigenvector rows onto sphere
        approxEV = approxEV/sqrt(rowSums(approxEV^2))
    }

    kmeansResult = bigkmeans(approxEVP, nBlocks, nstart = nIter)

    if(verbose == T) {
        return( list(cluster = kmeansResult$cluster[invMapping],
                     wcss = kmeansResult$tot.withinss,
                     eigenVals = subEigen$values[1:nBlocks],
                     approxEigenVecs = approxEV[invMapping,]) )
    } else {
        return(kmeansResult$cluster[invMapping])
    }
}
