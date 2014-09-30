#' Performs spectral clustering on an adjacency matrix.
#'
#' @param adjMat A graph adjacency matrix.
#' @param nBlocks The number of clusters in the graph.
#' @param method The form of the ajacency matrix to be used (default is
#' regLaplacian).
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
specClust <- function(adjMat, nBlocks, method = "regLaplacian",
                             rowNorm = T, nIter = 10, verbose = F) {
    
    # eigs does not support dsCMatrix type at this time
    # Matrix has Namespace problems when using dsCMatrix
    adjMat = as(adjMat, "dgCMatrix")

    similarityMat = getSimilarityMat(adjMat, method)

    eigsDecomp = eigs(similarityMat, nBlocks + 1)

    if(rowNorm == T) {
        eigsDecomp$vectors[,1:nBlocks] = eigsDecomp$vectors[,1:nBlocks] /
            sqrt(rowSums(eigsDecomp$vectors[,1:nBlocks]^2))
    }
    
    kmeansResult = bigkmeans(eigsDecomp$vectors[,1:nBlocks], nBlocks,
        nstart = nIter)

    if(verbose == T) {
        return( list(cluster = kmeansResult$cluster,
                     wcss = kmeansResult$tot.withinss,
                     eigenVals = eigsDecomp$values,
                     eigenVecs = eigsDecomp$vectors) )
    } else {
        return(kmeansResult$cluster)
    }
}

# ---------------------------------------------------------------------
# Helper methods
# ---------------------------------------------------------------------

# ---------------------------------------------------------------------
# Returns the graph similarity matrix corresponding to the given method
# ---------------------------------------------------------------------
getSimilarityMat <- function(adjMat, method) {
    if(method == "regLaplacian") {
        rSums = Matrix::rowSums(adjMat)
        tau = mean(rSums)
        normMat = Diagonal(length(rSums), 1/sqrt(rSums + tau))
        return(normMat %*% adjMat %*% normMat)
    }
    else if(method == "laplacian") {
        rSums = Matrix::rowSums(adjMat)
        normMat = Diagonal(length(rSums), 1/sqrt(rSums))
        return(normMat %*% adjMat %*% normMat)
    }
    else if(method == "adjacency"){
        return(adjMat)
    }
    else {
        stop(paste("Error: method =", method, "Not valid"))
    }
}
