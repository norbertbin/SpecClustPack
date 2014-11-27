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

    eigsDecomp = eigs(similarityMat, nBlocks + 3)

    if(rowNorm == T) {
        eigsDecomp$vectors[,1:nBlocks] = eigsDecomp$vectors[,1:nBlocks] /
            sqrt(rowSums(eigsDecomp$vectors[,1:nBlocks]^2))

        # if there were rows of zeros need to handle NaN's
        eigsDecomp$vectors[is.nan(eigsDecomp$vectors)] = 0
    }
    
    kmeansResult = bigkmeans(eigsDecomp$vectors[,1:nBlocks], nBlocks,
        nstart = nIter)

    if(verbose == T) {
        return( list(cluster = kmeansResult$cluster,
                     wcss = kmeansResult$tot.withinss,
                     eigenVals = eigsDecomp$values[1:(nBlocks+1)],
                     eigenVecs = eigsDecomp$vectors[,1:(nBlocks+1)]) )
    } else {
        return(kmeansResult$cluster)
    }
}

#' Spectral clustering of a covariate matrix (low rank).
#'
#' @param covMat A matrix of covariates with each column corresponding to
#' a covariate.
#' @param nBlocks The number of clusters.
#' @param nIter Number of kmeans iterations (default is 10).
#'
#' @export
#' @return Returns a vector of row cluster assignments.
#'
#' @examples
#' covProbMat = matrix(c(.6, .2, .2, .6), nrow = 2)
#' nMembers = c(50, 50)
#' covMat = simBernCovar(covProbMat, nMembers)
#' specClustCov(covMat, nBlocks = 2)
specClustCov <- function(covMat, nBlocks, nIter = 10) {

    #center and normalize columns
    covMat = scale(covMat, center = T,
        scale = sqrt(Matrix::colSums(covMat^2)))    
    
    svdDecomp = svd(covMat)

    kmeansResult = bigkmeans(svdDecomp$u[,1:nBlocks], nBlocks,
        nstart = nIter)
    
    return(kmeansResult$cluster)
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
