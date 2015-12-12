# casc
#' Covariate Assisted Spectral Clustering
#' 
#' @param adjMat An adjacency matrix
#' @param covMat A covariate matrix
#' @param nBlocks The number of clusters
#' @param nPoints Number of iterations to find the optimal tuning
#' parameter.
#' @param method The form of the adjacency matrix to be used.
#' @param rowNorm True if row normalization should be
#' done before running kmeans.
#' @param center A boolean indicating if the covariate matrix columns
#' should be centered.
#' @param verbose A boolean indicating if casc output should include eigendecomposition.
#' @param assortative A boolean indicating if the assortative version of casc should be used.
#' @param randStarts Number of random restarts for kmeans.
#'
#' @export
#' @return A list with node cluster assignments, the
#' the value of the tuning parameter used, the within
#' cluster sum of squares, and the eigengap.
#'
#' @keywords spectral clustering
casc <- function(adjMat, covMat, nBlocks, nPoints = 100,
                 method = "regLaplacian", rowNorm = F,
                 center = F, verbose = F,
                 assortative = F, randStarts = 20) {

    # Matrix has Namespace problems when using dsCMatrix
    adjMat = as(adjMat, "dgCMatrix")
    
    adjMat <- getGraphMatrix(adjMat, method)
    covMat <- scale(covMat, center = center,
                    scale = sqrt(Matrix::colSums(covMat^2)))

    # return
    getCascClusters(adjMat, covMat, nBlocks, nPoints,
                            rowNorm, verbose,
                            assortative, randStarts)    
}

# CCA
#' A modified CCA approach for spectral clustering.
#' 
#' Uses the eigenvectors of adjMat %*% covMat for spectral clustering.
#' 
#' @param adjMat An adjacency matrix
#' @param covMat A covariate matrix
#' @param nBlocks The number of clusters
#' @param method The form of the adjacency matrix to be used.
#' @param rowNorm True if row normalization should be
#' done before running kmeans.
#' @param center A boolean indicating if the covariate matrix columns
#' should be centered.
#' @param verbose A boolean indicating if casc output should include eigendecomposition.
#' @param randStarts Number of random restarts for kmeans.
#'
#' @export
#' @return A list with node cluster assignments, the
#' the value of the tuning parameter used, the within
#' cluster sum of squares, and the eigengap.
#'
#' @keywords spectral clustering
cca <- function(adjMat, covMat, nBlocks,
                 method = "regLaplacian", rowNorm = F,
                 center = F, verbose = F,
                 randStarts = 10) {

    nCov = ncol(covMat)

    # Matrix has Namespace problems when using dsCMatrix
    adjMat = as(adjMat, "dgCMatrix")
    
    adjMat <- getGraphMatrix(adjMat, method)
    covMat <- scale(covMat, center = center,
                    scale = sqrt(Matrix::colSums(covMat^2)))

    ccaSvd = svd(adjMat %*% covMat, nu = min(nBlocks, nCov), nv = 0)

    if(rowNorm == T) {
        ccaSvd$u = ccaSvd$u/sqrt(colSums(ccaSvd$u^2))
    }

    kmeansResults = kmeans(ccaSvd$u, nBlocks, nstart = randStarts)

    if(verbose == F) {
        return(list(cluster = kmeansResults$cluster,
            wcss = kmeansResults$tot.withinss
            ))
    } else {
        return(list(cluster = kmeansResults$cluster,
            wcss = kmeansResults$tot.withinss,
            ccaSvd = ccaSvd
            ))
    }
}


# ---------------------------------------------------------------------
# Helper methods
# ---------------------------------------------------------------------

# ---------------------------------------------------------------------
# returns the graph matrix corresponding to the given method
# ---------------------------------------------------------------------
getGraphMatrix = function(adjacencyMat, method) {

    if(method == "regLaplacian") {
        rSums = Matrix::rowSums(adjacencyMat)
        tau = mean(rSums)
        normMat = Diagonal(length(rSums), 1/sqrt(rSums + tau))
        return(normMat %*% adjacencyMat %*% normMat)
    }
    else if(method == "laplacian") {
        rSums = Matrix::rowSums(adjacencyMat)
        normMat = Diagonal(length(rSums), 1/sqrt(rSums))
        return(normMat %*% adjacencyMat %*% normMat)
    }
    else if(method == "adjacency"){
        return(adjacencyMat)
    }
    else {
        stop(paste("Error: method =", method, "Not valid"))
    }
}

# ---------------------------------------------------------------------
# returns CASC optimal h tuning parameter SVD
# ---------------------------------------------------------------------
getCascClusters = function(graphMat, covariates, nBlocks,
    nPoints, rowNorm, verbose, assortative,
    randStarts) {
    
    rangehTuning = getTuningRange(graphMat, covariates, nBlocks, 
        assortative)

    hTuningSeq = seq(rangehTuning$hmin, rangehTuning$hmax,
        length.out = nPoints)
    wcssVec = vector(length = nPoints)
    gapVec = vector(length = nPoints)
    orthoX = vector(length = nPoints)
    orthoL = vector(length = nPoints)
    
    for(i in 1:nPoints) {
        cascResults = getCascResults(graphMat, covariates, hTuningSeq[i],
            nBlocks, rowNorm, F, assortative, randStarts)
        orthoX[i] = cascResults$orthoX
        orthoL[i] = cascResults$orthoL
        wcssVec[i] = cascResults$wcss
        gapVec[i] = cascResults$singGap
    }

    hOpt = hTuningSeq[which.min(wcssVec)]

    hOptResults = getCascResults(graphMat, covariates, hOpt, nBlocks, rowNorm, 
        verbose, assortative, randStarts)
    
    if(verbose == F) {
        return( list(cluster = hOptResults$cluster,
                 h = hOpt,
                 wcss = hOptResults$wcss,
                 eigenGap = hOptResults$eigenGap) )
    } else {
        return( list(cluster = hOptResults$cluster,
                 h = hOpt,
                 wcss = hOptResults$wcss,
                 eigenGap = hOptResults$eigenGap,
                 cascSvd = hOptResults$cascSvd) )
    }
}

# ---------------------------------------------------------------------
# returns cluster memberships for CASC based clustering takes graphMat
# ---------------------------------------------------------------------
getCascResults = function(graphMat, covariates, hTuningParam,
    nBlocks, rowNorm, verbose, assortative, randStarts) {
    
    cascSvd = getCascSvd(graphMat, covariates, hTuningParam, nBlocks, assortative)

    ortho = getOrtho(graphMat, covariates, cascSvd$eVec, cascSvd$eVal,
        hTuningParam, nBlocks)

    if(rowNorm == T) {
        cascSvd$eVec = cascSvd$eVec/sqrt(colSums(cascSvd$eVec^2))
    }
    
    kmeansResults = kmeans(cascSvd$eVec, nBlocks, nstart = randStarts)
    
    if(verbose == F) {
        return( list(cluster = kmeansResults$cluster,
            wcss = kmeansResults$tot.withinss,
            singGap = cascSvd$eVal[nBlocks] -
                cascSvd$eVal[nBlocks + 1],
            orthoL = ortho$orthoL,
            orthoX = ortho$orthoX) )
    } else {
         return( list(cluster = kmeansResults$cluster,
            wcss = kmeansResults$tot.withinss,
            eigenGap = cascSvd$eVal[nBlocks] -
                cascSvd$eVal[nBlocks + 1],
            orthoL = ortho$orthoL,
            orthoX = ortho$orthoX,
            cascSvd = cascSvd) )
    }
}

# ---------------------------------------------------------------------
# returns left singular vectors and values for CASC based clustering
# ---------------------------------------------------------------------
getCascSvd = function(graphMat, covariates, hTuningParam, nBlocks, assortative) {

    # # define custom matrix vector multiply function
    # if(assortative == T) {
    #     matMult = function(x, args) {
    #         as.numeric(args$graphMat %*% x + 
    #             args$h * args$cov %*% crossprod(args$cov, x))
    #     }
    # } else {
    #     matMult = function(x, args) {
    #         as.numeric(args$graphMat %*% (args$graphMat %*% x) + 
    #             args$h * args$cov %*% crossprod(args$cov, x))
    #     }
    # }

    # eDecomp = eigs(matMult, nBlocks + 2, n = nrow(graphMat), which = "LR", 
    #     args = list(graphMat = graphMat, cov = covariates, h = hTuningParam))

    # # return
    # list(eVec = eDecomp$vectors[, 1:nBlocks],
    #         eVal = eDecomp$values[1:(nBlocks+1)],
    #         eVecKPlus = eDecomp$vectors[, nBlocks+1])

    # using irlba
    # define custom matrix vector multiply function
    if(assortative == T) {
        matMult = function(A, x, transpose) {
            as.numeric(graphMat %*% x + 
                hTuningParam * covariates %*% crossprod(covariates, x))
        }
    } else {
        matMult = function(A, x, transpose) {
            as.numeric(graphMat %*% (graphMat %*% x) + 
                hTuningParam * covariates %*% crossprod(covariates, x))
        }
    }

    sDecomp = irlba(A = graphMat, nu = nBlocks + 1, nv = 0,
        matmul = matMult)

    # return
    list(eVec = sDecomp$u[, 1:nBlocks],
            eVal = sDecomp$d[1:(nBlocks+1)],
            eVecKPlus = sDecomp$u[, nBlocks+1])
}


# ---------------------------------------------------------------------
# gets a good range for the tuning parameter in CASC
# ---------------------------------------------------------------------
getTuningRange = function(graphMatrix, covariates, nBlocks,
    assortative) {

    nCov = ncol(covariates)
    
    singValCov = svd(covariates, nu = min(nBlocks, nCov))$d

    if(assortative == T) {
        # eigenValGraph = eigs(graphMatrix, nBlocks + 2, which = "LR",
        #     opts = list(retvec = F))$values #eigenvalues only

        eigenValGraph = irlba(graphMatrix, nu = nBlocks + 1, nv = 0, 
            m_b = max(20, 2*nBlocks))$d

        if(nCov > nBlocks) {
            hmax = eigenValGraph[1]/(singValCov[nBlocks]^2 - singValCov[nBlocks+1]^2) 
        } else {
            hmax = eigenValGraph[1]/singValCov[nCov]^2 
        }
        hmin = (eigenValGraph[nBlocks] - eigenValGraph[nBlocks + 1])/singValCov[1]^2
    } else {
        # eigenValGraph = eigs(graphMatrix, nBlocks + 2,
        #     opts = list(retvec = F))$values #eigenvalues only
        # eigenValGraph = sort(eigenValGraph^2, decreasing=T)

        eigenValGraph = (irlba(graphMatrix, nu = nBlocks + 1, nv = 0, 
            m_b = max(20, 2*nBlocks))$d)^2

        if(nCov > nBlocks) {
            hmax = eigenValGraph[1]/(singValCov[nBlocks]^2 - singValCov[nBlocks+1]^2) 
        } else {
            hmax = eigenValGraph[1]/singValCov[nCov]^2 
        }
        hmin = (eigenValGraph[nBlocks] - eigenValGraph[nBlocks + 1])/singValCov[1]^2
    }

    # return
    list( hmax = hmax, hmin = hmin )
}

# ---------------------------------------------------------------------
# returns the proportion of the eigenvalues due to X in the top eigenspace
# ---------------------------------------------------------------------
getOrtho <- function(graphMat, covariates, cascSvdEVec, cascSvdEVal,
                     h, nBlocks) {
    orthoL <- as.numeric((t(cascSvdEVec[, nBlocks])%*%graphMat%*%
           cascSvdEVec[, nBlocks])/cascSvdEVal[nBlocks])
    orthoX <- as.numeric(h*(t(cascSvdEVec[, nBlocks])%*%covariates%*%
                          t(covariates)%*%cascSvdEVec[, nBlocks])/
                       cascSvdEVal[nBlocks])
    return( list(orthoL = orthoL/(orthoL + orthoX),
                 orthoX = orthoX/(orthoL + orthoX)) )
}
