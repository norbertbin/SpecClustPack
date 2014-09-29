#' Simulates a matrix from the Stochastic Blockmodel.
#'
#' @param blockPMat The block probability matrix.
#' @param nMembers A vector specifying the number of nodes in each block.
#'
#' @export
#' @return The sparse representation of an adjacency matrix simulated
#' according to the Stochastic Blockmodel.
#'
#' @keywords stochastic blockmodel
simSBM <- function(blockPMat, nMembers) {

    nBlocks = length(nMembers)
    nNodes = sum(nMembers)
    adjMat = NULL

    # simulate the matrix block by block and bind it together
    for(j in nBlocks:1) {
        adjCol = NULL
        for(i in nBlocks:1) {
            if(i <= j) {
                # to prevent overflow problems
                    nRowxnCol = round(as.numeric(nMembers[i])*
                        as.numeric(nMembers[j]))
                adjTemp = Matrix(simBernSparseVec(nRowxnCol, blockPMat[i,j]),
                    nrow = nMembers[i], ncol = nMembers[j])
            }  
            else {
                adjTemp = Matrix(0, nrow = nMembers[i],
                    ncol = nMembers[j])
            }
            adjCol = rBind(adjTemp, adjCol)
        }        
        adjMat = cBind(adjCol, adjMat)        
    }
    
    return( forceSymmetric(triu(adjMat, k=1)) )
}

#' Simulates a Bernoulli random vector.
#'
#' @param nElem The number of elements in the vector.
#' @param p The probability of an element being equal to one.
#'
#' export
#' @return The sparse representation of a simulated Bernoulli vector.
#'
#' @keywords simulate sparse Bernoulli
simBernSparseVec <- function(nElem, p) {

    if(p == 0) {
        return( sparseVector(0, 1, length = nElem) )
    }

    # get expectation and standard deviation of ones in the vector
    expNumOnes = nElem*p
    sdNumOnes = sqrt(nElem*p*(1-p))

    # vector with intervals at which ones occur in the simulated vector
    oneIntervals = rnbinom(expNumOnes + round(3*sdNumOnes), 1, p) + 1

    # take cumulative sum to get the index values for the ones
    oneIndices = cumsum(oneIntervals)

    # loop to ensure the indices cover the entire simulated vector
    while (max(oneIndices) < nElem) {
        oneIndices = c(oneIndices, max(oneIndices) +
            rnbinom(1, 1, p) + 1)
    }

    # truncate the vector to remove access values
    oneIndices = oneIndices[oneIndices <= nElem]

    return( sparseVector(1, oneIndices, nElem) )
}
