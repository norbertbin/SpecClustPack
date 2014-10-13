#' Simulates a block mixture of Bernoulli covariates.
#'
#' @param covProbMat A matrix with each row corresponding to a block and
#' each column corresponding to a covariate. The values are Bernoulli
#' probabilities.
#' @param nMembers The number of nodes in each block.
#'
#' @export
#' @return A matrix of Bernoulli covariates with a block probability structure.
#'
#' @examples
#' covProbMat = matrix(c(.6, .2, .2, .6), nrow = 2)
#' nMembers = c(50, 50)
#' simBernCovar(covProbMat, nMembers)
simBernCovar <- function(covProbMat, nMembers) {
    nBlocks = dim(covProbMat)[1]
    nCov = dim(covProbMat)[2]
    covMat = NULL
    
    for(j in nCov:1) {
        covCol = NULL
        for(i in nBlocks:1) {
            covTemp = Matrix( 
                simBernSparseVec(nMembers[i], covProbMat[i,j]),
                nrow = nMembers[i], ncol = 1)
          
            covCol = rBind(covTemp, covCol)
        }        
        covMat = cBind(covCol, covMat)        
    }
    return( covMat )
}
