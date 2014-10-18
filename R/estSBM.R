#' Estimates the Stochastic Blockmodel based on the adjacency matrix and given
#' block membership.
#'
#' @param adjMat The adjacency matrix.
#' @param blockAssign A vector of node block assignments.
#'
#' @export
#' @return The estimated block probability matrix from the Stochastic Blockmodel.
#' @keywords stochastic blockmodel
estSBM <- function(adjMat, blockAssign) {

   # make sure the matrix is stored in nonsymmetric form
   # for this implementation to work correctly 
   adjMat = as(adjMat, "dgCMatrix") 
    
   nMembers = tabulate(blockAssign)
    
   nBlocks = length(unique(blockAssign))
   blockPMat = matrix(0, nrow = nBlocks, ncol = nBlocks)  

   # only consider nonzero entries of adjMat
   adjMat = as.matrix(Matrix::summary(adjMat))

   # get node block assignments for each edge
   bi = blockAssign[adjMat[,1]]
   bj = blockAssign[adjMat[,2]]

   for(i in 1:nBlocks) {
       for(j in i:nBlocks) {
           if(i == j){
               blockPMat[i, i] = sum((bi == i & bj == i)) /
                   (nMembers[i]^2 - nMembers[i])
           }
           else {
               blockPMat[i, j] = sum((bi == i & bj == j)) /
                   (nMembers[i]*nMembers[j])
               blockPMat[j, i] = blockPMat[i, j]
           }
       }
   }

	return(blockPMat)
} 
