#' Plots the adjacency matrix.
#'
#' @param adjMat The adjacency matrix to be plotted.
#' @param membership The node membership used to reorder the adjacency matrix
#' before plotting.
#'
#' @export
#'  
plotAdj <- function(adjMat, membership = NULL) {
    
    if(!is.null(membership)) {
        adjMat = adjMat[membership, membership]
    }

   image(adjMat, scales = list(draw=F), xlab = NULL, ylab = NULL, lwd = .1,
         sub = NULL)
}

#' Plots the Stochastic Blockmodel.
#'
#' @param blockPMat The block probability matrix.
#' @param nMembers A vector of the number of nodes in each block.
#'
#' @export 
#' 
plotSBM <- function(blockPMat, nMembers = NULL) {

    nBlocks = dim(blockPMat)[1]
    # if nMembers not given assume equal number of nodes in each block
    if(is.null(nMembers)) {
        nMembers = rep(1, nBlocks)
    }
    
    #adjust the block end points based on number of members
    blockGrid = c(0, cumsum(nMembers/sum(nMembers)))

    zlimMax = min(1, 1.2*max(blockPMat))
    zlimMin = max(0, .8*min(blockPMat))
    
    image(blockGrid, blockGrid, blockPMat[,nBlocks:1],
          zlim = c(zlimMin, zlimMax), xlim = c(0,1), ylim = c(0,1),
          col = gray(1000:0/1000), axes = F, xlab = "", ylab = "", lwd = .1)
}
