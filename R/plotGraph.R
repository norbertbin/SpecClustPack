#' Plots the adjacency matrix.
#'
#' @param adjMat The adjacency matrix to be plotted.
#' @param membership The node membership used to reorder the adjacency matrix
#' before plotting.
#' @param lines If true, plot lines separating blocks when membership given.
#'
#' @export
#'  
plotAdj <- function(adjMat, membership = NULL, lines = T) {

    adjMat = as(adjMat, "dgCMatrix")
    df = data.frame(Matrix::summary(adjMat))
    
    if(!is.null(membership)) {
        newOrder = match(1:length(membership),
            sort(membership, index.return = T)$ix)
        df$i = newOrder[df$i]
        df$j = newOrder[df$j]
    }

p1 = ggplot(df, aes(i, j, fill = x)) +
     geom_raster(hjust=0, vjust=0) +
     theme(axis.line=element_blank(),
           axis.text.x=element_blank(),
           axis.text.y=element_blank(),
           axis.ticks=element_blank(),
           axis.title.x=element_blank(),
           axis.title.y=element_blank(),
           legend.position="none",
           panel.background=element_blank(),
           panel.border=element_blank(),
           panel.grid.major=element_blank(),
           panel.grid.minor=element_blank(),
           plot.background=element_blank()) +
     coord_cartesian(xlim = c(0,max(df$i)),
                     ylim = c(0,max(df$j))) +
     scale_fill_gradient(low = 'gray40',
                         high = 'black',
                         limits = c(quantile(df$x, probs=.05),
                                    quantile(df$x, probs=.95)))

    if(!is.null(membership) & lines) {
        border = c(0, cumsum(table(membership)))
        p1 = p1 + geom_hline(yintercept = border) +
            geom_vline(xintercept = border)
    }

    print(p1)
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
