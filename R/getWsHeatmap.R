##' This function provide a useful function to visualize gene expression matrix.
##' ratio adjusted gene-wise normalization will be utilized.
##' Cluster assignment and non-zero weight genes must be provided.
##'
##' nothing at this moment
##' @title A function for visualization gene expression matrix
##' @param astudy p by n matrix.
##' p represents number of features and n represents number of samples
##' @param aCs cluster assignment labels. Total length should be n.
##' @param ws Logical vector indicating whether some genes should show up in the heatmap.
##' Length of ws should be p
##' @param Rowv If TRUE, genes will be clustered using hierarchical clustering algorithm.
##' If NA, gene will use the order from argument geneOrder.
##' @param geneOrder If Rowv is NA, then the genes will be order by argument of geneOrder.
##' Length of geneOrder should be the number of non-zero weight.
##' @param \dots Other parameters inherited from function heatmap function.
##' @return Return heatmap object. Refer to ?heatmap
##' @author Zhiguang Huo
##' @export
##' @examples
##' ######################################
##' ## generate data
##' set.seed(15213)
##'
##' G = 1000
##' n11 = 100
##' n12 = 100
##' n13 = 150
##' label1 = c(rep(1,n11),rep(2,n12),rep(3,n13))
##'
##' P0 = 0.6
##' P1 = 0.1
##' P2 = 0.1
##' P3 = 0.1
##' P4 = 0.1
##' sd = 0.5
##'
##' G0 = G*P0  # nonDE genes
##' G1 = G*P1  # DE H-L
##' G2 = G*P2\t# DE L-H
##' G3 = G*P3
##' G4 = G*P4
##'
##'
##' mu111 = runif(G1,-0.25,0.25)
##' mu112 = runif(G1,0.5,1)
##' mu113 = runif(G1,-1,-0.5)
##'
##' mu121 = runif(G2,-1,-0.5)
##' mu122 = runif(G2,-0.25,0.25)
##' mu123 = runif(G2,0.5,1)
##'
##' mu131 = runif(G3,-1,-0.5)
##' mu132 = runif(G3,-0.25,0.25)
##' mu133 = runif(G3,0.5,1)
##'
##' mu14 = runif(G4,-0.25,0.25)
##' mu10 = runif(G0,-0.25,0.25)
##'
##' Data111 = matrix(rnorm(n11*G1,mu111,sd^2),nrow=G1)
##' Data112 = matrix(rnorm(n12*G1,mu112,sd^2),nrow=G1)
##' Data113 = matrix(rnorm(n13*G1,mu113,sd^2),nrow=G1)
##' Data11 = cbind(Data111,Data112,Data113)
##'
##' Data121 = matrix(rnorm(n11*G2,mu121,sd^2),nrow=G2)
##' Data122 = matrix(rnorm(n12*G2,mu122,sd^2),nrow=G2)
##' Data123 = matrix(rnorm(n13*G2,mu123,sd^2),nrow=G2)
##' Data12 = cbind(Data121,Data122,Data123)
##'
##' Data131 = matrix(rnorm(n11*G3,mu131,sd^2),nrow=G3)
##' Data132 = matrix(rnorm(n12*G3,mu132,sd^2),nrow=G3)
##' Data133 = matrix(rnorm(n13*G3,mu133,sd^2),nrow=G3)
##' Data13 = cbind(Data131,Data132,Data133)
##'
##' Data14 = matrix(rnorm((n11+n12+n13)*G4,mu14,sd^2),nrow=G4)
##'
##' Data10 = matrix(rnorm((n11+n12+n13)*G0,mu10,sd^2),nrow=G0)
##'
##' S1 = rbind(Data10,Data11,Data12,Data13,Data14)
##'
##' getWsHeatmap(S1,label1,main='two study before
##' metaSparseKMeans, S1')
##'
getWsHeatmap <- function(astudy, aCs, ws = NULL, Rowv = NA, geneOrder = NULL, maxIndex = 3, ...) {
    if (is.null(ws)) 
        ws = rep(1, nrow(astudy))
    astudy = astudy[ws != 0, ]
    coll <- NULL
    resHeatmap <- NULL
    label <- NULL
    uniCs = unique(aCs)
    uniorderCs = sort(uniCs)
    
    for (i in uniorderCs) {
        if (is.null(resHeatmap)) {
            resHeatmap = astudy[, i == aCs]
            label = rep(i, sum(i == aCs))
        } else {
            resHeatmap = cbind(resHeatmap, astudy[, i == aCs])
            label = c(label, rep(i, sum(i == aCs)))
        }
    }
    for (alabel in unique(label)) coll[which(label == alabel)] = palette()[alabel]
    finalRes = t(scale(t(resHeatmap)))
    if (!is.null(geneOrder)) {
        finalRes = finalRes[geneOrder, ]
    }
    outIndex = abs(finalRes) > maxIndex
    finalRes[outIndex] = maxIndex * sign(finalRes)[outIndex]
    B = 16
    # return(gplots::heatmap.2(finalRes, col='greenred' ,trace='none',Rowv=Rowv,ColSideColors=coll,
    # Colv=NA,keysize=1.3,...)  )
    return(heatmap(finalRes, Rowv = Rowv, ColSideColors = coll, col = rgb(c(rep(0, B), (0:B)/B), c((B:0)/16, 
        rep(0, B)), rep(0, 2 * B + 1)), scale = "none", Colv = NA, ...))
}
