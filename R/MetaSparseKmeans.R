##' Implementation for Zhiguang Huo, Ying Ding, Shuchang Liu, Steffi Oesterreic
##' and George Tseng. (2015) A Sparse K-means Meta-analysis framework combining
##' multiple transcriptomic studies for disease subtype discovery. Journal of the
##' American Statistical Association. accepted.
##'
##' Here is the instruction about the input
##' @title Main function to perform MetaSparseKmeans.
##' @param x A list for several microarray studies.
##' Each element of the list should be a p*n matrix.
##' p is number of features and n is number of samples.
##' Clustering is performed on sample level.
##' p has to be the same for each study.
##' Missing value should be set to be NA.
##' Current version won't support missing value, it will be allowed in the next ##' version.
##' @param K K specifies number of clusters.
##' We assume the number of clusters to be the same in each study.
##' @param wbounds wbounds is the tuning parameter that controls number of selected features.
##' Larger tuning parameter yield more selected features.
##' wbounds could be a number or a vector.
##' wbounds is suggested to be selected using prior information
##' (e.g. which tuning parameter generate best survival difference.)
##' @param nstart In the MetaSparseKmeans algorithm,
##' there are multiple places in which we will use Kmeans and weighted Kmeans.
##' nstart specify the number of starting point for each of these Kmeans and weighted Kmeans.
##' @param ntrial Since for high dimensional data,
##' it is likely to stuck into a local minimum.
##' ntrial specifies how many times we would repeat the algorithm.
##' The result with the best objective score will be used.
##' @param maxiter The algorithm iteratively update ws (feature weight), cs (clulster assignment) and matching.
##' maxiter specifies the max number of iteration for MetaSparseKmeans
##' @param lambda A tuning parameter controlling the balance between separation ability (BCSS/TSS) and matching function.
##' lambda is set to be 1/2 by default.
##' 1, exhaustive. 2, linear. 3, MCMC.
##' exhaustive will perform exhaustive search and it has 100 percent accurate but very time consuming if the configuration space is too large.
##' linear will perform stepwise search. The searching is very efficient but yield less accurate if the data is noisy.
##' MCMC combine stepwise search and simulated anealling.
##' It utilizes the stepwise search result as initial and perform simulated anealling to find the best matching.
##' MCMC is suggested when the configuration space is too large.In the manuscript, we suggested 14,400 to be the cutoff.
##' If the total number of configurations is greater than 14,400, MCMC is suggested and otherwise exhaustive is suggested.
##' @param sampleSizeAdjust logical argument,
##' controlling whether to adjust for sample size.
##' If true, that means study with larger sample size will have a larger impact.
##' If false, each study has equal contribution.
##' Without prior information, sampleSizeAdjust=FALSE is suggested since we are not sure about data quality.
##' @param wsPre If there is prior knowledge which genes are important, we could specify the initialization of the gene weight.
##' @param silence Logical parameter whether we should print out details.
##' @return The returning result could be a list or a list vector, depending whether the input wbounds is a number or a vector.
##' \item{ws }{Weight for a feature}
##' \item{Cs }{Resulting clustering assignment}
##' \item{wbound }{the used tuning parameter}
##' \item{score }{objective score, the larger the better}
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
##' G2 = G*P2  # DE L-H
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
##'
##' G = 1000
##' n21 = 150
##' n22 = 100
##' n23 = 100
##'
##' label2 = c(rep(1,n21),rep(2,n22),rep(3,n23))
##'
##' P0 = 0.6
##' P1 = 0.1 #common features
##' P2 = 0.1 #common features
##' P3 = 0.1 #noise in S1
##' P4 = 0.1 #noise in S2
##' sd = 0.5
##'
##' G0 = G*P0  # nonDE genes
##' G1 = G*P1  # DE H-L
##' G2 = G*P2  # DE L-H
##' G3 = G*P3  #noise in S1
##' G4 = G*P4  #noise in S2
##'
##' mu211 = runif(G1,-0.25,0.25)
##' mu212 = runif(G1,0.5,1)
##' mu213 = runif(G1,-1,-0.5)
##'
##' mu221 = runif(G2,-1,-0.5)
##' mu222 = runif(G2,-0.25,0.25)
##' mu223 = runif(G2,0.5,1)
##'
##' mu23 = runif(G3,-0.25,0.25)
##'
##' mu241 = runif(G4,-1,-0.5)
##' mu242 = runif(G4,-0.25,0.25)
##' mu243 = runif(G4,0.5,1)
##'
##' mu20 = runif(G0,-0.25,0.25)
##'
##' Data211 = matrix(rnorm(n21*G1,mu211,sd^2),nrow=G1)
##' Data212 = matrix(rnorm(n22*G1,mu212,sd^2),nrow=G1)
##' Data213 = matrix(rnorm(n23*G1,mu213,sd^2),nrow=G1)
##' Data21 = cbind(Data211,Data212,Data213)
##'
##' Data221 = matrix(rnorm(n21*G2,mu221,sd^2),nrow=G2)
##' Data222 = matrix(rnorm(n22*G2,mu222,sd^2),nrow=G2)
##' Data223 = matrix(rnorm(n23*G2,mu223,sd^2),nrow=G2)
##' Data22 = cbind(Data221,Data222,Data223)
##'
##' Data23 = matrix(rnorm((n21+n22+n23)*G3,mu23,sd^2),nrow=G3)
##'
##' Data241 = matrix(rnorm(n21*G4,mu241,sd^2),nrow=G4)
##' Data242 = matrix(rnorm(n22*G4,mu242,sd^2),nrow=G4)
##' Data243 = matrix(rnorm(n23*G4,mu243,sd^2),nrow=G4)
##' Data24 = cbind(Data241,Data242,Data243)
##'
##'
##' Data20 = matrix(rnorm((n21+n22+n23)*G0,mu20,sd^2),nrow=G0)
##'
##' S2 = rbind(Data20,Data21,Data22,Data23,Data24)
##'
##'
##' ## visualize the data
##'
##' S = list(t(S1),t(S2))
##' getWsHeatmap(t(S[[1]]),label1,main='two study before
##' metaSparseKMeans, S1')
##'
##' getWsHeatmap(t(S[[2]]),label2,main='two study before metaSparseKMeans, S2')
##'
##'
##' ## perform meta sparse Kmeans
##'
##' res = MetaSparseKmeans(x=S,K=3,wbounds=10,lambda=2)
##'
##' ## visualize the result
##'
##'
##' getWsHeatmap(t(S[[1]]),res$Cs[[1]],res$ws,main='two study after metaSparseKMeans, S1')
##'
##' getWsHeatmap(t(S[[2]]),res$Cs[[2]],res$ws,main='two study after metaSparseKMeans, S2')
##'
##' plot(res$ws,main='metaSparseKmeans weight dist',xlab='geneIndex')
##'
MetaSparseKmeans <- function(x, K = NULL, wbounds = NULL, nstart = 20, ntrial = 1, maxiter = 20, 
    lambda = 1/2, sampleSizeAdjust = FALSE, wsPre = NULL, silence = FALSE) {
    # x is list of data, (nxp) for each study wbounds is a vector of L1 constraints on w,
    # of the form sum(abs(w))<=wbounds[i] nstart: initial Kmeans searching space.
    # maxiter: maximum number of iterations.  lambda: a tuning parameter keeping the
    # balance between separation and wsPre: if specified, the initial weight.  silence: if
    # print some details
    
    ## check input
    
    ## get some basic information
    numStudies = length(x)
    
    ## create result
    bestOut <- vector("list", length(wbounds))
    out <- vector("list", length(wbounds))
    Cs0 <- list()
    tss.x <- list()
    
    for (i in 1:numStudies) {
        ## total sum of square for each study
        tss.x[[i]] <- apply(scale(x[[i]], center = TRUE, scale = FALSE)^2, 2, sum)
    }
    
    for (atrail in 1:ntrial) {
        # initialize initialize cluster by KMeans initialize w
        if (is.null(wsPre)) {
            for (i in 1:numStudies) {
                Cs0[[i]] <- kmeans(x[[i]], centers = K, nstart = nstart)$cluster
            }
        } else {
            if (length(wsPre) != ncol(x[[1]])) 
                stop("length of wsPre differs from number of genes")
            if (is.null(names(wsPre))) 
                stop("there is no name for wsPre")
            if (any(names(wsPre) != colnames(x[[1]]))) 
                stop("name of wsPre differs from gene name")
            for (i in 1:numStudies) Cs0[[i]] <- weightedKMeans(x = t(x[[i]]), K = K, ws = wsPre, 
                tss.x = tss.x[[i]])
        }
        # 
        for (w in 1:length(wbounds)) {
            awbound = wbounds[w]
            if (!is.null(wsPre)) {
                ws = wsPre
            } else {
                ws <- rep(1, ncol(x[[1]]))/(ncol(x[[1]])) * awbound  # Start with equal weights on each feature
            }
            
            ws.old <- rnorm(ncol(x[[1]]))
            store.ratio <- NULL
            niter <- 0
            trace <- list()
            Cs = Cs0
            
            ## iterate until converge or maxiteration
            while ((sum(abs(ws - ws.old))/sum(abs(ws.old))) > 1e-04 && niter < maxiter) {
                niter <- niter + 1
                ws.old <- ws
                if (niter > 1) 
                  Cs <- UpdateCs(x, K, ws, Cs, tss.x, nstart = nstart)  # if niter=1, no need to update!!
                
                fmatch = patternMatch(x, Cs, ws, silence = silence)
                # patternMatch_old(x, Cs, ws, silence = silence)
                ratio = GetRatio(x, Cs, tss.x, sampleSizeAdjust = sampleSizeAdjust)
                ws <- UpdateWs(x, Cs, awbound, ratio, lambda * (fmatch$perEng + 1)/2)
                store.ratio <- c(store.ratio, sum(ratio * ws))
                
                if (!silence) {
                  cat("iteration:")
                  cat(niter)
                  cat("\n")
                  cat("convergence criteria: ")
                  cat(sum(abs(ws - ws.old))/sum(abs(ws.old)))
                  cat("\n")
                }
                
            }
            score = sum((ratio + lambda * (fmatch$perEng + 1)/2) * ws)
            names(ws) <- colnames(x[[1]])
            out[[w]] <- list(ws = ws, Cs = fmatch$matchCs, wbound = awbound, score = score)
            if (is.null(bestOut[[w]])) {
                bestOut[[w]] = out[[w]]
            } else {
                if (bestOut[[w]]$score < out[[w]]$score) 
                  bestOut[[w]] = out[[w]]
            }
        }
        
    }
    if (length(bestOut) == 1) {
        bestOut = bestOut[[1]]
    }
    class(bestOut) <- "metaSparseKmeans"
    return(bestOut)
}
