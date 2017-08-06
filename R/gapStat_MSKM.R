##' This function could calculate original objective score and permuted objective score.
##'
##' nothing at this moment
##' @title Calculate Gap statistics
##' @param S Input multiple studies.
##' Same format as the input for MetaSparseKmeans fucntion.
##' @param K Number of clusters
##' @param B Permutation number.
##' @param wbounds Tuning parameter. Again is could be a number or a vector.
##' @param nval length of wbounds if no input for wbounds
##' @return Objective score
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
##' G2 = G*P2 # DE L-H
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
##'
##' B <- 2
##'
##' gapStatResult <- gapStat_MSKM(S,K=3,B=2)
##' plot(gapStatResult$wbounds,gapStatResult$gapStat,type='b',xlab='mu',ylab='gapStat') 
##' arrows(gapStatResult$wbounds, gapStatResult$gapStat-gapStatResult$se.score, gapStatResult$wbounds, gapStatResult$gapStat+gapStatResult$se.score, length=0.05, angle=90, code=3)
##'
##'
gapStat_MSKM <- function(S, K = 3, B = 10, wbounds = NULL, nvals = 10, silence = FALSE) {
    if (B != (B. <- as.integer(B)) || (B <- B.) <= 0) 
        stop("'B' has to be a positive integer")
    if (is.null(wbounds)) 
        wbounds <- exp(seq(log(1.2), log(sqrt(ncol(S[[1]])) * 0.9), len = nvals))
    if (min(wbounds) <= 1) 
        stop("Wbounds should be greater than 1, since otherwise only one weight will be nonzero.")
    if (length(wbounds) < 2) 
        stop("Wbounds should be a vector of at least two elements.")
    
    ## get true objective score
    set.seed(15213)
    tmpres <- MetaSparseKmeans(x = S, K = K, wbounds = wbounds, silence = TRUE)
    score <- sapply(tmpres, function(x) x$score)
    
    E.score.full <- NULL
    
    if (!silence) 
        cat("calculating gap stat, b = 1,2,..., B (= ", B, ")  [one \".\" per sample]:\n", sep = "")
    
    for (b in 1:B) {
        cat(".", if (b%%50 == 0) 
            paste(b, "\n"))
        
        set.seed(15213 + b)
        ax = permuteX(S)
        tmpres = MetaSparseKmeans(x = ax, K = K, wbounds = wbounds, silence = TRUE)
        ascore = sapply(tmpres, function(x) x$score)
        E.score.full <- rbind(E.score.full, ascore)
    }
    
    if (!silence && (B%%50 != 0)) 
        cat("", B, "\n")
    
    E.score <- colMeans(E.score.full)
    se.score <- apply(E.score.full, 2, sd)/nrow(E.score.full)
    
    gapStat <- score - E.score
    
    result <- data.frame(wbounds, score, E.score, gapStat, se.score)
    return(result)
}
