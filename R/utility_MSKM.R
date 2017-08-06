permuteX <- function(x) {
    permx <- list()
    for (i in 1:length(x)) {
        permx[[i]] <- matrix(NA, nrow = nrow(x[[i]]), ncol = ncol(x[[i]]))
        for (j in 1:ncol(x[[i]])) permx[[i]][, j] <- sample(x[[i]][, j])
    }
    return(permx)
}

BinarySearch <- function(argu, sumabs) {
    if (l2n(argu) == 0 || sum(abs(argu/l2n(argu))) <= sumabs) 
        return(0)
    lam1 <- 0
    lam2 <- max(argu) - 1e-05
    iter <- 1
    while (iter <= 15 && (lam2 - lam1) > (1e-04)) {
        su <- soft(argu, (lam1 + lam2)/2)
        if (sum(abs(su/l2n(su))) < sumabs) {
            lam2 <- (lam1 + lam2)/2
        } else {
            lam1 <- (lam1 + lam2)/2
        }
        iter <- iter + 1
    }
    return((lam1 + lam2)/2)
}


encode <- function(avec, n) {
    res <- 0
    for (i in seq_along(avec)) {
        res <- res + avec[i] * n^(i - 1)
    }
    res
}


prepareMCC <- function(x, Cs) {
    numS = length(Cs)
    p = ncol(x[[1]])
    res = list()
    for (i in 1:numS) {
        reduCs = unique(Cs[[i]])
        m = length(reduCs)
        reduMeanTot = colMeans(x[[i]])
        MeanInd = matrix(NA, nrow = m, ncol = p)
        VarInd = matrix(NA, nrow = m, ncol = p)
        MeanTot = colMeans(x[[i]])
        for (j in 1:m) {
            tmpVar = x[[i]][reduCs[j] == Cs[[i]], ]
            if (is.null(dim(tmpVar))) {
                MeanInd[j, ] = tmpVar
                VarInd[j, ] = 0
            } else {
                MeanInd[j, ] = colMeans(tmpVar)
                VarInd[j, ] = apply(tmpVar, 2, var)
            }
            zeroIndex = VarInd[j, ] == 0
            # VarInd[j,zeroIndex] = rep(0.001,sum(zeroIndex))
            VarInd[j, ][zeroIndex] = 1e-04
        }
        denInd = colMeans(VarInd) + colMeans((MeanInd - matrix(rep(MeanTot, each = m), 
            nrow = m))^2)
        res[[i]] = list(MeanInd = MeanInd, denInd = denInd)
    }
    return(res)
}


l2n <- function(vec) {
    return(sqrt(sum(vec^2)))
}


listLength <- function(alist) {
    if (is.list(alist)) 
        return(length(alist))
    return(1)
}


eng_MCC_pair <- function(corS1, corS2, reduCs1, reduCs2) {
    K1 = length(reduCs1)
    K2 = length(reduCs2)
    p = ncol(corS1$MeanInd)
    m = K1
    EXY = numeric(p)
    EX = numeric(p)
    EY = numeric(p)
    
    for (i in 1:m) {
        EXY = EXY + corS1$MeanInd[reduCs1 == i, ] * corS2$MeanInd[reduCs2 == i, ]
        EX = EX + corS1$MeanInd[reduCs1 == i, ]
        EY = EY + corS2$MeanInd[reduCs2 == i, ]
    }
    MCC_num = EXY/m - EX * EY/m^2
    MCC_den = sqrt(corS1$denInd * corS2$denInd)
    MCC_pair_per = MCC_num/MCC_den
    return(MCC_pair_per)
}


eng_cor_per <- function(corPre, reduCs) {
    p = ncol(corPre[[1]]$MeanInd)
    perEng = vector("numeric", p)
    numS = length(corPre)
    reduXComb = as.matrix(combn(numS, 2))
    
    for (i in 1:ncol(reduXComb)) {
        index1 = reduXComb[, i][1]
        index2 = reduXComb[, i][2]
        perEng = perEng + eng_MCC_pair(corPre[[index1]], corPre[[index2]], reduCs[[index1]], 
            reduCs[[index2]])
    }
    ## scale to comparable with bcss/tss
    return(perEng/choose(length(corPre), 2))
}


reorderLabel <- function(alabel, aorder) {
    uniLab = unique(alabel)
    resLabel = NULL
    for (i in 1:length(uniLab)) {
        resLabel[uniLab[i] == alabel] = aorder[i]
    }
    return(resLabel)
}


weightedKMeans <- function(x, K, ws, tss.x = NULL) {
    if (is.null(tss.x)) {
        tss.x <- apply(scale(t(x), center = TRUE, scale = FALSE)^2, 2, sum)
    }
    
    commonNonZeroNames = intersect(rownames(x), names(ws)[ws != 0])
    x <- x[commonNonZeroNames, ]
    tss.x = tss.x[commonNonZeroNames]
    z <- sweep(x, 1, sqrt(ws[commonNonZeroNames]/tss.x), "*")
    km <- kmeans(t(z), centers = K, nstart = 50)
    newCs <- km$cluster
    return(newCs)
}

UpdateCs <- function(x, K, ws, Cs, tss.x, nstart = nstart) {
    newCs <- list()
    for (i in 1:length(x)) {
        x[[i]] <- x[[i]][, ws != 0]
        z <- sweep(x[[i]], 2, sqrt((ws/tss.x[[i]])[ws != 0]), "*")
        nrowz <- nrow(z)
        mus <- NULL
        if (!is.null(Cs[[i]])) {
            for (k in unique(Cs[[i]])) {
                if (sum(Cs[[i]] == k) > 1) 
                  mus <- rbind(mus, apply(z[Cs[[i]] == k, ], 2, mean))
                if (sum(Cs[[i]] == k) == 1) 
                  mus <- rbind(mus, z[Cs[[i]] == k, ])
            }
        }
        if (is.null(mus)) {
            km <- kmeans(z, centers = K, nstart = nstart)
        } else {
            distmat <- as.matrix(dist(rbind(z, mus)))[1:nrowz, (nrowz + 1):(nrowz + K)]
            nearest <- apply(distmat, 1, which.min)
            if (length(unique(nearest)) == K) {
                km <- kmeans(z, centers = mus)
            } else {
                km <- kmeans(z, centers = K, nstart = nstart)
            }
        }
        newCs[[i]] <- km$cluster
    }
    return(newCs)
}


UpdateWs <- function(x, Cs, l1bound, ratio, permatch) {
    a = ratio + permatch
    lam <- BinarySearch(a, l1bound)
    ws.unscaled <- soft(a, lam)
    return(ws.unscaled/l2n(ws.unscaled))
}


soft <- function(x, d) {
    return(pmax(0, x - d))
}


GetRatio <- function(x, Cs, tss.x, sampleSizeAdjust = FALSE) {
    ratio.bcss_tss <- numeric(ncol(x[[1]]))
    for (i in 1:length(x)) {
        wcss.perfeature <- numeric(ncol(x[[i]]))
        tss.perfeature <- tss.x[[i]]
        for (k in unique(Cs[[i]])) {
            whichers <- (Cs[[i]] == k)
            if (sum(whichers) > 1) 
                wcss.perfeature <- wcss.perfeature + apply(scale(x[[i]][whichers, ], center = TRUE, 
                  scale = FALSE)^2, 2, sum)
        }
        aratio <- numeric(ncol(x[[1]]))
        bcss.perfeature = tss.perfeature - wcss.perfeature
        tssNonZeroIndex <- tss.perfeature != 0
        aratio[tssNonZeroIndex] <- bcss.perfeature[tssNonZeroIndex]/tss.perfeature[tssNonZeroIndex]
        
        if (sampleSizeAdjust) {
            thisNsamples <- sapply(x, function(x) nrow(x))
            ratio.bcss_tss = ratio.bcss_tss + aratio * thisNsamples[i]/sum(thisNsamples)
        } else {
            ratio.bcss_tss = ratio.bcss_tss + aratio/length(x)
        }
    }
    ratio.bcss_tss
}


patternMatch <- function(x, Cs, ws, silence = FALSE) {
    numS = length(Cs)
    
    ws2 <- ws[ws != 0]
    
    corPre = prepareMCC(x, Cs)
    corPre2 <- corPre
    for (s in 1:length(corPre2)) {
        corPre2[[s]]$MeanInd <- corPre2[[s]]$MeanInd[, ws != 0]
        corPre2[[s]]$denInd <- corPre2[[s]]$denInd[ws != 0]
    }
    
    S <- length(corPre)
    G <- length(ws)
    
    CombS <- as.data.frame(combn(S, 2))
    encodeS <- sapply(CombS, encode, S)
    
    K <- length(unique(Cs[[1]]))
    
    permK <- permn(K)
    
    energyS <- replicate(length(CombS), list())
    for (i in seq_along(CombS)) {
        index_a <- CombS[[i]][1]
        index_b <- CombS[[i]][2]
        
        corPre_a <- corPre2[[index_a]]
        corPre_b <- corPre2[[index_b]]
        
        energyK <- replicate(length(permK)^2, list())
        encodeK0 <- replicate(length(permK)^2, list())
        for (j1 in seq_along(permK)) {
            aRank <- permK[[j1]]
            for (j2 in seq_along(permK)) {
                bRank <- permK[[j2]]
                energyK[[j2 + (j1 - 1) * length(permK)]] <- sum(eng_MCC_pair(corPre_a, 
                  corPre_b, aRank, bRank) * ws2)
                encodeK0[[j2 + (j1 - 1) * length(permK)]] <- c(aRank, bRank)
            }
        }
        encodeK <- sapply(encodeK0, encode, K)
        
        ahashK <- hash(encodeK, energyK)
        energyS[[i]] <- ahashK
    }
    if (length(energyS) == 1) {
        energyS <- energyS[[1]]
    }
    
    hashS <- hash(encodeS, energyS)
    
    resCs <- as.data.frame(replicate(S, 1:K))
    
    ### exhaustive search start here
    lenCs <- rep(K, S)
    permRule = lapply(lenCs, permn)
    permFlag = rep(1, S)
    permEndFlag = sapply(permRule, listLength)
    
    tmpCs = resCs
    
    iniEnergy <- numeric(S)
    for (s in 1:S) {
        if (s > 1) {
            interEnergy <- 0
            for (ps in 1:(s - 1)) {
                aSencode <- as.character(encode(c(ps, s), S))
                aKencode <- as.character(encode(c(tmpCs[[ps]], tmpCs[[s]]), K))
                interEnergy <- interEnergy + hashS[[aSencode]][[aKencode]]
            }
            iniEnergy[s] <- iniEnergy[s - 1] + interEnergy
        }
    }
    
    highEng <- iniEnergy[S]
    
    while (permFlag[1] == 1) {
        
        tmpEng <- iniEnergy[S]
        print(tmpEng/choose(S, 2))
        if (tmpEng > highEng) {
            highEng = tmpEng
            resCs = tmpCs
        }
        
        # eng_cor_total(corPre, reduCs = tmpCs, ws = ws) eng_cor_total(corPre2, reduCs =
        # tmpCs, ws = ws2)
        
        permFlag[[S]] = permFlag[[S]] + 1
        
        updateFlag <- rep(0, S)
        updateFlag[S] <- 1
        for (s in S:2) {
            if (permFlag[[s]] > permEndFlag[[s]]) {
                permFlag[[s]] = 1
                permFlag[[s - 1]] = permFlag[[s - 1]] + 1
                updateFlag[s] <- 1
                updateFlag[s - 1] <- 1
                
            }
        }
        
        for (s in 2:S) {
            if (updateFlag[s] == 1) {
                tmpCs[[s]] <- permRule[[s]][[permFlag[[s]]]]
                
                interEnergy <- 0
                for (ps in 1:(s - 1)) {
                  aSencode <- as.character(encode(c(ps, s), S))
                  aKencode <- as.character(encode(c(tmpCs[[ps]], tmpCs[[s]]), K))
                  interEnergy <- interEnergy + hashS[[aSencode]][[aKencode]]
                }
                
                iniEnergy[s:S] <- iniEnergy[s:S] + interEnergy - iniEnergy[s] + iniEnergy[s - 
                  1]
            }
        }
        # print(iniEnergy/3)
    }
    
    
    #################### stop here
    perEng = eng_cor_per(corPre = corPre, reduCs = resCs)
    resumeCs = Cs
    for (s in 1:numS) {
        resumeCs[[s]] = reorderLabel(Cs[[s]], resCs[[s]])
    }
    ## return resumed matching Cs, high energy, energy per gene, trace
    return(list(matchCs = resumeCs, highEng = highEng, perEng = perEng))
}

if (F) {
    
    patternMatch_old <- function(x, Cs, ws, silence = FALSE) {
        numS = length(Cs)
        corPre = prepareMCC(x, Cs)
        lenCs = vector("numeric", numS)
        uniCs = vector("list", numS)
        
        for (s in 1:numS) {
            uniCs[[s]] = unique(Cs[[s]])
            lenCs[s] = length(uniCs[[s]])
        }
        resCs = lapply(uniCs, sort)  ## result Cluster label, start form 12345, then exhausive search
        
        ## how to calculate energy, resCs specify matching rule
        highEng = eng_cor_total(corPre, reduCs = resCs, ws = ws)
        
        ### exhaustive search start here
        stdMatch = 1:lenCs[1]
        combRule = lapply(lenCs, function(x) as.matrix(combn(stdMatch, x)))
        permRule = lapply(lenCs, permn)
        combFlag = rep(1, numS)
        permFlag = rep(1, numS)
        combEndFlag = sapply(combRule, ncol)
        permEndFlag = sapply(permRule, listLength)
        tmpCombFlag = combFlag
        tmpPermFlag = permFlag
        tmpReduCs = resCs
        trace = highEng
        while (tmpPermFlag[1] == 1) {
            combFlag = tmpCombFlag
            permFlag = tmpPermFlag
            tmpEng = eng_cor_total(corPre = corPre, reduCs = tmpReduCs, ws = ws)
            print(tmpEng)
            
            if (tmpEng > highEng) {
                highEng = tmpEng
                resCs = tmpReduCs
            }
            if (!silence) 
                # print(highEng)
            trace = c(trace, highEng)
            tmpPermFlag[[length(Cs)]] = tmpPermFlag[[length(Cs)]] + 1
            for (s in length(Cs):2) {
                if (tmpPermFlag[[s]] > permEndFlag[[s]]) {
                  tmpPermFlag[[s]] = 1
                  tmpCombFlag[[s]] = tmpCombFlag[[s]] + 1
                  if (tmpCombFlag[[s]] > combEndFlag[[s]]) {
                    tmpCombFlag[[s]] = 1
                    tmpPermFlag[[s - 1]] = tmpPermFlag[[s - 1]] + 1
                  }
                }
                tmpOrder = combRule[[s]][, tmpCombFlag[[s]]][permRule[[s]][[tmpPermFlag[s]]]]
                tmpReduCs[[s]] = reorderLabel(resCs[[s]], tmpOrder)
            }
        }
        
        #################### stop here
        perEng = eng_cor_per(corPre = corPre, reduCs = resCs)
        resumeCs = Cs
        for (s in 1:numS) {
            resumeCs[[s]] = reorderLabel(Cs[[s]], resCs[[s]])
        }
        ## return resumed matching Cs, high energy, energy per gene, trace
        return(list(matchCs = resumeCs, highEng = highEng, perEng = perEng))
    }
    
    
    eng_cor_total <- function(corPre, reduCs, ws) {
        non0ws = ws != 0
        perEng = vector("numeric", sum(non0ws))
        
        numS = length(corPre)
        for (i in 1:numS) {
            corPre[[i]]$MeanInd = corPre[[i]]$MeanInd[, non0ws]
            corPre[[i]]$denInd = corPre[[i]]$denInd[non0ws]
        }
        reduXComb = as.matrix(combn(numS, 2))
        
        for (i in 1:ncol(reduXComb)) {
            index1 = reduXComb[, i][1]
            index2 = reduXComb[, i][2]
            perEng = perEng + eng_MCC_pair(corPre2[[index1]], corPre2[[index2]], reduCs[[index1]], 
                reduCs[[index2]])
        }
        ## scale to comparable with bcss/tss
        perEng = perEng/choose(length(corPre), 2)
        return(sum(perEng * ws[non0ws]))
    }
    
}
