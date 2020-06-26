# data generation for cfa & ccfa with cross loadings & err. covariance library(MASS)
# N=1000;K=3;ipf=8;cpf=2;lam=.7;lac=.3;phi=.3;misp=.1
# necb=1;necw=1;ecb=.3;ecw=.3;cati=NULL;noc=c(4);rseed = 333
sim_lvm <- function(N = 1000, K = 3, ipf = 8, cpf = 2, lam = 0.7, lac = 0.3, phi = 0.3, ecb = 0.3, ecw = 0.3, 
    ome_out = F, LD = 1, cati = NULL, noc = c(4), misp = 0.1, rseed = 333) {
    
    set.seed <- rseed
    options(digits = 4)
    PHI <- matrix(phi, K, K)
    diag(PHI) <- 1
    lam0 <- c(rep(lam, ipf), rep(0, ipf - cpf), rep(lac, cpf), rep(0, (K - 2) * ipf))
    lam1 <- matrix(lam0, ipf, K)
    
    mla <- c()
    for (i in K:1) {
        # i <- 1
        ind <- (c(1:K) + i)%%K
        ind[ind == 0] <- K
        mla <- rbind(mla, lam1[, ind])
    }
    
    J <- K * ipf
    iecb <- iecw <- NULL
    li <- c(1:K) * ipf - cpf
    for (i in 1:K) {
        iecw <- rbind(iecw, c(li[i], li[i] - 1))
        i1 <- i%%K + 1
        iecb <- rbind(iecb, c(li[i] - 2, li[i1] - 3))
    }
    
    ecm <- matrix(0, J, J)
    
    if (LD == 1) {
        ecm[iecw] <- ecm[iecw[, 2:1]] <- ecw
        ecm[iecb] <- ecm[iecb[, 2:1]] <- ecb
    }
    
    evr <- rep(0, J)
    for (j in 1:J) {
        evr[j] = 1 - t(mla[j, ]) %*% PHI %*% mla[j, ]
    }  #end j
    
    diag(ecm) <- evr
    er <- mvrnorm(N, rep(0, J), ecm)
    eta <- mvrnorm(N, rep(0, K), PHI)
    y <- t(mla %*% t(eta)) + er
    scale <- apply(y, 2, sd)
    eigen <- diag(crossprod(mla))
    
    out <- list(N = N, PHI = PHI, MLA = mla, Eigen = eigen, PSX = ecm, scale = scale)
    if (ome_out) 
        out$OME <- eta
    
    pos <- lower.tri(ecm)
    ind <- which(pos, arr.ind = T)
    rind <- which(ecm[pos] > 0)
    out$ofd_ind <- ind[rind, ]
    
    Jp <- length(cati)
    if (Jp > 0) 
        {
            if (Jp == 1 && cati == -1) {
                cati <- c(1:J)
                Jp <- J
            }
            yc <- y[, cati]
            len <- length(noc)
            if (Jp%%len != 0) 
                stop("cati is not divisable by length(noc).", call. = F)
            for (i in 1:len) {
                M <- noc[i]  #categories
                UL <- 3
                LL <- -3
                stp <- LL + c(1:(M - 1)) * (UL - LL)/M  #M-1 step points
                i0 <- (i - 1) * Jp/len + 1
                i1 <- i * Jp/len
                for (j in i0:i1) {
                  tmp <- yc[, j] > matrix(stp, N, (M - 1), byrow = T)
                  yc[, j] <- rowSums(tmp) + 1  # value starting from 1
                }
            }  #end len
            y[, cati] <- yc
            out$noc <- noc
            out$cati <- cati
        }  #end Jp
    
    mind <- matrix(rbinom(N * J, 1, misp), N, J)
    y[mind == 1] <- NA
    out$misp <- misp
    out$dat <- y
    
    return(out)
}
