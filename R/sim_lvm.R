#' @title Simulating data with Latent Variable Modeling
#'
#' @description \code{sim_lvm} can simulate data based on factor analysis or
#' item response models with different response formats (continuous or categorical),
#' loading patterns and residual covariance (local dependence) structures.
#'
#' @name sim_lvm
#'
#' @param N Sample size.
#'
#' @param K Number of factors.
#'
#' @param ipf Items per factor.
#'
#' @param cpf Cross-loadings per factor.
#'
#' @param lam Number of formal iterations for posterior sampling.
#'
#' @param lac Number of iterations to update the sampling information.
#'
#' @param phi Homogeneous correlations between any two factors.
#'
#' @param ph1 Correlation between factor 1 and 2 (if it's different from \code{phi}.
#'
#' @param ecr Residual correlation (local dependence).
#'
#' @param ome_out Output factor score or not.
#'
#' @param cati The set of categorical (polytomous) items in sequence number (i.e., 1 to \eqn{J});
#' \code{NULL} for no and -1 for all (default is \code{NULL}).
#'
#' @param noc Number of categories for categorical items
#'
#' @param misp Proportion of missingness.
#'
#' @param rseed An integer for the random seed.
#'
#' @param necb Number of between-factor local dependence.
#'
#' @param necw Number of within-factor local dependence.
#'
#' @param add_ind (Additional) minor factor with cross-loadings.
#'
#' @param add_la Value of cross-loadings on (Additional) minor factor.
#'
#' @param add_phi Correlations between (Additional) minor factor and other factors.
#'
#' @param zero_it Surplus items with zero loading.
#'
#' @param digits Number of significant digits to print when printing numeric values.
#'
#' @return An object of class \code{list} containing the data, loading, and factorial correlation matrix.
#'
#'
#' @importFrom MASS mvrnorm
#'
#' @export
#'
#' @examples
#'
#' # for continuous data with cross-loadings and local dependence effect .3
#' out <- sim_lvm(N=1000,K=3,ipf=6,lam = .7, lac=.3,ecr=.3)
#' summary(out$dat)
#' out$MLA
#' out$ofd_ind
#'
#' # for categorical data with cross-loadings .4 and 10% missingness
#' out <- sim_lvm(N=1000,K=3,ipf=6,lam = .7, lac=.4,cati=-1,noc=4,misp=.1)
#' summary(out$dat)
#' out$MLA
#' out$ofd_ind
#'
sim_lvm <- function(N = 1000, K = 3, ipf = 8, cpf = 2, lam = 0.7, lac = 0.3, phi = 0.5, ph1 = -1,
        ecr = .0, ome_out = FALSE,cati = NULL, noc = c(4), misp = 0, rseed = 333,
        necw=K,necb=K,add_ind=c(),add_la=.5,add_phi=0,zero_it=0,digits = 4) {

    if (exists(".Random.seed", .GlobalEnv))
        oldseed <- .GlobalEnv$.Random.seed else oldseed <- NULL
        set.seed(rseed)
    set.seed <- rseed

    oo <- options()       # code line i
    on.exit(options(oo))  # code line i+1
    # old_digits <- getOption("digits")
    options(digits = digits)

    PHI <- matrix(phi, K, K)
    diag(PHI) <- 1
    if (ph1 > -1 && ph1 < 1) PHI[1,2] <- PHI[2,1]<-ph1

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
    K1<-K
    if(!is.null(add_ind)){
        # len<-length(add_ind)
        K1<-K+1
        mla<-cbind(mla,0)
        mla[add_ind,K1]<-add_la
        PHI<-cbind(PHI,add_phi)
        PHI<-rbind(PHI,add_phi)
        PHI[K1,K1]<-1
    }

    if (zero_it>0){
        J <- J + zero_it
        mla <- rbind(mla,matrix(0,zero_it,K1))
    }

    ecm <- matrix(0, J, J)
    if (ecr > 0) {
        iecb <- iecw <- NULL
        li <- c(1:K) * ipf - cpf
        for (i in 1:necw) {
            if(necw>0) ecm[li[i], li[i] - 1] <- ecm[li[i] - 1,li[i]] <- ecr
        }
        for (i in 1:necb) {
            i1 <- i%%K + 1
            if(necb>0) ecm[li[i] - 2, li[i1] - 3] <- ecm[li[i1] - 3,li[i] - 2]<- ecr
        }
    }

    # iecb <- iecw <- NULL
    # li <- c(1:K) * ipf - cpf
    # for (i in 1:necw) {
    #     if(necw>0) iecw <- rbind(iecw, c(li[i], li[i] - 1))
    #     # i1 <- i%%K + 1
    #     # iecb <- rbind(iecb, c(li[i] - 2, li[i1] - 3))
    # }
    # for (i in 1:necb) {
    #     # iecw <- rbind(iecw, c(li[i], li[i] - 1))
    #     i1 <- i%%K + 1
    #     if(necb>0) iecb <- rbind(iecb, c(li[i] - 2, li[i1] - 3))
    # }

    # ecm <- matrix(0, J, J)
    #
    # if (ecr > 0) {
    #     ecm[iecw] <- ecm[iecw[, 2:1]] <- ecr
    #     ecm[iecb] <- ecm[iecb[, 2:1]] <- ecr
    # }

    evr <- rep(0, J)
    for (j in 1:J) {
        evr[j] = 1 - t(mla[j, ]) %*% PHI %*% mla[j, ]
    }  #end j

    diag(ecm) <- evr
    er <- mvrnorm(N, rep(0, J), ecm)
    eta <- mvrnorm(N, rep(0, K1), PHI)
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

    if (!is.null(oldseed))
        .GlobalEnv$.Random.seed <- oldseed else rm(".Random.seed", envir = .GlobalEnv)

    return(out)
}
