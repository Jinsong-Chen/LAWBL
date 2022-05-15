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
#' @param mla Population loading matrix.
#'
#' @param K Number of factors (if \code{mla=NULL}).
#'
#' @param J Number of items (if \code{mla=NULL}).
#'
#' @param cpf Number of cross-loadings per factor (if \code{mla=NULL}).
#'
#' @param lam Number of formal iterations for posterior sampling.
#'
#' @param lac Number of iterations to update the sampling information.
#'
#' @param phi Homogeneous correlations between any two factors.
#'
#' @param ph12 Correlation between factor 1 and 2 (if it's different from \code{phi}.
#'
#' @param ecr Residual correlation (local dependence).
#'
#' @param ome_out Output factor score or not.
#'
#' @param P Number of observable predictors (for MIMIC model).
#'
#' @param b Coefficients of observable predictors (for MIMIC model).
#'
#' @param K1 Number of latent predictors (for MIMIC model).
#'
#' @param ph1 Correlations between latent predictors (for MIMIC model).
#'
#' @param b1 Coefficients of latent predictors (for MIMIC model).
#'
#' @param ilvl Specified levels of all items (i.e., need to specify Item 1 to \eqn{J+P});
#'  Any value smaller than 2 is considered as continuous item.
#'
#' @param cati The set of polytomous items in sequence number (i.e., can be any number set
#'  in between 1 and \eqn{J+P}); \code{NULL} for no and -1 for all (if \code{ilvl=NULL}).
#'
#' @param noc Number of levels for polytomous items.
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
#' out <- sim_lvm(N=1000,K=3,J=18,lam = .7, lac=.3,ecr=.3)
#' summary(out$dat)
#' out$MLA
#' out$ofd_ind
#'
#' # for categorical data with cross-loadings .4 and 10% missingness
#' out <- sim_lvm(N=1000,K=3,J=18,lam = .7, lac=.4,cati=-1,noc=4,misp=.1)
#' summary(out$dat)
#' out$MLA
#' out$ofd_ind
#'
sim_lvm <- function(N = 1000, mla=NULL, K = 3, J = 18, cpf = 0, lam = 0.7, lac = 0.3, phi = .3, ph12 = -1,
        ecr = .0,P = 0, b = .3,K1 = 0,ph1 = .2, b1 = .3, ilvl=NULL,cati = NULL, noc = c(4), misp = 0, ome_out = FALSE,
        necw=K,necb=K,add_ind=c(),add_la=.5,add_phi=0,zero_it=0, rseed = 333,digits = 4) {

    if (is.null(mla) && J%%K != 0)
        stop("J should be a multiple of K.", call. = FALSE)

    if (exists(".Random.seed", .GlobalEnv))
        oldseed <- .GlobalEnv$.Random.seed else oldseed <- NULL
        set.seed(rseed)
    set.seed <- rseed

    oo <- options()       # code line i
    on.exit(options(oo))  # code line i+1
    # old_digits <- getOption("digits")
    options(digits = digits)

    if (!is.null(mla)) {
        mla<-as.matrix(mla)
        K <- ncol(mla)
        J <- nrow(mla)
        cpf <- 0
    }


    if (P == 0 && K1 == 0){
        PHI <- matrix(phi, K, K)
        diag(PHI) <- 1
        if (ph12 > -1 && ph12 < 1) PHI[1,2] <- PHI[2,1]<-ph12

        PHX <- NULL
        eta <- mvrnorm(N, rep(0, K), PHI)
    }else{
        mb <- mb1 <- PHX <- NULL
        tmp <- eta <- 0
        if (P>0){
             PHX <- matrix(phi, P, P)
            diag(PHX) <- 1
            x <- mvrnorm(N, rep(0, P), PHX)
            if (length(b)==1){
                mb <- matrix(b,K-K1,P)
            }else{
                mb <- b
            }
             tmp <- mb %*% PHX %*% t(mb)
             eta <- t(mb %*% t(x))
        }


        if (K1 > 0){
            PH1 <- matrix(ph1, K1, K1)
            diag(PH1) <- 1
            ome <- mvrnorm(N, rep(0, K1), PH1)
            if (length(b1)==1){
                mb1 <- matrix(b1,K-K1,K1)
            }else{
                mb1 <- b1
            }
            tmp <- tmp + mb1 %*% PH1 %*% t(mb1)
            eta <- eta + t(mb1 %*% t(ome))
        }


        ber = 1 - diag(tmp) # factor variance standardized
        if (K > 1) ber <- diag(ber)
        mber <- mvrnorm(N, rep(0, K-K1), ber)
        eta <- eta + mber
        # diag(tmp)<-1

        if (K1 > 0){
            # eta <- eta + t(mb1 %*% t(ome))
            # PHI <- matrix(0,K,K)
            # PHI[1:(K-K1),1:(K-K1)]<-tmp
            # diag(PHI)<-1
            eta <- cbind(eta,ome)
        }

    } #end P == 0 && K1 == 0

    # J <- K * ipf
    ipf <- J / K

    if (is.null(mla)) {
            if (K > 1){
                lam0 <- c(rep(lam, ipf), rep(0, ipf - cpf), rep(lac, cpf), rep(0, (K - 2) * ipf))
            }else{
                lam0<-lam
            }

            lam1 <- matrix(lam0, ipf, K)

            if(K == 1){
                mla <- lam1
            }else{
                mla <- c()
                for (i in K:1) {
                    # i <- 1
                    ind <- (c(1:K) + i)%%K
                    ind[ind == 0] <- K
                    mla <- rbind(mla, lam1[, ind])
                }
            }

            if (P == 0 && K1 == 0){
                Kp1<-K
                if(!is.null(add_ind)){
                    # len<-length(add_ind)
                    Kp1<-K+1
                    mla<-cbind(mla,0)
                    mla[add_ind,Kp1]<-add_la
                    PHI<-cbind(PHI,add_phi)
                    PHI<-rbind(PHI,add_phi)
                    PHI[Kp1,Kp1]<-1
                }
                if (zero_it>0){
                    J <- J + zero_it
                    mla <- rbind(mla,matrix(0,zero_it,Kp1))
                }
                eta <- mvrnorm(N, rep(0, Kp1), PHI)
            }
        } #end mla


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

    PHI <- cor(eta)
    evr<-1 - diag(mla %*% PHI %*% t(mla))

    # evr <- rep(0, J)
    # for (j in 1:J) {
    #     evr[j] = 1 - t(mla[j, ]) %*% PHI %*% mla[j, ]
    # }  #end j

    diag(ecm) <- evr
    er <- mvrnorm(N, rep(0, J), ecm)

    # eta <- mvrnorm(N, rep(0, K1), PHI)

    y <- t(mla %*% t(eta)) + er
    if ( P >0){
        y <- cbind(y,x)
    }

    scale <- apply(y, 2, sd)
    eigen <- diag(crossprod(mla))

    out <- list(N = N, PHI = PHI, MLA = mla, Eigen = eigen, PSX = ecm, scale = scale)
    if (ome_out)
        out$OME <- eta

    pos <- lower.tri(ecm)
    ind <- which(pos, arr.ind = T)
    rind <- which(ecm[pos] > 0)
    out$ofd_ind <- ind[rind, ]

    UL <- 3;LL <- -3
    if (is.null(ilvl)) {
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

    }else{
        if (length(ilvl)!=J+P)
            stop("number of item levels should be J + P.", call. = F)

        for (i in 1:(J+P)) {
            M <- ilvl[i]  #categories

            if (M >=2) {
                M <- round(M)
                stp <- LL + c(1:(M - 1)) * (UL - LL)/M  #M-1 step points
                tmp <- y[, i] > matrix(stp, N, (M - 1), byrow = T)
                y[, i] <- rowSums(tmp) + 1  # value starting from 1
            }
        }  #end J+P
    } #end ilvl


    mind <- matrix(rbinom(N * J, 1, misp), N, J)
    y[mind == 1] <- NA
    out$misp <- misp
    out$dat <- y

    if (P != 0 || K1 != 0){
        out$mb <- mb
        out$mb1 <- mb1
        out$PHX <- PHX
    }


    if (!is.null(oldseed))
        .GlobalEnv$.Random.seed <- oldseed else rm(".Random.seed", envir = .GlobalEnv)

    return(out)
}
