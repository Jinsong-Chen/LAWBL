#' @title Simulating data with Latent Variable Modeling
#'
#' @description \code{sim_lvm} can simulate data with continuous latent variables (factors)
#' and continuous or categorical observed variables, plus a MIMIC-type structure.
#' One can also include an error covariance (local dependence) structure. Categorical
#' observed variables are generated with latent continuous responses normally distributed
#' and equally spaced within [-3,3].
#'
#' @name sim_lvm
#'
#' @param N Sample size.
#'
#' @param lam Loading value (for major loadings) or matrix (\eqn{J \times K}).
#'
#' @param K Number of factors (if \code{lam} is a value).
#'
#' @param J Number of items (if \code{lam} is a value).
#'
#' @param cpf Number of cross-loadings per factor (if \code{lam} is a value).
#'
#' @param lac Cross-loading value (if \code{lam} is a value).
#'
#' @param phi Factor correlation scalar or matrix, or error correlations (for MIMIC-type model).
#'
#' @param ecr Error covariance (local dependence) value.
#'
#' @param fac_score Output factor score or not.
#'
#' @param P Number of observable predictors (for MIMIC-type model).
#'
#' @param phix Observable predictor correlation value or matrix (for MIMIC-type model).
#'
#' @param b Coefficients of observable predictors (for MIMIC-type model), value or \eqn{K \times P}.
#'
#' @param lam1 Loading value (for major loadings) or matrix (\eqn{J1 \times K1}) for latent predictors (for MIMIC-type model).
#'
#' @param K1 Number of latent predictors (if \code{lam} is a value, for MIMIC-type model).
#'
#' @param J1 Number of items latent predictors (if \code{lam} is a value, for MIMIC-type model).
#'
#' @param phi1 Latent predictor correlation scalar or matrix (for MIMIC-type model).
#'
#' @param b1 Coefficients of latent predictors (for MIMIC-type model), value or \eqn{K \times K1}
#'
#' @param ilvl Specified levels of all items (i.e., need to specify Item 1 to \eqn{J+P+J1});
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
#' @param digits Number of significant digits to print when printing numeric values.
#'
#' @return An object of class \code{list} containing the data, loadings, factor correlations,
#' local dependence, and other information. The data consists of J items for the factors,
#' P items for observable predictors, and J1 items for latent predictors.
#'
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
#' out$lam
#' out$loc_dep
#'
#' # for categorical data with cross-loadings .4 and 10% missingness
#' out <- sim_lvm(N=1000,K=3,J=18,lam = .7, lac=.4,cati=-1,noc=4,misp=.1)
#' summary(out$dat)
#' out$lam
#' out$loc_dep
#'
sim_lvm <- function(N = 1000, lam = 0.7, K = 3, J = 18, cpf = 0, lac = 0.3, phi = .3,ecr = .0,
                    necw=K, necb=K, P = 0, phix = 0, b = 0, lam1 = 0, K1 = 0, J1 = 0, b1 = 0, phi1 = 0,
       ilvl=NULL, cati = NULL, noc = c(4), misp = 0, fac_score = FALSE, rseed = 333,digits = 4) {

    if (is.scalar(lam) && (J%%K != 0 || J<1) )
        stop("J should be a multiple of K when lam is a value.", call. = FALSE)

    if (exists(".Random.seed", .GlobalEnv))
        oldseed <- .GlobalEnv$.Random.seed else oldseed <- NULL
        set.seed(rseed)
    set.seed <- rseed

    oo <- options()       # code line i
    on.exit(options(oo))  # code line i+1
    # old_digits <- getOption("digits")
    options(digits = digits)

    mla <- function(la0,J,K,lac=0,cpf=0){
      ipf <- round(J / K)
      if (K > 1){
        lam0 <- c(rep(la0, ipf), rep(0, ipf - cpf), rep(lac, cpf), rep(0, (K - 2) * ipf))
      }else{
        lam0<-la0
      }
      mla1 <- matrix(lam0, ipf, K)
      if(K == 1){
        mla <- mla1
      }else{
        mla <- c()
        for (i in K:1) {
          # i <- 1
          ind <- (c(1:K) + i)%%K
          ind[ind == 0] <- K
          mla <- rbind(mla, mla1[, ind])
        }
      }
      return(mla)
     }


    if (!is.scalar(lam1)) {
      lam1<-as.matrix(lam1)
      K1 <- ncol(lam1)
      J1 <- nrow(lam1)
    }else{
      if (K1 > 0){
        if (J1%%K1 != 0  || J1 < 1)
          stop("J1 should be a multiple of K1 when lam1 is a value.", call. = FALSE)
        lam1<-mla(lam1,J1,K1)
      }else{
        # lam1 <- NULL
        J1 <- 0
      }
    }

    PHI <- fac <- 0
    x <- y1 <- fac1 <- eigen1 <- mb <- mb1 <- evr1<- NULL
    if (P>0){
      PHIX <- cm_check(phix,P)
      if(any(PHIX<=-1))
        stop("phix should be a correlation scalar or matrix (P*P & PD).", call. = FALSE)

      x <- mvrnorm(N, rep(0, P), PHIX)
      if (is.scalar(b)){
        mb <- matrix(b,K,P)
      }else{
        mb <- b #K x P
      }

      PHI <- PHI + mb %*% PHIX %*% t(mb)
      fac <- fac + t(mb %*% t(x))
    }

   if (K1 > 0){
     PHI1 <- cm_check(phi1,K1)
     if(any(PHI1<=-1))
       stop("phi1 should be a correlation scalar or matrix (K1*K1 & PD).", call. = FALSE)

     fac1 <- mvrnorm(N, rep(0, K1), PHI1)
     evr1<-1 - diag(lam1 %*% PHI1 %*% t(lam1))
     err1 <- mvrnorm(N, rep(0, J1), diag(evr1))
     y1 <- t(lam1 %*% t(fac1)) + err1
     eigen1 <- diag(crossprod(lam1))

     if (is.scalar(b1)){
       mb1 <- matrix(b1,K,K1)
     }else{
       mb1 <- b1
     }
     PHI <- PHI + mb1 %*% PHI1 %*% t(mb1)
     fac <- fac + t(mb1 %*% t(fac1))
   }

    if (!is.scalar(lam)) {
      lam<-as.matrix(lam)
      K <- ncol(lam)
      J <- nrow(lam)
      cpf <- 0
    }else{
      lam <- mla(lam,J,K,lac,cpf)
    }

    PHI0 <- cm_check(phi,K)
    if(any(PHI0<=-1))
      stop("phi should be a correlation scalar or matrix (K*K & PD).", call. = FALSE)

   if (K > 1 && (P+K1)>0)
     diag(PHI0) <- 1 - diag(PHI) # factor variance standardized

    fac <- fac + mvrnorm(N, rep(0, K), PHI0)
    PHI <- PHI + PHI0

    ipf <- round(J / K)
    ecm <- matrix(0, J, J)
    if (ecr > 0) {
        # iecb <- iecw <- NULL
        li <- c(1:K) * ipf - cpf
        for (i in 1:necw) {
            if(necw>0) ecm[li[i], li[i] - 1] <- ecm[li[i] - 1,li[i]] <- ecr
        }
        for (i in 1:necb) {
            i1 <- i%%K + 1
            if(necb>0) ecm[li[i] - 2, li[i1] - 3] <- ecm[li[i1] - 3,li[i] - 2]<- ecr
        }
    }

    # PHI <- cor(eta)
    evr0<-1 - diag(lam %*% PHI %*% t(lam))
    diag(ecm) <- evr0
    err <- mvrnorm(N, rep(0, J), ecm)
    y0 <- t(lam %*% t(fac)) + err

    y<-cbind(y0,x,y1)

    scale <- apply(y, 2, sd)
    eigen <- c(diag(crossprod(lam)),eigen1)
    evr <- c(evr0,evr1)

    out <- list(lam=lam, PHI = PHI, Eigen = eigen, scale = scale, var_ie=evr,PSX = ecm)
    if (fac_score)
        out$fac <- cbind(fac,fac1)

    pos <- lower.tri(ecm)
    ind <- which(pos, arr.ind = T)
    rind <- which(ecm[pos] > 0)
    out$loc_dep <- ind[rind, ]

    JA <- J + P + J1
    UL <- 3;LL <- -3
    if (is.null(ilvl)) {
        Jp <- length(cati)
        if (Jp > 0)
            {
                if (Jp == 1 && cati == -1) {
                    cati <- c(1:JA)
                    Jp <- JA
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
        if (length(ilvl)!=JA)
            stop("Number of item levels should be J + P + J1.", call. = F)

        for (i in 1:(JA)) {
            M <- ilvl[i]  #categories

            if (M >=2) {
                M <- round(M)
                stp <- LL + c(1:(M - 1)) * (UL - LL)/M  #M-1 step points
                tmp <- y[, i] > matrix(stp, N, (M - 1), byrow = T)
                y[, i] <- rowSums(tmp) + 1  # value starting from 1
            }
        }  #end J+P
    } #end ilvl

  if (misp > 0){
    mind <- matrix(rbinom(N * J, 1, misp), N, J)
    y0<-y[,1:J]
    y0[mind == 1] <- NA
    y[,1:J] <- y0
  }
  colnames(y)<-paste0("y",c(1:ncol(y)))
  out$dat <- y

    if (P != 0 || K1 != 0){
        out$mb <- mb
        out$mb1 <- mb1
        out$var_fe <- PHI0
        out$PHIX <- PHIX
        out$lam1 <- lam1
    }

    if (!is.null(oldseed))
        .GlobalEnv$.Random.seed <- oldseed else rm(".Random.seed", envir = .GlobalEnv)

    return(out)
}
