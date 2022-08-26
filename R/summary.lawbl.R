#' @title Summary method for \code{lawbl} objects
#'
#' @description Provide summaries of posterior information for a \code{lawbl} object, .
#'
#' @name summary.lawbl
#'
#' @param object A \code{lawbl} object
#'
#' @param what A list of options for what to summarize.
#'
#' \itemize{
#'     \item \code{basic}: Basic information about the model and posteriors.
#'     \item \code{lambda}: Loading estimates.
#'     \item \code{qlambda}: Loading estimates in pattern/Q-matrix format.
#'     \item \code{eigen}: Factorial eigen value.
#'     \item \code{dpsx}: Diagonal elements in the residual covariance matrix \code{PSX}.
#'     \item \code{offpsx}: Off-diagonal elements in \code{PSX}; local dependence terms.
#'     \item \code{phi}: Factorial correlations.
#'     \item \code{thd}: Threshold estimates.
#'     \item \code{int}: Intercept estimates (for  \code{\link{pcirm}} only).
#'     \item \code{shrink}: (Ave) shrinkage for each factor's loadings and LD (if \code{LD} in \code{pcfa} = T).
#'     \item \code{factor}: Are the factors true of not (for EFA).
#'     \item \code{all}: All above information.
#'  }
#'
#' @param med logical; if the posterior median (\code{TRUE}) or mean (\code{FALSE}) is used as the estimate.
#'
#' @param SL Significance level for interval estimate. The default is .05.
#'
#' @param detail logical; if only significant (\code{FALSE}) or all (\code{TRUE}) estimates are presented.
#'
#' @param digits Number of significant digits to print when printing numeric values.
#'
#' @param istart Starting point of the Markov chain for summary.
#'
#' @param iend Ending point of the Markov chain for summary; -1 for the actual final point.
#'
#' @param ... additional arguments
#'
#' @return A list or matrix containing the summarized information based on the option \code{what}.
#'
#' @export
#'
#' @examples
#' \donttest{
#' dat <- sim18cfa0$dat
#' J <- ncol(dat)
#' K <- 3
#' Q<-matrix(-1,J,K);
#' Q[1:2,1]<-Q[7:8,2]<-Q[13:14,3]<-1
#'
#' m0 <- pcfa(dat = dat, Q = Q, LD = FALSE,burn = 1000, iter = 1000)
#' summary(m0) # summarize basic information
#' summary(m0, what = 'lambda') #summarize significant loadings
#' summary(m0, what = 'qlambda') #summarize significant loadings in pattern/Q-matrix format
#' summary(m0, what = 'offpsx') #summarize significant LD terms
#' }

summary.lawbl <- function(object, what = "basic", med = FALSE, SL = 0.05, detail = FALSE, digits = 4, istart = 1, iend = -1, ...) {

    Q <- object$Q
    J <- nrow(Q)
    K <- ncol(Q)
    LD <- object$LD
    ELA <- object$LA
    iter <- object$iter
    N <- dim(object$Omega)[2]
    # TF_ind = object$TF_ind
    # if(is.null(TF_ind)) TF_ind <- rep(TRUE, K)

    # oo <- options()       # code line i
    # on.exit(options(oo))  # code line i+1
    # options(digits = digits)

    if (iend == -1 | iend > iter)
      iend <- iter
    ELA <- ELA[istart:iend,]

    eig_eps <- object$eig_eps
    if(is.null(eig_eps)) eig_eps <- .1
    TF_ind<-(colMeans(object$Eigen[istart:iend,])>eig_eps)

    # Loading med=TRUE;SL=.05
    tmp <- result(ELA, med, SL)
    nt <- dim(tmp)[2]
    sig <- tmp[, nt]
    MLA <- Q

    if (detail) {
        # all Loading
        ind <- which(Q != 0, arr.ind = TRUE)
        colnames(ind) <- c("Item", "F")
        LAM <- (cbind(ind, tmp))
        LAM <- LAM[ind[,2] %in% which(TF_ind),]
        MLA[Q != 0] <- tmp[, 1]
        # MLA  <- MLA[,TF_ind]
    } else {
        # Sig. Loading
        MLA[Q != 0] <- sig
        ind <- which(MLA > 0, arr.ind = TRUE)
        MLA[MLA > 0] <- tmp[sig > 0, 1]

        colnames(ind) <- c("Item", "F")
        LAM <- (cbind(ind, tmp[sig > 0, ]))
    }
    MLA  <- MLA[,TF_ind]
    colnames(MLA)<-which(TF_ind)
    NSLA <- sum(sig > 0)
    row.names(MLA) <- paste0("I", 1:J)

    # eigenvalue
    # poq = which(Q != 0, arr.ind = TRUE)
    # eig_arr <- array(0, dim = c(iter, K))
    # for (k in 1:K) {
    #     ind1 <- (poq[, 2] == k)
    #     eig_arr[, k] <- rowSums(ELA[, ind1]^2)
    # }
    # colnames(eig_arr) <- paste0("F", c(1:K))
    #
    # # if (detail){
    # eigen <- result(eig_arr, med, SL)
    # # }else{ eigen <- apply(mcmc(eig_arr), 2, mean) }

    eigen <- result(object$Eigen[istart:iend,], med, SL)
    row.names(eigen) <- paste0("F", c(1:K))

    PSX <- matrix(0,J,J)
    # Sig. PSX
    if (LD) {
        tmp <- result(object$PSX[istart:iend,], med, SL)
        pos <- lower.tri(PSX, diag = TRUE)
        ind <- which(pos, arr.ind = TRUE)
        PSX[ind]<- tmp[,1]
        tPSX <- t(PSX)
        PSX[upper.tri(PSX)]<-tPSX[upper.tri(tPSX)]

        dpsx <- tmp[ind[, 1] == ind[, 2], ]
        ofind <- ind[, 1] != ind[, 2]
        offpsx <- cbind(ind[ofind, ], tmp[ofind, ])
        nt <- dim(offpsx)[2]
        no_ofd <- sum(offpsx[, nt] > 0)
        if (!detail)
            offpsx <- offpsx[offpsx[, nt] > 0, ]
    } else {
        dpsx <- result(object$PSX[istart:iend,], med, SL)
        diag(PSX) <- dpsx[,1]
        offpsx <- NULL
    }
    row.names(dpsx) <- paste0("I", 1:J)

    # sign_chg<-object$chg_count
    # row.names(sign_chg) <- c("Burn-in", "Iteration")

    out0 <- list(NJK = c(N, J, K), `Miss%` = object$Nmis/J/N * 100, `LD Allowed` = LD,
                 `Burn in` = object$burn+istart-1,Iteration = iend-istart+1, `No. of sig lambda` = NSLA,
                Selected = TF_ind,'Auto, NCONV, MCONV'=object$auto_conv)


    if (!detail) eigen <- eigen[TF_ind,]

    # APSR = object$APSR
    # if(!is.null(nrow(APSR))) row.names(APSR) <- paste0("F", c(ind))
    EPSR <- schain.grd(object$Eigen[istart:iend,TF_ind])
    out0$EPSR <- round(EPSR, digits)

    if (K > 1){
      tmp <- result(object$PHI[istart:iend,], med, SL)
      pos <- lower.tri(matrix(0, K, K))
      ind0 <- which(pos, arr.ind = TRUE)
      tmp0 <- cbind(ind0, tmp)
      ind <- which(TF_ind)
      if (detail){
        phi<-tmp0

      }else{
      sind <- NULL
      for (i in 1:dim(tmp0)[1]){
        if(all(ind0[i,] %in% ind)) sind <- c(sind, i)
      }
      phi <- tmp0[sind,]
      }
    }else{
      phi<-1
    }

    Qb <- object$Qb
    if (is.null(Qb)){
      gammal <- result(object$gammal[istart:iend,], med, SL)
      row.names(gammal) <- paste0("F", 1:K)
      gammal <- round(gammal,digits)
      gammab<-coef<-qcoef<-coef.er<-NULL
    }else{
      qcoef <- Qb
      tmp <- result(object$B[istart:iend,], med, SL)
      nt <- dim(tmp)[2]
      sig <- tmp[, nt]
      EPSRb <- schain.grd(object$R2[istart:iend,])
      out0$EPSRb <- round(EPSRb, digits)
      out0$`No. of sig coef` <- sum(sig > 0)

      if (detail) {
        # all Loading
        ind <- which(Qb != 0, arr.ind = TRUE)
        colnames(ind) <- c("X", "F")
        coef <- (cbind(ind, tmp))
        qcoef[Qb != 0] <- tmp[, 1]
      } else {
        # Sig. Loading
        qcoef[Qb != 0] <- sig
        ind <- which(qcoef > 0, arr.ind = TRUE)
        colnames(ind) <- c("X", "F")
        coef <- (cbind(ind, tmp[sig > 0, ]))
        qcoef[qcoef > 0] <- tmp[sig > 0, 1]
      }
      coef <- round(coef,digits)
      rownames(coef) <- NULL
      qcoef <- round(qcoef,digits)

      coef.er<-result(object$PSXb[istart:iend,],med,SL)
      row.names(coef.er) <- paste0("F", 1:ncol(Qb))
      coef.er <- round(coef.er,digits)

      R2<-result(object$R2[istart:iend,],med,SL)
      row.names(R2) <- paste0("F", 1:ncol(Qb))
      R2 <- round(R2,digits)
      if (sum(Qb == -1)>0){
        gammab <- result(object$gammab[istart:iend,], med, SL)
        # ind <- which(Qb == -1, arr.ind = TRUE)
        # colnames(ind) <- c("F", "X")
        # gammab <- (cbind(ind, tmp))
        row.names(gammab) <- paste0("F", 1:K)
        gammab <- round(gammab,digits)
      }else{
        gammab <- NULL
      }
      gammal<-NULL

    }

    if (LD){
        out0$"No. of sig LD terms" = no_ofd
        # tgam <- cbind(object$gammal, object$gammas)
        gammas <- result(as.matrix(object$gammas[istart:iend], med, SL))
        # row.names(gammas) <- c(paste0("F", 1:K), "PSX")
        gammas <- round(gammas,digits)
        rownames(gammas) <- NULL
    }else{
      gammas <- NULL
    }
    # else{
    #     tgam <- (object$gammal)
    #     allgam <- result(tgam, med, SL)
    #     row.names(allgam) <- paste0("F", 1:K)
    # }

    Jp <- length(object$cati)
    if (Jp > 0) {
        nthd <- object$mnoc - 1
        # if (nthd == 1) {
        #     Mthd <- result(object$THD[, , 1], med, SL)
        #     row.names(Mthd) <- paste0("I", 1:Jp)
        # } else {

            if (!detail) {
                Mthd <- matrix(0, Jp, nthd)
                row.names(Mthd) <- paste0("I", 1:Jp)
            } else {
                Mthd <- NULL
            }
            for (thd in 1:(nthd)) {
                if (!detail) {
                  Mthd[, thd] <- colMeans(object$THD[istart:iend, , thd])  # mean estimates only

                } else {
                  tmpt <- cbind(thd, result(object$THD[istart:iend, , thd], med, SL))
                  Mthd <- rbind(Mthd, tmpt)
                }
            }
        # } # end nthd
        out0$"Cat Items" <- object$cati
        out0$"max No. of categories" = nthd + 1
        Mthd <- round(Mthd,digits)
        rownames(Mthd) <- NULL
        # out$THD = Mthd
    } else {
        Mthd <- NULL
    }

    if (!is.null(object$PPP))
        out0$PPP <- mean(object$PPP[istart:iend])

    if(!is.null(object$MU)){
        MU <- result(object$MU[istart:iend,], med, SL)
        MU <- round(MU,digits)
        row.names(MU) <- paste0("I", 1:J)
    }else{
        MU <- NULL
    }

    lpry<- object$lpry
    if (!is.null(lpry)){
      Yc <- object$Y - MLA %*% object$Omega[TF_ind,]  # J*N
      tmp<-(t(Yc) %*%chol(chol2inv(chol(PSX))))^2
      lhat <- sum(tmp)+N*(log(det(PSX))+log(2*pi))
      # D_hat <- sum(tmp)+N*(log(det(PSX)))
      # out0$DIC <- 2*lsum-lhat

      npar<-sum(Q!=0)+K*(K-1)/2+J+LD*J*(J-1)/2
      DIC <- lpry + 2*npar
      BIC <- lhat + npar*log(N)
      AIC <- lhat +2*npar
      out0$"DIC, BIC, AIC" <- c(DIC, BIC, AIC)
    }
    out0$Time <- object$time

    LAM <- round(LAM,digits)
    rownames(LAM) <- NULL
    MLA <- round(MLA,digits)
    eigen <- round(eigen,digits)
    dpsx <- round(dpsx,digits)
    if (!is.null(offpsx)) {
      offpsx <- round(offpsx,digits)
      rownames(offpsx) <- NULL
    }
    phi <- round(phi,digits)
    rownames(phi) <- NULL
    if (!is.null(gammal)) gammal <- round(gammal,digits)
    if (!is.null(gammas)) gammas <- round(gammas,digits)

    out <- switch(what, basic = out0, lambda = LAM, qlambda = MLA, eigen = eigen,
                  dpsx = dpsx, offpsx = offpsx,phi = phi, gammal = gammal,gammas = gammas,
                  thd = Mthd, int = MU, factor = TF_ind,R2=R2,
                  coef = coef, qcoef=qcoef,coef.er = coef.er, gammab = gammab, all = {
            out1 <- out0
            out1$lambda <- LAM
            out1$qlambda <- MLA
            out1$eigen <- eigen
            out1$dpsx <- dpsx
            out1$offpsx <- offpsx
            out1$phi <- phi
            out1$gammal <- gammal
            out1$gammas <- gammas
            out1$thd <- Mthd
            out1$int <- MU
            # out1$factor <- TF_ind
            out1$coef <- coef
            out1$coef.er <- coef.er
            out1$gammab <- gammab
            out1
        }, stop(sprintf("Cannot show element '%s'", what), call. = FALSE))

    # options(digits = old_digits)
    return(out)
}

######## end of Summary #################################################
