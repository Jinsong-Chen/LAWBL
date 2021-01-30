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

summary.lawbl <- function(object, what = "basic", med = FALSE, SL = 0.05, detail = FALSE, digits = 4, ...) {

    Q <- object$Q
    J <- nrow(Q)
    K <- ncol(Q)
    LD <- object$LD
    ELA <- object$LA
    iter <- dim(ELA)[1]
    N <- dim(object$Omega)[2]

    oo <- options()       # code line i
    on.exit(options(oo))  # code line i+1

    # old_digits <- getOption("digits")
    options(digits = digits)

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
        MLA[Q != 0] <- tmp[, 1]
    } else {
        # Sig. Loading
        MLA[Q != 0] <- sig
        ind <- which(MLA > 0, arr.ind = TRUE)
        colnames(ind) <- c("Item", "F")
        LAM <- (cbind(ind, tmp[sig > 0, ]))
        MLA[MLA > 0] <- tmp[sig > 0, 1]
    }

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

    eigen <- result(object$Eigen, med, SL)
    row.names(eigen) <- paste0("F", c(1:K))

    # Sig. PSX
    if (LD) {
        tmp <- result(object$PSX, med, SL)
        pos <- lower.tri(matrix(0, J, J), diag = TRUE)
        ind <- which(pos, arr.ind = TRUE)
        dpsx <- tmp[ind[, 1] == ind[, 2], ]
        ofind <- ind[, 1] != ind[, 2]
        offpsx <- cbind(ind[ofind, ], tmp[ofind, ])
        nt <- dim(offpsx)[2]
        no_ofd <- sum(offpsx[, nt] > 0)
        if (!detail)
            offpsx <- offpsx[offpsx[, nt] > 0, ]
    } else {
        dpsx <- result(object$PSX, med, SL)
        offpsx <- NULL
    }
    row.names(dpsx) <- paste0("I", 1:J)

    sign_chg<-object$chg_count
    row.names(sign_chg) <- c("Burn-in", "Iteration")
    out0 <- list(N = N, J = J, K = K, `Miss%` = object$Nmis/J/N * 100, `LD enabled` = LD, `Burn in` = object$burn,
                 Iteration = object$iter, 'Sign Change'= sign_chg, `No. of sig lambda` = NSLA)

    TF_ind = object$TF_ind
    if(is.null(TF_ind)) TF_ind <- rep(TRUE, K)
    eigen <- eigen[TF_ind,]
    KE <- sum(TF_ind)
    out0$'True Factor' = TF_ind
    ind <- which(TF_ind)
    APSR = object$APSR
    row.names(APSR) <- paste0("F", c(ind))
    out0$'Adj. PSR' = APSR

    tmp <- result(object$PHI, med, SL)
    pos <- lower.tri(matrix(0, K, K))
    ind0 <- which(pos, arr.ind = TRUE)
    tmp0 <- cbind(ind0, tmp)

    sind <- NULL
    for (i in 1:dim(tmp0)[1]){
        if(all(ind0[i,] %in% ind)) sind <- c(sind, i)
    }
    phi <- tmp0[sind,]

    if (LD){
        out0$"No. of sig LD terms" = no_ofd
        tgam <- cbind(object$gammal, object$gammas)
        allgam <- result(tgam, med, SL)
        row.names(allgam) <- c(paste0("F", 1:K), "PSX")
    }else{
        tgam <- (object$gammal)
        allgam <- result(tgam, med, SL)
        row.names(allgam) <- paste0("F", 1:K)
    }

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
                  Mthd[, thd] <- colMeans(object$THD[, , thd])  # mean estimates only

                } else {
                  tmpt <- cbind(thd, result(object$THD[, , thd], med, SL))
                  Mthd <- rbind(Mthd, tmpt)
                }
            }
        # } # end nthd
        out0$"Cat Items" <- object$cati
        out0$"max No. of categories" = nthd + 1
        # out$THD = Mthd
    } else {
        Mthd <- NULL
    }

    if (!is.null(object$PPP))
        out0$PPP <- mean(object$PPP)

    if(!is.null(object$MU)){
        MU <- result(object$MU, med, SL)
        row.names(MU) <- paste0("I", 1:J)
    }else{
        MU <- NULL
    }

    out <- switch(what, basic = out0, lambda = LAM, qlambda = MLA, eigen = eigen, dpsx = dpsx, offpsx = offpsx,
        phi = phi, shrink = allgam, thd = Mthd, int = MU, factor = TF_ind, all = {
            out1 <- out0
            out1$lambda <- LAM
            out1$qlambda <- MLA
            out1$eigen <- eigen
            out1$dpsx <- dpsx
            out1$offpsx <- offpsx
            out1$phi <- phi
            out1$shrink <- allgam
            out1$thd <- Mthd
            out1$int <- MU
            out1$factor <- TF_ind
            out1
        }, stop(sprintf("Can not show element '%s'", what), call. = FALSE))

    # options(digits = old_digits)
    return(out)
}

######## end of Summary #################################################
