#' @title Summary method for \code{pcfa} objects
#'
#' @description Provide basic information for an \code{PCFA} object, and summarize various posteriors.
#'
#' @name summary.pcfa
#'
#' @param obj A \code{pcfa} object
#'
#' @param what A list of options for what to summarize.
#'
#' \itemize{
#'      \item \code{basic}: Basic information about the model and posteriors.
#'     \item \code{lambda}: Loading estimates.
#'     \item \code{qlambda}: Loading estimates in pattern/Q-matrix format.
#'     \item \code{eigen}: Factorial eigen value.
#'     \item \code{dpsx}: Diagonal elements in the residual covariance matrix \code{PSX}.
#'     \item \code{offpsx}: Off-diagonal elements in \code{PSX}; local dependence terms.
#'     \item \code{phi}: Factorial correlations.
#'     \item \code{shrink}: Shrinkage parameters for the loadings and LD (if \code{LD} in \code{pcfa} = T).
#'     \item \code{all}: All above information.
#'  }
#'
#' @param med logical; if the posterior median (\code{T}) or mean (\code{F}) is used as the estimate.
#'
#' @param SL Significance level for interval estimate. The default is .05.
#'
#' @param details logical; if only significant (\code{F}) or all (\code{T}) estimates are presented.
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
#' \dontrun{
#' dat <- sim18cfa0$dat
#' J <- ncol(dat)
#' K <- 3
#' Q<-matrix(-1,J,K);
#' Q[1:2,1]<-Q[9:10,2]<-Q[13:14,3]<-1
#'
#' mod0 <- pcfa(dat = dat, Q = Q, LD = F,burn = 2000, iter = 2000)
#' summary(mod0) # summarize basic information
#' summary(mod0, what = 'lambda') #summarize significant loadings
#' summary(mod0, what = 'qlambda') #summarize significant loadings in pattern/Q-matrix format
#' summary(mod0, what = 'offpsx') #summarize significant LD terms
#' summary(mod0, what = 'all') #summarize all information
#' }

summary.pcfa <- function(obj, what = "basic", med = F, SL = 0.05, detail = F, digits = 4, ...) {
    
    Q <- obj$Q
    J <- nrow(Q)
    K <- ncol(Q)
    LD <- obj$LD
    ELA <- obj$LA
    iter <- dim(ELA)[1]
    N <- dim(obj$Omega)[2]
    
    old_digits <- getOption("digits")
    options(digits = digits)
    
    # Loading med=T;SL=.05
    tmp <- result(ELA, med, SL)
    nt <- dim(tmp)[2]
    sig <- tmp[, nt]
    MLA <- Q
    
    if (detail) {
        # all Loading
        ind <- which(Q != 0, arr.ind = T)
        colnames(ind) <- c("Item", "F")
        LAM <- (cbind(ind, tmp))
        MLA[Q != 0] <- tmp[, 1]
    } else {
        # Sig. Loading
        MLA[Q != 0] <- sig
        ind <- which(MLA > 0, arr.ind = T)
        colnames(ind) <- c("Item", "F")
        LAM <- (cbind(ind, tmp[sig > 0, ]))
        MLA[MLA > 0] <- tmp[sig > 0, 1]
    }
    
    NSLA <- sum(sig > 0)
    row.names(MLA) <- paste0("I", 1:J)
    
    # eigenvalue
    poq = which(Q != 0, arr.ind = T)
    eig_arr <- array(0, dim = c(iter, K))
    for (k in 1:K) {
        ind1 <- (poq[, 2] == k)
        eig_arr[, k] <- rowSums(ELA[, ind1]^2)
    }
    colnames(eig_arr) <- paste0("F", c(1:K))
    
    # if (detail){
    eigen <- result(eig_arr, med, SL)
    # }else{ eigen <- apply(mcmc(eig_arr), 2, mean) }
    
    # Sig. PSX
    if (LD) {
        tmp <- result(obj$PSX, med, SL)
        pos <- lower.tri(matrix(0, J, J), diag = T)
        ind <- which(pos, arr.ind = T)
        dpsx <- tmp[ind[, 1] == ind[, 2], ]
        ofind <- ind[, 1] != ind[, 2]
        offpsx <- cbind(ind[ofind, ], tmp[ofind, ])
        nt <- dim(offpsx)[2]
        no_ofd <- sum(offpsx[, nt] > 0)
        if (!detail) 
            offpsx <- offpsx[offpsx[, nt] > 0, ]
    } else {
        dpsx <- result(obj$PSX, med, SL)
        offpsx <- NULL
    }
    row.names(dpsx) <- paste0("I", 1:J)
    
    tmp <- result(obj$PHI, med, SL)
    pos <- lower.tri(matrix(0, K, K))
    ind <- which(pos, arr.ind = T)
    phi <- cbind(ind, tmp)
    
    
    tlam <- cbind(obj$gammal, obj$gammas)
    allgam <- result(tlam, med, SL)
    row.names(allgam) <- c(paste0("F", 1:K), "PSX")
    # # GRD <- schain.grd(ELA) GRD_max<-obj$GRD_max names(GRD_max) <- names(obj$GRD_mean)
    
    out0 <- list(N = N, J = J, K = K, `Miss%` = obj$Nmis/J/N * 100, `LD enabled` = LD, `Burn in` = obj$burn, 
        Iteration = obj$iter, `No. of sig lambda` = NSLA)
    # 'Auto Conv' = obj$conv, 'mean GRD' = obj$GRD_mean, 'max GRD' = GRD_max,
    
    if (LD) {
        out0$"No. of sig LD terms" = no_ofd
        # out$ofd_PSX = offpsx
    }
    
    Jp <- length(obj$cati)
    if (Jp > 0) {
        nthd <- obj$mnoc - 1
        if (nthd == 1) {
            Mthd <- result(obj$THD, med, SL)
            row.names(Mthd) <- paste0("I", 1:Jp)
        } else {
            
            if (!detail) {
                Mthd <- matrix(0, Jp, nthd)
            } else {
                Mthd <- NULL
            }
            for (thd in 1:(nthd)) {
                if (!detail) {
                  Mthd[, thd] <- colMeans(obj$THD[, , thd])  # mean estimates only
                  
                } else {
                  tmpt <- cbind(thd, result(obj$THD[, , thd], med, SL))
                  Mthd <- rbind(Mthd, tmpt)
                }
            }
        }
        out0$"Cat Items" <- obj$cati
        out0$"max No. of categories" = nthd + 1
        # out$THD = Mthd
    } else {
        Mthd <- NULL
    }
    
    if (!is.null(obj$PPP)) 
        out0$PPP <- mean(obj$PPP)
    
    out <- switch(what, basic = out0, lambda = LAM, qlambda = MLA, eigen = eigen, dpsx = dpsx, offpsx = offpsx, 
        phi = phi, shrink = allgam, threshold = Mthd, all = {
            out1 <- out0
            out1$lambda <- LAM
            out1$qlambda <- MLA
            out1$eigen <- eigen
            out1$dpsx <- dpsx
            out1$offpsx <- offpsx
            out1$phi <- phi
            out1$shrink <- allgam
            out1$threshold <- Mthd
            out1
        }, stop(sprintf("Can not show element '%s'", what), call. = FALSE))
    
    options(digits = old_digits)
    return(out)
}


######## end of Summary #################################################
