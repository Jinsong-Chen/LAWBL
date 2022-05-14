#' @title Partially Confirmatory Item Response Model
#'
#' @description \code{pcirm} is a partially confirmatory approach to item response models (Chen, 2020),
#' which estimates the intercept for continuous and dichotomous data. Similar to PCFA and GPCFA,
#'  there are two major model variants with different constraints for identification. One assumes local
#'  independence (LI) with a more exploratory tendency, which can be also called the E-step.
#'  The other allows local dependence (LD) with a more confirmatory tendency, which can be also
#'  called the C-step. Parameters are obtained by sampling from the posterior distributions with
#'  the Markov chain Monte Carlo (MCMC) techniques. Different Bayesian Lasso methods are used to
#'  regularize the loading pattern and LD. The estimation results can be summarized with \code{\link{summary.lawbl}}
#'  and the factorial eigenvalue can be plotted with \code{\link{plot_lawbl}}.
#'
#' @name pcirm
#'
#' @param dat A \eqn{N \times J} data \emph{matrix} or \emph{data.frame} consisting of the
#'     responses of \eqn{N} individuals to \eqn{J} items. Only continuous and dichotomous data are supported.
#'
#' @param Q A \eqn{J \times K} design matrix for the loading pattern with \eqn{K} factors and \eqn{J} items.
#' Elements are 1, -1, and 0 for specified, unspecified, and zero-fixed loadings, respectively. For models with
#' LI or the E-step, one can specify a few (e.g., 2) loadings per factor. For models with LD or the C-step, the
#' sufficient condition of one specified loading per item is suggested, although there can be a few items
#' without any specified loading. See \code{Examples}.
#'
#' @param LD logical; \code{TRUE} for allowing LD (model with LD or C-step).
#'
#' @param cati The set of dichotomous items in sequence number (i.e., 1 to \eqn{J});
#' \code{NULL} for no and -1 for all items (default is \code{NULL}).
#'
#' @param PPMC logical; \code{TRUE} for conducting posterior predictive model checking.
#'
#' @param burn Number of burn-in iterations before posterior sampling.
#'
#' @param iter Number of formal iterations for posterior sampling (> 0).
#'
#' @param update Number of iterations to update the sampling information.
#'
#' @param missing Value for missing data (default is \code{NA}).
#'
#' @param alas logical; for adaptive Lasso or not. The default is \code{FALSE}.
#'
#' @param rseed An integer for the random seed.
#'
#' @param digits Number of significant digits to print when printing numeric values.
#'
#' @param verbose logical; to display the sampling information every \code{update} or not.
#' \itemize{
#'     \item \code{Feigen}: Eigenvalue for each factor.
#'     \item \code{NLA_le3}: Number of Loading estimates >= .3 for each factor.
#'     \item \code{Shrink}: Shrinkage (or ave. shrinkage for each factor for adaptive Lasso).
#'     \item \code{EPSR & NCOV}: EPSR for each factor & # of convergence.
#'     \item \code{Ave. Int.}: Ave. item intercept.
#'     \item \code{LD>.2 >.1 LD>.2 >.1}: # of LD terms larger than .2 and .1, and LD's shrinkage parameter.
#'     \item \code{#Sign_sw}: Number of sign switch for each factor.
#' }
#'
#' @param sign_check logical; \code{TRUE} for checking sign switch of loading vector.
#'
#' @param sign_eps minimum value for switch sign of loading vector (if \code{sign_check=TRUE}).
#'
#' @param auto_stop logical; \code{TRUE} for enabling auto stop based on \code{EPSR<1.1}.
#'
#' @param max_conv maximum consecutive number of convergence for auto stop.
#'
#' @return \code{pcirm} returns an object of class \code{lawbl} with item intercepts. It contains a lot of information about
#' the posteriors that can be summarized using \code{\link{summary.lawbl}}.
#'
#' @references
#'
#' Chen, J. (2020). A partially confirmatory approach to the multidimensional item response theory with the Bayesian Lasso.
#' \emph{Psychometrika}. 85(3), 738-774. DOI:10.1007/s11336-020-09724-3.
#'
#' @importFrom MASS mvrnorm
#'
#' @export
#'
#' @examples
#' \donttest{
#'####################################
#'#  Example 1: Estimation with LD   #
#'####################################
#'
#' dat <- sim24ccfa21$dat
#' J <- ncol(dat)
#' K <- 3
#' Q<-matrix(-1,J,K);
#' Q[1:8,1]<-Q[9:16,2]<-Q[17:24,3]<-1
#'
#' m0 <- pcirm(dat = dat, Q = Q, LD = TRUE, cati = -1, burn = 2000,iter = 2000)
#' summary(m0) # summarize basic information
#' summary(m0, what = 'qlambda') #summarize significant loadings in pattern/Q-matrix format
#' summary(m0, what = 'offpsx') #summarize significant LD terms
#'
#'####################################
#'#  Example 2: Estimation with LD   #
#'####################################
#'
#' Q<-cbind(Q,-1);
#' Q[15:16,4]<-1
#'
#' m1 <- pcirm(dat = dat, Q = Q, LD = FALSE, cati = -1, burn = 2000,iter = 2000)
#' summary(m1) # summarize basic information
#' summary(m1, what = 'qlambda') #summarize significant loadings in pattern/Q-matrix format
#' summary(m1, what = 'offpsx') #summarize significant LD terms
#' }
pcirm <- function(dat, Q, LD = TRUE,cati = NULL, PPMC = FALSE, burn = 5000, iter = 5000, update = 1000, missing = NA, rseed = 12345,
                  sign_check = FALSE, sign_eps = -.5, auto_stop=FALSE,max_conv=10,digits = 4, alas = FALSE, verbose = FALSE) {

    # cand_thd = 0.2
    # conv = 0

    if (nrow(Q) != ncol(dat))
        stop("The numbers of items in data and Q are unequal.", call. = FALSE)

    if (iter == 0)
        stop("Parameter iter must be larger than zero.", call. = FALSE)

    if (exists(".Random.seed", .GlobalEnv))
        oldseed <- .GlobalEnv$.Random.seed else oldseed <- NULL
    set.seed(rseed)

    oo <- options()       # code line i
    on.exit(options(oo))  # code line i+1
    # old_digits <- getOption("digits")
    options(digits = digits)

    Y <- t(dat)
    Y[which(Y == missing)] <- NA
    N <- ncol(Y)
    J <- nrow(Y)
    int<-T #intercept retained or not
    Q <- as.matrix(Q)
    K <- ncol(Q)
    Jp <- length(cati)
    if (Jp == 1 && cati == -1) {
        cati <- c(1:J)
        Jp <- J
    }

    Nmis <- sum(is.na(Y))
    mind <- which(is.na(Y), arr.ind = TRUE)

    const <- list(N = N, J = J, K = K, Q = Q, cati = cati, Jp = Jp, Nmis = Nmis, int = int)

    ######## Init ########################################################
    miter <- iter + burn
    mOmega <- array(0, dim = c(K, N))  # mean latent variable omega, K*N
    NLA <- sum(Q != 0)  #number of all lambda need to be estimated
    ELA <- array(0, dim = c(iter, NLA))  #Store retained trace of Lambda
    EMU <- array(0, dim = c(iter, J))  #Store retained trace of MU
    # EPSX <- array(0, dim = c(iter, J, J)) #Store retained trace of PSX EPHI <- array(0, dim = c(iter,
    # K, K)) #Store retained trace of PHI
    if (LD) {
        EPSX <- array(0, dim = c(iter, J * (J + 1)/2))
    } else {
        EPSX <- array(0, dim = c(iter, J))
    }
    # EPSX <- ifelse(LD,array(0, dim = c(iter, J*(J+1)/2)),array(0, dim = c(iter, J))) Store retained
    # trace of PSX
    EPHI <- array(0, dim = c(iter, K * (K - 1)/2))  #Store retained trace of PHI
    Egammas <- array(0, dim = c(iter, 1))  #Store retained trace of shrink par gammas
    Egammal <- array(0, dim = c(iter, K))  #Store retained trace of shrink par gammal (per factor)
    # Delta<<-array(0,dim=c(iter,J)) #Store retained trace of Ys scale
    Eppmc <- NULL
    if (PPMC)
        Eppmc <- array(0, dim = c(iter))

    init <- init(y = Y, const = const)
    Y <- init$y
    const <- init$const
    prior <- init$prior
    PSX <- init$PSX
    inv.PSX <- chol2inv(chol(PSX))
    PHI <- init$PHI
    LA <- init$LA
    THD <- init$THD
    gammas <- init$gammas
    gammal_sq <- init$gammal_sq

    if (Jp > 0) {
        mnoc <- const$mnoc
        if (mnoc>2)
            stop("PCIRM is only for continuous or dichotomous data.", call. = FALSE)
        Etd <- array(0, dim = c(iter, Jp, mnoc - 1))
    }


    accrate <- 0
    if (Nmis > 0)
        Y[mind] <- rnorm(Nmis)

    # OME <- t(mvrnorm(N,mu=rep(0,K),Sigma=diag(1,K))) # J*N
    sign_sw <- rep(0, K)
    sign_eps <- -.5

    Eigen <- array(0, dim = c(iter, K))  #Store retained trace of Eigen
    tmp <- which(Q!=0,arr.ind=TRUE)
    pof <- matrix(0,NLA,K) #pos of est lam for each factor
    for (k in 1:K){
        ind<-(tmp[,2]==k)
        pof[ind,k]<-1
    }
    no_conv <- 0

    ######## end of Init #################################################

    ptm <- proc.time()
    for (ii in 1:miter) {
        # ii <- 1
        g = ii - burn
        OME <- Gibbs_Omega(y = Y, la = LA, phi = PHI, inv.psx = inv.PSX, N = N, K = K)
        # PHI<-MH_PHI(ph0=PHI,ome=Omega)
        MU<-Gibbs_MU(y=Y,ome=OME,la=LA,inv.psx=inv.PSX, N = N, J = J,prior = prior)

        if (LD) {
            tmp <- Gibbs_PSX(y = Y, ome = OME, la = LA, psx = PSX, inv.psx = inv.PSX, const = const,
                prior = prior)
            PSX <- tmp$obj
            inv.PSX <- tmp$inv
            gammas <- tmp$gammas
            LAY <- Gibbs_LA_IYC(y = Y, mu = MU, ome = OME, la = LA, psx = PSX, gammal_sq = gammal_sq, thd = THD,
                const = const, prior = prior, alas = alas)
        } else {
            LAY <- Gibbs_LA_IYE(y = Y, mu = MU, ome = OME, la = LA, psx = PSX, gammal_sq = gammal_sq, thd = THD,
                const = const, prior = prior, alas = alas)
            PSX <- LAY$psx
            inv.PSX <- chol2inv(chol(PSX))
        }

        LA <- LAY$la # OME <- LAY$ome

        if (sign_check) {
            chg <- (colSums(LA)<= sign_eps)
            # chg <- (colSums(LA)<= LA_eps)
            if (any(chg)) {
                sign <- diag(1 - 2 * chg)
                sign_sw <- sign_sw + chg
                # if(g<0){chg0_count <- chg0_count + chg}else{chg_count <- chg_count + chg}
                LA <- LA %*% sign
                OME <- t(t(OME) %*% sign)
                print(c("ii=", ii), quote = FALSE)
                cat(sign_sw, fill = TRUE, labels = "#Sign switch:")
            }
        } #end if

        gammal_sq <- LAY$gammal_sq
        # OME <- Gibbs_Omega(y = Y, la = LA, phi = PHI, inv.psx = inv.PSX, N = N, K = K)
        PHI <- MH_PHI(phi = PHI, ome = OME, N = N, K = K, s0 = prior$s_PHI)
        # Omega<-Gibbs_Omega(y=Y,la=LA,phi=PHI,inv.psx=inv.PSX)
        if (Jp > 0) {
            Y[cati, ] <- LAY$ys
            # THD <- LAY$thd
            accrate <- accrate + LAY$accr
        }
        if (Nmis > 0)
            Y[mind] <- LAY$ysm[mind]

        # Save results
        if ((g > 0)) {
            mOmega <- mOmega + OME
            ELA[g, ] <- LA[Q != 0]
            Eigen[g,] <- (ELA[g,]^2)%*%pof
            # EPSX[g, , ] <- PSX EPSX[g, ] <- ifelse(LD,PSX[lower.tri(PSX,diag=TRUE)],diag(PSX))
            if (LD) {
                EPSX[g, ] <- PSX[lower.tri(PSX, diag = TRUE)]
            } else {
                EPSX[g, ] <- diag(PSX)
            }
            Egammas[g, ] <- gammas
            Egammal[g, ] <- colMeans(sqrt(gammal_sq))
            # EPHI[g, , ] <- PHI[, ]
            EPHI[g, ] <- PHI[lower.tri(PHI)]
            EMU[g,]<-MU
            # if (Jp > 0)
            #     Etd[g, , ] <- THD[, 2:mnoc]
            if (PPMC)
                Eppmc[g] <- post_pp(y = Y, ome = OME, la = LA, psx = PSX, inv.psx = inv.PSX, N = N,
                  J = J)
        }

        if (ii%%update == 0){
                if (g > 0) {

                    APSR <- schain.grd(Eigen[1:g,])
                    if (auto_stop) {
                        if (max(APSR[,1]) < 1.1) {
                            no_conv <- no_conv + 1
                        } else{
                            no_conv <- 0
                        }
                    } # end auto_stop
                } # end g


            if(verbose){
                # print(proc.time() - ptm)
                Shrink <- colMeans(sqrt(gammal_sq))
                Feigen <- diag(crossprod(LA))
                NLA_le3 <- colSums(abs(LA) >= 0.3)
                # Meigen <- colMeans(Eigen)
                # Mlambda<-colMeans(LA)

                cat(ii, fill = TRUE, labels = "\nTot. Iter =")
                print(rbind(Feigen, NLA_le3, Shrink))
                cat(mean(MU), fill = TRUE, labels = "Ave. Int.:")
                if (g > 0) cat(c(t(APSR[,1]),no_conv), fill = TRUE, labels = "EPSR & NCONV")

                if (LD) {
                    opsx <- abs(PSX[row(PSX) != col(PSX)])
                    tmp <- c(sum(opsx > 0.2)/2, sum(opsx > 0.1)/2, gammas)
                    print(c("LD>.2", ">.1", "Shrink"), quote = FALSE)
                    print(tmp)
                }

                # if (Jp > 0) {
                #     cat(colMeans(THD), fill = TRUE, labels = "Ave. Thd:")
                #     cat(accrate/update, fill = TRUE, labels = "Acc Rate:")
                #     accrate <- 0
                # }
                # print(chg_count) cat(chg_count, fill = TRUE, labels = '#Sign change:')

            }#end verbose

            if (auto_stop * no_conv >= max_conv) break
            }  # end update

    }  #end of g MCMAX

    if(verbose){
        # cat(chg0_count,chg_count, fill = TRUE, labels = "\n#Sign change:")
        print(proc.time()-ptm)
    }


    if (auto_stop * no_conv >= max_conv) {
        ELA <- ELA[1:g, ]
        Eigen <- Eigen[1:g, ]
        EPSX <- EPSX[1:g, ]
        EPHI <- EPHI[1:g, ]
        Egammal <- Egammal[1:g, ]
        Eppmc <- Eppmc[1:g]
        Egammas <- Egammas[1:g, ]
        # if (Jp > 0)
        #     Etd <- Etd[1:g, , ]
        iter <- g
    }


    # if (conv != 0) {
    #     st <- g/2 + 1
    #     ELA <- ELA[(st):g, ]
    #     # EPSX <- EPSX[(st):g, , ] EPHI <- EPHI[(st):g, , ]
    #     EPSX <- EPSX[(st):g, ]
    #     EPHI <- EPHI[(st):g, ]
    #     Egammal <- Egammal[(st):g, ]
    #     Egammas <- Egammas[(st):g, ]
    #     EMU <- EMU[(st):g,]
    #     # if (Jp > 0)
    #     #     Etd <- Etd[(st):g, , ]
    #     Eppmc <- Eppmc[(st):g]
    #     mOmega <- mOmega/2
    #     iter <- burn <- g/2
    # }

    # chg1_count<-rbind(chg0_count,chg_count)
    out <- list(Q = Q, LD = LD, LA = ELA, Omega = mOmega/iter, PSX = EPSX, iter = iter, burn = burn,
                PHI = EPHI, gammal = Egammal, gammas = Egammas, Nmis = Nmis, PPP = Eppmc,MU = EMU,
                auto_conv = c(auto_stop, no_conv, max_conv),Eigen = Eigen, time = (proc.time()-ptm))

    if (Jp > 0) {
        out$cati = cati
        out$THD = Etd
        out$mnoc = mnoc
    }

    class(out) <- c("lawbl")

    if (!is.null(oldseed))
        .GlobalEnv$.Random.seed <- oldseed else rm(".Random.seed", envir = .GlobalEnv)

    # options(digits = old_digits)

    return(out)
}

