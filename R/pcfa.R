#' @title Partially Confirmatory Factor Analysis
#'
#' @description \code{PCFA} is a partially confirmatory approach covering a wide range of
#' the exploratory-confirmatory continuum in factor analytic models (Chen, Guo, Zhang, & Pan, 2020).
#'  There are two major model variants with different constraints for identification. One assumes local
#'  independence (LI) with a more exploratory tendency, which can be also called the E-step.
#'  The other allows local dependence (LD) with a more confirmatory tendency, which can be also
#'  called the C-step. Parameters are obtained by sampling from the posterior distributions with
#'  the Markov chain Monte Carlo (MCMC) techniques. Different Bayesian Lasso methods are used to
#'  regularize the loading pattern and LD. The estimation results can be summarized with \code{\link{summary.lawbl}}
#'  and the factorial eigenvalue can be plotted with \code{\link{plot_eigen}}.
#'
#' @name pcfa
#'
#' @param dat A \eqn{N \times J} \code{matrix} or \code{data.frame} consisting of the
#'     responses of \eqn{N} individuals to \eqn{J} items. Only continuous data are supported currently.
#'
#' @param Q A \eqn{J \times K} design matrix for the loading pattern with \eqn{K} factors and \eqn{J} items.
#' Elements are 1, -1, and 0 for specified, unspecified, and zero-fixed loadings, respectively. For models with
#' LI or the E-step, one can specify a few (e.g., 2) loadings per factor. For models with LD or the C-step, the
#' sufficient condition of one specified loading per item is suggested, although there can be a few items
#' without any specified loading (Chen, Guo, Zhang, & Pan, 2020). See \code{Examples}.
#'
#' @param LD logical; \code{TRUE} for allowing LD (model with LD or C-step).
#'
#' @param PPMC logical; \code{TRUE} for conducting posterior predictive model checking.
#'
#' @param burn Number of burn-in iterations before posterior sampling.
#'
#' @param iter Number of formal iterations for posterior sampling.
#'
#' @param update Number of iterations to update the sampling information.
#'
#' @param missing Value for missing data (default is \code{NA}) (under development).
#'
#' @param alas logical; for adaptive Lasso or not. The default is \code{FALSE}, which seems slightly stabler.
#'
#' @param rseed An integer for the random seed.
#'
#' @param digits Number of significant digits to print when printing numeric values.
#'
#' @return \code{pcfa} returns an object of class \code{pcfa}. It contains a lot of information about
#' the posteriors that can be summarized using \code{\link{summary.lawbl}}. The factorial eigenvalue
#'  can be plotted with \code{\link{plot_eigen}}.
#'
#' @references
#'
#' Chen, J., Guo, Z., Zhang, L., & Pan, J. (2020). A partially confirmatory approach to scale development
#'  with the Bayesian Lasso. \emph{Psychological Methods}.
#'  Advance online publication. DOI:10.1037/met0000293
#'
#' @importFrom MASS mvrnorm
#'
#' @export
#'
#' @examples
#' \donttest{
#'####################################
#'#  Example 1: Estimation with LI   #
#'####################################
#'
#' dat <- sim18cfa1$dat
#' J <- ncol(dat)
#' K <- 3
#' Q<-matrix(-1,J,K);
#' Q[1:2,1]<-Q[9:10,2]<-Q[13:14,3]<-1
#' Q
#' m0 <- pcfa(dat = dat, Q = Q, LD = FALSE,burn = 1000, iter = 1000)
#' summary(m0) # summarize basic information
#' summary(m0, what = 'lambda') #summarize significant loadings
#' summary(m0, what = 'qlambda') #summarize significant loadings in pattern/Q-matrix format
#' summary(m0, what = 'offpsx') #summarize significant LD terms
#'
#'####################################
#'#  Example 2: Estimation with LD   #
#'####################################
#'
#' Q<-matrix(-1,J,K);
#' Q[1:6,1]<-Q[7:12,2]<-Q[13:18,3]<-1
#' Q
#' m1 <- pcfa(dat = dat, Q = Q,burn = 1000, iter = 1000)
#' summary(m1) # summarize basic information
#' summary(m1, what = 'qlambda') #summarize significant loadings in pattern/Q-matrix format
#' summary(m1, what = 'offpsx') #summarize significant LD terms
#' }
pcfa <- function(dat, Q, LD = TRUE, PPMC = FALSE, burn = 5000, iter = 5000, update = 1000, missing = NA, rseed = 12345,
    digits = 4, alas = FALSE) {
    cati = NULL
    conv = 0
    cand_thd = 0

    if (nrow(Q) != ncol(dat))
        stop("The numbers of items in data and Q are unequal.", call. = FALSE)

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
    Q <- as.matrix(Q)
    K <- ncol(Q)
    Jp <- length(cati)

    Nmis <- sum(is.na(Y))
    mind <- which(is.na(Y), arr.ind = TRUE)

    const <- list(N = N, J = J, K = K, Q = Q, cati = cati, Jp = Jp, Nmis = Nmis, cand_thd = cand_thd)

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


    Eigen <- array(0, dim = c(iter, K))  #Store retained trace of Eigen
    tmp <- which(Q!=0,arr.ind=TRUE)
    pof <- matrix(0,NLA,K) #pos of est lam for each factor
    for (k in 1:K){
        ind<-(tmp[,2]==k)
        pof[ind,k]<-1
    }

    ######## end of Init #################################################

    ptm <- proc.time()
    for (ii in 1:miter) {
        # i <- 1
        g = ii - burn
        OME <- Gibbs_Omega(y = Y, la = LA, phi = PHI, inv.psx = inv.PSX, N = N, K = K)
        # PHI<-MH_PHI(ph0=PHI,ome=Omega)

        if (LD) {
            tmp <- Gibbs_PSX(y = Y, ome = OME, la = LA, psx = PSX, inv.psx = inv.PSX, const = const,
                prior = prior)
            PSX <- tmp$obj
            inv.PSX <- tmp$inv
            gammas <- tmp$gammas
            LAY <- GwMH_LA_MYC0(y = Y, ome = OME, la = LA, psx = PSX, gammal_sq = gammal_sq, thd = THD,
                const = const, prior = prior, alas = alas)
        } else {
            LAY <- GwMH_LA_MYE0(y = Y, ome = OME, la = LA, psx = PSX, gammal_sq = gammal_sq, thd = THD,
                const = const, prior = prior, alas = alas)
            PSX <- LAY$psx
            inv.PSX <- chol2inv(chol(PSX))
        }
        LA <- LAY$la

        gammal_sq <- LAY$gammal_sq
        # OME <- Gibbs_Omega(y = Y, la = LA, phi = PHI, inv.psx = inv.PSX, N = N, K = K)
        PHI <- MH_PHI(phi = PHI, ome = OME, N = N, K = K, prior = prior)
        # Omega<-Gibbs_Omega(y=Y,la=LA,phi=PHI,inv.psx=inv.PSX)

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
            # EMU[g,]<-MU
            # if (Jp > 0)
            #     Etd[g, , ] <- THD[, 2:mnoc]
            if (PPMC)
                Eppmc[g] <- post_pp(y = Y, ome = OME, la = LA, psx = PSX, inv.psx = inv.PSX, N = N,
                  J = J)
        }

        if (ii%%update == 0)
            {
                # cat(ii, fill = TRUE, labels = "\nTot. Iter =")
                # print(proc.time() - ptm)
                Shrink <- colMeans(sqrt(gammal_sq))
                Feigen <- diag(crossprod(LA))
                NLA_lg3 <- colSums(abs(LA) > 0.3)
                # print(rbind(Feigen, NLA_lg3, Shrink))

                if (g > 0) {
                      APSR <- schain.grd(Eigen[1:g,])
                      # cat(t(APSR[,1]), fill = TRUE, labels = "Adj PSR")

                }

                if (LD) {
                    opsx <- abs(PSX[row(PSX) != col(PSX)])
                    tmp <- c(sum(opsx > 0.2)/2, sum(opsx > 0.1)/2, gammas)
                    # print(c("LD>.2", ">.1", "Shrink"), quote = FALSE)
                    # print(tmp)
                }

            }  # end update

    }  #end of g MCMAX


    out <- list(Q = Q, LD = LD, LA = ELA, Omega = mOmega/iter, PSX = EPSX, iter = iter, burn = burn,
                PHI = EPHI, gammal = Egammal, gammas = Egammas, Nmis = Nmis, PPP = Eppmc, conv = conv,
                Eigen = Eigen, APSR = APSR)


    class(out) <- c("lawbl")

    if (!is.null(oldseed))
        .GlobalEnv$.Random.seed <- oldseed else rm(".Random.seed", envir = .GlobalEnv)

    # options(digits = old_digits)

    return(out)
}

