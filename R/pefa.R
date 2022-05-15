#' @title Partially Exploratory Factor Analysis
#'
#' @description \code{PEFA} is a partially exploratory approach to factor analysis, which can incorporate
#' partial knowledge together with unknown number of factors, using bi-level Bayesian regularization.
#' When partial knowledge is not needed, it reduces to the fully exploratory factor analysis (\code{FEFA}; Chen, 2021).
#' A large number of factors can be imposed for selection where true factors will be identified against spurious factors.
#' The loading vector is reparameterized to tackle model sparsity at the factor and loading levels
#' with the multivariate spike and slab priors. Parameters are obtained by sampling from the posterior
#' distributions with the Markov chain Monte Carlo (MCMC) techniques. The estimation results can be summarized
#' with \code{\link{summary.lawbl}} and the trace or density of the posterior can be plotted with \code{\link{plot_lawbl}}.
#'
#' @name pefa
#'
#' @param dat A \eqn{N \times J} data \emph{matrix} or \emph{data.frame} consisting of the
#'     responses of \eqn{N} individuals to \eqn{J} items.
#'
#' @param Q A \eqn{J \times K} design matrix for the loading pattern with \eqn{K} factors and \eqn{J} items for \code{PEFA}.
#' Elements are 1, -1, and 0 for specified, unspecified, and zero-fixed loadings, respectively. It's not needed
#' for \code{FEFA}, which is the default. See \code{Examples}.
#'
#' @param K Maximum number of factors for selection under \code{FEFA}. Not used for \code{PEFA}.
#'
#' @param mjf Minimum number of items per factor.
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
#' @param rseed An integer for the random seed.
#'
#' @param digits Number of significant digits to print when printing numeric values.
#'
#' @param verbose logical; to display the sampling information every \code{update} or not.
#' \itemize{
#'     \item \code{Feigen}: Eigenvalue for each factor.
#'     \item \code{NLA_lg0}: Number of Loading magnitudes > 0 for each factor.
#'     \item \code{iShrink}: Inverted shrinkage parameter for each factor.
#'     \item \code{True Fa}: Is the factor identified as true or not.
#'     \item \code{EPSR & NCOV}: EPSR for each factor & # of convergence.
#'     \item \code{ROW: LA overflow,sign switch,bk=0, <eig_eps}: Loading overflow, sign switch,
#'      vector bk=0 and eigenvalue<eig_eps.
#' }
#'
#' @param rfit logical; \code{TRUE} for providing relative fit (DIC, BIC, AIC).
#'
#' @param eig_eps minimum eigenvalue for factor extraction.
#'
#' @param sign_eps minimum value for switch sign of loading vector.
#'
#' @param rs logical; \code{TRUE} for enabling recommendation system.
#'
#' @param auto_stop logical; \code{TRUE} for enabling auto stop based on EPSR.
#'
#' @param max_conv maximum consecutive number of convergence for auto stop.
#'
#'
#' @return \code{pcfa} returns an object of class \code{lawbl} without item intercepts. It contains a lot of information about
#' the posteriors that can be summarized using \code{\link{summary.lawbl}}.
#'
#' @references
#'
#' Chen, J. (2021). A Bayesian regularized approach to exploratory factor analysis in one step.
#'  \emph{Structural Equation Modeling: A Multidisciplinary Journal},
#'   28(4), 518-528. DOI: 10.1080/10705511.2020.1854763.
#'
#'	Chen, J. (In Press). Fully and partially exploratory factor analysis with bi-level Bayesian regularization.
#'	 \emph{Behavior Research Methods}.
#'
#' @importFrom MASS mvrnorm
#'
#' @export
#'
#' @examples
#' \donttest{
#'#####################################################
#'#  Example 1: Fully EFA                             #
#'#####################################################
#'
#' dat <- sim18cfa0$dat
#'
#' m0 <- pefa(dat = dat, K=5, burn = 2000, iter = 2000,verbose = TRUE)
#' summary(m0) # summarize basic information
#' summary(m0, what = 'qlambda') #summarize significant loadings in pattern/Q-matrix format
#' summary(m0, what = 'phi') #summarize factorial correlations
#' summary(m0, what = 'eigen') #summarize factorial eigenvalue
#'
#'##########################################################
#'#  Example 2: PEFA with two factors partially specified  #
#'##########################################################
#'
#' J <- ncol(dat)
#' K <- 5
#' Q<-matrix(-1,J,K);
#' Q[1:2,1]<-Q[7:8,2]<-1
#' Q
#'
#' m1 <- pefa(dat = dat, Q = Q,burn = 2000, iter = 2000,verbose = TRUE)
#' summary(m1)
#' summary(m1, what = 'qlambda')
#' summary(m1, what = 'phi')
#' summary(m1,what='eigen')
#' }
pefa <- function(dat, Q=NULL, K=8, mjf=3, PPMC = FALSE, burn = 5000, iter = 5000,missing = NA,
                  eig_eps = 1,sign_eps = 0, rfit = TRUE, rs = FALSE, update = 1000, rseed = 12345, verbose = FALSE,
                 auto_stop=FALSE,max_conv=10,digits = 4)  #ini.burn = 100
{

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
    ysig<-apply(Y, 1, sd, na.rm=TRUE)
    ybar<-apply(Y, 1, mean, na.rm=TRUE)

    Nmis <- sum(is.na(Y))
    mind <- which(is.na(Y), arr.ind = TRUE)
    # }

    N <- ncol(Y)
    J <- nrow(Y)
    if (is.null(Q)){
        Q <- matrix(-1,J,K)
    }else{
        Q <- as.matrix(Q)
        K <- ncol(Q)
    }

    if (nrow(Q) != ncol(dat))
        stop("The numbers of items in data and Q are unequal.", call. = FALSE)

    int<-F #intercept retained or not

    t_num<-1000

    const <- list(N = N, J = J, K = K, Q = Q, cati = NULL, Jp = 0, int = int,
                  t_num=t_num, mjf=mjf)

    ######## Init ########################################################
    miter <- iter + burn
    mOmega <- array(0, dim = c(K, N))  # mean latent variable omega, K*N
    NLA <- sum(Q != 0)  #number of all lambda need to be estimated
    ELA <- array(0, dim = c(iter, NLA))  #Store retained trace of Lambda
    EMU <- array(0, dim = c(iter, J))  #Store retained trace of MU
    EPSX <- array(0, dim = c(iter, J)) #Store retained trace of PSX
    # EPSX <- array(0, dim = c(iter, J * (J + 1)/2))

    EPHI <- array(0, dim = c(iter, K * (K - 1)/2))  #Store retained trace of PHI
    # Egammas <- array(0, dim = c(iter, 1))  #Store retained trace of shrink par gammas
    Egammal <- array(0, dim = c(iter, K))  #Store retained trace of shrink par gammal (per factor)
    # Delta<<-array(0,dim=c(iter,J)) #Store retained trace of Ys scale
    Eppmc <- NULL
    if (PPMC)
        Eppmc <- array(0, dim = c(iter))

    ini <- init(y = Y, const = const)
    Y <- ini$y
    const <- ini$const
    prior <- ini$prior

    PSX <- ini$PSX
    inv.PSX <- chol2inv(chol(PSX))
    PHI <- ini$PHI
    LA <- ini$LA
    # LA[Q == -1] <- .1+(runif(sum(Q==-1))-.5)/100
    # THD <- init$THD
    # gammas <- init$gammas
    # gammal_sq <- init$gammal_sq

    # # accrate <- 0
    if (Nmis > 0)
        Y[mind] <- rnorm(Nmis)

    yps<-0

    Eigen <- array(0, dim = c(iter, K))  #Store retained trace of Eigen
    tmp <- which(Q!=0,arr.ind=TRUE)
    pof <- matrix(0,NLA,K) #pos of est lam for each factor
    for (k in 1:K){
        ind<-(tmp[,2]==k)
        pof[ind,k]<-1
    }

    # LA[Q != 0] <- .1 #+(runif(sum(Q==-1))-.5)/100
    # bk_pr<-rep(0,K)
    sksq<-rep(1,K)
    sksq_t<-sksq
    FIS <-const$FIS<-!colSums(Q==1) # factor involved in selection

    tausq<-array(1,dim=c(J,K))
    tausq[Q==0]<-0
    LA_OF <- 1.1 #.999 #overflow value
    PHI0 <- diag(K)
    count <- matrix(0,K,4) #LA overflow,sign switch, bk = 0, < eigen_eps
    # sign_eps <- 0 #-.01

    # K0 <- sum(!FIS)
    # TK <- K
    # ini.burn <- 100
    OME <- Gibbs_Omega(y = Y, la = LA, phi = PHI, inv.psx = inv.PSX, N = N, K = K)
    # OME <- t(mvrnorm(N,mu=rep(0,K),Sigma=diag(1,K))) # J*N
    lsum <- 0
    no_conv <- 0

    ######## end of Init #################################################

    ptm <- proc.time()
    for (ii in 1:miter) {
        # ii <- 1
        g = ii - burn

        main<-Gibbs_BLR_SSP(y=Y,ome=OME,ly=LA,prior=prior,tausq=tausq,
                            sksq=sksq,sksq_t=sksq_t,const=const)
        dPSX<-main$dpsx
        inv.PSX<- diag(1/dPSX)
        # bk<-main$bk
        tausq<-main$tausq
        sksq<-main$sksq
        sksq_t<-main$sksq_t
        # sksq<-pmax(main$sksq,1e-200)

        bk_pr<-main$bk_pr
        count[,3]<-count[,3]+!bk_pr

        LA<-main$ly
        # tausq1<-main$tausq
        tmp1<-abs(LA)>LA_OF
        if (any(tmp1)){
        # if (any(abs(LA)>.99)){

            # la_overf<-la_overf+(colSums(tmp1)>0)
            count[,1] <- count[,1]+(colSums(tmp1)>0)
            # LA1[tmp1]<-LA[tmp1]
            LA[LA>LA_OF]<-LA_OF
            LA[LA< -LA_OF]<--LA_OF
            # tausq1[tmp1]<-tausq[tmp1]
        }


        Feigen <- diag(crossprod(LA))
        # Feigen <- colSums(LA^2)
        ind <- !(FIS & (Feigen < eig_eps))
        # Feigen[!ind] <- 0
        # count[,4]<-count[,4]+!ind-!bk_pr
        chg <- (colSums(LA)< sign_eps)
        # chg <- colSums(LA * LA1) < LA_eps
        if (any(chg)) {
            sign <- diag(1 - 2 * chg)
            count[,2] <- count[,2]+chg
            LA <- LA %*% sign
            OME <- t(t(OME) %*% sign)
            # count[ind,2] <- count[ind,2]+chg[ind]
            # LA[,ind] <- LA[,ind] %*% sign
            # OME[ind,] <- t(t(OME[ind,]) %*% sign)
          }

        # if(ii > ini.burn) LA[,!ind] <- 0
        PHI <- MH_PHI(phi = PHI, ome = OME, N = N, K = K, s0 = prior$s_PHI)
        OME <- Gibbs_Omega(y = Y, la = LA, phi = PHI, inv.psx = inv.PSX, N = N, K = K)
        # if(ii > ini.burn) LA[,!ind] <- 0
        LA[,!ind] <- 0
        Feigen[!ind] <- 0
        count[,4]<-count[,4]+!ind-!bk_pr

        # TK <- sum(ind)
        # if(TK == K || ii < ini.burn){
        #     # PHI <- MH_PHI(phi = PHI, ome = OME, N = N, K = K, prior = prior)
        #     OME <- Gibbs_Omega(y = Y, la = LA, phi = PHI, inv.psx = inv.PSX, N = N, K = K)
        # }else{
        #     # ind <- Feigen > 0
        #     phi0<-PHI[ind,ind]
        #     # PHI <- PHI0
        #     # PHI[ind,ind] <- MH_PHI1(phi = phi0, ome = OME[ind,], N = N, K = TK, s0 = prior$s_PHI[ind,ind])
        #     OME[ind,] <- Gibbs_Omega(y = Y, la = LA[,ind], phi = phi0, inv.psx = inv.PSX, N = N, K = TK)
        #     LA[,!ind] <- 0
        #     # PHI[!ind,!ind]<-0
        #     count[,4]<-count[,4]+!ind-!bk_pr
        # }

        if (Nmis > 0){
            ytmp <- matrix(rnorm(N * J), J, N) * sqrt(diag(PSX))+ LA %*% OME
            Y[mind] <- ytmp[mind]
        }

        # Save results
        if ((g > 0)) {
            mOmega <- mOmega + OME
            ELA[g, ] <- LA[Q != 0]
            # Eigen[g,] <- (ELA[g,]^2)%*%pof
            Eigen[g,] <- Feigen
            EPSX[g, ] <- dPSX
            # EPSX[g, ] <- PSX[lower.tri(PSX, diag = TRUE)]
            Egammal[g, ] <- 1/sqrt(sksq)
            EPHI[g, ] <- PHI[lower.tri(PHI)]
            # Epig[g,]<-pig

            if (PPMC)
                Eppmc[g] <- post_pp(y = Y, ome = OME, la = LA, psx = PSX, inv.psx = inv.PSX, N = N,
                  J = J)
            if (rfit){
              Yc <- Y - LA %*% OME  # J*N
              tmp<-(t(Yc) %*%chol(inv.PSX))^2
              # tmp<-(t(Yc) %*%chol(chol2inv(chol(PSX))))^2
              # lsum<-lsum+sum(tmp)+N*(log(det(PSX)))
              lsum<-lsum+sum(tmp)+N*(sum(log(dPSX)))
            } #end dic

            if (rs) {
                if (Nmis == 0)
                    # ytmp <- matrix(rnorm(N * J), J, N) + LA %*% OME/sqrt(diag(PSX))
                    ytmp <- matrix(rnorm(N * J), J, N) * sqrt(diag(PSX))+ LA %*% OME
                    # ytmp <- ytmp/apply(ytmp, 1, sd)
                yps<-yps + ytmp
            }
        }

        if (ii%%update == 0){
                if (g > 0) {
                    TF_ind<-(colMeans(Eigen[1:g,])>eig_eps)
                    APSR <- schain.grd(Eigen[1:g,TF_ind])
                  # if (auto_stop) {
                      if (max(APSR[,1]) <= 1.1) {
                        no_conv <- no_conv + 1
                      } else{
                        no_conv <- 0
                      }
                    # } # end auto_stop
                } # end g

            if(verbose){
                iShrink <- sqrt(sksq)
                NLA_lg0 <- colSums(abs(LA) > 0)

                cat(ii, fill = TRUE, labels = "\nTot. Iter =")
                print(rbind(Feigen,NLA_lg0,iShrink))
                # cat(chg_count, fill = TRUE, labels = '#Sign change:')
                if (g > 0) {
                    cat(t(TF_ind+0), fill = TRUE, labels = "Tru Fac")
                    cat(c(t(APSR[,1]),no_conv), fill = TRUE, labels = "EPSR & NCONV")
                }
                print("ROW: LA overflow, sign switch, bk=0, <eig_eps")
                print(t(count))
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
        iter <- g
    }

    if (rfit){
      lpry <- lsum/iter+N*log(2*pi)
      # lsum <- lsum/iter
    }
    # chg1_count<-rbind(chg0_count,chg_count)
    out <- list(Q = Q, LD = FALSE, LA = ELA,PSX = EPSX, Omega = mOmega/iter, iter = iter, burn = burn,
                PHI = EPHI, gammal = Egammal,  Nmis = Nmis, PPP = Eppmc, Eigen = Eigen, APSR = APSR,Y = Y,
                auto_conv = c(auto_stop, no_conv, max_conv), eig_eps = eig_eps,lpry=lpry, time = (proc.time()-ptm))

    if (rs) {
        yp<-yps/iter*ysig+ybar
        out$yp = t(yp)
    }
    class(out) <- c("lawbl")

    if (!is.null(oldseed))
        .GlobalEnv$.Random.seed <- oldseed else rm(".Random.seed", envir = .GlobalEnv)

    # options(digits = old_digits)

    return(out)
}
