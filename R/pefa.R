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
#' @param rseed An integer for the random seed.
#'
#' @param digits Number of significant digits to print when printing numeric values.
#'
#' @param verbose logical; to display the sampling information every \code{update} or not.
#' \itemize{
#'     \item \code{Feigen}: Eigenvalue for each factor.
#'     \item \code{NLA_le3}: Number of Loading estimates >= .3 for each factor.
#'     \item \code{Shrink}: Shrinkage parameter for each factor.
#'     \item \code{True Fa}: Is the factor identified as true or not.
#'     \item \code{Adj PSR}: Adjusted PSR for each factor.
#'     \item \code{ROW: LA overflow,sign switch,MJF total, current, reset}: Loading overflow, sign switch,
#'      min item per factor reached and related diagnostic information.
#' }
#'
#' @return \code{pcfa} returns an object of class \code{lawbl} without item intercepts. It contains a lot of information about
#' the posteriors that can be summarized using \code{\link{summary.lawbl}}.
#'
#' @references
#'
#' Chen, J. (2021). A Bayesian regularized approach to exploratory factor analysis in one step.
#'  \emph{Structural Equation Modeling: A Multidisciplinary Journal}. DOI: 10.1080/10705511.2020.1854763.
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
pefa <- function(dat, Q=NULL, K=8, mjf=3, PPMC = FALSE, burn = 5000, iter = 5000,
                 update = 1000, rseed = 12345, digits = 4, verbose = FALSE) {

    conv = 0

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
    # Y[which(Y == missing)] <- NA
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

    # Nmis <- sum(is.na(Y))
    # mind <- which(is.na(Y), arr.ind = TRUE)
    Nmis <- 0

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

    EPHI <- array(0, dim = c(iter, K * (K - 1)/2))  #Store retained trace of PHI
    # Egammas <- array(0, dim = c(iter, 1))  #Store retained trace of shrink par gammas
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
    # THD <- init$THD
    # gammas <- init$gammas
    # gammal_sq <- init$gammal_sq

    # accrate <- 0
    # if (Nmis > 0)
    #     Y[mind] <- rnorm(Nmis)



    # # OME <- t(mvrnorm(N,mu=rep(0,K),Sigma=diag(1,K))) # J*N
    # chg_count <- rep(0, K)
    # chg0_count <- rep(0, K)
    # # chg1_count <- rep(0, K)
    # # Jest <- colSums(Q != 0)
    # LA_eps <- -1
    # la_overf <- rep(0, K)

    Eigen <- array(0, dim = c(iter, K))  #Store retained trace of Eigen
    tmp <- which(Q!=0,arr.ind=TRUE)
    pof <- matrix(0,NLA,K) #pos of est lam for each factor
    for (k in 1:K){
        ind<-(tmp[,2]==k)
        pof[ind,k]<-1
    }

    Eigen_eps <- .1
    pig<-rep(1,K)
    tmp<-colSums(Q==1)
    pig[tmp>0]<--1
    # Epig<-array(0,dim=c(iter,K))
    tausq<-array(1,dim=c(J,K))
    tausq[Q==0]<-0
    ilamsq<-rep(1,K)
    ilamsq_t<-ilamsq
    # Elamsq<-array(0,dim=c(iter,K))
    const$indg<-(pig!=-1)

    tmp <- max(round(J*N/500),10)
    const$max_var <- min(10^tmp,10^60) #10^10 to 10^60
    # const$max_var <- 10^24
    const$no_mjf <- N*(-9/800)+49/4 # 10 for 200; 1 for 1000
    # const$no_mjf <- 0

    count <- matrix(0,K,5) #LA overflow,sign switch,MJF total, current, reset
    sign_eps <- -.5
    # la0<-.2
    # i <- 1
    # for (k in 1:K){
    #     if(const$indg[k]) {LA[,k]<-la0; la0 <- la0 - .05/i; i <- i +1}
    # }
    ######## end of Init #################################################

    ptm <- proc.time()
    for (ii in 1:miter) {
        # i <- 1
        g = ii - burn
        OME <- Gibbs_Omega(y = Y, la = LA, phi = PHI, inv.psx = inv.PSX, N = N, K = K)
        # PHI<-MH_PHI(ph0=PHI,ome=Omega)


        # OME<-Gibbs_Omega(y=Y,ly=LA,phi=PHI,inv.psx=inv.PSX)
        main<-Gibbs_pefa_main(y=Y,ome=OME,ly=LA,psx=PSX,tausq=tausq,pig=pig,prior=prior,
                              ilamsq=ilamsq,ilamsq_t=ilamsq_t,const=const,count=count)
        PSX<-main$psx
        inv.PSX<-chol2inv(chol(PSX))
        count <- main$count

        # LA1<-main$ly
        # tmp1<-(abs(LA)>=.3)*(abs(LA1)>=.3)
        # chg0 <- colSums(LA * LA1 * tmp1 < 0)
        # if (any(chg0>1)) { # possible sign change
        #   # chg0_count <- chg0_count + chg0
        #   tmp2<-colSums(tmp1)
        #   # chg <- (chg0 > (chg2/2))*(chg0!=0) #all non-zero la change sign
        #   chg <- (chg0 > (tmp2/2))
        #   if (any(chg)){
        #     sign <- diag(1 - 2 * chg)
        #     # chg_count <- chg_count + chg
        #     if(g<0){chg0_count <- chg0_count + chg}else{chg_count <- chg_count + chg}
        #     LA1 <- LA1 %*% sign
        #     OME <- t(t(OME) %*% sign)
        #     print(c("ii=", ii), quote = FALSE)
        #     # cat(chg_count, fill = TRUE, labels = "#Sign change:")
        #     print(rbind(tmp2,chg0))
        #   }
        # }
        # LA <- LA1


        LA1<-main$ly
        # tausq1<-main$tausq
        tmp1<-abs(LA1)>.99
        # LA1[tmp1]<-LA[tmp1]
        if (any(tmp1)){

            # la_overf<-la_overf+(colSums(tmp1)>0)
            count[,1] <- count[,1]+(colSums(tmp1)>0)
            LA1[tmp1]<-LA[tmp1]
            # tausq1[tmp1]<-tausq[tmp1]
        }

        chg <- (colSums(LA)<= sign_eps)
        # chg <- colSums(LA * LA1) < LA_eps
        if (any(chg)) {
            # print(c("ii=", ii), quote = FALSE)
            # print(colSums(LA * LA1))
            # # print(cbind(LA, LA1))
            sign <- diag(1 - 2 * chg)
            count[,2] <- count[,2]+chg
            # if(g<0){chg0_count <- chg0_count + chg}else{chg_count <- chg_count + chg}
            LA1 <- LA1 %*% sign
            OME <- t(t(OME) %*% sign)
        }

        LA <- LA1
        # tausq<-tausq1
        tausq<-main$tausq

        # OME<-main$ome
        ilamsq<-main$ilamsq
        ilamsq_t<-main$ilamsq_t
        # ilamsq_add<-main$ilamsq_add

        pig<-main$pig

        # gammal_sq <- LAY$gammal_sq
        # OME <- Gibbs_Omega(y = Y, la = LA, phi = PHI, inv.psx = inv.PSX, N = N, K = K)
        PHI <- MH_PHI(phi = PHI, ome = OME, N = N, K = K, prior = prior)
        # Omega<-Gibbs_Omega(y=Y,la=LA,phi=PHI,inv.psx=inv.PSX)


        # Save results
        if ((g > 0)) {
            mOmega <- mOmega + OME
            ELA[g, ] <- LA[Q != 0]
            Eigen[g,] <- (ELA[g,]^2)%*%pof
            EPSX[g, ] <- diag(PSX)

            Egammal[g, ] <- 1/sqrt(ilamsq)
            # Egammal[g, ] <- colMeans(sqrt(gammal_sq))

            # EPHI[g, , ] <- PHI[, ]
            EPHI[g, ] <- PHI[lower.tri(PHI)]
            # # EMU[g,]<-MU

            # Epig[g,]<-pig


            if (PPMC)
                Eppmc[g] <- post_pp(y = Y, ome = OME, la = LA, psx = PSX, inv.psx = inv.PSX, N = N,
                  J = J)
        }

        if (ii%%update == 0)
            {
                if (g > 0) {
                    TF_ind<-(colMeans(Eigen[1:g,])>Eigen_eps)
                  if (conv == 0) {
                      APSR <- schain.grd(Eigen[1:g,TF_ind])
                      # cat(t(APSR[,1]), fill = TRUE, labels = "Adj PSR")

                    # tmp <- schain.grd(ELA[1:g, ])
                    # GRD_mean <- colMeans(tmp)
                    # GRD_max <- c(max(tmp[, 1]), max(tmp[, 2]))
                    # GRD <- c(GRD_mean[1], GRD_max[1])
                    # cat(GRD, fill = TRUE, labels = "PGR mean & max:")
                  } else {
                      APSR <- schain.grd(Eigen[1:g,TF_ind], auto = TRUE)
                      # cat(t(APSR[,1]), fill = TRUE, labels = "Adj PSR")
                      if (mean(APSR[,1])<1.1) break

                    # tmp <- schain.grd(ELA[1:g, ], auto = TRUE)
                    # GRD_mean <- colMeans(tmp)
                    # GRD_max <- c(max(tmp[, 1]), max(tmp[, 2]))
                    # # cat(GRD_mean, fill = TRUE, labels = 'GR Mean & UL:') cat(GRD_max, fill = TRUE, labels = 'GR Max & UL:')
                    # GRD <- c(GRD_mean[1], GRD_max[1])
                    # cat(GRD, fill = TRUE, labels = "PGR mean & max:")
                    # if (conv == 1 && GRD[2] < 1.1)
                    #   break
                    # if (conv == 2 && GRD[1] < 1.1 && GRD[2] < 1.2)
                    #   break
                  }
                }


            if(verbose){
                # print(proc.time() - ptm)
                iShrink <- sqrt(ilamsq)
                Feigen <- diag(crossprod(LA))
                NLA_le3 <- colSums(abs(LA) >= 0.3)
                # Meigen <- colMeans(Eigen)
                # NLA <- (colSums(LA!=0))

                cat(ii, fill = TRUE, labels = "\nTot. Iter =")
                print(rbind(Feigen,NLA_le3, iShrink))
                # cat(chg_count, fill = TRUE, labels = '#Sign change:')
                if (g > 0) {
                    cat(t(TF_ind+0), fill = TRUE, labels = "True Fa")
                    cat(t(APSR[,1]), fill = TRUE, labels = "Adj PSR")
                }
                # cat(t(la_overf), fill = TRUE, labels = "LA Overf")
                # cat(t(count_mjf), fill = TRUE, labels = "MJF")
                # cat(t(ilamsq_t), fill = TRUE, labels = "ilamsq_t")
                print("ROW: LA overflow,sign switch,MJF total, current, reset")
                print(t(count))
            }#end verbose

            }  # end update

    }  #end of g MCMAX

    if(verbose){
        # cat(chg0_count,chg_count, fill = TRUE, labels = "\n#Sign change:")
        print(proc.time()-ptm)
    }

    if (conv != 0) {
        st <- g/2 + 1
        ELA <- ELA[(st):g, ]
        # EPSX <- EPSX[(st):g, , ] EPHI <- EPHI[(st):g, , ]
        # EPSX <- EPSX[(st):g, ]
        EPHI <- EPHI[(st):g, ]
        Egammal <- Egammal[(st):g, ]
        # Egammas <- Egammas[(st):g, ]
        # if (Jp > 0)
        #     Etd <- Etd[(st):g, , ]
        Eppmc <- Eppmc[(st):g]
        mOmega <- mOmega/2
        iter <- burn <- g/2
    }

    # chg1_count<-rbind(chg0_count,chg_count)
    out <- list(Q = Q, LD = FALSE, LA = ELA,PSX = EPSX, Omega = mOmega/iter, iter = iter, burn = burn,
                PHI = EPHI, gammal = Egammal,  Nmis = Nmis, PPP = Eppmc, conv = conv,
                Eigen = Eigen, APSR = APSR,TF_ind=TF_ind)

    class(out) <- c("lawbl")

    if (!is.null(oldseed))
        .GlobalEnv$.Random.seed <- oldseed else rm(".Random.seed", envir = .GlobalEnv)

    # options(digits = old_digits)

    return(out)
}

