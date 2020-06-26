################### update PSX #################################################################
Gibbs_PSX <- function(y, mu = 0, ome, la, psx, inv.psx, const, prior) {
    # y=Y;ome=Omega;la=LA;psx=PSX;inv.psx=inv.PSX;mu=0
    J <- const$J
    N <- const$N
    upperind <- const$upperind
    lowerind <- const$lowerind
    ind_nod <- const$ind_nod
    a_gams <- prior$a_gams
    b_gams <- prior$b_gams
    temp <- y - mu - la %*% ome  # J*N
    S <- temp %*% t(temp)  # J*J
    
    # sample gammas
    a_gams <- 1
    b_gams <- 0.1
    apost <- a_gams + J * (J + 1)/2
    bpost <- b_gams + sum(abs(inv.psx))/2  # C is the presicion matrix
    gammas <- rgamma(1, shape = apost, rate = bpost)
    
    # sample tau off-diagonal
    Cadj <- pmax(abs(inv.psx[upperind]), 10^(-6))
    mu_p <- pmin(gammas/Cadj, 10^12)
    gammas_p <- gammas^2
    len <- length(mu_p)
    taus_tmp <- 1/rinvgauss1(len, mean = mu_p, dispersion = 1/gammas_p)
    
    taus <- matrix(0, J, J)
    # taus[upperind]<-taus_tmp
    taus[upperind] <- taus[lowerind] <- taus_tmp
    
    # sample PSX and inv(PSX)
    for (i in 1:J) {
        ind_noi <- ind_nod[, i]
        # tau_temp1<-tau[ind_noi,i]
        Sig11 <- psx[ind_noi, ind_noi]
        Sig12 <- psx[ind_noi, i]
        invC11 <- Sig11 - Sig12 %*% t(Sig12)/psx[i, i]
        # Ci<-(S[i,i]+gammas)*invC11+diag(1/tau_temp1)
        Ci <- (S[i, i] + gammas) * invC11 + diag(1/taus[ind_noi, i])
        Sigma <- chol2inv(chol(Ci))
        mu_i <- -Sigma %*% S[ind_noi, i]
        beta <- mvrnorm(1, mu_i, Sigma)
        inv.psx[ind_noi, i] <- beta
        inv.psx[i, ind_noi] <- beta
        gam <- rgamma(1, shape = N/2 + 1, rate = (S[i, i] + gammas)/2)
        inv.psx[i, i] <- gam + t(beta) %*% invC11 %*% beta
        
        # below updating covariance matrix according to one-column change of precision matrix
        invC11beta <- invC11 %*% beta
        psx[ind_noi, ind_noi] <- invC11 + invC11beta %*% t(invC11beta)/gam
        Sig12 <- -invC11beta/gam
        psx[ind_noi, i] <- Sig12
        psx[i, ind_noi] <- t(Sig12)
        psx[i, i] <- 1/gam
        
    }  # end of i, sample Sig and C=inv(Sig)
    # inv.sqrt.psx<-chol(inv.PSX)
    return(list(obj = psx, inv = inv.psx, gammas = gammas))
}

################## end of update PSX ##########################################################


