is.scalar <- function(x) is.atomic(x) && length(x) == 1L

# x: iter * nvar
schain.grd <- function(x0, auto = F) {
    x <- as.matrix(x0)
    len <- dim(x)[1]
    if (auto) {
        len <- len/2
        x <- x[(1 + len):(2 * len), ]
    }
    x1 <- mcmc(x[1:(len/2), ])
    x2 <- mcmc(x[(len/2 + 1):len, ])
    xx <- mcmc.list(x1, x2)
    mpsr <- gelman.diag(xx, autoburnin = F)
    return(mpsr$psrf)
}

result <- function(dat, med = F, SL = 0.05) {
    if (!med) {
        est = apply(dat, 2, mean)
    } else {
        est = apply(dat, 2, median)
    }

    sd = apply(dat, 2, sd)
    hpd = HPDinterval(mcmc(dat), prob = 1 - SL)
    # hpd.90=HPDinterval(mcmc(dat),prob=.90)
    len = dim(dat)[2]
    sig <- array(F, dim = c(len))
    k <- 1
    for (i in 1:len) {
        if (!(hpd[i, 1] <= 0 && hpd[i, 2] >= 0))
            sig[k] <- T
        # if(!(hpd.95[i,1]<=0 && hpd.95[i,2]>=0)) sig[k]<-2;
        k <- k + 1
    }
    # res=cbind(mean,sd,hpd.95,hpd.90,sig)
    res = cbind(est, sd, hpd, sig)
    return(res)
}

# Wishart dist.
rwish1 <- function(v, S) {
    p <- nrow(S)
    CC <- chol(S)
    Z <- matrix(0, p, p)
    diag(Z) <- sqrt(rchisq(p, v:(v - p + 1)))
    if (p > 1) {
        pseq <- 1:(p - 1)
        Z[rep(p * pseq, pseq) + unlist(lapply(pseq, seq))] <- rnorm(p * (p - 1)/2)
    }
    return(crossprod(Z %*% CC))
}

# slightly adjusted rinvgauss function
rinvgauss1 <- function(n, mean = 1, dispersion = 1) {

    mu <- rep_len(mean, n)
    phi <- rep_len(dispersion, n)
    r <- rep_len(0, n)

    Y <- rnorm(n)^2
    # Yphi <- Y * phi[i] * mu[i]
    Yphi <- Y * phi * mu
    bigphi <- (Yphi > 5e+05)
    if (any(bigphi)) {
        X1 <- Y
        X1[bigphi] <- 1/Yphi[bigphi]
        X1[!bigphi] <- 1 + Yphi[!bigphi]/2 * (1 - sqrt(1 + 4/Yphi[!bigphi]))
    } else {
        X1 <- 1 + Yphi/2 * (1 - sqrt(1 + 4/Yphi))
    }
    firstroot <- (runif(n) < 1/(1 + X1))

    r[][firstroot] <- X1[firstroot]
    r[][!firstroot] <- 1/X1[!firstroot]
    r <- mu * r
    r
}

post_pp <- function(y, mu = 0, ome, la, psx, inv.psx, N, J) {

    Ycen <- y - mu - la %*% ome  # NY*N
    tmp1 <- diag(t(Ycen) %*% inv.psx %*% Ycen)
    # O.tmp<-t(mvrnorm(N,mu=rep(0,K),Sigma=PHI))

    # Ym<-MU+LD%*%Omega # NY*N Y.tmp<-t(mvrnorm(N,mu=rep(0,J),Sigma=PSX)) Y.cen<-Y.tmp+Ym-MU-LD%*%Omega

    # Ycen<-t(mvrnorm(N,mu=rep(0,J),Sigma=psx))

    cpsx <- chol(psx)
    Ycen <- matrix(rnorm(N * J), J, N)
    Ycen <- t(t(Ycen) %*% cpsx)

    tmp2 <- diag(t(Ycen) %*% inv.psx %*% Ycen)

    return(sum(tmp2 - tmp1) > 0)

}


cor_check<-function(phi,K){
  if (length(phi)==1){
    phi<-matrix(phi, K, K)
    diag(phi) <- 1
  }
  c0 <- any(dim(phi)!=K) #not K x K
  if (!c0) c1 <- any(eigen(phi)$values<=0) #not positive definite
  c2<- any(diag(phi)!=1)
  c3<- any(phi!=t(phi))
  if (c0 | c1 | c2 | c3)
    stop("phi should be a scalar in (small negative value,1) or K*K PD correlation matrix.", call. = FALSE)
  return(phi)
}

# check for corr or cov matrix (K*K & PD)
# scalar only for corr; type = 1 for corr, 2 for cov
cm_check<-function(cm,K,diag=1,type=1){
  if (is.scalar(cm)){
    cm<-matrix(cm, K, K)
    diag(cm) <- diag
  }
  c0 <- any(dim(cm)!=K) #not K x K
  c1<- any(cm!=t(cm)) #not symmetric
  if (type ==1){
    c2<- any(diag(cm)!=1)
  }else{
    c2 <- any(diag(cm)<=0)
  }
  if (c0 || c1 || c2) return(-1)

  if (any(eigen(cm)$values<=0)) return(-2) #not positive definite

  return(cm)
}



