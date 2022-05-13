################## update Omega ##############################################################
Gibbs_Omega <- function(y, mu = 0, la, phi, inv.psx, N, K) {
    # y=Y;ly=LD;psx=PSX;mu=0;inv.psx=inv.PSX
    inv.sqrt.psx <- chol(inv.psx)
    inv.phi <- chol2inv(chol(phi))
    ISG <- crossprod(inv.sqrt.psx %*% la) + inv.phi
    SIG <- chol2inv(chol(ISG))
    Ycen <- y - mu
    mean <- SIG %*% t(la) %*% inv.psx %*% Ycen

    CSIG <- chol(SIG)
    ome <- matrix(rnorm(N * K), K, N)
    ome <- t(t(ome) %*% CSIG) + mean

    return(ome)
}
################## end of update Omega #######################################################

######## update PHI ########################################################
MH_PHI <- function(phi, ome, N, K, s0) {
    # ome=Omega;ph0=PHI
    # v0 <- prior$v_PHI
    # s0 <- prior$s_PHI
    inv.cov <- rwish1(K + 2 + N, solve(tcrossprod(ome) + s0))
    cov <- chol2inv(chol(inv.cov))

    tmp <- chol2inv(chol(sqrt(diag(diag(cov)))))
    phi1 <- tmp %*% cov %*% tmp
    acc <- exp((K + 1)/2 * (log(det(phi1)) - log(det(phi))))
    acid <- acc > runif(1)
    phi <- phi1 * acid + phi * (1 - acid)
    # inv.CPH<-chol2inv(chol(CPH))
    return(phi)
}
######## end of update PHI #################################################

######## update PHI ########################################################
# old version; all prior needed
MH_PHI1 <- function(phi, ome, N, K, prior) {
    # ome=Omega;ph0=PHI
    v0 <- prior$v_PHI
    s0 <- prior$s_PHI
    S<-solve(tcrossprod(ome) + s0)
    inv.cov <- rwish1(v0 + N, S)
    cov <- chol2inv(chol(inv.cov))

    tmp <- chol2inv(chol(sqrt(diag(diag(cov)))))
    phi1 <- tmp %*% cov %*% tmp
    acc <- exp((K + 1)/2 * (log(det(phi1)) - log(det(phi))))
    acid <- acc > runif(1)
    phi <- phi1 * acid + phi * (1 - acid)
    # inv.CPH<-chol2inv(chol(CPH))
    return(phi)
}
######## end of update PHI #################################################

################## update Omega ##############################################################
# allow factor mean to be specified (i.e., no zero)
Gibbs_Omega1 <- function(y, mu = 0, la, phi, inv.psx, N, K,int) {
    # y=Y;la=LA;psx=PSX;mu=0;inv.psx=inv.PSX
    inv.sqrt.psx <- chol(inv.psx)
    # inv.phi <- solve(phi)
    inv.phi <- chol2inv(chol(phi))
    ISG <- crossprod(inv.sqrt.psx %*% la) + inv.phi
    SIG <- chol2inv(chol(ISG))
    # SIG <- solve(ISG)
    Ycen <- y - mu
    mean <- SIG %*% (t(la) %*% inv.psx %*% Ycen + inv.phi %*% int)

    CSIG <- chol(SIG)
    ome <- matrix(rnorm(N * K), K, N)
    ome <- t(t(ome) %*% CSIG) + mean

    return(ome)
}
################## end of update Omega #######################################################

################## update Omega ##############################################################
# allow factor mean and variance to be specified
Gibbs_Omega2 <- function(y, mu = 0, la, inv.psx, N, K,m.add,s.add) {
    # y=Y;la=LA;psx=PSX;mu=0;inv.psx=inv.PSX
    inv.sqrt.psx <- chol(inv.psx)
    # inv.phi <- solve(phi)
    # inv.phi <- chol2inv(chol(phi))
    ISG <- crossprod(inv.sqrt.psx %*% la) + s.add
    SIG <- chol2inv(chol(ISG))
    # SIG <- solve(ISG)
    Ycen <- y - mu
    mean <- SIG %*% (t(la) %*% inv.psx %*% Ycen + m.add)

    CSIG <- chol(SIG)
    ome <- matrix(rnorm(N * K), K, N)
    ome <- t(t(ome) %*% CSIG) + mean

    return(ome)
}
################## end of update Omega #######################################################

