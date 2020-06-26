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
MH_PHI <- function(phi, ome, N, K, prior) {
    # ome=Omega;ph0=PHI
    v0 <- prior$v_PHI
    s0 <- prior$s_PHI
    inv.cov <- rwish1(v0 + N, solve(tcrossprod(ome) + s0))
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

