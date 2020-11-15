######## Init ########################################################
init <- function(y, const) {

    N <- const$N
    J <- const$J
    int <- const$int
    K <- const$K
    Q <- const$Q
    pind <- const$cati
    Jp <- const$Jp
    const$sub_sl <- Q == 1
    const$len_sl <- rowSums(Q == 1)
    const$sub_ul <- Q == -1
    const$len_ul <- rowSums(Q == -1)

    indmx <- matrix(1:J^2, nrow = J, ncol = J)
    temp <- indmx[upper.tri(indmx)]
    upperind <- temp[temp > 0]
    indmx_t <- t(indmx)
    temp <- indmx_t[upper.tri(indmx_t)]
    lowerind <- temp[temp > 0]
    const$upperind <- upperind
    const$lowerind <- lowerind

    ind_nod <- array(0, dim = c(J - 1, J))
    for (j in 1:J) {
        if (j == 1) {
            tmp <- 2:J
        } else if (j == J) {
            tmp <- 1:(J - 1)
        } else tmp <- c(1:(j - 1), (j + 1):J)
        ind_nod[, j] <- tmp
    }  # end of j
    const$ind_nod <- ind_nod

    # Prior mean of MU & Loading
    prior <- list()
    prior$m_LA <- 0
    prior$m_MU <- rep(0, J)
    # Inverse prior variance of MU & loading
    prior$s_LA <- 0.25
    prior$s_MU <- 0.25

    # hyperparameters of Wishart distribution
    prior$v_PHI <- K + 2
    prior$s_PHI <- matrix(0.1, K, K)
    diag(prior$s_PHI) <- 1

    # #Hyperparameters of Gamma distribution for the shrinkage parameter
    prior$a_gams <- 1
    prior$b_gams <- .1
    prior$a_gaml_sq <- 1
    prior$b_gaml_sq <- .1

    # initial values
    LA <- matrix(0, J, K)
    LA[Q == 1] <- .7
    LA[Q == -1] <- .1
    # LA <- matrix(runif(J * K), J, K)
    # LA[Q == 0] <- 0

    PHI <- matrix(0.1, K, K)
    diag(PHI[, ]) <- 1

    PSX <- matrix(0, J, J)
    diag(PSX) <- 0.3
    # inv.PSX <-chol2inv(chol(PSX))

    # shrinkage for PSX
    gammas <- 0

    # shrinkage for loading
    gammal_sq <- array(1, dim = c(J, K))
    gammal_sq[Q != -1] <- 0

    if (Jp > 0) {
        Z <- y[pind, ]
        # inf<-.Machine$double.xmax
        inf <- 1e+200
        # val <- NULL
        noc <- rep(0, Jp)
        zind <- Z  #index start from 1
        for (j in 1:Jp) {
            tmp <- sort(unique(na.omit(Z[j, ])))
            # val <- c(val, tmp)
            noc[j] <- length(tmp)
            for (m in 1:noc[j]) {
                zind[j, Z[j, ] == tmp[m]] <- m
            }
        }
        mnoc <- max(noc)
        THD <- matrix(0, Jp, mnoc + 1)
        for (j in 1:Jp) {
            # j=1
            THD[j, 1:(noc[j] + 1)] <- c(-inf, 0:(noc[j] - 2), inf)
            if (noc[j] < mnoc)
                THD[j, (noc[j] + 2):(mnoc + 1)] <- inf
        }
        ys <- t(mvrnorm(N, mu = rep(0, Jp), Sigma = diag(1, Jp)))  # J*N
        y[pind, ] <- ys
        const$cand_std <- const$cand_thd/mean(noc)  #cand sd for td
        const$inf <- inf
        const$mnoc <- mnoc
        const$zind <- zind
        ycs <- y[-pind, ]
        y[-pind, ] <- t(scale(t(ycs), center = T))  #J * N
    } else {
        y <- t(scale(t(y), center = !int))
        THD <- NULL
    }  #end Jp

    out <- list(y = y, const = const, prior = prior, PSX = PSX, PHI = PHI, LA = LA, THD = THD, gammas = gammas,
        gammal_sq = gammal_sq)

}
######## end of Init #################################################
