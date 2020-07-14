################## update Loading ##############################################################
GwMH_LA_MYE0 <- function(y, mu = 0, ome, la, psx, gammal_sq, thd, const, prior, alas) {
    # y=Y;ome=OME;la=LA;psx=PSX;mu=0;thd=THD
    Q <- const$Q
    J <- const$J
    N <- const$N
    K <- const$K
    Jp <- const$Jp
    Nmis <- const$Nmis
    a_gamma <- prior$a_gaml_sq
    b_gamma <- prior$b_gaml_sq
    Pmean <- prior$m_LA
    Sigla <- prior$s_LA
    sub_sl <- const$sub_sl
    len_sl <- const$len_sl
    sub_ul <- const$sub_ul
    len_ul <- const$len_ul
    taul_sq <- gammal_sq
    a_gams <- prior$a_gams
    b_gams <- prior$b_gams

    temp <- y - mu - la %*% ome  # NY*N
    S <- temp %*% t(temp)  # NY*NY

    for (j in 1:J) {
        # j=1 for specified loadings
        subs <- sub_sl[j, ]
        len <- len_sl[j]
        if (len > 0)
            {
                # yj<-(y[j,]-mst[j,])-matrix(la[j,(!subs)],nrow=1)%*%matrix(ome[(!subs),],ncol=N)
                yj <- y[j, ] - matrix(la[j, (!subs)], nrow = 1) %*% matrix(ome[(!subs), ], ncol = N)
                yj <- as.vector(yj)  # vector
                if (len == 1) {
                  omesub <- matrix(ome[subs, ], nrow = 1)
                } else {
                  omesub <- ome[subs, ]
                }

                PSiginv <- diag(len) * Sigla
                # diag(PSiginv)<-rep(Sigly,len) Pmean<-PLA[j,subs]
                vtmp <- chol2inv(chol(tcrossprod(omesub)/psx[j, j] + PSiginv))
                mtmp <- (omesub %*% yj/psx[j, j] + PSiginv %*% rep(Pmean,len))
                la[j, subs] <- mvrnorm(1, vtmp %*% mtmp, Sigma = vtmp)
            }  # end len>0

        # subs<-(Q[j,]==-1) len<-sum(subs) for specified loadings
        subs <- sub_ul[j, ]
        len <- len_ul[j]
        if (len > 0) {
            # yj<-(y[j,]-mst[j,])-matrix(la[j,(!subs)],nrow=1)%*%matrix(ome[(!subs),],ncol=N)
            yj <- y[j, ] - matrix(la[j, (!subs)], nrow = 1) %*% matrix(ome[(!subs), ], ncol = N)
            yj <- as.vector(yj)  # vector
            Cadj <- pmax((la[j, subs])^2, 10^(-6))
            mu_p <- pmin(sqrt(gammal_sq[j, subs]/Cadj), 10^12)
            taul_sq[j, subs] <- 1/rinvgauss1(len, mean = mu_p, dispersion = 1/gammal_sq[j, subs])
            if(alas) gammal_sq[j, subs] <- rgamma(len, shape = a_gamma + 1, rate = b_gamma + taul_sq[j, subs]/2)

            if (len == 1) {
                omesub <- matrix(ome[subs, ], nrow = 1)
                invD_tau <- 1/taul_sq[j, subs]
            } else {
                omesub <- ome[subs, ]
                invD_tau <- diag(1/taul_sq[j, subs])
            }
            vtmp <- chol2inv(chol(tcrossprod(omesub)/psx[j, j] + invD_tau))
            mtmp <- (omesub %*% yj/psx[j, j])
            la[j, subs] <- mvrnorm(1, vtmp %*% mtmp, Sigma = vtmp)

            tmp <- t(la[j, subs]) %*% invD_tau %*% la[j, subs]
            psx[j, j] <- 1/rgamma(1, shape = a_gams + (N + len)/2 - 1, rate = b_gams + (S[j, j] +
                tmp)/2)
        } else {
            psx[j, j] <- 1/rgamma(1, shape = a_gams + (N - 1)/2, rate = b_gams + (S[j, j])/2)
        }  # end len>0
    }  # end of J

    if(!alas)
        gammal_sq[Q==-1]<- rgamma(1, shape=a_gamma+sum(Q==-1), rate=b_gamma + sum(taul_sq)/2)


        out <- list(la = la, gammal_sq = gammal_sq, psx = psx)

    return(out)
}
