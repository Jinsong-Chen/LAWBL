################## update Loading ##############################################################
Gibbs_LA_IYC <- function(y, mu, ome, la, psx, gammal_sq, thd, const, prior, alas) {
    # y=Y;ome=OME;la=LA;psx=PSX;mu=0;thd=THD;const=const;prior=prior
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

    mst <- matrix(0, J, N)
    psxjs <- rep(0, J)
    # taul_sq <- gammal_sq Pmean <- 0
    for (j in 1:J) {
        # j=1
        psx_tmp <- psx[j, -j] %*% chol2inv(chol(psx[-j, -j]))
        psxjs[j] <- psx[j, j] - psx_tmp %*% psx[-j, j]
        mst[j, ] <- psx_tmp %*% (y[-j, ] - la[-j, ] %*% ome)  # 1*N
        # subs<-(Q[j,]==1) len<-sum(subs) for specified loadings
        subs <- sub_sl[j, ]
        len <- len_sl[j]
        if (len > 0)
            {
                yj <- (y[j, ] - mst[j, ]) - matrix(la[j, (!subs)], nrow = 1) %*% matrix(ome[(!subs),
                  ], ncol = N)
                yj <- as.vector(yj)  # vector
                if (len == 1) {
                  omesub <- matrix(ome[subs, ], nrow = 1)
                } else {
                  omesub <- ome[subs, ]
                }

                PSiginv <- diag(len) * Sigla
                # diag(PSiginv)<-rep(Sigly,len) Pmean<-PLA[j,subs]
                vtmp <- chol2inv(chol(tcrossprod(omesub)/psxjs[j] + PSiginv))
                mtmp <- (omesub %*% yj/psxjs[j] + PSiginv %*% rep(Pmean,len))
                la[j, subs] <- mvrnorm(1, vtmp %*% mtmp, Sigma = vtmp)
            }  # end len>0

        # subs<-(Q[j,]==-1) len<-sum(subs) for specified loadings
        subs <- sub_ul[j, ]
        len <- len_ul[j]
        if (len > 0)
            {
                yj <- (y[j, ] - mst[j, ]) - matrix(la[j, (!subs)], nrow = 1) %*% matrix(ome[(!subs),
                  ], ncol = N)
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
                vtmp <- chol2inv(chol(tcrossprod(omesub)/psxjs[j] + invD_tau))
                mtmp <- (omesub %*% yj/psxjs[j])
                la[j, subs] <- mvrnorm(1, vtmp %*% mtmp, Sigma = vtmp)
            }  # end len>0
    }  # end of J

    if(!alas)
        gammal_sq[Q==-1]<- rgamma(1, shape=a_gamma+sum(Q==-1), rate=b_gamma + sum(taul_sq)/2)

    if (Nmis > 0 || Jp > 0) {
        # yst<-la[pind,]%*%ome+mst[pind,]
        ysta <- la %*% ome + mst
        spsxa <- sqrt(psxjs)
        ysa <- matrix(rnorm(N * J), J, N) + ysta/spsxa
        ysa <- ysa/apply(ysa, 1, sd)

        if (Jp > 0) {

            pind <- const$cati
            # inf <- const$inf
            zind <- const$zind

            ys <- ysa[pind, ]
            acc <- ((ys > 0) ==(zind>1))
            ys <- ys * acc + (1 - acc) * y[pind, ]
            # accr<-c(mean(accind),mean(acc1),mean(sel),mean(ptd0),mean(ptd1))
            accr <- c( mean(acc, na.rm = T))

            out <- list(la = la, gammal_sq = gammal_sq, ys = ys, thd = thd, accr = accr, psx = psx,
                ysm = ysa)
        } else {
            out <- list(la = la, gammal_sq = gammal_sq, psx = psx, ysm = ysa)
        }
    } else {
        out <- list(la = la, gammal_sq = gammal_sq, psx = psx)
    }  #end Nmis || Jp

    # return(list(la=la,lamsq=lamsq,ome=ome,tausq=tausq,ys=ys,sdy=sdy,sup=sup))
    return(out)
}

################## end of update Loading ########################################################
