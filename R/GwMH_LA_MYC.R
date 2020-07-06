################## update Loading ##############################################################
GwMH_LA_MYC <- function(y, mu = 0, ome, la, psx, gammal_sq, thd, const, prior, alas) {
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

            td1 <- thd
            mnoc <- const$mnoc
            cand_std <- const$cand_std
            pind <- const$cati
            inf <- const$inf
            zind <- const$zind
            # sel <- matrix(0,Jp,mnoc-1)
            for (m in 2:mnoc) {
                # m<-5 td1[,m]<-rtnorm(J,mean=td0[,m],sd=cand.std,lower=td1[,m-1],upper=td0[,m+1])
                tmp <- rnorm(Jp, mean = thd[, m], sd = cand_std)
                # sel[,m-1]<-(tmp>td1[,m-1]&tmp<=thd[,m+1]) td1[,m]<-tmp*sel[,m-1]+thd[,m]*(1-sel[,m-1])
                sel <- (tmp > td1[, m - 1] & tmp <= thd[, m + 1])
                td1[, m] <- tmp * sel + thd[, m] * (1 - sel)

                # tmp0<-pnorm((td0[,m+1]-td0[,m])/cand_std)-pnorm((td1[,m-1]-td0[,m])/cand_std) tmp0[tmp0<=0]<-1
                # tmp1<-pnorm((td1[,m+1]-td1[,m])/cand_std)-pnorm((td0[,m-1]-td1[,m])/cand_std) tmp1[tmp1<=0]<-1
                # ptd0[,m-1]<-log(tmp0) ptd1[,m-1]<-log(tmp1)
            }

            tmp00 <- matrix(thd[1:Jp + (zind - 1) * Jp], Jp, N)
            tmp01 <- matrix(thd[1:Jp + (zind) * Jp], Jp, N)
            tmp10 <- matrix(td1[1:Jp + (zind - 1) * Jp], Jp, N)
            tmp11 <- matrix(td1[1:Jp + (zind) * Jp], Jp, N)

            # yst<-la[pind,]%*%ome+mst[pind,] spsx<-sqrt(psxjs[pind])
            yst <- ysta[pind, ]
            spsx <- spsxa[pind]
            tmp1 <- log(pnorm((tmp11 - yst)/spsx) - pnorm((tmp10 - yst)/spsx))
            tmp1[tmp1 < (-inf)] <- -inf
            tmp0 <- log(pnorm((tmp01 - yst)/spsx) - pnorm((tmp00 - yst)/spsx))
            tmp0[tmp0 < (-inf)] <- -inf

            # acc<-exp(rowSums(tmp1-tmp0)+rowSums(ptd0-ptd1))
            acc <- exp(rowSums(tmp1 - tmp0, na.rm = T))
            accind <- (acc > runif(Jp))
            thd[accind, ] <- td1[accind, ]
            # accrate<-mean(accind)

            tmp00 <- matrix(thd[1:Jp + (zind - 1) * Jp], Jp, N)
            tmp01 <- matrix(thd[1:Jp + (zind) * Jp], Jp, N)

            # ys<-matrix(rnorm(N*Jp),Jp,N)+yst/spsx # sdy<-apply(ys,1,sd) # ys<-ys/sdy #using with
            # #ys[j,]<-rnorm(N,tmp,sqrt(convar[j]))# ys<-ys/apply(ys,1,sd) # ys <- t(scale(t(ys), center = F,
            # scale = apply(ys, 1, sd, na.rm = T)))

            ys <- ysa[pind, ]
            acc1 <- ((ys > tmp00) & (ys <= tmp01))
            ys <- ys * acc1 + (1 - acc1) * y[pind, ]
            # accr<-c(mean(accind),mean(acc1),mean(sel),mean(ptd0),mean(ptd1))
            accr <- c(mean(accind), mean(acc1, na.rm = T))

            # out <- list(la = la, gammal_sq = gammal_sq, ome = ome, ys = ys, thd = thd, accr = accr, psx = psx,
            # ysm = ysa)
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
