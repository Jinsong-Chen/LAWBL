##################  update pefa main##############################################################
# Gibbs sampler with bilevel loading regularization and Spike and Slab priors
Gibbs_BLR_SSP<-function(y,mu=0,ome,ly,tausq,prior,sksq,sksq_t,const){
  #y=Y;ome=OME;ly=LA;mu=0
  Q <- const$Q
  J <- const$J
  N <- const$N
  K <- const$K
  FIS<-const$FIS
  t_num <- const$t_num
  a_gamma <- prior$a_gaml_sq
  b_gamma <- prior$b_gaml_sq
  Pmean <- prior$m_LA
  Sigla <- prior$s_LA
  sub_sl <- const$sub_sl
  len_sl <- const$len_sl

  yc <-y-mu
  bk <- tausq
  bk_pr<-rep(1,K)
  tmp <- yc - ly %*% ome  # J*N
  S <- tcrossprod(tmp)  # J*J
  dpsx<-1/rgamma(J, shape=a_gamma+(N-1)/2, rate=b_gamma+diag(S)/2)
  # psx <- diag(dpsx)

  if (any(!FIS)){
    for(j in 1:J){
    # j =1
    ind1 <- sub_sl[j, ]
    len <- len_sl[j]
    if(len>0){

      yj <- as.vector(yc[j,]) # vector
      # yj <- as.vector(yc[j, ] - matrix(ly[j, (!ind1)], nrow = 1) %*% matrix(ome[(!ind1), ], ncol = N))
      if(len==1){
        # omesub<-matrix(ome[subs,],nrow=1)
        omesub<-t(ome[ind1,])
      }else{omesub<-ome[ind1,]}
      PSiginv<-Sigla*diag(len)
      # Pmean<-PLY[j,ind1]

      vtmp<-chol2inv(chol(tcrossprod(omesub)/dpsx[j]+PSiginv))
      # mtmp<-(omesub%*%as.vector(yc[j,])/dpsx[j]+PSiginv%*%rep(Pmean,len))
      mtmp<-(omesub%*%yj/dpsx[j]+PSiginv%*%rep(Pmean,len))
      ly[j,ind1]<-mvrnorm(1,vtmp%*%mtmp,Sigma = vtmp)
    } # end len>0
  }# end of j
  }#end if


  for(k in 1:K){
    #k=1

    Vh<-diag(sqrt(tausq[,k]))
    Sig<-chol2inv(chol(Vh%*%Vh*sum(ome[k,]^2)/dpsx+diag(J)))

    if(FIS[k]){

        # Dk <- Vh%*%as.vector(ome[k,]%*%t(yc-ly[,!FIS]%*%ome[!FIS,]))/dpsx
        Dk <- Vh%*%as.vector(ome[k,]%*%t(yc))/dpsx
        # Dk <-Vh%*%as.vector(ome[k,]%*%t(yc-ly[,-k]%*%ome[-k,]))/dpsx

        tmp1<-determinant(Sig,logarithm = TRUE)$modulus
        # tmp2<-t(Dk)%*%Sig%*%Dk
        d1<-as.vector(ome[k,]%*%t(yc-ly[,-k]%*%ome[-k,]))/dpsx
        D1 <- Vh%*%d1
        tmp2<-t(D1)%*%Sig%*%D1
        tmp<-exp(.5*(tmp1+tmp2))

        # pi0<-rbeta(1,sum(1-pig[indg])+1,sum(pig[indg])+1)
        pi0<-.5
        tmp0<-pi0/(pi0+(1-pi0)*tmp)
        bk_pr[k]<-1-(tmp0>runif(1))
        if (bk_pr[k]){
          bk[,k]<- mvrnorm(1,Sig%*%Dk,Sigma=(Sig))
        }else{
          bk[,k]<- 0
        }

        vjk<-1/(1/(sksq[k]+1e-200)+sum(ome[k,]^2)*bk[,k]^2/dpsx)
        # ujk<-(yc-ly[,-k]%*%ome[-k,])%*%ome[k,]*vjk*bk[,k]/dpsx
        ujk<-d1*vjk*bk[,k]

        tmp<-.5*(log(vjk)-log(sksq[k])+ujk^2/vjk)+pnorm(ujk/sqrt(vjk),log.p=TRUE)
        tmp<-exp(tmp)

        # tmp1<-sum(tausq==0)
        # # pi1<-rbeta(1,J*K-tmp1-sum(Q!=-1)+1,tmp1+1)
        # pi1<-rbeta(1,tmp1+1,J*K-tmp1+1)

        # tmpt<-sum(tausq[,k]==0)
        # pi1<-rbeta(1,tmpt+1,J-tmpt+1)

        pi1<-.5

        tmp1 <-pi1/(pi1+2*(1-pi1)*tmp)
        tau_pr<-1-(tmp1>runif(J))

        if (sum(tau_pr) < const$mjf){
          # count[k,3]<-count[k,3]+1
          tausq[,k] <- 0
        }else{
          tausq[,k]<-rnorm(J,ujk,sqrt(vjk))^2*tau_pr #truncated normal square
          tausq[,k]<-pmin(tausq[,k],1e200)
        }

        sk.tmp<- rgamma(t_num, shape=1+sum(tausq[,k]!=0)/2, rate=sksq_t[k]+sum(tausq[,k])/2)
        sksq_t[k]<-1/mean(sk.tmp)
        sksq[k]<-1/rgamma(1, shape=1+sum(tausq[,k]!=0)/2, rate=sksq_t[k]+sum(tausq[,k])/2)

        ly[,k]<-bk[,k]*sqrt(tausq[,k])
    # } #end if FIS

      }else{
        ind<-(Q[,k]==-1)
        # Dk <-Vh%*%as.vector(ome[k,]%*%t(yc-ly[,-k]%*%ome[-k,]))/dpsx
        # Dk<-Vh%*%as.vector(ome[k,]%*%t(yc-ly[,FIS]%*%ome[FIS,]))/dpsx
        Dk<-Vh%*%as.vector(ome[k,]%*%t(yc))/dpsx
        bk[,k]<- mvrnorm(1,Sig%*%Dk,Sigma=Sig)
        # bk[ind,k]<- mvrnorm(1,Sig[ind,ind]%*%Dk[ind],Sigma=(Sig[ind,ind]))
        vjk<-1/(1/sksq[k]+sum(ome[k,]^2)*bk[,k]^2/dpsx)
        ujk<-(yc-ly[,-k]%*%ome[-k,])%*%ome[k,]*vjk*bk[,k]/dpsx

        tmp<-.5*(log(vjk)-log(sksq[k])+ujk^2/vjk)+pnorm(ujk/sqrt(vjk),log.p=TRUE)
        tmp<-exp(tmp)

        # tmp1<-sum(tausq==0)
        # # pi1<-rbeta(1,J*K-tmp1-sum(Q!=-1)+1,tmp1+1)
        # pi1<-rbeta(1,tmp1+1,J*K-tmp1+1)

        # tmpt<-sum(tausq[,k]==0)
        # pi1<-rbeta(1,tmpt+1,J-tmpt+1)

        pi1<-.5
        tmp1 <-pi1/(pi1+2*(1-pi1)*tmp)

        tau_pr<-1-(tmp1>runif(J))
        tmp <-rnorm(J,ujk,sqrt(vjk))^2*tau_pr
        tausq[ind,k]<-tmp[ind]

        # len<-sum(ind)
        # tau_pr<-1-(tmp1[ind]>runif(len))
        # tausq[ind,k]<-rnorm(len,ujk[ind],sqrt(vjk[ind]))^2*tau_pr

        ly[ind,k]<-bk[ind,k]*sqrt(tausq[ind,k])

        sk.tmp<- rgamma(t_num, shape=1+sum(tausq[,k]!=0)/2, rate=sksq_t[k]+sum(tausq[,k])/2)
        sksq_t[k]<-1/mean(sk.tmp)
        sksq[k]<-1/rgamma(1, shape=1+sum(tausq[,k]!=0)/2, rate=sksq_t[k]+sum(tausq[,k])/2)

        # sk.tmp<- rgamma(t_num, shape=1+sum(tausq[ind,k]!=0)/2, rate=sksq_t[k]+sum(tausq[ind,k])/2)
        # sksq_t[k]<-1/mean(sk.tmp)
        # sksq[k]<-1/rgamma(1, shape=1+sum(tausq[ind,k]!=0)/2, rate=sksq_t[k]+sum(tausq[ind,k])/2)

      } #end if FIS
  } # end of k

# if(any(!FIS)){
#   sk.tmp<- rgamma(t_num, shape=1+sum(tausq[,!FIS]!=0)/2, rate=mean(sksq_t[!FIS])+sum(tausq[,!FIS])/2)
#   sksq_t[!FIS] <- tmp_t <- 1/mean(sk.tmp)
#   sksq[!FIS]<-1/rgamma(1, shape=1+sum(tausq[,!FIS]!=0)/2, rate=tmp_t+sum(tausq[,!FIS])/2)
# }

    # sk.tmp<- rgamma(t_num, shape=1+sum(tausq[,FIS]!=0)/2, rate=mean(sksq_t[FIS])+sum(tausq[,FIS])/2)
    # sksq_t[FIS] <- tmp_t <- 1/mean(sk.tmp)
    # sksq[FIS]<-1/rgamma(1, shape=1+sum(tausq[,FIS]!=0)/2, rate=tmp_t+sum(tausq[,FIS])/2)

  # ilamsq<-pmin(ilamsq,1e200)

# return(list(ly=ly,lamsq=lamsq,ome=ome,pig=pig,psx=psx,tausq=tausq,lamsq_add=lamsq_add,lamsq_t=lamsq_t))
return(list(ly=ly,bk=bk,bk_pr=bk_pr,dpsx=dpsx,sksq=sksq,tausq=tausq,sksq_t=sksq_t))
}

##################  end of update bsfa main ########################################################
