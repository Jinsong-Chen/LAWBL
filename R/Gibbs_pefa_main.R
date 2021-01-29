##################  update bsfa main##############################################################
Gibbs_pefa_main<-function(y,mu=0,ome,ly,psx, tausq,pig,prior,ilamsq,ilamsq_t,const){
  Q <- const$Q
  J <- const$J
  N <- const$N
  K <- const$K
  t_num <- const$t_num
  mjf <- const$mjf

  Ycen<-y-mu
  indg<-(pig!=-1)
  lyb<-tausq
  lyb[tausq!=0]<-ly[tausq!=0]/sqrt(tausq[tausq!=0])
  a_gamma <- prior$a_gaml_sq
  b_gamma <- prior$b_gaml_sq

  Pmean <- prior$m_LA
  Sigla <- prior$s_LA
  sub_sl <- const$sub_sl
  len_sl <- const$len_sl
  # sub_ul <- const$sub_ul
  # len_ul <- const$len_ul

  temp <- y - mu - ly %*% ome  # J*N
  S <- temp %*% t(temp)  # J*J

  for(j in 1:J){

    # tmp<-Ycen[j,]-ly[j,]%*%ome
    # S<-tmp%*%t(tmp)

    psx[j,j]<-1/rgamma(1, shape=a_gamma+(N-1)/2, rate=b_gamma+(S[j,j])/2)

    # ind1<-(Q[j,]==1)
    # len<-sum(ind1)
    ind1 <- sub_sl[j, ]
    len <- len_sl[j]
    if(len>0){

      yj<-as.vector(Ycen[j,]) # vector
      if(len==1){
        # omesub<-matrix(ome[subs,],nrow=1)
        omesub<-t(ome[ind1,])
      }else{omesub<-ome[ind1,]}
      PSiginv<-Sigla*diag(len)
      # Pmean<-PLY[j,ind1]

      vtmp<-chol2inv(chol(tcrossprod(omesub)/psx[j,j]+PSiginv))
      mtmp<-(omesub%*%yj/psx[j,j]+PSiginv%*%rep(Pmean,len))
      # LYnpsx<-calsmnpsx%*%temp
      ly[j,ind1]<-mvrnorm(1,vtmp%*%mtmp,Sig=vtmp)

    } # end len>0

  }# end of j

  for(k in 1:K){
    #k=4
    # if(indg[k]){

        vgh<-diag(sqrt(tausq[,k]))
        sigg<-chol2inv(chol(vgh%*%vgh*sum(ome[k,]^2)/diag(psx)+diag(J)))
        mug<-sigg%*%vgh%*%as.vector(ome[k,]%*%t(Ycen))/diag(psx)

        lyb[,k]<- mvrnorm(1,mug,Sig=(sigg))

      if(indg[k]){

        m1<-(sigg)%*%vgh%*%as.vector(ome[k,]%*%t(Ycen-ly[,!indg]%*%ome[!indg,]))/diag(psx)
        tmp<-log(det(sigg))/2+t(m1)%*%chol2inv(chol(sigg))%*%m1/2

        tmp<-exp(tmp)

        p0<-rbeta(1,sum(1-pig[indg])+1,sum(pig[indg])+1)
        # p0<-.5

        pg<-p0/(p0+(1-p0)*tmp)

        pig[k]<-(1-pg>runif(1))
        # } #end if eigen
        lyb[,k]<- pig[k]*lyb[,k]

        vg2<-pmin(1/(1/ilamsq[k]+sum(ome[k,]^2)*lyb[,k]^2/diag(psx)),10^12)
        # vg2<-1/(1/ilamsq[,k]+sum(ome[k,]^2)*lyb[,k]^2/diag(psx))
        ug<-(Ycen-ly[,-k]%*%ome[-k,])%*%ome[k,]*vg2*lyb[,k]/diag(psx)
        # ug<-(Ycen)%*%ome[k,]*vg2*lyb[,k]/diag(psx)

        tmp<-(log(vg2)-log(ilamsq[k])+ug^2/vg2+2*pnorm(ug/sqrt(vg2)))/2
        # tmp<-(log(vg21)-log(ilamsq[k])+ug1^2/vg21+2*pnorm(ug1/sqrt(vg21)))/2
        tmp<-exp(tmp)

        tmp1<-sum(tausq==0)
        # pi1<-rbeta(1,J*K-tmp1-sum(Q!=-1)+1,tmp1+1)
        pi1<-rbeta(1,tmp1+1,J*K-tmp1+1)
        # pi1<-.5

        qg<-pi1/(pi1+2*(1-pi1)*tmp)
        pij<-1-(qg>runif(J))

        # if (!pig[k] || sum(pij) <= 3){
        if (sum(pij) < mjf){
          tausq[,k] <- 0
        }else{
          # tausq[,k]<-(rnorm(J,ug,vg2))^2*pij #truncated normal
          tausq[,k]<-(rtnorm(J,ug,vg2,lower=0))^2*pij #truncated normal
          # tausq[,k]<-(rtnorm(J,ug,diag(psx)*vg2,lower=0))^2*pij #truncated normal
        }


        ly[,k]<-lyb[,k]*sqrt(tausq[,k])
        # lamsq[,k]<-rgamma(1, shape=a_lamsq+sum(tausq[,k]!=0), rate=b_lamsq+sum(tausq[,k])/2)

        lam.tmp<- rgamma(t_num, shape=1+sum(tausq[,k]!=0)/2, rate=ilamsq_t[k]+sum(tausq[,k])/2)
        ilamsq_t[k]<-1/mean(lam.tmp)
        ilamsq[k]<-1/rgamma(1, shape=1+sum(tausq[,k]!=0)/2, rate=ilamsq_t[k]+sum(tausq[,k])/2)


      }else{

          vg2<-pmin(1/(1/ilamsq[k]+sum(ome[k,]^2)*lyb[,k]^2/diag(psx)),10^12)
          # vg2<-1/(1/ilamsq[,k]+sum(ome[k,]^2)*lyb[,k]^2/diag(psx))
          ug<-(Ycen-ly[,-k]%*%ome[-k,])%*%ome[k,]*vg2*lyb[,k]/diag(psx)

        tmp<-(log(vg2)-log(ilamsq[k])+ug^2/vg2+2*pnorm(ug/sqrt(vg2)))/2
        # tmp<-(log(vg2)-log(ilamsq[,k])+ug^2/vg2+2*pnorm(ug/sqrt(vg2)))/2
        tmp<-exp(tmp)

        tmp1<-sum(tausq==0)
        # pi1<-rbeta(1,J*K-tmp1-sum(Q!=-1)+1,tmp1+1)
        pi1<-rbeta(1,tmp1+1,J*K-tmp1+1)
        # pi1<-.5

        qg<-pi1/(pi1+2*(1-pi1)*tmp)
        pij<-1-(qg>runif(J))
        tmp2<-(rtnorm(J,ug,vg2,lower=0))^2*pij #truncated normal
        # tmp2<-(rtnorm(J,ug,diag(psx)*vg2,lower=0))^2*pij #truncated normal

        ind<-(Q[,k]==-1)
        len<-sum(ind)
        tausq[ind,k]<-tmp2[ind]
        ly[ind,k]<-lyb[ind,k]*sqrt(tausq[ind,k])
        # ly[,k]<-lyb[,k]*sqrt(tausq[,k])

      } #end if indg

  } # end of k

ly[ly>.999]<-.999
ly[ly< -.999]<--.999
# ly[abs(ly)< 10^(-12)]<-0


if(sum(!indg)>0){
  lam.tmp<- rgamma(t_num, shape=1+sum(tausq[,!indg]!=0)/2, rate=mean(ilamsq_t[!indg])+sum(tausq[,!indg])/2)
  ilamsq_t[!indg] <- tmp_t <- 1/mean(lam.tmp)
  ilamsq[!indg]<-1/rgamma(1, shape=1+sum(tausq[,!indg]!=0)/2, rate=tmp_t+sum(tausq[,!indg])/2)
}

# return(list(ly=ly,lamsq=lamsq,ome=ome,pig=pig,psx=psx,tausq=tausq,lamsq_add=lamsq_add,lamsq_t=lamsq_t))
return(list(ly=ly,ilamsq=ilamsq,pig=pig,psx=psx,tausq=tausq,ilamsq_t=ilamsq_t))
}

##################  end of update bsfa main ########################################################
