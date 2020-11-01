###################    update MU #################################################################
Gibbs_MU<-function(y,ome,la,inv.psx, N, J, prior){
   # y=Ys;z=Z;ome=Omega;ly=LD;inv.psx=inv.PSX;
   PMU <- prior$m_MU
   Sigmu <- prior$s_MU
   calsm<-chol2inv(chol(N*inv.psx+diag(rep(Sigmu,J)))) # inv[sigma0^(-1)+N*inv.PSX]
   Ycen<-y-la%*%ome
   temp<-rowSums(Ycen)
   mumu<-calsm%*%(inv.psx%*%temp+rep(Sigmu,J)*PMU)
   mu<-mvrnorm(1,mumu,Sigma=calsm)
   # mu<-mvrnorm(1,mu=rep(0,J),Sig=calsm)
   # mu<-as.vector(mu+mumu)
	 return(mu)
}
###################    end of update MU  ##########################################################
