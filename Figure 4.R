library(mnormt)
library(Matrix)
library(ggplot2)
#========================================
S<-function(p,rho){
  M<-matrix(nrow=p,ncol=p)
  for (i in 1:p){
    for(j in 1:p){
      M[i,j]<-rho^abs(i-j)
    }
  } 
  return(M)
}
#------------------------
Covmat<-function(gamma,Sigmaxx,Sigmaxy,pvec){
  p<-sum(pvec)
  k<-length(pvec)
  rmat<-kronecker(diag(k), Sigmaxx)+  gamma^2*kronecker(matrix(rep(1,k^2),nrow = k),Sigmaxy)
  return(rmat)
}
#===================Dette 2020===========
LRT<-function(n,pvec,Data){#Data should be as a nxp matrix
  k<-length(pvec)
  pv<-cumsum(c(0,pvec))
  Sjj<-vector("list",length = k)
  S<-var(Data)
  for(i in 1:k){Sjj[[i]]<-var(Data[,(pv[i]+1):pv[i+1]])}
  A0<-(det(S)/(prod(as.numeric(lapply(Sjj, det)))))^(n/2)
  c1<-(sum(pvec)^3-sum(pvec^3))
  c2<-(sum(pvec)^2-sum(pvec^2))
  eta0<-1-(2*c1+9*c2)/(6*n*c2)
  return(-2*eta0*log(A0))
}
#------------------------------------
Tdette<-function(n,pvec,Data){#Data should be as a nxp matrix
  data<-t(Data)
  k<-length(pvec)
  y<-pvec/n
  c<-sum(y)
  sigma2<-2*log(prod(1-y)/(1-c))
  A1<-sum((n-pvec-1)*log(1-pvec/(n-1)))
  A2<-(n-sum(pvec)-1)*log(1-sum(pvec)/(n-1))
  mun<-A1-A2-0.25*sigma2
  pv<-cumsum(c(0,pvec))
  Sjj<-vector("list",length = k)
  S<-var(data)
  for(i in 1:k){Sjj[[i]]<-var(data[,(pv[i]+1):pv[i+1]])}
  A0<-(det(S)/(prod(as.numeric(lapply(Sjj, det))))) 
  Y<- -( log(A0)-mun)/sqrt(sigma2)
  return(Y)
}
#===================Hydo 2020==========
tr<- function(A) {
  n <- dim(A)[1] # get dimension of matrix
  tr <- 0 # 
  for (k in 1:n) {
    l <- A[k,k]
    tr <- tr + l
  }
  return(tr)
}
#------------------------------------------------
Sgh<-function(g,h,nvec,pvec,Data){
  k<-length(nvec)
  nv<-(nvec)
  pv<-cumsum(c(0,pvec))
  Xhj<-vector("list",length = k)
  for(i in 1:k){
    Xhj[[i]]<-Data[(pv[i]+1):pv[i+1],1:nv[i]]#---k p0x1 of size 1:ni
  }
  t<-min(g,h)
  r<-max(g,h)
  if(t==r){
    S<-(Xhj[[t]]-rowMeans(Xhj[[t]]))%*%t(Xhj[[t]]-rowMeans(Xhj[[t]]))/(nvec[t]-1)
  }else{
    S<-(Xhj[[t]]-rowMeans(Xhj[[t]]))%*%t((Xhj[[r]][,1:nv[t]]-rowMeans(Xhj[[r]][,1:nv[t]])))/(nvec[t]-1)
  }
  return(S)
}
#-------------------------------------- 
SIGMAhatF2<-function(g,h,nvec,pvec,Data){
  t<-min(g,h)
  r<-max(g,h)
  w<-(nvec[t]-1)^2/((nvec[t]-2)*(nvec[t]+1))*(norm(Sgh(t,r,nvec,pvec,Data), "F")^2- tr(Sgh(t,t,nvec,pvec,Data))*tr(Sgh(r,r,nvec,pvec,Data))/(nvec[t]-1))  
  return(w)
}
#---------------------------------------------------------------
HRVgh<-function(g,h,nvec,pvec,Data){SIGMAhatF2(g,h,nvec,pvec,Data)/(sqrt(SIGMAhatF2(g,g,nvec,pvec,Data))*sqrt(SIGMAhatF2(h,h,nvec,pvec,Data)))}
#-------------------------------------------------------------------------------------
THydo<-function(nvec,pvec,Data){
  k<-length(nvec)
  Tstat<-0
  sigma2<-0
  for(g in 1:(k-1)){
    for(h in (g+1):k){Tstat<-Tstat+HRVgh(g,h,nvec,pvec,Data)}
    sigma2<-sigma2+2*(k-g)/nvec[g]^2
  }
  return(Tstat/sqrt(sigma2))
}
#================Dariush 2021===================
Dt<-function(U){if(is.matrix(U)){det(U)}else{abs(U)}}
#----------------------------
Par<-function(n,pvec){
  p<-sum(pvec)
  S4<-p^4-sum(pvec^4)
  S3<-p^3-sum(pvec^3)
  S2<-p^2-sum(pvec^2)
  d1<-2*S3+9*S2
  d2<-6*n*S2
  eta<-1-(d1/d2)
  f<-0.5*S2
  nu<-((1/48)*S4-(5/96)*S2-(S3^2)/(72*S2))/(n^2*eta^2)
  return(c(eta,f,nu))
}
#-----------------------
Rmat<-function(c,pvec){
  k<-length(pvec)
  M<-vector("list",length = k)
  E<-mapply(rnorm, n=c*pvec)
  for(i in 1:k){M[[i]]<-matrix(E[,i],ncol = c)}
  rmat<-.bdiag(M)
  return(rmat)#======px(kc)=pxr matrix
}
#-----------------------------
Tnt<-function(n,k,eta,c,pvec,Data){
  #data is a pxn matrix
  R<-Rmat(c,pvec)
  Y<-t(as.matrix(t(R)%*%Data))#nxr dim(Y)
  pv<-cumsum(c(0,rep(c,k)))
  Yj<-vector("list",length = k)
  for(i in 1:k){Yj[[i]]<-Y[,(pv[i]+1):pv[i+1] ]}#nxc dim(Yj[[1]])
  SR<-Dt(var(Y))
  SRD<-prod(as.numeric(lapply(lapply(Yj,var),Dt)))
  A1<-(SR/SRD)
  out<- -eta*n*log(A1)
  return(out)
}
#-----------------------------------
Tnm<-function(n,k,eta,c,m,pvec,Data){
  x<-sapply(1:m, function(x){Tnt(n,k,eta,c,pvec,Data)})
  return(x)
}
#------------Bao et al(2017)---------------------------
#------------Bao et al(2017)---------------------------
#------------Bao et al(2017)---------------------------
NegRootMat<-function(A){
  e <- eigen(A)
  V <- e$vectors
  B <- V %*%diag((e$values)^(-0.5))%*%t(V)
  return(B)
}
#------------------------------------------------------
TBao<-function(n,k,pvec,Data){
  p<-sum(pvec)
  X=Data-matrix(rep(rowMeans(Data),n),nrow=p,byrow=F)
  pv<-c(0,cumsum(pvec))
  Xi<-vector("list",length = k)
  for(i in 1:k){Xi[[i]]<-X[(1+pv[i]):pv[i+1], ]}
  M<-lapply(1:k,function(i){NegRootMat(Xi[[i]]%*%t(Xi[[i]]))})
  D<-.bdiag(M)
  B<-D%*%(X%*%t(X))%*%D
  SB<-0.5*sum(diag(B%*%B))-0.5*p
  Pmat<-pvec%*%t(pvec)
  an<-2*sum(Pmat[lower.tri(Pmat,diag = F)])/(2*(n-1))
  PNmat<-(pvec*(n-1-pvec))%*%t(pvec*(n-1-pvec))
  bn<- 2*sum(PNmat[lower.tri(PNmat,diag = F)])/((n-1)^4)
  TB<-(SB-an)/sqrt(bn)
  return(TB)
}
#----------------------------------------------
Sim<-function(rep,alpha,n,p,k,gamma){
  pvec<-rep(p/k,k)
  #SigmaMat<-Covmat(gamma=0.3,Sigmaxx=diag(p/k),Sigmaxy=diag(p/k),pvec)
  SigmaMat<-Covmat(gamma,Sigmaxx=S(p/k,rho=-0.95),Sigmaxy=S(p/k,rho=0.95),pvec)
  p<-sum(pvec)
  k<-length(pvec)
  #c<-floor(min(n/(2*k),pvec/k))
  c<-1
  nvec<-cumsum(rep(n/k,k))
  AsPar<-Par(n,rep(c,k))
  eta<-AsPar[1]
  f<-AsPar[2]
  nu<-abs(AsPar[3])
  z<-qnorm(1-(alpha),lower.tail = TRUE)
  Tmndist<-function(y,m){nu*pchisq(y,f+4)+(1-nu)*pchisq(y,f)-(1-alpha)^(1/m)}
  SimvalHyd<-SimvalDet<-SimvalBao<-numeric(rep)
  SimvalDar5<-SimvalDar25<-SimvalDar50<-SimvalDar100<-SimvalDar150<-numeric(rep)
  Quant5<-uniroot(Tmndist, c(0,10*f*(1+nu)),m=5,extendInt="yes")$root
  Quant25<-uniroot(Tmndist, c(0,10*f*(1+nu)),m=25,extendInt="yes")$root
  Quant50<-uniroot(Tmndist, c(0,10*f*(1+nu)),m=50,extendInt="yes")$root
  Quant100<-uniroot(Tmndist, c(0,10*f*(1+nu)),m=100,extendInt="yes")$root
  Quant150<-uniroot(Tmndist, c(0,10*f*(1+nu)),m=150,extendInt="yes")$root
  if(prod(n>pvec)==1 & n>p){
    for(j in 1:rep){
      Data<-t(rmnorm(n,numeric(p),SigmaMat))
      TntVals5<-Tnm(n,k,eta,c,m=5,pvec,Data)
      TntVals25<-Tnm(n,k,eta,c,m=25,pvec,Data)
      TntVals50<-Tnm(n,k,eta,c,m=50,pvec,Data)
      TntVals100<-Tnm(n,k,eta,c,m=100,pvec,Data)
      TntVals150<-Tnm(n,k,eta,c,m=150,pvec,Data)
      SimvalDar5[j]<-as.numeric(max(TntVals5)>=Quant5)
      SimvalDar25[j]<-as.numeric(max(TntVals25)>=Quant25)
      SimvalDar50[j]<-as.numeric(max(TntVals50)>=Quant50)
      SimvalDar100[j]<-as.numeric(max(TntVals100)>=Quant100)
      SimvalDar150[j]<-as.numeric(max(TntVals150)>=Quant150)
      h<-THydo(nvec,pvec,Data)
      SimvalHyd[j]<-as.numeric(h>z)
      d<-Tdette(n,pvec,Data)
      SimvalDet[j]<-as.numeric(d>z)
      Bao<-TBao(n,k,pvec,Data)
      SimvalBao[j]<-as.numeric(Bao>z)
    }
    out<-c(mean(SimvalDar5),mean(SimvalDar25),mean(SimvalDar50),mean(SimvalDar100),mean(SimvalDar150),mean(SimvalDet),mean(SimvalHyd),mean(SimvalBao))
  }
  if(prod(n>pvec)==1 & n<p){
    for(j in 1:rep){
      Data<-t(rmnorm(n,numeric(p),SigmaMat))
      TntVals5<-Tnm(n,k,eta,c,m=5,pvec,Data)
      TntVals25<-Tnm(n,k,eta,c,m=25,pvec,Data)
      TntVals50<-Tnm(n,k,eta,c,m=50,pvec,Data)
      TntVals100<-Tnm(n,k,eta,c,m=100,pvec,Data)
      TntVals150<-Tnm(n,k,eta,c,m=150,pvec,Data)
      SimvalDar5[j]<-as.numeric(max(TntVals5)>=Quant5)
      SimvalDar25[j]<-as.numeric(max(TntVals25)>=Quant25)
      SimvalDar50[j]<-as.numeric(max(TntVals50)>=Quant50)
      SimvalDar100[j]<-as.numeric(max(TntVals100)>=Quant100)
      SimvalDar150[j]<-as.numeric(max(TntVals150)>=Quant150)
      h<-THydo(nvec,pvec,Data)
      SimvalHyd[j]<-as.numeric(h>z)
      #d<-Tdette(n,pvec,Data)
      #SimvalDet[j]<-as.numeric(d>z)
      Bao<-TBao(n,k,pvec,Data)
      SimvalBao[j]<-as.numeric(Bao>z)
    }
    out<-c(mean(SimvalDar5),mean(SimvalDar25),mean(SimvalDar50),mean(SimvalDar100),mean(SimvalDar150),NA,mean(SimvalHyd),mean(SimvalBao))
  }
  if(prod(n>pvec)==0){
    for(j in 1:rep){
      Data<-t(rmnorm(n,numeric(p),SigmaMat))
      TntVals5<-Tnm(n,k,eta,c,m=5,pvec,Data)
      TntVals25<-Tnm(n,k,eta,c,m=25,pvec,Data)
      TntVals50<-Tnm(n,k,eta,c,m=50,pvec,Data)
      TntVals100<-Tnm(n,k,eta,c,m=100,pvec,Data)
      TntVals150<-Tnm(n,k,eta,c,m=150,pvec,Data)
      SimvalDar5[j]<-as.numeric(max(TntVals5)>=Quant5)
      SimvalDar25[j]<-as.numeric(max(TntVals25)>=Quant25)
      SimvalDar50[j]<-as.numeric(max(TntVals50)>=Quant50)
      SimvalDar100[j]<-as.numeric(max(TntVals100)>=Quant100)
      SimvalDar150[j]<-as.numeric(max(TntVals150)>=Quant150)
      h<-THydo(nvec,pvec,Data)
      SimvalHyd[j]<-as.numeric(h>z)
      #d<-Tdette(n,pvec,Data)
      #SimvalDet[j]<-as.numeric(d>z)
      #Bao<-TBao(n,k,pvec,Data)
      #SimvalBao[j]<-as.numeric(Bao>z)
    }
    out<-c(mean(SimvalDar5),mean(SimvalDar25),mean(SimvalDar50),mean(SimvalDar100),mean(SimvalDar150),NA,mean(SimvalHyd),NA)
  }
  names(out)<-c("T5","T25","T50","T100","T150","TD","TH","TB")
  return(out)
}
#------------------------------------------------
#===================================================
gamma=seq(0,2,0.05)
RES<-mapply(Sim, rep=5000,alpha=0.05,n=180,p=150,k=2,gamma)
RES
#save.image(file = "C:/Users/Dariu/OneDrive/Desktop/Figure4.RData")

mm<-length(gamma)
G<-rep(gamma,8)
Tests<-c(rep("m=5",mm),rep("m=25",mm),rep("m=50",mm),rep("m=100",mm),rep("m=150",mm),rep("Dette",mm),rep("Hyodo",mm),rep("Bao",mm))
Power<-as.vector(t(RES))
Simdata<-data.frame(G,Tests,Power)


gg.mc10 <- ggplot(Simdata, aes(x = G, y = Power, colour = Tests,lty= Tests)) + geom_line(size=1) 
gg.mc20<-gg.mc10+labs(x=expression(~gamma),y=expression(Power)) 
gg.mc20
gg.mc30<-gg.mc20+   theme_bw( base_size = 14)+ theme(legend.title = element_blank(),
                                                     legend.spacing.y = unit(0, "mm"), 
                                                     panel.border = element_rect(colour = "black", fill=NA),
                                                     axis.text = element_text(colour = 1, size = 12),
                                                     legend.background = element_blank(),
                                                     legend.box.background = element_rect(colour = "black"))
gg.mc30
#scale_colour_manual(values=1:7,labels = expression(hat(phi)[T[D]],hat(phi)[T[H]],hat(phi)[100],hat(phi)[150],hat(phi)[25 ],hat(phi)[5],hat(phi)[50]))

library(ggforce)
library(tidyverse)
theme_set(theme_bw(16))
gg.mc30+facet_zoom(G <0.25,Power>0.04 & Power< 0.08,horizontal = F,show.area=F,zoom.size =0.4,shrink=T)
#facet_zoom(ylim = c(0.9,1),zoom.size = 0.5)















