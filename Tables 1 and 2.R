#=========================================
.bdiag<-function(M){
  mats<-vector("list",length=length(M))
  for (i in 1:length(M)) {
    mats[[i]]<-matrix(0,nrow =dim(M[[i]])[1] ,ncol=dim(M[[i]])[2])
  }
  if(length(M)==2){
    
    mats[[1]]<-matrix(0,nrow =dim(M[[1]])[1] ,ncol=dim(M[[2]])[2])
    mats[[2]]<-matrix(0,nrow =dim(M[[2]])[1] ,ncol=dim(M[[1]])[2])
    Q<-rbind(do.call(cbind, list(M[[1]],mats[[1]])),
             do.call(cbind, list(mats[[2]],M[[2]])))
  }
  if(length(M)==3){
    mats[[1]]<-matrix(0,nrow =dim(M[[1]])[1] ,ncol=sum(dim(M[[2]])[2]+dim(M[[3]])[2]))
    mats[[2]]<-matrix(0,nrow =dim(M[[2]])[1] ,ncol=dim(M[[1]])[2])
    mats[[3]]<-matrix(0,nrow =dim(M[[2]])[1] ,ncol=dim(M[[3]])[2])
    mats[[4]]<-matrix(0,nrow =dim(M[[3]])[1] ,ncol=sum(dim(M[[1]])[2]+dim(M[[2]])[2]))
    
    Q<-rbind(do.call(cbind, list(M[[1]],mats[[1]])),
             do.call(cbind, list(mats[[2]],M[[2]],mats[[3]])),
             do.call(cbind, list(mats[[4]],M[[3]])))
  }
  stopifnot(length(M)!=2|length(M)!=3)
  return(Q)
}
#mat1 <- matrix(1:4, nrow = 2)      # 2x2 matrix
#mat2 <- matrix(7:15, nrow = 3)     # 3x3 matrix
#mat3 <- matrix(16:27, nrow = 4)    # 3x3 matrix
#.bdiag(M<-list(mat1,mat2,mat3))
#=========================================
rmnorm<-function(n,mu,Sigma) {
  # generate n random vectors from MVN(mu, Sigma)
  # dimension is inferred from mu and Sigma
  d <- length(mu)
  Q <- chol(Sigma) # Choleski factorization of Sigma
  Z <- matrix(rnorm(n*d), nrow=n, ncol=d)
  X <- Z %*% Q + matrix(mu, n, d, byrow=TRUE)
  return(X)
}
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
#=============Bao et al(2017)===========================
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
#=============Jiang(2013)===========================
tr<-function(M){sum(diag(M))}
#---------------------------------------------------
HatMat<-function(M){t(M)%*%solve(M%*%t(M))%*%M}
#---------------------------------------------------
Jiang<-function(n,k,pvec,Data){
  p<-sum(pvec)
  X<-Data
  pv<-c(0,cumsum(pvec))
  k<-length(pvec)
  Yimad<-Pi<-list()
  for(i in 1:k){Yimad[[i]]<-X[(1+pv[i]):pv[k+1], ]
                Pi[[i]]<-HatMat(X[(1+pv[i]):pv[i+1], ])
                }
  Fival<-Mival<-list()
  egi<-rn1i<-rn2i<-vgi<-NULL
  for(i in 2:k){
    pbasti<-sum(pvec[1:(i-1)])
    Part1<-Yimad[[i]]%*%Pi[[i-1]]%*%t(Yimad[[i]])
    Part2<-solve(Yimad[[i]]%*%(diag(n)-Pi[[i-1]])%*%t(Yimad[[i]]))
    Fival[[i]]<- ((n-pbasti)/pbasti)*Part1%*%Part2
    Mival[[i]]<-Fival[[i]]%*%solve(Fival[[i]]+(n-pbasti)/(pbasti)*diag(p-pbasti))
    
    rn1i[i]<-pvec[i]/pbasti
    rn2i[i]<-pvec[i]/(n-pbasti)
    egi[i]<-rn2i[i]/(rn1i[i]+rn2i[i])
    
    h<-sqrt(rn1i[i]+rn2i[i]-rn1i[i]*rn2i[i])
    vgi[i]<-(2*(h*rn1i[i]*rn2i[i])^2)/((rn1i[i]+rn2i[i])^4)
  }
  LNQ<-(sum(unlist(lapply(Mival,tr)))-sum(pvec*egi,na.rm = T))/sqrt(sum(vgi,na.rm = T))
  return(LNQ)
}
#----------------------------------------------
EmpSize<-function(reppower,repempinull,rovec,n,p,k,gamma){
  pvec<-rep(p/k,k)
  SigmaMat0<-Covmat(gamma=0,Sigmaxx=S(p/k,rho=rovec[1]),Sigmaxy=S(p/k,rho=rovec[2]),pvec)
  SigmaMat<-Covmat(gamma,Sigmaxx=S(p/k,rho=rovec[1]),Sigmaxy=S(p/k,rho=rovec[2]),pvec)
  p<-sum(pvec)
  k<-length(pvec)
  c<-1
  AsPar<-Par(n,rep(c,k))
  eta<-AsPar[1]
  f<-AsPar[2]
  nu<-abs(AsPar[3])
  
  nvec<-cumsum(rep(n/k,k))
  SimvalJia<-SimvalHyd<-SimvalDet<-SimvalBao<-numeric(max(reppower,repempinull))
  SimvalDar5<-SimvalDar25<-SimvalDar50<-SimvalDar100<-SimvalDar150<-numeric(max(reppower,repempinull))
  EmpQuant5=EmpQuant25=EmpQuant50=EmpQuant100=EmpQuant150=EmpQuantHyd=EmpQuantDet=EmpQuantBao=EmpQuantJia=NA
#-------------------------------------------------------------------------------
  z<-qnorm(0.95,lower.tail = TRUE)
  Tmndist<-function(y,m){nu*pchisq(y,f+4)+(1-nu)*pchisq(y,f)-(1-0.05)^(1/m)}
  Quant5<-uniroot(Tmndist, c(0,10*f*(1+nu)),m=5,extendInt="yes")$root
  Quant25<-uniroot(Tmndist, c(0,10*f*(1+nu)),m=25,extendInt="yes")$root
  Quant50<-uniroot(Tmndist, c(0,10*f*(1+nu)),m=50,extendInt="yes")$root
  Quant100<-uniroot(Tmndist, c(0,10*f*(1+nu)),m=100,extendInt="yes")$root
  Quant150<-uniroot(Tmndist, c(0,10*f*(1+nu)),m=150,extendInt="yes")$root
#-------------------------------------------------------------------------------
  SizeHyd<-SizeDet<-SizeBao<-SizeJia<-NA
  SizeDar5<-SizeDar25<-SizeDar50<-SizeDar100<-SizeDar150<-NA
#--------------------------------------------------------------------------------
  if(prod(n>pvec)==1 & n>p){
    for(j in 1:repempinull){
      Data<-t(rmnorm(n,numeric(p),SigmaMat0))
      SimvalDar5[j]<-max(Tnm(n,k,eta,c=1,m=5,pvec,Data))
      SimvalDar25[j]<-max(Tnm(n,k,eta,c=1,m=25,pvec,Data))
      SimvalDar50[j]<-max(Tnm(n,k,eta,c=1,m=50,pvec,Data))
      SimvalDar100[j]<-max(Tnm(n,k,eta,c=1,m=100,pvec,Data))
      SimvalDar150[j]<-max(Tnm(n,k,eta,c=1,m=150,pvec,Data))
      SimvalHyd[j]<-THydo(nvec,pvec,Data)
      SimvalDet[j]<-Tdette(n,pvec,Data)
      SimvalBao[j]<-TBao(n,k,pvec,Data)
      SimvalJia[j]<-Jiang(n,k,pvec,Data)
    }
    #--------------------------------------------
    SizeDar5<-mean(SimvalDar5>Quant5)
    SizeDar25<-mean(SimvalDar25>Quant25)
    SizeDar50<-mean(SimvalDar50>Quant50)
    SizeDar100<-mean(SimvalDar100>Quant100)
    SizeDar150<-mean(SimvalDar150>Quant150)
  
    SizeHyd<-mean(SimvalHyd>z)
    SizeDet<-mean(SimvalDet>z)
    SizeBao<-mean(SimvalBao>z)
    SizeJia<-mean(SimvalJia>z)
    
    out0<-c(SizeDar5,SizeDar25,SizeDar50,SizeDar100,SizeDar150,SizeDet,SizeHyd,SizeBao,SizeJia)
    names(out0)<-c("T5","T25","T50","T100","T150","TD","TH","TB","TJ")
  }
  #------------------------------------------------------------------------------------
  if(prod(n>pvec)==1 & n<p){
    for(j in 1:repempinull){
      Data<-t(rmnorm(n,numeric(p),SigmaMat0))
      SimvalDar5[j]<-max(Tnm(n,k,eta,c=1,m=5,pvec,Data))
      SimvalDar25[j]<-max(Tnm(n,k,eta,c=1,m=25,pvec,Data))
      SimvalDar50[j]<-max(Tnm(n,k,eta,c=1,m=50,pvec,Data))
      SimvalDar100[j]<-max(Tnm(n,k,eta,c=1,m=100,pvec,Data))
      SimvalDar150[j]<-max(Tnm(n,k,eta,c=1,m=150,pvec,Data))
      SimvalHyd[j]<-THydo(nvec,pvec,Data)
      SimvalBao[j]<-TBao(n,k,pvec,Data)
    }
    #--------------------------------------------
    SizeDar5<-mean(SimvalDar5>Quant5)
    SizeDar25<-mean(SimvalDar25>Quant25)
    SizeDar50<-mean(SimvalDar50>Quant50)
    SizeDar100<-mean(SimvalDar100>Quant100)
    SizeDar150<-mean(SimvalDar150>Quant150)
 
    SizeHyd<-mean(SimvalHyd>z)
    SizeBao<-mean(SimvalBao>z)
    #-----------------------------------------------
    out0<-c(SizeDar5,SizeDar25,SizeDar50,SizeDar100,SizeDar150,NA,SizeHyd,SizeBao,NA)
    names(out0)<-c("T5","T25","T50","T100","T150","TD","TH","TB","TJ")
  }
  #------------------------------------------------------------------------------------
  if(prod(n>pvec)==0){
    for(j in 1:repempinull){
      Data<-t(rmnorm(n,numeric(p),SigmaMat0))
      SimvalDar5[j]<-max(Tnm(n,k,eta,c=1,m=5,pvec,Data))
      SimvalDar25[j]<-max(Tnm(n,k,eta,c=1,m=25,pvec,Data))
      SimvalDar50[j]<-max(Tnm(n,k,eta,c=1,m=50,pvec,Data))
      SimvalDar100[j]<-max(Tnm(n,k,eta,c=1,m=100,pvec,Data))
      SimvalDar150[j]<-max(Tnm(n,k,eta,c=1,m=150,pvec,Data))
      SimvalHyd[j]<-THydo(nvec,pvec,Data)
    }
    #--------------------------------------------
    SizeDar5<-mean(SimvalDar5>Quant5)
    SizeDar25<-mean(SimvalDar25>Quant25)
    SizeDar50<-mean(SimvalDar50>Quant50)
    SizeDar100<-mean(SimvalDar100>Quant100)
    SizeDar150<-mean(SimvalDar150>Quant150)

    SizeHyd<-mean(SimvalHyd>z)
    #-----------------------------------------------
    out0<-c(SizeDar5,SizeDar25,SizeDar50,SizeDar100,SizeDar150,NA,SizeHyd,NA,NA)
    names(out0)<-c("T5","T25","T50","T100","T150","TD","TH","TB","TJ")
  }
  #------------------------------------------------------------------------------------
   y<-c(gamma=gamma,n=n,p=p,round(out0,digits=3))
  return(y)
}
#----------------------------------------------
Alpha_AdjPower<-function(reppower,repempinull,rovec,n,p,k,gamma){
  pvec<-rep(p/k,k)
  SigmaMat0<-Covmat(gamma=0,Sigmaxx=S(p/k,rho=rovec[1]),Sigmaxy=S(p/k,rho=rovec[2]),pvec)
  SigmaMat<-Covmat(gamma,Sigmaxx=S(p/k,rho=rovec[1]),Sigmaxy=S(p/k,rho=rovec[2]),pvec)
  p<-sum(pvec)
  k<-length(pvec)
  c<-1
  AsPar<-Par(n,rep(c,k))
  eta<-AsPar[1]
  f<-AsPar[2]
  nu<-abs(AsPar[3])
  
  nvec<-cumsum(rep(n/k,k))
  SimvalJia<-SimvalHyd<-SimvalDet<-SimvalBao<-numeric(max(reppower,repempinull))
  SimvalDar5<-SimvalDar25<-SimvalDar50<-SimvalDar100<-SimvalDar150<-numeric(max(reppower,repempinull))
  EmpQuant5=EmpQuant25=EmpQuant50=EmpQuant100=EmpQuant150=EmpQuantHyd=EmpQuantDet=EmpQuantBao=EmpQuantJia=NA
  #--------------------------------------------------------------------------------
  if(prod(n>pvec)==1 & n>p){
    for(j in 1:repempinull){
      Data<-t(rmnorm(n,numeric(p),SigmaMat0))
      SimvalDar5[j]<-max(Tnm(n,k,eta,c=1,m=5,pvec,Data))
      SimvalDar25[j]<-max(Tnm(n,k,eta,c=1,m=25,pvec,Data))
      SimvalDar50[j]<-max(Tnm(n,k,eta,c=1,m=50,pvec,Data))
      SimvalDar100[j]<-max(Tnm(n,k,eta,c=1,m=100,pvec,Data))
      SimvalDar150[j]<-max(Tnm(n,k,eta,c=1,m=150,pvec,Data))
      SimvalHyd[j]<-THydo(nvec,pvec,Data)
      SimvalDet[j]<-Tdette(n,pvec,Data)
      SimvalBao[j]<-TBao(n,k,pvec,Data)
      SimvalJia[j]<-Jiang(n,k,pvec,Data)
    }
    EmpQuant5<-quantile(SimvalDar5,probs = 0.95)
    EmpQuant25<-quantile(SimvalDar25,probs = 0.95)
    EmpQuant50<-quantile(SimvalDar50,probs = 0.95)
    EmpQuant100<-quantile(SimvalDar100,probs = 0.95)
    EmpQuant150<-quantile(SimvalDar150,probs = 0.95)
    EmpQuantHyd<-quantile(SimvalHyd,probs = 0.95)
    EmpQuantDet<-quantile(SimvalDet,probs = 0.95)
    EmpQuantBao<-quantile(SimvalBao,probs = 0.95)
    EmpQuantJia<-quantile(SimvalJia,probs = 0.95)
    #-----------------------------------------------
    for(j in 1:reppower){
      Data<-t(rmnorm(n,numeric(p),SigmaMat))
      TntVals5<-max(Tnm(n,k,eta,c=1,m=5,pvec,Data))
      TntVals25<-max(Tnm(n,k,eta,c=1,m=25,pvec,Data))
      TntVals50<-max(Tnm(n,k,eta,c=1,m=50,pvec,Data))
      TntVals100<-max(Tnm(n,k,eta,c=1,m=100,pvec,Data))
      TntVals150<-max(Tnm(n,k,eta,c=1,m=150,pvec,Data))
      SimvalDar5[j]<-as.numeric(TntVals5>=EmpQuant5)
      SimvalDar25[j]<-as.numeric(TntVals25>=EmpQuant25)
      SimvalDar50[j]<-as.numeric(TntVals50>=EmpQuant50)
      SimvalDar100[j]<-as.numeric(TntVals100>=EmpQuant100)
      SimvalDar150[j]<-as.numeric(TntVals150>=EmpQuant150)
      h<-THydo(nvec,pvec,Data)
      SimvalHyd[j]<-as.numeric(h>EmpQuantHyd)
      d<-Tdette(n,pvec,Data)
      SimvalDet[j]<-as.numeric(d>EmpQuantDet)
      Bao<-TBao(n,k,pvec,Data)
      SimvalBao[j]<-as.numeric(Bao>EmpQuantBao)
      J<-Jiang(n,k,pvec,Data)
      SimvalJia[j]<-as.numeric(J>EmpQuantJia)
    }
    out<-c(mean(SimvalDar5),mean(SimvalDar25),mean(SimvalDar50),mean(SimvalDar100),mean(SimvalDar150),mean(SimvalDet),mean(SimvalHyd),mean(SimvalBao),mean(SimvalJia))
    names(out)<-c("T5","T25","T50","T100","T150","TD","TH","TB","TJ")
  }
  #------------------------------------------------------------------------------------
  if(prod(n>pvec)==1 & n<p){
    for(j in 1:repempinull){
      Data<-t(rmnorm(n,numeric(p),SigmaMat0))
      SimvalDar5[j]<-max(Tnm(n,k,eta,c=1,m=5,pvec,Data))
      SimvalDar25[j]<-max(Tnm(n,k,eta,c=1,m=25,pvec,Data))
      SimvalDar50[j]<-max(Tnm(n,k,eta,c=1,m=50,pvec,Data))
      SimvalDar100[j]<-max(Tnm(n,k,eta,c=1,m=100,pvec,Data))
      SimvalDar150[j]<-max(Tnm(n,k,eta,c=1,m=150,pvec,Data))
      SimvalHyd[j]<-THydo(nvec,pvec,Data)
      SimvalBao[j]<-TBao(n,k,pvec,Data)
    }
    EmpQuant5<-quantile(SimvalDar5,probs = 0.95)
    EmpQuant25<-quantile(SimvalDar25,probs = 0.95)
    EmpQuant50<-quantile(SimvalDar50,probs = 0.95)
    EmpQuant100<-quantile(SimvalDar100,probs = 0.95)
    EmpQuant150<-quantile(SimvalDar150,probs = 0.95)
    EmpQuantHyd<-quantile(SimvalHyd,probs = 0.95)
    EmpQuantBao<-quantile(SimvalBao,probs = 0.95)
    #--------------------------------------------
    for(j in 1:reppower){
      Data<-t(rmnorm(n,numeric(p),SigmaMat))
      TntVals5<-max(Tnm(n,k,eta,c=1,m=5,pvec,Data))
      TntVals25<-max(Tnm(n,k,eta,c=1,m=25,pvec,Data))
      TntVals50<-max(Tnm(n,k,eta,c=1,m=50,pvec,Data))
      TntVals100<-max(Tnm(n,k,eta,c=1,m=100,pvec,Data))
      TntVals150<-max(Tnm(n,k,eta,c=1,m=150,pvec,Data))
      SimvalDar5[j]<-as.numeric(TntVals5>=EmpQuant5)
      SimvalDar25[j]<-as.numeric(TntVals25>=EmpQuant25)
      SimvalDar50[j]<-as.numeric(TntVals50>=EmpQuant50)
      SimvalDar100[j]<-as.numeric(TntVals100>=EmpQuant100)
      SimvalDar150[j]<-as.numeric(TntVals150>=EmpQuant150)
      h<-THydo(nvec,pvec,Data)
      SimvalHyd[j]<-as.numeric(h>EmpQuantHyd)
      Bao<-TBao(n,k,pvec,Data)
      SimvalBao[j]<-as.numeric(Bao>EmpQuantBao)
    }
    out<-c(mean(SimvalDar5),mean(SimvalDar25),mean(SimvalDar50),mean(SimvalDar100),mean(SimvalDar150),NA,mean(SimvalHyd),mean(SimvalBao),NA)
    names(out)<-c("T5","T25","T50","T100","T150","TD","TH","TB","TJ")
  }
  #------------------------------------------------------------------------------------
  if(prod(n>pvec)==0){
    for(j in 1:repempinull){
      Data<-t(rmnorm(n,numeric(p),SigmaMat0))
      SimvalDar5[j]<-max(Tnm(n,k,eta,c=1,m=5,pvec,Data))
      SimvalDar25[j]<-max(Tnm(n,k,eta,c=1,m=25,pvec,Data))
      SimvalDar50[j]<-max(Tnm(n,k,eta,c=1,m=50,pvec,Data))
      SimvalDar100[j]<-max(Tnm(n,k,eta,c=1,m=100,pvec,Data))
      SimvalDar150[j]<-max(Tnm(n,k,eta,c=1,m=150,pvec,Data))
      SimvalHyd[j]<-THydo(nvec,pvec,Data)
    }
    EmpQuant5<-quantile(SimvalDar5,probs = 0.95)
    EmpQuant25<-quantile(SimvalDar25,probs = 0.95)
    EmpQuant50<-quantile(SimvalDar50,probs = 0.95)
    EmpQuant100<-quantile(SimvalDar100,probs = 0.95)
    EmpQuant150<-quantile(SimvalDar150,probs = 0.95)
    EmpQuantHyd<-quantile(SimvalHyd,probs = 0.95)
    #--------------------------------------------
    for(j in 1:reppower){
      Data<-t(rmnorm(n,numeric(p),SigmaMat))
      TntVals5<-max(Tnm(n,k,eta,c=1,m=5,pvec,Data))
      TntVals25<-max(Tnm(n,k,eta,c=1,m=25,pvec,Data))
      TntVals50<-max(Tnm(n,k,eta,c=1,m=50,pvec,Data))
      TntVals100<-max(Tnm(n,k,eta,c=1,m=100,pvec,Data))
      TntVals150<-max(Tnm(n,k,eta,c=1,m=150,pvec,Data))
      SimvalDar5[j]<-as.numeric(TntVals5>=EmpQuant5)
      SimvalDar25[j]<-as.numeric(TntVals25>=EmpQuant25)
      SimvalDar50[j]<-as.numeric(TntVals50>=EmpQuant50)
      SimvalDar100[j]<-as.numeric(TntVals100>=EmpQuant100)
      SimvalDar150[j]<-as.numeric(TntVals150>=EmpQuant150)
      h<-THydo(nvec,pvec,Data)
      SimvalHyd[j]<-as.numeric(h>EmpQuantHyd)
    }
    out<-c(mean(SimvalDar5),mean(SimvalDar25),mean(SimvalDar50),mean(SimvalDar100),mean(SimvalDar150),NA,mean(SimvalHyd),NA,NA)
    names(out)<-c("T5","T25","T50","T100","T150","TD","TH","TB","TJ")
  }
  #------------------------------------------------------------------------------------
  y<-c(gamma=gamma,n=n,p=p,round(out,digits=3))
  return(y)
}
#================table 1:k=2 ============================
#----------------gamma=0-----------------------
Tab011<-EmpSize(reppower=5000,repempinull=5000,rovec=c(0.95,0.95),n=60,p=30,k=2,gamma=0)
Tab012<-EmpSize(reppower=5000,repempinull=5000,rovec=c(0.95,0.95),n=120,p=90,k=2,gamma=0)
Tab013<-EmpSize(reppower=5000,repempinull=5000,rovec=c(0.95,0.95),n=180,p=150,k=2,gamma=0)
#----------------gamma=0.3-----------------------
Tab11<-Alpha_AdjPower(reppower=5000,repempinull=5000,rovec=c(0.95,0.95),n=60,p=30,k=2,gamma=0.3)
Tab12<-Alpha_AdjPower(reppower=5000,repempinull=5000,rovec=c(0.95,0.95),n=120,p=90,k=2,gamma=0.3)
Tab13<-Alpha_AdjPower(reppower=5000,repempinull=5000,rovec=c(0.95,0.95),n=180,p=150,k=2,gamma=0.3)
#----------------gamma=0.8-----------------------
Tab14<-Alpha_AdjPower(reppower=5000,repempinull=5000,rovec=c(0.95,0.95),n=60,p=30,k=2,gamma=0.8)
Tab15<-Alpha_AdjPower(reppower=5000,repempinull=5000,rovec=c(0.95,0.95),n=120,p=90,k=2,gamma=0.8)
Tab16<-Alpha_AdjPower(reppower=5000,repempinull=5000,rovec=c(0.95,0.95),n=180,p=150,k=2,gamma=0.8)
#---------------------------------------------
TABLE1<-rbind(Tab011,Tab012,Tab013,Tab11,Tab12,Tab13,Tab14,Tab15,Tab16)
#================table 2:k=3 ============================
#----------------gamma=0-----------------------
Tab021<-EmpSize(reppower=5000,repempinull=5000,rovec=c(0.95,0.95),n=60,p=30,k=3,gamma=0)
Tab022<-EmpSize(reppower=5000,repempinull=5000,rovec=c(0.95,0.95),n=120,p=90,k=3,gamma=0)
Tab023<-EmpSize(reppower=5000,repempinull=5000,rovec=c(0.95,0.95),n=180,p=150,k=3,gamma=0)
#----------------gamma=0.3-----------------------
Tab21<-Alpha_AdjPower(reppower=5000,repempinull=5000,rovec=c(0.95,0.95),n=60,p=30,k=3,gamma=0.3)
Tab22<-Alpha_AdjPower(reppower=5000,repempinull=5000,rovec=c(0.95,0.95),n=120,p=90,k=3,gamma=0.3)
Tab23<-Alpha_AdjPower(reppower=5000,repempinull=5000,rovec=c(0.95,0.95),n=180,p=150,k=3,gamma=0.3)
#----------------gamma=0.8-----------------------
Tab24<-Alpha_AdjPower(reppower=5000,repempinull=5000,rovec=c(0.95,0.95),n=60,p=30,k=3,gamma=0.8)
Tab25<-Alpha_AdjPower(reppower=5000,repempinull=5000,rovec=c(0.95,0.95),n=120,p=90,k=3,gamma=0.8)
Tab26<-Alpha_AdjPower(reppower=5000,repempinull=5000,rovec=c(0.95,0.95),n=180,p=150,k=3,gamma=0.8)
#-----------------------------------------------
TABLE2<-rbind(Tab021,Tab022,Tab023,Tab21,Tab22,Tab23,Tab24,Tab25,Tab26)
 


