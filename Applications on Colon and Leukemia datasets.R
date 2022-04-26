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
  M<-vector("list",length=k)
  for(i in 1:k){
    M[[i]]<-rnorm(n=c*pvec[i])
  }
  rmat<-.bdiag(M)
  return(rmat)#======px(kc)=pxr matrix
}
#-----------------------------
Tnt<-function(n,k,eta,c,pvec,Data){
    R<-Rmat(c,pvec)#dim(R)#dim(Data)
    Y<-(as.matrix(Data%*%R))#nxr dim(Y)
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
#-----------------------------------------------
ProjectionTest<-function(Data,pvec,alpha,m){
  p<-dim(Data)[1]
  n<-dim(Data)[2]
  k<-length(pvec)
    c<-1
    nvec<-cumsum(rep(n/k,k))
    AsPar<-Par(n,rep(c,k))
    eta<-AsPar[1]
    f<-AsPar[2]
    nu<-abs(AsPar[3])
    Tmndist<-function(y,m){nu*pchisq(y,f+4)+(1-nu)*pchisq(y,f)-(1-alpha)^(1/m)}
    Quant<-uniroot(Tmndist, c(0,10*f*(1+nu)),m,extendInt="yes")$root
    TntVals<-Tnm(n,k,eta,c,m,pvec,Data)
    return(list(TestStatistic=max(TntVals),f1=f,nu1=nu,CriticalValue=Quant))
}
#-----------------------------------------------
library(Matrix)
library(ShrinkCovMat)#data colon
data(colon)
head(colon)

Data<-t(colon)#data should be a n x p matrix, here 2000 x 62
dim(Data)
k=5
ProjectionTest(Data,pvec=rep(2000/k,k),alpha=0.05,m=1000)
k=10
ProjectionTest(Data,pvec=rep(2000/k,k),alpha=0.05,m=1000)
k=50
ProjectionTest(Data,pvec=rep(2000/k,k),alpha=0.05,m=1000)
k=100
ProjectionTest(Data,pvec=rep(2000/k,k),alpha=0.05,m=1000)




if (!requireNamespace("BiocManager", quietly = TRUE))install.packages("BiocManager")
BiocManager::install("Biobase")
library(propOverlap)
data(leukaemia)
dim(leukaemia)n
head(leukaemia)

Data<-as.matrix(t(leukaemia[1:7000,]))#data should be a n x p matrix 
dim(Data)
k=5
ProjectionTest(Data,pvec=rep(7000/k,k),alpha=0.05,m=1000)
k=10
ProjectionTest(Data,pvec=rep(7000/k,k),alpha=0.05,m=1000)
k=35
ProjectionTest(Data,pvec=rep(7000/k,k),alpha=0.05,m=1000)
k=100
ProjectionTest(Data,pvec=rep(7000/k,k),alpha=0.05,m=1000)





citation(package = "ShrinkCovMat")





