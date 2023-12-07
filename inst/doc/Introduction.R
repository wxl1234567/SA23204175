## ----eval=FALSE---------------------------------------------------------------
#  PMLE.Iter.R<-function(x,y,lambda0,max.it=5000){
#    n<-nrow(x);p<-ncol(x)
#    beta.f<-rep(0,p)
#    beta<-beta.f;sigma<-0.1
#    for(i in 1:max.it){
#      sigma.hat<-sqrt(sum(y*(y-x%*%beta))/n)
#      fit.hat<-lars(x,y,type="lasso",intercept=F,normalize=F,use.Gram=F)
#      beta.hat<-predict.lars(fit.hat,s=n*lambda0*sigma.hat,type="coefficients",mode="lambda")$coef
#      if(abs(sigma-sigma.hat)<1e-6) break
#      beta<-beta.hat;sigma<-sigma.hat
#    }
#    list(beta.pmle=beta.hat,sigma.pmle=sigma.hat)
#  }

## ----eval=FALSE---------------------------------------------------------------
#  BC.APMLE.R<-function(x,y,lambda0,beta.pmle,sigma.pmle){
#    n<-nrow(x);p<-ncol(x)
#    f<-lars(x,y,type="lasso",intercept=F,normalize=F,use.Gram=F)
#    y.hat<-predict.lars(f,x,s=n*lambda0*sigma.pmle,type="fit",mode="lambda")$fit
#    sigma.bc<-sqrt(mean((y-y.hat)^2))
#    beta.bc<-predict.lars(f,s=n*lambda0*sigma.bc,type="coefficients",mode="lambda")$coef
#    list(beta.bc=beta.bc,sigma.bc=sigma.bc)
#  }

## -----------------------------------------------------------------------------
library(SA23204175)
library(scalreg)
library(lars)

## -----------------------------------------------------------------------------
#Initialize parameters
n<-200;p<-2000;r0<-0.5
lambda1<-sqrt(log(p)/n)
lambda2<-sqrt(2*log(p)/n)
lambda3<-sqrt(4*log(p)/n)

#Generate covariance matrix Sigma of Gaussian distribution
#library(pracma) #create basic matrices like eye(n,m=n), ones(n,m=n) and zeros(n,m=n)
Sigma1<-pracma::eye(p);Sigma2<-0.5*pracma::eye(p)+0.5*pracma::ones(p) 

#Generate design matrix X
library(MASS);set.seed(123)
X1<-mvrnorm(n,rep(0,p),Sigma1)
X2<-mvrnorm(n,rep(0,p),Sigma2)

#Initialize parameters
sigma<-1;beta.star<-matrix(0,nrow=p,ncol=1);beta.star[1,1]<-1/sqrt(3)
beta.star[2,1]<-1/sqrt(3);beta.star[3,1]<-1/sqrt(3)

#Generate response vector y
set.seed(123)
y1<-as.vector(X1%*%beta.star)+rnorm(n,0,sigma^2)
y2<-as.vector(X2%*%beta.star)+rnorm(n,0,sigma^2)

## -----------------------------------------------------------------------------
#Method-PMLE Dataset-(X1,y1) LSE-FALSE

PMLE.LSEF.1<-list(l1=PMLE.Iter.R(X1,y1,lambda0=lambda1),l2=PMLE.Iter.R(X1,y1,lambda0=lambda2),l3=PMLE.Iter.R(X1,y1,lambda0=lambda3))

#Method-BC Dataset-(X1,y1) LSE-FALSE 

BC.LSEF.1<-list(l1=BC.APMLE.R(X1,y1,lambda0=lambda1,PMLE.LSEF$l1$beta.pmle,PMLE.LSEF.1$l1$sigma.pmle),
                l2=BC.APMLE.R(X1,y1,lambda0=lambda2,PMLE.LSEF$l2$beta.pmle,PMLE.LSEF.1$l2$sigma.pmle),
                l3=BC.APMLE.R(X1,y1,lambda0=lambda3,PMLE.LSEF$l3$beta.pmle,PMLE.LSEF.1$l3$sigma.pmle))

#Method-SLASSO Dataset-(X1,y1) LSE-FALSE 

SL.LSEF.1<-list(l1=scalreg(X1,y1,lam0=lambda1,LSE=FALSE),
                l2=scalreg(X1,y1,lam0=lambda2,LSE=FALSE),
                l3=scalreg(X1,y1,lam0=lambda3,LSE=FALSE))

#Method-PMLE Dataset-(X2,y2) LSE-FALSE

PMLE.LSEF.2<-list(l1=PMLE.Iter.R(X2,y2,lambda0=lambda1),l2=PMLE.Iter.R(X2,y2,lambda0=lambda2),l3=PMLE.Iter.R(X2,y2,lambda0=lambda3))

#Method-BC Dataset-(X2,y2) LSE-FALSE 

BC.LSEF.2<-list(l1=BC.APMLE.R(X2,y2,lambda0=lambda1,PMLE.LSEF$l1$beta.pmle,PMLE.LSEF.2$l1$sigma.pmle),
                l2=BC.APMLE.R(X2,y2,lambda0=lambda2,PMLE.LSEF$l2$beta.pmle,PMLE.LSEF.2$l2$sigma.pmle),
                l3=BC.APMLE.R(X2,y2,lambda0=lambda3,PMLE.LSEF$l3$beta.pmle,PMLE.LSEF.2$l3$sigma.pmle))

#Method-SLASSO Dataset-(X2,y2) LSE-FALSE 

SL.LSEF.2<-list(l1=scalreg(X2,y2,lam0=lambda1,LSE=FALSE),
                l2=scalreg(X2,y2,lam0=lambda2,LSE=FALSE),
                l3=scalreg(X2,y2,lam0=lambda3,LSE=FALSE))

## -----------------------------------------------------------------------------
result1<-matrix(0,nrow=9,ncol=2)
rownames(result1)<-c("PMLE/Lambda1","PMLE/Lambda2","PMLE/Lambda3",
                     "BC/Lambda1","BC/Lambda2","BC/Lambda3",
                     "SLASSO/Lambda1","SLASSO/Lambda2","SLASSO/Lambda3")
colnames(result1)<-c("sigma.hat/sigma*","AMS")
result1[1,1]<-PMLE.LSEF.1$l1$sigma.pmle
result1[2,1]<-PMLE.LSEF.1$l2$sigma.pmle
result1[3,1]<-PMLE.LSEF.1$l3$sigma.pmle
result1[4,1]<-BC.LSEF.1$l1$sigma.bc
result1[5,1]<-BC.LSEF.1$l2$sigma.bc
result1[6,1]<-BC.LSEF.1$l3$sigma.bc
result1[7,1]<-SL.LSEF.1$l1$hsigma
result1[8,1]<-SL.LSEF.1$l2$hsigma
result1[9,1]<-SL.LSEF.1$l3$hsigma

result1[1,2]<-sum(PMLE.LSEF.1$l1$beta.pmle!=0)
result1[2,2]<-sum(PMLE.LSEF.1$l2$beta.pmle!=0)
result1[3,2]<-sum(PMLE.LSEF.1$l3$beta.pmle!=0)
result1[4,2]<-sum(BC.LSEF.1$l1$beta.bc!=0)
result1[5,2]<-sum(BC.LSEF.1$l2$beta.bc!=0)
result1[6,2]<-sum(BC.LSEF.1$l3$beta.bc!=0)
result1[7,2]<-sum(SL.LSEF.1$l1$coefficients!=0)
result1[8,2]<-sum(SL.LSEF.1$l2$coefficients!=0)
result1[9,2]<-sum(SL.LSEF.1$l3$coefficients!=0)

result2<-matrix(0,nrow=9,ncol=2)
rownames(result2)<-c("PMLE/Lambda1","PMLE/Lambda2","PMLE/Lambda3",
                     "BC/Lambda1","BC/Lambda2","BC/Lambda3",
                     "SLASSO/Lambda1","SLASSO/Lambda2","SLASSO/Lambda3")
colnames(result2)<-c("sigma.hat/sigma*","AMS")
result2[1,1]<-PMLE.LSEF.2$l1$sigma.pmle
result2[2,1]<-PMLE.LSEF.2$l2$sigma.pmle
result2[3,1]<-PMLE.LSEF.2$l3$sigma.pmle
result2[4,1]<-BC.LSEF.2$l1$sigma.bc
result2[5,1]<-BC.LSEF.2$l2$sigma.bc
result2[6,1]<-BC.LSEF.2$l3$sigma.bc
result2[7,1]<-SL.LSEF.2$l1$hsigma
result2[8,1]<-SL.LSEF.2$l2$hsigma
result2[9,1]<-SL.LSEF.2$l3$hsigma

result2[1,2]<-sum(PMLE.LSEF.2$l1$beta.pmle!=0)
result2[2,2]<-sum(PMLE.LSEF.2$l2$beta.pmle!=0)
result2[3,2]<-sum(PMLE.LSEF.2$l3$beta.pmle!=0)
result2[4,2]<-sum(BC.LSEF.2$l1$beta.bc!=0)
result2[5,2]<-sum(BC.LSEF.2$l2$beta.bc!=0)
result2[6,2]<-sum(BC.LSEF.2$l3$beta.bc!=0)
result2[7,2]<-sum(SL.LSEF.2$l1$coefficients!=0)
result2[8,2]<-sum(SL.LSEF.2$l2$coefficients!=0)
result2[9,2]<-sum(SL.LSEF.2$l3$coefficients!=0)

## -----------------------------------------------------------------------------
cat("result1 with (X1,y1)(r0=0)\n");result1
cat("\nresult2 with (X2,y2)(r0=0.5)\n");result2

