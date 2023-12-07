#' @import microbenchmark
#' @import coda
#' @import DAAG
#' @import boot
#' @import lars
#' @import MASS
#' @import bootstrap
#' @import scalreg
#' @importFrom pracma eye ones 
#' @importFrom Rcpp evalCpp
#' @importFrom stats rbeta rbinom
#' @useDynLib SA23204175
NULL

gibbsR<-function(a,b,n,N){
  mat<-matrix(nrow=N,ncol=2)
  mat[1,1]<-3;mat[1,2]<-0.5
  for(i in 2:N){
    y<-mat[i-1,2]
    x<-rbinom(1,n,y)
    y<-rbeta(1,x+a,n-x+b)
    mat[i,]<-c(x,y)
  }
  mat
}

#' @title A illustration dataset with r0=0
#' @name d1
#' @description A dataset with r0=0 used to illustrate the performance of pmle and bc of versions R and Rcpp
#' @examples
#' \dontrun{
#'     data(d1)
#'     attach(d1)
#' }
d1<-NULL

#' @title A illustration dataset with r0=0.5
#' @name d2
#' @description A dataset with r0=0.5 used to illustrate the performance of pmle and bc of versions R and Rcpp
#' @examples
#' \dontrun{
#'     data(d2)
#'     attach(d2)
#' }
d2<-NULL

#' @title Compute The L1 Penalized Maximum Likelihood Estimator
#' @description More details refer to scaled sparse linear regression
#' @param x the design matrix
#' @param y the response vector
#' @param lambda0 the pre-fixed penalty level
#' @param max.it the maximum number of iterations
#' @return beta.pmle and sigma.pmle
#' @examples 
#' \dontrun{
#'     PMLE.Iter.R(d1[,-c(1,2)],d1[,1001],lambda0=d1[1,1002])
#'     PMLE.Iter.R(d2[,-c(1,2)],d2[,1001],lambda0=d2[1,1002])
#' }
#' @export
PMLE.Iter.R<-function(x,y,lambda0,max.it=5000){
  n<-nrow(x);p<-ncol(x)
  beta.f<-rep(0,p)
  beta<-beta.f;sigma<-0.1
  for(i in 1:max.it){
    sigma.hat<-sqrt(sum(y*(y-x%*%beta))/n)
    fit.hat<-lars(x,y,type="lasso",intercept=F,normalize=F,use.Gram=F)
    beta.hat<-predict.lars(fit.hat,s=n*lambda0*sigma.hat,type="coefficients",mode="lambda")$coef
    if(abs(sigma-sigma.hat)<1e-6) break
    beta<-beta.hat;sigma<-sigma.hat
  }  
  list(beta.pmle=beta.hat,sigma.pmle=sigma.hat)
}

#' @title Compute The Bias Correction After The L1 Penalized Maximum Likelihood Estimator
#' @description More details refer to scaled sparse linear regression
#' @param x the design matrix
#' @param y the response vector
#' @param lambda0 the pre-fixed penalty level
#' @param beta.pmle penalized maximum likelihood estimator of beta
#' @param sigma.pmle penalized maximum likelihood estimator of sigma
#' @return beta.bc and sigma.bc
#' @examples
#' \dontrun{
#'     r1.pmle<-PMLE.Iter.R(d1[,-c(1,2)],d1[,1001],lambda0=d1[1,1002])
#'     BC.APMLE.R(
#'       d1[,-c(1,2)],
#'       d1[,1001],
#'       lambda0=d1[1,1002],
#'       r1.pmle$beta.pmle,
#'       r1.pmle$sigma.pmle
#'     )
#' }
#' @export
BC.APMLE.R<-function(x,y,lambda0,beta.pmle,sigma.pmle){
  n<-nrow(x);p<-ncol(x)
  f<-lars(x,y,type="lasso",intercept=F,normalize=F,use.Gram=F)
  y.hat<-predict.lars(f,x,s=n*lambda0*sigma.pmle,type="fit",mode="lambda")$fit
  sigma.bc<-sqrt(mean((y-y.hat)^2))
  beta.bc<-predict.lars(f,s=n*lambda0*sigma.bc,type="coefficients",mode="lambda")$coef
  list(beta.bc=beta.bc,sigma.bc=sigma.bc) 
}