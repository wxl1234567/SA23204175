## -----------------------------------------------------------------------------
#If n=2 and \bar{y}=0.5, then the corresponding likelihood could be draw below.
Likelihood_n2_ybar0.5<-function(mu) mu*exp(-2*mu) ## likelihood function
curve(Likelihood_n2_ybar0.5,0,10,xname="mu") ## curve of likelihood with respect to mu

## -----------------------------------------------------------------------------
plot(function(x) dnorm(x),-4,4,main="Normal Density")

## -----------------------------------------------------------------------------
mysample<-function(x,size=NULL,prob=NULL){
  if(length(x)==1){
    x<-seq(1,x) #first transform integer x into 1:x
  }
  n<-length(x) #then derive the length of x
  if(is.null(prob)){
    prob<-rep(1,n)/n #second, for "prob=NULL", put rep(1,n) on prob
  }
  if(is.null(size)){
    size<-length(x) #third, for "size=NULL", directly set size=length(x)
  }
  cum_prob<-cumsum(prob)
  cum_prob<-c(0,cum_prob) #forth, obtain the cumulative sums of prob
  
##key point: generate the uniform distributed variables and sample x depending on corresponding cumulative sums
  u<-runif(size);y<-rep(0,size) 
  for(i in 1:size){
    j<-1
    while(u[i]>=cum_prob[j]){j<-j+1}
    y[i]=x[j-1]
  }
  return(y)
}
#we could also use "findInterval" to find intervals  

## -----------------------------------------------------------------------------
table(mysample(0:1,size=1000,prob=rep(0.5,2)))
table(mysample(1:3,size=10000,prob=c(0.2,0.3,0.5))) #list two examples without/with weights in large sample cases

table(mysample(c("b","c","a"),size=100,prob=c(0.2,0.3,0.5))) #check for the situation that x is a vector of letters

table(mysample(4,size=100,prob=rep(0.25,4))) #check the situation that x is an integer

table(mysample(1:10,prob=rep(0.1,10))) #check the funciton of "size=NULL"

table(mysample(0:1,size=100)) #check the function of "prob=NULL"

## -----------------------------------------------------------------------------
n<-1000
u<-runif(n)
x<-ifelse(u<=0.5,log(2*u),-log(2*(1-u)))
  #generate a random sample of size 1000 from the standard Laplace distribution

hist(x,prob=TRUE,xlim=c(-2,2),breaks=100,main=expression(f(x)==frac(1,2)*e^(-abs(x))))
y<-seq(-2,2,.01)
lines(y,0.5*exp(-abs(y))) #compare the generated sample with the target distribution using histogram

## ----echo=FALSE---------------------------------------------------------------
#This solution only considers the condition that a and b are positive integers. Otherwise, more advanced functions for integration computation or maximum searching are required to solution.

#The pdf of Beta(a,b) is $f(x)=\frac1{B(a,b)}x^{a-1}(1-x)^{b-1},0\le x\le1$, where $B(a,b)=\frac{(a-1)!(b-1)!}{(a+b-1)!}$ with positive integers a&b.

#The logarithm of this pdf is $\log f(x)\propto(a-1)\log x+(b-1)\log(1-x)$.

#Differentiate this logarithm: $[\log f(x)]'=\frac{a-1}{x}-\frac{b-1}{1-x}$.

#Set this formula to zero and derive  $x=\frac{a-1}{a+b-2}$.

#Thus, the maximum of the pdf is $f(\frac{a-1}{a+b-2})=\frac{(a-1)^{a-1}(b-1)^{b-1}}{(a-1)!(b-1)!}\frac{(a+b-1)!}{(a+b-2)^{a+b-2}}$.

## -----------------------------------------------------------------------------
mybeta<-function(n,a,b){ #see the beta(2,2) algorithm for reference
  k<-0
  while(k<n){
    u<-runif(1)
    x<-runif(1)
    if(x^(a-1)*(1-x)^(b-1)>u){ #accept x
      k<-k+1 
      y[k]<-x
    }
  }
  return(y)
}
  #write a function to generate samples from beta(a,b) by acceptance-rejection method

n_sample<-1000;a_sample<-3;b_sample<-2
beta_sample<-mybeta(n_sample,a_sample,b_sample)
  #give the parameters (n,a,b)=(1000,3,2) and then generate 1000 samples

hist(beta_sample,prob=TRUE,xlim=c(0,1),breaks=50,main=expression(f(x)==frac(1,B(a,b))*x^(a-1)*(1-x)^(b-1)))
y<-seq(0,1,.01)
lines(y,y^(a_sample-1)*(1-y)^(b_sample-1)/beta(a_sample,b_sample)) #compare the generated sample with Beta(3,2) distribution using histogram

## -----------------------------------------------------------------------------
myepa<-function(n){
  u_seq<-runif(3*n,min=-1,max=1)
  u_mat<-matrix(u_seq,ncol=3)
  epa<-ifelse(abs(u_mat[,3])>=abs(u_mat[,2])&abs(u_mat[,3])>=abs(u_mat[,1]),u_mat[,2],u_mat[,3])
  return(epa)
} 
  #write a function to generate variables with the algorithm in Exercise 3.9

n_sample<-2000;epa_sample<-myepa(n_sample) 
  #give a concrete set up and generate the corresponding variables 

hist(epa_sample,prob=TRUE,xlim=c(-1,1),breaks=50,main=expression(f[e](x)==frac(3,4)*(1-x^2)))
y<-seq(-1,1,.01)
lines(y,0.75*(1-y^2)) #compare the generated sample with the target distribution using histogram

## -----------------------------------------------------------------------------
set.seed(123)
m<-1e6;K<-100 #number of repeated simulations
n_rho0.5<-rbinom(K,m,2*0.5/pi)
n_rho0.8<-rbinom(K,m,2*0.8/pi)
n_rho1<-rbinom(K,m,2*1/pi)
  #generate binomial variables with parameters (n=100,size=1e6,prob=2*[0.5/0.8/1]/pi), where prob=p=2*rho/pi
pi_rho0.5<-2*0.5/(n_rho0.5/m)
pi_rho0.8<-2*0.8/(n_rho0.8/m)
pi_rho1<-2*1/(n_rho1/m)
  #generate corresponding estimates of pi with rho=0.5/0.8/1
cat("","variance of hat_pi with rho=0.5:",var(pi_rho0.5)/m,"\n","variance of hat_pi with rho=0.8:",var(pi_rho0.8)/m,"\n","variance of hat_pi with rho=1:",var(pi_rho1)/m,"\n")

## -----------------------------------------------------------------------------
-exp(2)+3*exp(1)-1 #give the approximate value of the covariance

## -----------------------------------------------------------------------------
-3*exp(2)+10*exp(1)-5 #give the approximate value of the variance

## -----------------------------------------------------------------------------
(2*exp(2)-6*exp(1)+2)/(-exp(2)+4*exp(1)-3) #give the approximate value of the percent reduction and find the antithetic variate approach works well

## -----------------------------------------------------------------------------
set.seed(111)
m<-1000000
u<-runif(m) #generate uniform variables
test_simpleMC<-exp(u) #simple MC
test_antithetic_variate<-(exp(u[1:m/2])+exp(1-u[1:m/2]))/2 #antithetic variate
cat("","theta value:",exp(1)-1,"\n","theta estimate of simple MC:",mean(test_simpleMC),"\n","theta estimate of antithetic variate:",mean(test_antithetic_variate),"\n")

## -----------------------------------------------------------------------------
cat("","theoretical value of percent reduction:",(2*exp(2)-6*exp(1)+2)/(-exp(2)+4*exp(1)-3),"\n","empirical estimate of percent reduction:",1-var(test_antithetic_variate)/var(test_simpleMC),"\n")

## -----------------------------------------------------------------------------
FWER.B<-FWER.BH<-numeric(1000)
FDR.B<-FDR.BH<-numeric(1000)
TPR.B<-TPR.BH<-numeric(1000)
set.seed(123)
for(m in 1:1000){
  p<-c(runif(950),rbeta(50,0.1,1)) #generate pvalues
  p.b<-p.adjust(p,method="bonferroni") #Bonferroni adjustment
  p.bh<-p.adjust(p,method="BH") #BH adjustment
  FWER.B[m]<-ifelse(any(p.b[1:950]<0.1),1,0) 
  FWER.BH[m]<-ifelse(any(p.bh[1:950]<0.1),1,0)
    #check number of [true H_0 and reject H_0] is not less than 1
  FDR.B[m]<-sum(ifelse((p.b[1:950]<0.1),1,0))/sum(ifelse(p.b,1,0))
  FDR.BH[m]<-sum(ifelse((p.bh[1:950]<0.1),1,0))/sum(ifelse((p.bh<0.1),1,0))
    #number of [true H_o and reject H_0]/number of [reject H_0]
  TPR.B[m]<-sum(ifelse((p.b[951:1000]<0.1),1,0))/50
  TPR.BH[m]<-sum(ifelse((p.bh[951:1000]<0.1),1,0))/50
    #number of [true H_a and reject H_0]/number of [true H_a]
}
m<-matrix(c(mean(FWER.B),mean(FDR.B),mean(TPR.B),mean(FWER.BH),mean(FDR.BH),mean(TPR.BH)),nrow=2,byrow=TRUE)
rownames(m)<-c("Bonferroni","BH")
colnames(m)<-c("FWER","FDR","TPR")
m #print results with matrix

## -----------------------------------------------------------------------------
#sample size:n=5
set.seed(1)
la<-2;n<-5;x<-rexp(n,la)
lambdastar<-numeric(1000)
lambdabias<-numeric(1000)
lambdasd<-numeric(1000)
lambda<-1/mean(x)
for(m in 1:1000){
  for(b in 1:1000){
    xstar<-sample(x,replace=TRUE)
    lambdastar[b]<-1/mean(xstar)
  }  
  lambdabias[m]<-mean(lambdastar)-lambda
  lambdasd[m]<-sd(lambdastar)
}
round(c(bias.boot=mean(lambdabias),se.boot=mean(lambdasd)),4)
round(c(bias.thm=la/(n-1),se.thm=la*n/((n-1)*sqrt(n-2))),4)

## -----------------------------------------------------------------------------
#sample size:n=10
set.seed(1)
la<-2;n<-10;x<-rexp(n,la)
lambdastar<-numeric(1000)
lambdabias<-numeric(1000)
lambdasd<-numeric(1000)
lambda<-1/mean(x)
for(m in 1:1000){
  for(b in 1:1000){
    xstar<-sample(x,replace=TRUE)
    lambdastar[b]<-1/mean(xstar)
  }  
  lambdabias[m]<-mean(lambdastar)-lambda
  lambdasd[m]<-sd(lambdastar)
}
round(c(bias.boot=mean(lambdabias),se.boot=mean(lambdasd)),4)
round(c(bias.thm=la/(n-1),se.thm=la*n/((n-1)*sqrt(n-2))),4)

## -----------------------------------------------------------------------------
#sample size:n=20
set.seed(1)
la<-2;n<-20;x<-rexp(n,la)
lambdastar<-numeric(1000)
lambdabias<-numeric(1000)
lambdasd<-numeric(1000)
lambda<-1/mean(x)
for(m in 1:1000){
  for(b in 1:1000){
    xstar<-sample(x,replace=TRUE)
    lambdastar[b]<-1/mean(xstar)
  }  
  lambdabias[m]<-mean(lambdastar)-lambda
  lambdasd[m]<-sd(lambdastar)
}
round(c(bias.boot=mean(lambdabias),se.boot=mean(lambdasd)),4)
round(c(bias.thm=la/(n-1),se.thm=la*n/((n-1)*sqrt(n-2))),4)

## -----------------------------------------------------------------------------
#select alpha as 0.05
library(bootstrap)
set.seed(1)
theta.b<-numeric(1000)
theta.sd.b<-numeric(1000)
theta.bm<-numeric(100)
theta<-cor(law$LSAT,law$GPA)
for(b in 1:1000){
  index<-sample(1:15,size=15,replace=TRUE)
  LSAT<-law$LSAT[index]
  GPA<-law$GPA[index]
  theta.b[b]<-cor(LSAT,GPA)
  for(m in 1:100){
    indexx<-sample(index,size=15,replace=TRUE)
    LSAT<-law$LSAT[indexx]
    GPA<-law$GPA[indexx]
    theta.bm[m]<-cor(LSAT,GPA)
  }
  theta.sd.b[b]<-sd(theta.bm)
}
tstar1<-quantile((theta.b-theta)/theta.sd.b,0.025)
tstar2<-quantile((theta.b-theta)/theta.sd.b,0.975)
theta.sd<-sd(theta.b)
cat("bootstrap t confidence interval estimate:","\n","[",theta-tstar2*theta.sd,",",theta-tstar1*theta.sd,"]","\n")

## -----------------------------------------------------------------------------
#standard normal, basic, percentile and BCa
#not use function boot.ci
#select alpha as 0.05
library(boot)
set.seed(1)
mt<-mean(aircondit[,1])
mt.b<-numeric(1000)
for(b in 1:1000){
  index<-sample(12,replace=TRUE)
  mt.b[b]<-mean(aircondit[index,1])
}
mt.sd<-sd(mt.b)

z0.hat<-qnorm(mean(mt.b<mt))
a.hat<-sum((mt-mt.b)^3)/(6*sum(((mt-mt.b)^2)^(3/2))) 
  #still questioning
alpha1<-pnorm(z0.hat+(z0.hat+qnorm(0.025))/(1-a.hat*(z0.hat+qnorm(0.025))))
alpha2<-pnorm(z0.hat+(z0.hat+qnorm(0.975))/(1-a.hat*(z0.hat+qnorm(0.975))))

cat("standard normal bootstrap confidence interval estimate:","\n","[",mt-qnorm(0.975)*mt.sd,",",mt-qnorm(0.025)*mt.sd,"]","\n")
cat("basic bootstrap confidence interval estimate:","\n","[",2*mt-quantile(mt.b,0.975),",",2*mt-quantile(mt.b,0.025),"]","\n")
cat("percentile bootstrap confidence interval estimate:","\n","[",quantile(mt.b,0.025),",",quantile(mt.b,0.975),"]","\n")
cat("BCa bootstrap confidence interval estimate:","\n","[",quantile(mt.b,alpha1),",",quantile(mt.b,alpha2),"]","\n")

## -----------------------------------------------------------------------------
#standard normal, basic, percentile and BCa
#use function boot.ci
#select alpha as 0.05
set.seed(1)
boot.mean<-function(x,i) mean(x[i])
mt.boot<-boot(data=aircondit[,1],statistic=boot.mean,R=1000)
ci<-boot.ci(mt.boot,type=c("norm","basic","perc","bca"))
ci.norm<-ci$norm[2:3];ci.basic<-ci$basic[4:5]
ci.perc<-ci$percent[4:5];ci.bca<-ci$bca[4:5]
cat("standard normal bootstrap confidence interval estimate:","\n","[",ci.norm[1],",",ci.norm[2],"]","\n")
cat("basic bootstrap confidence interval estimate:","\n","[",ci.basic[1],",",ci.basic[2],"]","\n")
cat("percentile bootstrap confidence interval estimate:","\n","[",ci.perc[1],",",ci.perc[2],"]","\n")
cat("BCa bootstrap confidence interval estimate:","\n","[",ci.bca[1],",",ci.bca[2],"]","\n")

## -----------------------------------------------------------------------------
library(bootstrap)
n<-length(scor[,1])
eigenvalues<-eigen(var(scor))$values
theta.hat<-eigenvalues[1]/sum(eigenvalues)
theta.jack<-numeric(n)
for(i in 1:n){
  scor.jack<-scor[-i,]
  eigenvalues.jack<-eigen(var(scor.jack))$values
  theta.jack[i]<-eigenvalues.jack[1]/sum(eigenvalues.jack)
}
bias.jack<-(n-1)*(mean(theta.jack)-theta.hat)
se.jack<-sqrt((n-1)*mean((theta.jack-theta.hat)^2))
cat("estimate of theta:",theta.hat,"\n")
cat("jackknife estimate of bias",bias.jack,"\n")
cat("jackknife estimate of standard error",se.jack,"\n")

## -----------------------------------------------------------------------------
#refer to leave-one-out cross validation in Example 7.18
#use leave-two-out cross validation here
library(DAAG)
set.seed(1)
n<-length(ironslag[,1])
chemical<-ironslag[,1]
magnetic<-ironslag[,2]
e1<-e2<-e3<-e4<-numeric((n-1)*n/2);k<-0
for(i in 1:(n-1)){
  for(j in (i+1):n){
    k<-k+1
    x<-chemical[-c(i,j)]
    y<-magnetic[-c(i,j)]
  
    J1<-lm(y~x)
    yhat1<-J1$coef[1]+J1$coef[2]*chemical[c(i,j)]
    e1[k]<-sum((magnetic[c(i,j)]-yhat1)^2)
  
    J2<-lm(y~x+I(x^2))
    yhat2<-J2$coef[1]+J2$coef[2]*chemical[c(i,j)]+J2$coef[3]*chemical[c(i,j)]^2
    e2[k]<-sum((magnetic[c(i,j)]-yhat2)^2)
  
    J3<-lm(log(y)~x)
    logyhat3<-J3$coef[1]+J3$coef[2]*chemical[c(i,j)]
    yhat3<-exp(logyhat3)
    e3[k]<-sum((magnetic[c(i,j)]-yhat3)^2)

    J4<-lm(log(y)~log(x))
    logyhat4<-J4$coef[1]+J4$coef[2]*log(chemical[c(i,j)])
    yhat4<-exp(logyhat4)
    e4[k]<-sum((magnetic[c(i,j)]-yhat4)^2)  
  }
}
cat("prediction error estimates with leave-two-out:","\n",c(mean(e1),mean(e2),mean(e3),mean(e4)),"\n")

## -----------------------------------------------------------------------------
attach(chickwts)
x<-sort(as.vector(weight[feed=="soybean"]))
y<-sort(as.vector(chickwts$weight[feed=="linseed"]))
detach(chickwts)

## -----------------------------------------------------------------------------
#Apply Cramer-von Mises statistic to data in Examples 8.1 and 8.2.
#Refer to codes of Example 8.1.
set.seed(1)
R<-999 #number of replicates
z<-c(x,y) #pooled sample
K<-1:26;reps<-numeric(R) # storage for replicates
CM<-function(x,y){ #Cramer-von Mises
  n<-length(x);m<-length(y);W2<-0
  for(i in 1:n)
    W2<-W2+((sum(x<=x[i])+1)/(n+1)-(sum(y<=x[i])+1)/(m+1))^2
  for(j in 1:m)
    W2<-W2+((sum(x<=y[j])+1)/(n+1)-(sum(y<=y[j])+1)/(m+1))^2
  W2*m*n/(m+n)^2
}
CM0<-CM(x,y)

for(i in 1:R){
  k<-sample(K,size=14,replace=FALSE) #generate index for first sample
  xx<-z[k];yy<-z[-k] #complement of xx
  reps[i]<-CM(xx,yy)
}
p<-mean(c(CM0,reps)>=CM0);cat("Derive p-value with Cramer-von Mises statistic:",p,"\n")

## -----------------------------------------------------------------------------
#Refer to codes of Example 8.1.
hist(reps,main="",freq=FALSE,xlab="CM(p=0.393)",breaks="scott")
points(CM0,0,cex=1,pch=16) #observed CM

## -----------------------------------------------------------------------------
Count5<-function(x,y){ #Count 5 Test
  X<-x-mean(x);Y<-y-mean(y)
  outx<-sum(X>max(Y))+sum(X<min(Y))
  outy<-sum(Y>max(X))+sum(Y<min(X))
  as.integer(max(c(outx,outy))>5)
}

Count5.PermTest<-function(x,y,rep=999){ 
  #Permutation Test With Count 5
  Count5.0<-Count5(x,y)
  z<-c(x,y);K<-1:(length(z));C<-numeric(rep)
  for(i in 1:rep){
    k<-sample(length(z),size=length(x),replace=FALSE)
    xx<-z[k];yy<-z[-k]
    C[i]<-Count5(xx,yy)
  }
  mean(c(Count5.0,C))
}

#Simulation
#Empirical Type I error: reject true H0
#Same variance of normal distribution
#Set mu1=1,mu2=2 and (n1,n2)={(50,50),(50,60),(50,70)}
set.seed(1)
n1<-rep(50,3);n2<-c(50,60,70)
mu1<-1;mu2<-2
sigma1<-1;sigma2<-1
m<-100
tests<-replicate(m,expr={
  C5P<-numeric(3)
  for(i in 1:3){
    x<-numeric(n1[i])
    y<-numeric(n2[i])
    x<-rnorm(n1[i],mu1,sigma1)
    y<-rnorm(n2[i],mu2,sigma2)
    C5P[i]<-Count5.PermTest(x,y)
  }
  C5P
})
type1error.hat<-rowMeans(tests)
cat("Empirical Type I error of permutation test of Count 5","\n","with (n1,n2)={(50,50),(50,60),(50,70)}:",type1error.hat,"\n")

#Empirical Power: not reject false H0
#Differing variance of normal distribution
#Also set mu1=1,mu2=2 and (n1,n2)={(50,50),(50,60),(50,70)}
set.seed(1)
n1<-rep(50,3);n2<-c(50,60,70)
mu1<-1;mu2<-2
sigma1<-1;sigma2<-2
m<-100
tests<-replicate(m,expr={
  C5P<-numeric(3)
  for(i in 1:3){
    x<-numeric(n1[i])
    y<-numeric(n2[i])
    x<-rnorm(n1[i],mu1,sigma1)
    y<-rnorm(n2[i],mu2,sigma2)
    C5P[i]<-Count5.PermTest(x,y)
  }
  C5P
})
power.hat<-1-rowMeans(tests)
cat("Empirical Power of permutation test of Count 5:","\n","with (n1,n2)={(50,50),(50,60),(50,70)}",power.hat,"\n")

## -----------------------------------------------------------------------------
#a.Design a function with inputs N,b1,b2,b3,f0 to produce a 
a.generator<-function(N,b1,b2,b3,f0){
  x1<-rpois(N,lambda=1)
  x2<-rexp(N,rate=1)
  x3<-rbinom(N,1,0.5)
  p.alpha<-function(alpha){
    p<-1/(1+exp(-alpha-b1*x1-b2*x2-b3*x3))
    mean(p)-f0
  }
  solu<-uniroot(p.alpha,c(-20,20))
  round(unlist(solu),5)
}

#Set input values as N=10^6,b1=0,b2=1,b3=-1,f0=0.1,0.01,0.001,0.0001
#Derive corresponding alpha for each setting
N<-1000000
b1<-0;b2<-1;b3<--1
f0<-c(0.1,0.01,0.001,0.0001)
print("f0=0.1")
a.generator(N,b1,b2,b3,f0[1])
print("f0=0.01")
a.generator(N,b1,b2,b3,f0[2])
print("f0=0.001")
a.generator(N,b1,b2,b3,f0[3])
print("f0=0.0001")
a.generator(N,b1,b2,b3,f0[4])

## -----------------------------------------------------------------------------
#Plot -logf0 vs a
neglogf0.vs.a<-function(neglogf0){
  as.numeric(a.generator(N,b1,b2,b3,exp(-neglogf0))[1])
}
nlf0<-seq(2.3,9.2,0.1)
a<-numeric(70)
for(i in 1:70){
  a[i]<-neglogf0.vs.a(nlf0[i])
}
plot(nlf0,a,type="l",xlab="-logf0",ylab="a")

## -----------------------------------------------------------------------------
#Refer to codes of page 214 in Bayesian Analysis
rw.Metropolis<-function(sigma,x0,N){
  #sigma:standard variance of proposal distribution N(xt,sigma2)
  #x0:initial value
  #N:size of random numbers required
  x<-numeric(N)
  x[1]<-x0
  u<-runif(N)  
  k<-0
  for(i in 2:N){
    y<-rnorm(1,x[i-1],sigma)
    if(u[i]<=exp(-abs(y)+abs(x[i-1])))
      x[i]<-y
    else{
      x[i]<-x[i-1]
      k<-k+1
    }
  }
  return(list(x=x,k=k))
}

#Set initial values
sigma<-c(0.05,0.5,2.5,16)
N<-3000;x0<-25

#Derive corresponding chains with rw.Metropolis
rw1<-rw.Metropolis(sigma[1],x0,N)
rw2<-rw.Metropolis(sigma[2],x0,N)
rw3<-rw.Metropolis(sigma[3],x0,N)
rw4<-rw.Metropolis(sigma[4],x0,N)

## -----------------------------------------------------------------------------
#Draw graphs of chains
par(mar=c(1,1,1,1)) #display four graphs together
laplace.quantile<-c(-log(20),log(20))
refline<-laplace.quantile
rw<-cbind(rw1$x,rw2$x,rw3$x,rw4$x)
for(j in 1:4){
  plot(rw[,j],type="l",
       xlab=bquote(sigma==.(round(sigma[j],3))),
       ylab="X",ylim=range(rw[,j]))
  abline(h=refline)
}

## -----------------------------------------------------------------------------
#Compare rates of acceptance
cat("Rates of acceptance with sigma=0.05,0.5,2.5,16:","\n",1-c(rw1$k,rw2$k,rw3$k,rw4$k)/N,"\n")

## -----------------------------------------------------------------------------
#Refer to codes of page 227 in Bayesian Analysis
#Initialize constants and parameters
N<-5000 #length of chain
burn<-500 #length of burn-in
Z<-matrix(0,N,2) #storage of chain of bivariate sample
rho<-0.9 #correlation
mu1<-0;mu2<-0
sigma1<-1;sigma2<-1
s1<-sqrt(1-rho^2)*sigma1
s2<-sqrt(1-rho^2)*sigma2

#Generate corresponding chain
Z[1,]<-c(mu1,mu2) #initialize
for(i in 2:N){
  z2<-Z[i-1,2]
  m1<-mu1+rho*(z2-mu2)*sigma1/sigma2
  Z[i,1]<-rnorm(1,m1,s1)
  z1<-Z[i,1]
  m2<-mu2+rho*(z1-mu1)*sigma2/sigma1
  Z[i,2]<-rnorm(1,m2,s2)
}

#Plot generated sample after discarding a burn-in sample
b<-burn+1;z<-Z[b:N,]
plot(z,main="",cex=0.5,xlab=bquote(X),ylab=bquote(Y),ylim=range(z[,2]))

## -----------------------------------------------------------------------------
colnames(Z)<-c("x","y")
lm.gibbs<-lm(y~x,data=data.frame(Z))
cat("Coefficients of linear regression model:","\n")
lm.gibbs$coefficients
cat("Estimated residuals of linear regression model:",mean((Z[,2]-lm.gibbs$fitted.values)^2),"\n")

## -----------------------------------------------------------------------------
#Check the residuals of the model for normality and constant variance
res<-lm.gibbs$residuals
qqnorm(res);qqline(res)

## -----------------------------------------------------------------------------
#Refer to codes of page 206 in Bayesian Analysis
Gelman.Rubin<-function(psi){
  #psi[i,j] is the statistic psi(X[i,1:j])
  #for chain in i-th row of X
  psi<-as.matrix(psi)
  n<-ncol(psi)
  k<-nrow(psi)
  psi.means<-rowMeans(psi)
  B<-n*var(psi.means) #between variance estimate
  psi.w<-apply(psi,1,"var") #within variances
  W<-mean(psi.w) #within estimate
  v.hat<-W*(n-1)/n+(B/n) #upper variance estimate
  r.hat<-v.hat/W #Gelman Rubin statistic
  return(r.hat)
}

f<-function(x,sigma){
  if(any(x<0)) return(0)
  stopifnot(sigma>0)
  return((x/sigma^2)*exp(-x^2/(2*sigma^2)))
}

rayleigh.MH<-function(sigma,m,x0){
  x<-numeric(m)
  u<-runif(m)
  x[1]<-x0
  for(i in 2:m){
    xt<-x[i-1]
    y<-rchisq(1,df=xt)
    num<-f(y,sigma)*dchisq(xt,df=y)
    den<-f(xt,sigma)*dchisq(y,df=xt)
    if(u[i]<=num/den) x[i]<-y
    else
      x[i]<-xt
  }
  x
}

sigma<-4;m<-10000
k<-4;b<-1000

#Choose overdispersed initial values
#x0<-c(1,5,10,15)
x0<-c(0.5,1,5,10)

#Generate chains
X<-matrix(0,nrow=k,ncol=m)
for(i in 1:k)
  X[i,]<-rayleigh.MH(sigma,m,x0[i])

#Compute diagnostic statistics
psi<-t(apply(X,1,cumsum))
for(i in 1:nrow(psi))
  psi[i,]<-psi[i,]/(1:ncol(psi))
print(Gelman.Rubin(psi))

#Plot psi for four chains
par(mar=c(1,1,1,1))
for(i in 1:k)
  plot(psi[i,(b+1):m],type="l",
       xlab=i,ylab=bquote(psi))
par(mfrow=c(1,1))

## -----------------------------------------------------------------------------
#Plot sequence of R.hat statistics
rhat<-rep(0,m)
for(j in (b+1):m)
  rhat[j]<-Gelman.Rubin(psi[,1:j])
plot(rhat[(b+1):m],type="l",xlab="",ylab="R")
abline(h=1.1,lwd=1,col="blue",lty=2)
abline(h=1.2,lwd=1,col="blue",lty=1)

## -----------------------------------------------------------------------------
#Use coda
library(coda)
coda<-mcmc.list(as.mcmc(X[1,]),as.mcmc(X[2,]),as.mcmc(X[3,]),as.mcmc(X[4,]))
gelman.diag(coda)
gelman.plot(coda)

## -----------------------------------------------------------------------------
#Compute MLE maximizing log-likelihood of incomplete data
u<-c(11,8,27,13,16,0,23,10,24,2)
v<-c(12,9,28,14,17,1,24,11,25,3)
loglikelihood.incomplete<-function(lambda){
  lli<-0
  for(i in 1:10)
    lli<-lli+log(exp(-lambda*u[i])-exp(-lambda*v[i]))
  lli
}
curve(loglikelihood.incomplete,from=0.01,to=1.01)
optimize(loglikelihood.incomplete,lower=0.01,upper=0.201,maximum=TRUE)$maximum

## -----------------------------------------------------------------------------
#Compute MLE with EM  
u<-c(11,8,27,13,16,0,23,10,24,2)
v<-c(12,9,28,14,17,1,24,11,25,3)
EM.MLE<-function(u,v,lambda0=1,max.it=500,eps=1e-7){
  lambda<-c(lambda0,-1)
  for(i in 1:max.it){
    glambda<-lambda[1]*mean((u*exp(-lambda[1]*u)-v*exp(-lambda[1]*v))/(exp(-lambda[1]*u)-exp(-lambda[1]*v)))
    lambda[2]<-lambda[1]/(1+glambda)
    if(abs(lambda[1]-lambda[2])<eps) break
    lambda[1]<-lambda[2]
  }
  c(lambda[1],abs(lambda[1]-lambda[2]))
}
EM.MLE(u,v)[1]

## -----------------------------------------------------------------------------
#Refer to codes in section 11.9
solve.game<-function(A){
  #solve the two player zero-sum game by simplex method
  #optimize for player 1, then player 2
  #maximize v subject to...
  #let x strategies 1:m and put v as extra variable
  #A1 the<=constraints
  min.A<-min(A);A<-A-min.A 
  max.A<-max(A);A<-A/max.A
  m<-nrow(A);n<-ncol(A)
  it<-n^3;a<-c(rep(0,m),1) #objective function
  A1<--cbind(t(A),rep(-1,n)) #constraints <=
  b1<-rep(0,n)
  A3<-t(as.matrix(c(rep(1,m),0))) #constraints sum(x)=1
  b3<-1
  sx<-simplex(a=a,A1=A1,b1=b1,A3=A3,b3=b3,maxi=TRUE,n.iter=it)
  #the solution is [x1,x2..xm|value of game]
  #minimize v subject to...
  #let y strategies 1:n and put v as extra variable
  a<-c(rep(0,n),1) #objective function
  A1<-cbind(A,rep(-1,m)) #constraints <=
  b1<-rep(0,m)
  A3<-t(as.matrix(c(rep(1,n),0))) #constraints sum(y)=1
  b3<-1
  sy<-simplex(a=a,A1=A1,b1=b1,A3=A3,b3=b3,maxi=FALSE)
  soln<-list("A"=A*max.A+min.A,
             "x"=sx$soln[1:m],
             "y"=sy$soln[1:n],
             "v"=sx$soln[m+1]*max.A+min.A)
  soln
}
A<-matrix(c(0,-2,-2,3,0,0,4,0,0,
            2,0,0,0,-3,-3,4,0,0,
            2,0,0,3,0,0,0,-4,-4,
            -3,0,-3,0,4,0,0,5,0,
            0,3,0,-4,0,-4,0,5,0,
            0,3,0,0,4,0,-5,0,-5,
            -4,-4,0,0,0,5,0,0,6,
            0,0,4,-5,-5,0,0,0,6,
            0,0,4,0,0,5,-6,-6,0),9,9)
B<-A+2
library(boot)
sa<-solve.game(A)
sb<-solve.game(B)
sabx<-round(cbind(sa$x,sa$y,sb$x,sb$y),7)
colnames(sabx)<-c("Player1GameA","Player2GameA","Player1GameB","Player2GameB");sabx
cat("Value of game A:",sa$v,"\n")
cat("Value of game B:",sb$v,"\n")

## -----------------------------------------------------------------------------
scale01<-function(x){
  rng<-range(x,na.rm=TRUE)
  (x-rng[1])/(rng[2]-rng[1])
}

## -----------------------------------------------------------------------------
cat("list.example\n")
list(x=matrix(0,2,2),y=seq(4,1))
cat("\nunlist(list.example)\n")
unlist(list(x=matrix(0,2,2),y=seq(4,1)))
cat("\nas.vector(lise.example)\n")
as.vector(list(x=matrix(0,2,2),y=seq(4,1)))
cat("\nas.vector doesn't work\n")
identical(list(x=matrix(0,2,2),y=seq(4,1)),as.vector(list(x=matrix(0,2,2),y=seq(4,1))))

## -----------------------------------------------------------------------------
cat("vector.example\n")
seq(5,1)
cat("\ndim(vector.example)\n")
dim(seq(5,1))

## -----------------------------------------------------------------------------
cat("matrix.example\n")
matrix(seq(4,1),2,2)
cat("\nis.matrix(matrix.example)\n")
is.matrix(matrix(seq(4,1),2,2))
cat("\nis.array(matrix.example)\n")
is.array(matrix(seq(4,1),2,2))

## -----------------------------------------------------------------------------
cat("numeric.dataframe.example\n")
df1<-data.frame(x=1:3,y=seq(3,1),z=rep(1,3));df1
df1.scale<-as.data.frame(lapply(df1[,1:ncol(df1)],scale01));df1.scale
cat("\ndataframe.example\n")
df2<-data.frame(x=1:3,y=seq(3,1),z=c("a","b","c"));df2
df2.scale<-data.frame(as.data.frame(lapply(df2[,1:2],scale01)),z=df2[,3]);df2.scale

## -----------------------------------------------------------------------------
df1;vapply(df1[,1:3],sd,FUN.VALUE=c(1))

## -----------------------------------------------------------------------------
df2;vapply(df2[,1:2],sd,FUN.VALUE=c(1))

## ----eval=FALSE---------------------------------------------------------------
#  gibbsR<-function(a,b,n,N){
#    mat<-matrix(nrow=N,ncol=2)
#    mat[1,1]<-3;mat[1,2]<-0.5
#    for(i in 2:N){
#      y<-mat[i-1,2]
#      x<-rbinom(1,n,y)
#      y<-rbeta(1,x+a,n-x+b)
#      mat[i,]<-c(x,y)
#    }
#    mat
#  }

## ----eval=FALSE---------------------------------------------------------------
#  #include <Rcpp.h>
#  using namespace Rcpp;
#  
#  // [[Rcpp::export]]
#  NumericMatrix gibbsC(double a,double b,int n,int N){
#    NumericMatrix mat(N,2);
#    double x,y;
#    mat(0,0)=3;mat(0,1)=0.5;
#    for(int i=1;i<N;i++){
#    	y=mat(i-1,1);
#      mat(i,0)=rbinom(1,n,y)[0];
#      x=mat(i,0);
#      mat(i,1)=rbeta(1,x+a,n-x+b)[0];
#    }
#    return(mat);
#  }

## ----echo=FALSE,eval=FALSE----------------------------------------------------
#  library(Rcpp)
#  library(microbenchmark)
#  source(file="gibbsR.R")
#  sourceCpp(file="gibbsC.cpp")
#  ts<-microbenchmark(gibbR=gibbsR(3,4,15,100),
#                     gibbC=gibbsC(3,4,15,100))
#  summary(ts)[,c(1,3,5,6)]

