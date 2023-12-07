#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix gibbsC(double a,double b,int n,int N){
  NumericMatrix mat(N,2);
  double x,y;
  mat(0,0)=3;mat(0,1)=0.5;
  for(int i=1;i<N;i++){
    y=mat(i-1,1); 
    mat(i,0)=rbinom(1,n,y)[0];
    x=mat(i,0);
    mat(i,1)=rbeta(1,x+a,n-x+b)[0];
  }
  return(mat);
}

//' @title Compute sigma.hat/sigma.star for evaluation
//' @description More details refer to scaled sparse linear regression
//' @param x the design matrix
//' @param y the response vector
//' @param beta_star the true coefficient vector of beta
//' @param sigma_hat the estimator of sigma
//' @return sigma_hat/sigma_star
//' @examples 
//' \dontrun{
//'    compute_ratio_cpp(d1[,1:1000],as.vector(d1[,1001]),rep(0,1000),1)
//' }
//' @export
// [[Rcpp::export]]
double compute_ratio_cpp(NumericMatrix x, NumericVector y, NumericVector beta_star, double sigma_hat) {
  int n = y.length();
  NumericVector residuals = y - x * beta_star;
  double sigma_star = sqrt(sum(pow(residuals,2)) / n);
  return sigma_hat / sigma_star;
}
