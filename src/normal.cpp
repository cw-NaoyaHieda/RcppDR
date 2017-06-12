#include <Rcpp.h>
#include <RcppCommon.h>

using namespace Rcpp;


// [[Rcpp::export]]
NumericVector rnorm_cpp(const int N,double mean,double sd) {
  RNGScope scope;		// ensure RNG gets set/reset
  NumericVector X(N);
  X = rnorm(N,mean,sd);
  return X;
}

// [[Rcpp::export]]
NumericVector dnorm_cpp(std::vector<double> Q,double mean,double sd) {
  RNGScope scope;		// ensure RNG gets set/reset
  NumericVector q(Q.size());
  q = Q;
  NumericVector X(q.size());
  X = dnorm(q,mean,sd);
  return X;
}

// [[Rcpp::export]]
NumericVector g_DR_cpp(std::vector<double> dr,std::vector<double> pd,std::vector<double> RHO) {
  NumericVector DR(dr.size());
  NumericVector PD(pd.size());
  NumericVector rho(RHO.size());
  DR = dr;
  PD = pd;
  rho = RHO;
  NumericVector density = sqrt((1 - rho)/rho)*exp(0.5*(pow(qnorm(DR),2)-pow((sqrt(1-rho)*qnorm(DR)-qnorm(PD))/sqrt(rho),2)));
  return density;
}

// [[Rcpp::export]]
NumericMatrix rmvnorm_rcpp(int N, NumericVector mean, NumericVector sigma){
  RNGScope scope;		// ensure RNG gets set/reset
  NumericMatrix z(N,mean.size());
  for(int i = 0;i < mean.size(); i++){
    z(_,i) = rnorm(N,mean(i),sigma(i));
  }
  return(z);
}

// [[Rcpp::export]]
double dmvnorm_rcpp(NumericVector x, NumericVector mean, NumericVector sigma){
  double pi = 3.141592653589793238462643383280;
  double y = 1;
  for(int i = 0;i < mean.size(); i++){
    y = y * 1/sqrt(2*pi)*sigma(i) * exp(pow(x(i)-mean(i),2)/sqrt(2*pi*sigma(i)));
  }
  return(y);
}
// [[Rcpp::export]]
double tmp(double x){
  return(x);
}

// [[Rcpp::export]]
NumericVector vector_f_sum(List input, Function f) {
  int n = input.size();
  List out(n);
  NumericVector output(n);
  NumericVector tmp(n);
  for(int i = 0; i < n; i++) {
    out[i] = f(input[i]);
    tmp = out[i];
    output[i] = sum(tmp);
  }
  return output;
}


