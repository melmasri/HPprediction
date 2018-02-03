#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <RcppEigen.h>

// [[Rcpp::depends(BH)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppArmadillo, RcppEigen)]]
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
mat muliC_arma(mat & a, mat & Z){
  /* a function for matrix multiplication faster than R */
  return a * Z;
}

// [[Rcpp::export]]
NumericVector raffinity_MH_Rcpp(NumericVector & old,NumericVector & z,NumericVector & Ud, NumericVector & sig, NumericVector & hyper){
  /* MH sampler for affinity parameters */
  unsigned int size = old.size();
  RNGScope scope;
  NumericVector epsilon = sig * Rcpp::rnorm(size);
  NumericVector newpara = ifelse(old  + epsilon  <= 0, old - epsilon,old+ epsilon);
  NumericVector likeli = (z+hyper(0)-1)*log(newpara/old) -  (newpara - old)*(hyper(1) +Ud);
  NumericVector u = runif(size);
  return Rcpp::wrap(ifelse(u <= pmin(1, exp(likeli)), newpara, old));
}

// // [[Rcpp::export]]
// SEXP to_matrix_coph(SEXP eb, Function & f, IntegerVector & lind, IntegerVector & uind, int & n){
//   NumericVector M(n^2);
//   NumericVector res =  f(eb);
//   M[lind] = res;
//   // M.attr("dim") = Dimension(n,n);
//   return M;
// }

// // [[Rcpp::export]]
// IntegerVector rEta_C(const SEXP & eb, Function & f, const int nrow){
//
//   IntegerVector i = seq_len(nrow-1);
//   IntegerVector vec = sapply(i, g);
//   return vec;
// }


// // [[rcpp::export]]
// SEXP lowerTri(int nrow){
//   IntegerVector x = seq_len(nrow-1);
//   Rcout << x << std::endl;
//   IntegerVector vec = sapply(x,index_fun);
//   return 0; 
// }




