// this file is only for fcns that call other user defined functions
//  [ ] logLikNB - calculating the log-likelihood of the right-truncated negative binomial distr.
//  [ ] getr     - curve is the estimated epidemic curve (i.e., back-calculated counts & nowcasted counts)
#include <RcppArmadillo.h>
#include <fstream>
using namespace Rcpp;
using namespace std;

double logLikNB(Rcpp::NumericVector delay_vec, Rcpp::NumericMatrix x12, 
                Rcpp::NumericVector disp, Rcpp::NumericVector betaplus, int maxdelay);

Rcpp::NumericVector getr(Rcpp::NumericVector curve, Rcpp::NumericVector si, int size);
