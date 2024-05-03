// this file is only for fcns that Dont call other user defined functions
//  [x] findmiss - finds missing values in a vector
//  [x] dnb      - get the PMF value for a given x in a standard negative binomial distribution
//  [x] pnb      - get the CDF value for a given x in a standard negative binomial distribution
//  [x] rnb      - get a random draw from a negative binomial distribution
//  [x] get_mu_vec - function to get a vector of mu, line-specific estimates of mean of right-trunc neg.binom
//  [x] dummy    - create a matrix of indicator values for is_week_i and is_weekend
//  [ ] lambda   - create a lambda function to compute mean of poisson dist.
//  [ ] prop     - function to make proportion of counts from 0 to maxdelay
#include <RcppArmadillo.h>
#include <fstream>
using namespace Rcpp;
using namespace std;

Rcpp::IntegerVector findmiss(Rcpp::NumericVector x);
double dnb(double x, double r, double m);
double pnb(int x, double r, double m);
double rnb(double r, double p);
Rcpp::NumericVector get_mu_vec(Rcpp::NumericMatrix x12, Rcpp::NumericVector beta);
Rcpp::NumericMatrix dummy(Rcpp::IntegerVector week, Rcpp::IntegerVector weekend);
Rcpp::NumericVector prop(Rcpp::NumericVector x, Rcpp::NumericVector onset, int maxdelay, int cd);
Rcpp::NumericVector lambda(Rcpp::NumericVector curve, Rcpp::NumericVector si);
