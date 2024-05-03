#include <RcppArmadillo.h>
#include <fstream>
using namespace Rcpp;
using namespace std;
// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::depends(RcppArmadillo)]]

#include "util.h"

//' Compute Log-Likelihood for a Right-Truncated Negative Binomial Model
//'
//' This function calculates the log-likelihood of observing the given data under a 
//' right-truncated negative binomial distribution. It is used to assess the fit of 
//' a model involving delay times in reporting cases, where each case has a 
//' delay modeled by a negative binomial distribution that is truncated at a specified maximum delay.
//'
//' @param delay_vec NumericVector representing the observed delay times for each case.
//' @param x12 NumericMatrix of covariates used to model the mean of the negative binomial distribution.
//'        Each row corresponds to a case and columns correspond to covariates including time since onset and others.
//' @param disp NumericVector indicating the dispersion parameter for each case. This can affect the distribution 
//'        of counts and is used to model heterogeneity in reporting delays.
//' @param betaplus NumericVector containing current estimates of the model parameters, which may include 
//'        coefficients for covariates and additional parameters.
//' @param maxdelay Integer specifying the maximum reporting delay, truncating the distribution at this value.
//' @return Double representing the log-likelihood of the data given the model parameters. This value 
//'         measures how well the model with the current parameter estimates fits the observed data.
//' @details This function computes the log-likelihood by calculating the likelihood of each observed 
//'          delay under the specified model parameters, considering the truncation at `maxdelay`. 
//'          The parameters `disp` and `betaplus` allow for flexibility in modeling different types of 
//'          heterogeneity and covariate effects.
//' @examples
//' # Example usage with arbitrary data:
//' delay_vec <- rnorm(100, mean=10, sd=3)  # Simulated delay times
//' x12 <- matrix(rnorm(300), ncol=3)       # Simulated covariates
//' disp <- rep(1, 100)                     # Dispersion parameter, constant for simplicity
//' betaplus <- runif(4)                    # Simulated parameter estimates
//' maxdelay <- 15                          # Maximum delay for truncation
//' loglik_value <- logLikNB(delay_vec, x12, disp, betaplus, maxdelay)
//' @export
// [[Rcpp::export]]
double logLikNB(NumericVector delay_vec, NumericMatrix x12, 
            NumericVector disp, NumericVector betaplus, int maxdelay) {
  // Description:
  //   TL: faster version of like for large sample
  //   CM: this is the log-likelihood function for the right-truncated negative binomial distribution
  //
  // Arguments
  //   delay_vec:  vector, reporting delay vector
  //   x12:        matrix, cov, the indicator values for is_week_i & is_weekend
  //   disp:       vector, the flag for each outcome, =1 if cd is null, otherwise =1 and =2 if some other condition is met
  //                       I think this is the dispersion associated with changing cd? unclear
  //   betaplus:   vector, current row of parameter estimates, betas(gammas) and .... other params???
  //   maxdelay:   integer, maximum reporting delay
  //
  // Returns:
  //   the logLiklihood of the right-truncated negative binom distribution given the current param set
  
  // ----------------
  // get the number of columns of the indicator matrix
  int nbeta = x12.ncol();        // number of parameters estimated for is_week_i and is_weekend
  int n     = delay_vec.size();  // number of people in the line-list dataset
  NumericVector s_vec;       // vector of the size parameters
  double r1;                 // either the first index of the two sizes for the two dispersion or not
  
  // ----------------
  // splits betaplus into two vectors, from 0:(nbeta-1), and then nbeta:extra
  // beta = the beta coefficients for each week and is_weekend
  NumericVector beta = betaplus[seq(0, nbeta - 1)];   

  // the remainder are the other parameters being estimated which are ...?
  // TODO: The additional parameters that are being estimated, which are ...?
  // if cd is defined then this is just one extra
  // if cd is not defined then there are two extra
  NumericVector r    = betaplus[seq(nbeta, betaplus.size() - 1)]; // these are the additional params, either 1 or 2
  if (! ((r.size() == 1) || (r.size() == 2))) {
    stop("size of r is not in c(1, 2)");
  }

  // ----------------
  // and then mu_vec is the expected value for the mean for each person's reporting delay distribution
  // get mu, = exp(X1*beta + X2*gamma), one value for each person, E(delay)
  NumericVector mu_vec   = get_mu_vec(x12, beta); 

  // ----------------
  // TODO: ???
  // (basically, if this is type = 0, if cd is not defined)
  if (max(disp) == 1){ // originally this was max but we probably want all() right?
    // ### ERROR CHECKING ###
    // if(! (all(disp == 1))) {
    //   stop("all disp != 1");
    // }
    if(! (r.size() == 1)) {
      stop("r.size() must = 1");
    }
    // ####
    r1 = as<double>(r); // type conversion to a scalar
    s_vec = disp * r1;  // you know that disp = 1 everywhere, so this is just repeated r right?
  } else {
    // ### ERROR CHECKING ###
    // these is the case where cd is defined and where disp == 2
    if(! (r.size() == 2)) {
      stop("r.size() must = 2");
    }
    //if(! ((max(disp) == 2) && (min(disp) == 1))) {
    //  stop("all disp not in 1 or 2");
    //}
    // ####
    r1 = r[0];            // set r1 equal to the first dispersion parameter
    s_vec = rep(r1, n);   // right, why not do this above?
    LogicalVector select = (disp==2);
    double r2 = r[1];     // the two dispersions
    s_vec[select] = r2;
  }
  
  // ************
  // Now get the log-liklihood of the right-truncated negative binomial distribution
  // This is manifested as the ratio between P(x[i]) and sum{P(x[i] <= K)} for each combination of variables
  // (i.e., size = s[i] and mean = mu[i])

  // (1) First, the numerator
  //     for each reporting delay value (x[i] in delay_vec)
  //     calculate the P(x[i]), the probability of that delay
  //     given a negative binomial distribution defined by size = s[i] and mean = mu[i]
  NumericVector nbinom_pmf_vec = mapply(delay_vec, s_vec, mu_vec, dnb);

  // (2) Then, calculate the denominator(?)
  //     create a vector that looks like delay_vec but is instead the max reporting delay
  IntegerVector maxdelay_vec  = rep(maxdelay, n);

  // and then get the CDF value for max reporting delay for a distribution given size? = s, and mean = m1
  // this varies for every value because the size parameter and mu are different for each distribution
  NumericVector max_nbinom_cdf = mapply(maxdelay_vec, s_vec, mu_vec, pnb);

  // then, the result (log-likelihood of all the right-truncated negative binomial distribution^s^) 
  // is the difference of all summed numerators and denominators
  double result = sum(nbinom_pmf_vec) - sum(max_nbinom_cdf);
  // ************
  
  //
  return result;
}
