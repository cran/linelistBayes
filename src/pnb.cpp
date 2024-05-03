#include <RcppArmadillo.h>
#include <fstream>
using namespace Rcpp;
using namespace std;
// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::depends(RcppArmadillo)]]

//' Compute Cumulative Distribution Function for Negative Binomial Distribution
//'
//' This function calculates the cumulative distribution function (CDF) of the negative binomial distribution
//' given a number of successes, a dispersion parameter, and the mean. The negative binomial distribution 
//' is commonly used to model count data with overdispersion relative to a Poisson distribution.
//'
//' @param x Non-negative integer specifying the number of successes for which the CDF is computed.
//'          This value must be non-negative as it represents the number of successes in the distribution.
//' @param r Dispersion parameter of the distribution, a positive real number.
//'          Higher values of 'r' indicate a higher probability of counts clustering around the mean, reducing overdispersion.
//' @param m Mean of the distribution, a positive real number indicating the expected number of successes.
//'          The mean must be positive.
//' @return y The probability of observing up to 'x' successes in a negative binomial distribution
//'           parameterized by 'r' (dispersion) and 'm' (mean). This function returns the CDF value.
//' @details The negative binomial distribution can be parameterized by a dispersion parameter 'r' and a mean 'm',
//'          which together determine the shape of the distribution. This function is essential for modeling and
//'          probability calculations in various fields such as epidemiology and ecology where the negative binomial
//'          distribution is used to model count data.
//' @examples
//' pnb(10, 5, 10)
//' @export
// [[Rcpp::export]]
double pnb(int x, double r, double m){
  // Description:
  //   pnb function used by mapply for likelihood.
  //   These functions provide information about the 
  //   Negative binomial distribution with dispersion r and mean m
  //   R:: implementation chosen for speed, returns a scalar
  //
  // Arguments:
  //   x: a double random variable, 0 to inf
  //   r: the dispersion parameter, size. higher r means tighter around mu 
  //   m: the mean, mu > 0
  //
  // Returns: 
  //   the Cumulative density function value at x 
  //   of a distribution given by r and m
  
  // implied in the function:
  //   lower: true --> Calculate the probability of the region where 
  //                 the random variable is less than or equal to x
  //   log: true --> probabilities p are given as log(p)
  // 
  double y = R::pnbinom_mu(x, r, m, true, true);
  return y;
}
