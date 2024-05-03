#include <RcppArmadillo.h>
#include <fstream>
using namespace Rcpp;
using namespace std;
// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::depends(RcppArmadillo)]]

//' Calculate Log-Probability Density for Negative Binomial Distribution
//'
//' This function computes the log-probability density function of the negative
//' binomial distribution, given the number of successes, the dispersion parameter,
//' and the mean of the distribution. This is useful for probabilistic models
//' where negative binomial assumptions are applicable, such as certain types of
//' count data.
//'
//' @param x A numeric value representing the number of successes, which should be
//'          a non-negative integer.
//' @param r A numeric value representing the dispersion parameter of the negative
//'          binomial distribution. A higher value indicates a distribution more tightly
//'          concentrated around the mean.
//' @param m A numeric value representing the mean of the distribution.
//' @return A numeric value representing the log-probability density function value
//'         of observing 'x' successes given the mean 'm' and dispersion 'r'.
//'         This return value is given on the log scale to facilitate calculations
//'         that involve very small probabilities.
//' @examples
//' dnb(5, 2, 10);
//' @export
// [[Rcpp::export]]
double dnb(double x, double r, double m){
 // Implementation details:
 //   The function uses the R::dnbinom_mu function to calculate the log-probability
 //   density of a negative binomial distribution with parameters r (dispersion)
 //   and m (mean). The calculation is performed on the log scale for numerical
 //   stability when dealing with small probabilities.
 //
 // Arguments:
 //   x: a non-negative integer, the number of successes
 //   r: a positive real number, the dispersion parameter
 //   m: a positive real number, the mean of the distribution
 //
 // Returns: 
 //   The log of the probability density function value at x, which is a real number
 
 double y = R::dnbinom_mu(x, r, m, true);
 return y;
}
 
