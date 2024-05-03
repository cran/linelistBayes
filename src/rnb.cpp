#include <RcppArmadillo.h>
#include <fstream>
using namespace Rcpp;
using namespace std;

//' Generate Random Samples from a Negative Binomial Distribution
//'
//' This function generates random samples from a negative binomial distribution with the specified dispersion parameter (\code{r}) and success probability (\code{p}).
//'
//' @param r The number of failures before achieving a specified number of successes in a negative binomial experiment. It also serves as the dispersion parameter which controls the variance of the distribution.
//' @param p The probability of success on each independent Bernoulli trial within the negative binomial experiment.
//' @return A random value sampled from the negative binomial distribution with parameters \code{r} and \code{p}.
//'
//' @examples
//' r <- 2
//' p <- 0.3
//' rnb(r, p)
//' @export
// [[Rcpp::export]]
double rnb(double r, double p){
  // Description:
  // rnb function used by mapply for likelihood.
  // These functions provide information about the 
  // Negative binomial distribution with dispersion r and mean m
  // R:: implementation chosen for speed, returns a scalar
  //
  // Arguments:
  //   r: the dispersion parameter, size
  //   p: a probability, ...?
  //
  // Returns: 
  //   a single random value in this distribution, given p?
  
  double y = R::rnbinom(r, p);
  return y;
}
