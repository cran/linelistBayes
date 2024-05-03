#include <RcppArmadillo.h>
#include <fstream>
using namespace Rcpp;
using namespace std;
// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::depends(RcppArmadillo)]]

//' Calculate Exponential of Linear Combinations
//'
//' Computes the exponential of linear combinations of beta coefficients and a matrix
//' of predictors, typically used in Poisson or logistic regression models for estimating
//' rates or probabilities. This function specifically handles the exponential transformation,
//' which is commonly used to ensure that rates or probabilities are non-negative.
//'
//' @param x12 NumericMatrix representing a matrix of predictors, where each row corresponds
//'        to an observation and columns correspond to different predictor variables (e.g., weeks and weekends).
//' @param beta NumericVector of coefficients corresponding to the predictors in `x12`.
//'        This should include coefficients for both weekly effects and potentially an additional
//'        coefficient for weekends.
//' @return NumericVector where each element is the exponential of the linear combination
//'         of the predictors and coefficients for a given observation. This vector represents
//'         the model-estimated mean values for each observation.
//' @details The function multiplies the matrix `x12` by the vector `beta` to get the linear predictors,
//'          then applies the exponential function to convert these linear predictors to a scale
//'          suitable for models where the response variable is a count or probability.
//'          This is a critical step in generalized linear models where the link function is
//'          the natural logarithm.
//' @examples
//' # Assuming x12 is a matrix with 10 observations and 3 predictors
//' # and beta is a vector of 3 coefficients
//' x12 <- matrix(rnorm(30), ncol=3)
//' beta <- c(0.1, -0.2, 0.05)
//' get_mu_vec(x12, beta)
//' @export
// [[Rcpp::export]]
NumericVector get_mu_vec(NumericMatrix x12, NumericVector beta) {
  // Description:
  //   Function to compute mu=exp(X1*beta + X2*gamma), the parameter representing the mean value 
  //
  // Arguments:
  //   x12:  matrix of indicator variables for is_week_i and is_weekend, n x p
  //   beta: existing week-wise beta coefficients and a single weekend parameter, p x 1
  //         this includes what is called beta in the text
  //         as well as gamma, the single parameter for weekend
  //         but shortened to just beta for this script
  //
  // Returns: 
  //   a vector, exp(x12 %*% beta)
  
  // Conversions into armadillo syntax for ease of processing in 
  // Armadillo linear algebra
  // the a_ represents the same variable as above but in armadillo
  arma::mat a_x12  = as<arma::mat>(x12);  // dim(x1) is n observations x n weeks
  arma::vec a_beta = as<arma::vec>(beta); // beta are corresponding model coefficients

  // matrix multiplication, X1*B + X2*gamma
  arma::mat a_mu_mat = a_x12 * a_beta;

  // wrap() is conversion back into an R object
  NumericMatrix mu_mat = wrap(a_mu_mat);

  // getting all the rows (_) of the first column (0) of lc
  // it MUST have one 1 column
  if (mu_mat.ncol() != 1) {
        stop("ncol of mu_vector must = 1");
  }
  NumericVector mu_vec = mu_mat(_, 0); 

  // taking the exponent and returning
  NumericVector exp_mu_vec  = exp(mu_vec);

  return exp_mu_vec;
}
