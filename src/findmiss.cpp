#include <RcppArmadillo.h>
#include <fstream>
using namespace Rcpp;
using namespace std;
// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::depends(RcppArmadillo)]]

//' Locate Missing Values in a Numeric Vector
//'
//' This function identifies the indices of missing values (NA) in a given numeric vector.
//' It is useful for data cleaning and preprocessing steps where identification of missing
//' data points is required.
//'
//' @param x A numeric vector potentially containing NA values.
//'          The values can range from -infinity to infinity.
//' @return An integer vector containing the indices where NA values are found in `x`.
//'         These indices can be used directly to reference or manipulate elements in other vectors
//'         of the same length or for subsetting the original vector.
//' @examples
//' vec <- c(1, 2, NA, 4, NA, 6)
//' findmiss(vec)
//' @export
// [[Rcpp::export]]
IntegerVector findmiss(NumericVector x){
  // Description:
  //   function that returns the indices of NA values of a vector
  //
  // Arguments:
  //   x: a numeric vector, -inf to inf
  //
  // Returns: 
  //   the indices of x for which is_na(x) returns True
  IntegerVector y = seq(0, x.size() - 1);
  LogicalVector z = is_na(x);
  return y[z];
}
