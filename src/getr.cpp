#include <RcppArmadillo.h>
#include <fstream>
using namespace Rcpp;
using namespace std;
// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::depends(RcppArmadillo)]]

#include "util.h"

//' Calculate Time-Varying Reproduction Number R(t)
//'
//' This function estimates the time-varying reproduction number, R(t), based on the epidemic curve and
//' serial interval distribution. R(t) is calculated for each day using a moving window approach,
//' which involves taking a segment of the epidemic curve and applying a transformation based on the
//' serial interval to estimate how many subsequent cases are generated by cases within the window.
//'
//' @param curve NumericVector representing the estimated epidemic curve with daily counts.
//'              This curve can include both back-calculated and nowcasted counts of infections.
//' @param si NumericVector representing the serial interval distribution, expressed as a vector
//'           where each element corresponds to the probability of a delay of that many days between
//'           successive cases.
//' @param size Integer specifying the size of the moving window used to calculate the mean reproduction number.
//'             This window size determines how many days are included in the calculation of R(t) at each step.
//' @return NumericVector containing the estimated mean reproduction numbers (R(t)) for each day.
//'         The length of this vector will be the length of `curve` minus `size` minus one,
//'         reflecting the fact that the last few days do not have enough data to fill the window.
//' @examples
//' # Assume curve is a numeric vector of daily case counts and si is the serial interval distribution
//' curve <- rnorm(100, mean=10, sd=3)  # example epidemic curve
//' si <- rep(0.1, 10)  # example serial interval distribution
//' size <- 6  # example size of the moving window
//' getr(curve, si, size)
//' @export
// [[Rcpp::export]]
NumericVector getr(NumericVector curve, NumericVector si, int size){
  // Description:
  //    curve is the estimated epidemic curve (i.e., back-calculated counts & nowcasted counts)
  //    Compute the mean only
  //
  // Arguments
  //    curve: i.e., the current backcaulation, dim of (nd + maxdelay)   
  //    si: serial interval of 14?
  //    size: NB_size (default is 6)
  // Returns:
  //
  
  int n = curve.size();
  int nr = n - size - 1;
  NumericVector result (nr);
  double shape;
  double scale;
  NumericVector incid;
  NumericVector dem = lambda(curve, si);
  NumericVector dem1;
  for (int i = 0; i < nr; ++i){
    incid     = curve[seq(i + 1, i + size + 1)];  //daily counts for days i to i + size
    shape     = sum(incid) + 1;
    dem1      = dem[seq(i, i + size)];
    scale     = 1 / (sum(dem1) + 0.2);
    result[i] = shape * scale ; 
  }
  return result;
}
