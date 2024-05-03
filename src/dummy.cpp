#include <RcppArmadillo.h>
#include <fstream>
using namespace Rcpp;
using namespace std;
// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::depends(RcppArmadillo)]]

//' Generate Dummy Variables Matrix for Weeks and Weekends
//'
//' This function creates a matrix of dummy variables based on reported weeks and weekend indicators.
//' Each column in the resulting matrix corresponds to a specific week, except for the last column,
//' which indicates whether the date falls on a weekend. This matrix is typically used in regression
//' models where week-specific effects are to be adjusted along with the effect of weekends.
//'
//' @param week An integer vector representing the week number of each observation.
//'        Each element denotes the week during which a specific event occurred.
//' @param weekend A binary integer vector (elements being 0 or 1) indicating whether
//'        each observation corresponds to a weekend. Here, `1` indicates a weekend and `0` a weekday.
//' @return A numeric matrix where each row corresponds to an observation and each column
//'         to a week, with an additional final column for weekend indicators. The elements
//'         of the matrix are dummy variables (0 or 1); each row contains exactly one '1' in
//'         one of the first several columns corresponding to the week of the observation,
//'         and a '1' or '0' in the last column indicating weekend status.
//' @examples
//' week <- c(1, 1, 1, 2, 2, 3, 3, 3, 3, 3)
//' weekend <- c(0, 1, 0, 0, 1, 0, 0, 1, 1, 0)
//' dummy_vars <- dummy(week, weekend)
//' @export
// [[Rcpp::export]]
NumericMatrix dummy(IntegerVector week, IntegerVector weekend) {
  // Description:
  //   
  //
  // Arguments
  //   week:     d$week of the report, allows you to assume that things don't change within this window
  //   weekend:  d$weekend, 0 or 1, indicator of weekend (1)
  //
  // Returns:
  //   cov: a matrix of dummy variables where the first n columns are is_week_i
  //        and the last column is is_weekend. the rows are data for individuals
  //
  
  int nw = max(week);   // maximum number of weeks
  int nc = nw + 1;      // adds an extra column for weekends
  int n  = week.size(); // this is the rows of line-list data (so 1 per person) 
  NumericMatrix cov(n, nc); // the final matrix has a row for each person and a column for each week + 
                            // one additional column for is_weekend

  for (int i = 0; i < nc; ++i) {
    if (i < nw) {
      cov(_,i) = ifelse(week == (i+1), 1, 0); // see above
    } else {
      cov(_,i) = weekend; // see above
    }
  }
  
  return cov;
}
