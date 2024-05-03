#include <RcppArmadillo.h>
#include <fstream>
using namespace Rcpp;
using namespace std;

//' Calculate Proportions of Event Counts Within a Specified Time Range
//'
//' This function calculates the proportion of event counts that occur within each
//' unit of time from a specified starting point (onset time) up to a maximum delay.
//' The proportions are computed relative to the total number of events occurring within
//' the specified time range.
//'
//' @param x NumericVector representing the count of events at each time point.
//' @param onset NumericVector representing the onset times for each event count.
//'          Each element in the vector indicates the time at which an event was initiated.
//' @param maxdelay The maximum time delay (inclusive) for which the proportions are calculated.
//'          This parameter defines the upper bound of the time interval over which the proportions
//'          of event counts are evaluated.
//' @param cd The onset time (exclusive) from which to start calculating proportions.
//'          Events that start at times less than or equal to this value are excluded from the calculation.
//' @return NumericVector containing the proportions of the total events falling within each unit
//'         of time from the specified 'cd' up to 'maxdelay'. The proportions are cumulative,
//'         with each element representing the proportion of events that have occurred by that time point,
//'         starting from 'cd'.
//'
//' @examples
//' x <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
//' onset <- c(1, 1, 2, 2, 3, 3, 4, 4, 5, 5)
//' maxdelay <- 5
//' cd <- 1
//' prop(x, onset, maxdelay, cd)
//' @export
// [[Rcpp::export]]
NumericVector prop(NumericVector x, NumericVector onset, int maxdelay, int cd) {

  LogicalVector v = (x <= maxdelay) & (onset >= cd);
  NumericVector x1 = x[v];
  int dem = x1.size();
  NumericVector p1 (maxdelay);

  //
  for (int i = 0; i < maxdelay; ++i){
    p1[i] = sum(x1 == (maxdelay - i));
  }

  NumericVector p = p1 / dem;
  NumericVector p2 = cumsum(p);
  NumericVector result = 1 - p2;

  return result;
}
