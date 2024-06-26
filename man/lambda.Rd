% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{lambda}
\alias{lambda}
\title{Calculate Lambda Values for Poisson Distribution Mean}
\usage{
lambda(curve, si)
}
\arguments{
\item{curve}{NumericVector representing the daily counts of new cases, estimated from the epidemic data.
This vector should include both historically observed data and nowcasted estimates.}

\item{si}{NumericVector representing the serial interval distribution, a probability distribution
describing the time delay between successive cases in an infectious disease transmission chain.}
}
\value{
NumericVector containing the estimated mean values for the Poisson distribution
        at each day, excluding the last day as the serial interval cannot be fully applied.
}
\description{
This function computes the convolution of the epidemic curve with the serial interval distribution
to estimate the mean of the Poisson distribution for each day. This mean is crucial
in models where the number of new cases follows a Poisson process. The calculation is effectively
a weighted sum of past case counts, where weights are given by the serial interval distribution,
representing the expected number of new cases caused by an individual case on each subsequent day.
}
\details{
The function applies the serial interval to the epidemic curve via convolution, 
         essentially calculating the expected number of secondary cases generated by each primary case
         over the serial interval. The length of the resulting vector is one less than that of `curve`
         because the last day's value cannot be calculated without full serial interval data.
}
\examples{
curve <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
si <- c(0.5, 0.3, 0.2)  # Declining probability over three days
lambda_values <- lambda(curve, si)
}
