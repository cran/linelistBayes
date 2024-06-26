% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{findmiss}
\alias{findmiss}
\title{Locate Missing Values in a Numeric Vector}
\usage{
findmiss(x)
}
\arguments{
\item{x}{A numeric vector potentially containing NA values.
The values can range from -infinity to infinity.}
}
\value{
An integer vector containing the indices where NA values are found in `x`.
        These indices can be used directly to reference or manipulate elements in other vectors
        of the same length or for subsetting the original vector.
}
\description{
This function identifies the indices of missing values (NA) in a given numeric vector.
It is useful for data cleaning and preprocessing steps where identification of missing
data points is required.
}
\examples{
vec <- c(1, 2, NA, 4, NA, 6)
findmiss(vec)
}
