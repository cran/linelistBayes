#' Calculate Serial Interval Distribution for COVID-19
#'
#' This function computes the probability distribution function (PDF) of the serial interval for COVID-19
#' using a gamma distribution with specified shape and rate parameters. The serial interval is defined
#' as the time between successive cases in a chain of transmission. This implementation generates a discrete
#' PDF over a given number of days.
#'
#' The function uses the `pgamma` function to calculate cumulative probabilities for each day up to `ndays`
#' and then differences these to get daily probabilities. The resulting probabilities are normalized to sum to 1,
#' ensuring that they represent a valid probability distribution.
#'
#' @param ndays Integer, the number of days over which to calculate the serial interval distribution.
#' @param alpha Numeric, the shape parameter of the gamma distribution.
#' @param beta Numeric, the rate parameter of the gamma distribution.
#'
#' @return Numeric vector representing the serial interval probabilities for each of the first `ndays` days.
#'         The probabilities are normalized so that their sum is 1.
#'
#' @examples
#' sip <- si(14, 4.29, 1.18)
#' @references
#' Nishiura, H., Linton, N. M., & Akhmetzhanov, A. R. (2020). Serial interval of novel coronavirus (COVID-19) infections.
#' International Journal of Infectious Diseases, 93, 284-286.
#' [Link to the article](https://www.sciencedirect.com/science/article/pii/S1201971220306111)
#'
#' @export
si <- function(ndays, alpha, beta) {
  
  prob <- numeric(ndays) # creates a numeric vector
  
  for (i in 1:ndays){
    prob[i] <- pgamma(i, shape = alpha, rate = beta) -
      pgamma(i - 1, shape = alpha, rate = beta)
  }
  
  result <- prob/sum(prob) # normalizes the whole vector
  
  return(result)
}
