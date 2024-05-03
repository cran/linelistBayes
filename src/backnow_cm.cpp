// -------------------------------------------------------------------------------
// Author: Li Tenglong
// Annotations by CWM
// Date started: 12.20.2023
// Notes:
// * anytime you see something in double square brackets, 
//   its important for the compiler
//
// Functions:
//  findmiss - finds missing values in a vector
//  dnb      - get the PMF value for a given x in a standard negative binomial distribution
//  pnb      - get the CDF value for a given x in a standard negative binomial distribution
//  rnb      - get a random draw from a negative binomial distribution
//  get_mu_vec - function to get a vector of mu, line-specific estimates of mean of right-trunc neg.binom
//  dummy    - create a matrix of indicator values for is_week_i and is_weekend
//  logLikNB - calculating the log-likelihood of the right-truncated negative binomial distr.
//  prop     - function to make proportion of counts from 0 to maxdelay
//  lambda   - create a lambda function to compute mean of poisson dist.
//  getr     - curve is the estimated epidemic curve (i.e., back-calculated counts & nowcasted counts)
//  backnow  - main function to do Gibbs MCMC, impute missing delays, back-calculation, then now-casting
//
// -------------------------------------------------------------------------------
#include <RcppArmadillo.h>
#include <fstream>
using namespace Rcpp;
using namespace std;
// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::depends(RcppArmadillo)]]

#include "util.h"
#include "header.h"

// -------------------------------------------------------------------------------
//' Get Bayesian Back-calculation Estimates and Model Diagnostics
//' 
//' This function performs Bayesian back-calculation, imputation of missing delays, and nowcasting based on the provided data.
//' 
//' @param outcome Vector of outcomes; difference between report and onset times
//' @param days Vector of days when the report is given, aligned from the minimum report day
//' @param week Vector indicating the week of the report, assumes no change within the week
//' @param weekend Binary vector indicating if the outcome was reported during a weekend
//' @param iter Number of iterations for the Bayesian back-calculation algorithm
//' @param sigma The standard deviation for the normal distribution
//' @param maxdelay The maximum delay parameter for the negative binomial distribution
//' @param si Serial interval vector
//' @param size The size parameter for the negative binomial distribution
//' @param workerID Identifier for the parallel worker
//' @param printProgress Flag to print the progress information
//' @param cd second size parameter, unused
//' @return output A list object that contains the back-calculated estimates and model diagnostics
//' @examples
//' \donttest{
//' data("sample_onset_dates")
//' data("sample_report_dates")
//' line_list <- create_linelist(sample_report_dates, sample_onset_dates)
//' sip <- si(14, 4.29, 1.18)
//' results <- run_backnow(
//'  line_list, 
//'   MAX_ITER = as.integer(2000), 
//'   norm_sigma = 0.5, 
//'   sip = sip,
//'   NB_maxdelay = as.integer(20), 
//'   NB_size = as.integer(6), 
//'   workerID = 1, 
//'   printProgress = 1, 
//'  preCalcTime = TRUE)
//' }
//' @export
// [[Rcpp::export]]
List backnow_cm(NumericVector outcome, NumericVector days, 
             IntegerVector week, IntegerVector weekend, 
             int iter, double sigma, int maxdelay, 
             NumericVector si, int size, 
             int workerID,
             int printProgress,
             Nullable<int> cd = R_NilValue){
  // Description:
  //
  //
  // Arguments:
  //   outcome:  d$delay, which = report - onset
  //   days:     d$report, the day that the report is given on, center on min(report), so 1, ...., max()
  //   week:     d$week of the report, allows you to assume that things don't change within this window
  //   weekend:  d$weekend, 0 or 1, indicator of weekend (1)
  //   iter:     the number of iterations to do in the big loop
  //   sigma:    set to 0.2, the standard deviation of the normal distribution used for parameter draws
  //   maxdelay: the maxmimum reporting delay, 20 days
  //   si:       the PDF of the serial interval by day, stands in for the generation time between infection and onset
  //   size:     ?? 6, seems like tau, the sliding window size in days
  //   workerID: workerID number
  //   printProgress: prints progress using worker ID
  //   cd:       ?? I believe its a placeholder for incorrect reporting delays?
  //              Note: cd must be earlier than the first day of nowcasting period
  //
  // Returns:
  //   Back:
  //   R:
  // --------------------------------
  // VARIABLE DEFINITIONS
  int cday;                                 // ??? unclear, a placeholder for the incorrect reporting delay?
  int type;                                 // ??? a flag for something, either 0 or 1
  int nc;                                   // ??? n parameters to define, 1 for every week + either 2 or 3
  int mi;                                   // ???
  int dv;                                   // ???
  int nm;                                   // ???

  double ratio;                             // ??? what is ratio
  double decision;                          // ??? what is decision
  double r1;                                // ???
  double m;                                 // ???

  NumericMatrix cov;                        // ??? 
  NumericVector par0;                       // create a placeholder for initial parameter estimates
  NumericVector oldpar;                     // initialize the loop oldpar vector
  NumericVector newpar;                     // initialize the loop newpar vector
  NumericVector mm;                         // ???
  NumericVector par1;                       // ???
  NumericVector mean;                       // ???
  NumericVector prob;                       // ???
  NumericVector p;                          // ???
  NumericVector s;                          // ???

  // Variables that depend on inputs
  NumericVector outcome1 = clone(outcome);           // ??? creates a copy of the outcome: why?
  int n  = outcome1.size();                          // the length of the outcome vector
  NumericVector dpind = rep(1.0, n);                 // ??? dispersion param changes? specific flag for each outcome
  int nw = max(week);                                // total number of weeks reporting data exist over
  int nd = max(days);                                // total number of days reporting data exist over
  IntegerVector x = seq(0, maxdelay);                // ???
  NumericVector dayseq = as<NumericVector>(x);       // ???
  NumericMatrix rt (iter, nd + maxdelay - size - 1); // the r(t) variable at the end, minus size??
  
  // -------------------------------- 
  // INITIALIZE VARIABLE SIZES & FLAGS
  // Ok so `nc` is the number of parameters that are being estimated
  // there are one for every week, and 3 extra if c is specified, 2 otherwise
  // perhaps this is the weekly value of the reporting delay, which doesn't change much within week
  if (cd.isNotNull()){      // if cd (which is ...???) is specified as input
    cday = as <int> (cd);   // make sure that its an integer
    nc   = nw + 3;          // set nc to nw + 3, which is + 1 for is_weekend then + 2 for TODO:... reasons???
    LogicalVector kk = (days > cday);  // for any report that has a day > reporting day max, set kk = T, otherwise kk = F
    dpind[kk] = 2;          // for any delay where kk=2, set dpind flag to 2    
    type      = 1;          // this type of model is called type 1
  } else {
    cday = -maxdelay;       // otherwise if cd is not specified, set it equal to -1 * maxdelay
    nc   = nw + 2;          // and nc now equals nw + 2, so +1 for is_weekend then +1 for ....?
    type = 0;               // this type of model is called type 0
  }

  // So now you can define parameter and back
  // matrix of Bayesian parameter estimates, and how they evolve over iterations
  NumericMatrix parameter (iter, nc);       
  // matrix of back-calculated counts; and how they evolve over iterations
  NumericMatrix back (iter, nd + maxdelay); 

  // --------------------------------
  // FILL IN MISSING
  // find all the missing delays (aka when report - onset = NA, because onset was missing)
  IntegerVector missind = findmiss(outcome1);  
  // the number of missing indices
  int nmiss = missind.size();                  
  // randomly sample from 1 to maxdelay, with replacement
  NumericVector miss0 = as<NumericVector>(sample(maxdelay, nmiss, true)); 
  // fill these values to delay vector
  outcome1[missind] = miss0; 
  
  // --------------------------------
  // INITIALIZE BAYESIAN PARAMETER ESTIMATES
  par0 = runif(nc);        // create a vector of `nc` uniform random variables between 0 and 1
  parameter(0,_) = par0;   // set the first row of the parameter estimation to be this runif() variable

  // --------------------------------
  // CREATE INDICATOR MATRIX AND MISSING MATRIX
  // creates a matrix of indicator values, first n columns are is_week_i, last column is is_weekend
  cov = dummy(week, weekend); 

  // Creates a separate indicator value matrix for the subset of days that are/were missing data
  NumericMatrix misscov (nmiss, cov.ncol());
  for (int i=0;i < nmiss; ++i){
    mi = missind[i];
    misscov(i,_) = cov(mi,_);
  }
  
  // --------------------------------
  std::string name_base = "./tmp/w" + std::to_string(workerID);
  std::string old_name = name_base + "-" + "0" + ".txt";
  std::string new_name = old_name;
  std::ofstream file(old_name.c_str()); // Open the file
  file.close();

  // --------------------------------  
  ////////////////////////////////////////////////////////////////////
  /// The MCMC loop //////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////
  //Rcpp::Rcout << "ITER# ";
  for (int i = 1; i < iter; ++i) { 
    
    // Print for every 500
    if(printProgress == 1) {
      if(i % 500 == 0) {
        new_name = name_base + "-" + std::to_string(i) + ".txt";
        //Rcpp::Rcout << i << "\t";
        if (std::rename(old_name.c_str(), new_name.c_str()) != 0) {
          // If std::rename returns a non-zero value, the rename operation failed
          Rcpp::stop("Failed to rename the file.");
        }
        old_name = new_name;
      }
    }

    // initialize the starting old parameters
    oldpar = parameter(i - 1, _);

    // --------------------------------
    // LOOP PART 1: The "Gibbs" Sampler, although this is not Gibbs, its metropolis-within-gibbs
    // This runs once for each parameter. Future iterations could look at randomizing this or making
    // some similarity beween 
    for (int k = 0; k < nc; ++k) {

      // for the variable that is i, ...., nw, 1 for each week and then 1 for is_weekend
      // THIS IS THE WHOLE REASON WE HAVE TO DO GIBBS?! otherwise we could just sample independently for each week
      if (k < nw + 1) {
        
        // make a copy of the old parameter values
        newpar    = clone(oldpar);

        // make a draw from a normal distribution with mean=oldpar[k] and sd=sigma
        newpar[k] = R::rnorm(oldpar[k], sigma);

        // compute the exp of the difference of the log-likelihoods? aka the ratio?
        ratio     = exp(logLikNB(outcome1, cov, dpind, newpar, maxdelay) - 
                        logLikNB(outcome1, cov, dpind, oldpar, maxdelay));

        // same as decision = ifelse(ratio > 1, 1, ratio)
        decision  = (ratio>1)?1:ratio;

        // ok so then this first gets a random uniform variable R::runif(0,1)
        // and if this is less than decision, make the parameter value = newpark[k],
        // but if its >= decision, keep the oldparameter value
        parameter(i,k) = (R::runif(0,1) < decision)?newpar[k]:oldpar[k];

        // updated oldpar[k]
        oldpar[k] = parameter(i, k);
      
      // for the extra 1 or 2 variables, which are ......
      } else {
        newpar    = clone(oldpar);
        newpar[k] = exp(R::rnorm(log(oldpar[k]), sigma));
        
        if (newpar[k] < 100) {
          ratio = (newpar[k]/oldpar[k])* 
            exp(logLikNB(outcome1, cov, dpind, newpar, maxdelay) - 
                logLikNB(outcome1, cov, dpind, oldpar, maxdelay));
        } else {
          ratio = 0;
        }
        
        decision = (ratio>1)?1:ratio;
        parameter(i,k) = (R::runif(0,1)<decision)?newpar[k]:oldpar[k];
        oldpar[k] = parameter(i,k);
      }
    }
    
    // --------------------------------
    // LOOP Part 2: Impute the missing delays
    // Rprintf(">> Part 2: Imputation of missing delays\n");
    // TODO: Why would you not just do this on the last iteration?

    // TYPE == 1 means there are two types of reporting delay sizes, aka r1 and r2
    if (type == 1) {
      par1      = oldpar[seq(0, nc - 3)];
      r1        = oldpar[nc - 2];
      double r2 = oldpar[nc - 1];
      mm        = get_mu_vec(misscov, par1);
      mean      = unique(mm);
      dv        = mean.size();
      //
      for (int k = 0 ; k < dv; ++k) {
        m = mean[k];
        LogicalVector select1 = (mm==m)&(dpind==1);
        LogicalVector select2 = (mm==m)&(dpind==2);
        int nm1 = sum(select1);
        int nm2 = sum(select2);

        if (nm1 > 0){
          prob = Rcpp::dnbinom_mu(dayseq, r1, m);
          p = prob / sum(prob);
          s = sample(dayseq, nm1, true, p); 
          IntegerVector ind = missind[select1];
          outcome1[ind] = s;
        }

        if (nm2 > 0){
          prob = Rcpp::dnbinom_mu(dayseq, r2, m);
          p = prob / sum(prob);
          s = sample(dayseq, nm2, true, p); 
          IntegerVector ind = missind[select2];
          outcome1[ind] = s;
        }
      }
    } else { // type != 1, is simpler, means there is just one size
      par1 = oldpar[seq(0, nc - 2)];
      r1   = oldpar[nc - 1];
      mm   = get_mu_vec(misscov, par1);
      mean = unique(mm);
      dv   = mean.size();
      for (int k = 0; k < dv; ++k){
        m = mean[k];
        LogicalVector select = (mm == m);
        nm = sum(select);
        prob = Rcpp::dnbinom_mu(dayseq, r1, m);
        p = prob / sum(prob);
        s = sample(dayseq, nm, true, p); 
        IntegerVector ind = missind[select];
        outcome1[ind] = s;
      }
    }
    
    // --------------------------------
    // LOOP Part 3: Back-calculation
    // seems like the upper limit of this is nd + maxdelay
    NumericVector backc (nd + maxdelay);     // initialize an empty vector for backc
    NumericVector backd = days - outcome1;   // initialize a vector for backd
    for (int j = 0; j <  (nd + maxdelay); ++j){
      backc[j] = sum(backd == (j - maxdelay + 1));
    }
    
    // --------------------------------
    // LOOP Part 4: Nowcasting
    NumericVector weights = prop(outcome1, backd, maxdelay, cday);
    NumericVector back1 = backc[seq(nd, nd + maxdelay - 1)];
    NumericVector check0 (back1.size());
    LogicalVector l = (back1 == 0);
    check0[l] = 1;
    NumericVector back2 = back1 + check0;  
    NumericVector trunc = mapply(back2, weights, rnb);
    NumericVector now   = back1 + trunc;
    NumericVector now1  = now - check0;
    LogicalVector check = (now1 < 0);
    now1[check] = 0;
    backc[seq(nd, nd + maxdelay - 1)] = now1; // set the rest of the backc vector = to now1
    back(i,_) = backc;
    rt(i,_) = getr(backc, si, size);
  } ///
  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////

  // Final print
  if(printProgress == 1) {
    //if(i % 500 == 0) {
      new_name = name_base + "-" + std::to_string(iter) + ".txt";
      //Rcpp::Rcout << i << "\t";
      if (std::rename(old_name.c_str(), new_name.c_str()) != 0) {
        // If std::rename returns a non-zero value, the rename operation failed
        Rcpp::stop("Failed to rename the file.");
      }
      old_name = new_name;
    //}
  }
  
  // --------------------------------
  // OUTPUT
  List output = List::create(Named("Back") = back, Named("R") = rt);
  return output;
}

