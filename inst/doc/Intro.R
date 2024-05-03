## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(linelistBayes)

## ----example1-----------------------------------------------------------------
data("sample_dates")
data("sample_location")
data("sample_cases")

head(sample_dates)
head(sample_cases)
head(sample_location)

## ----casecounts---------------------------------------------------------------
caseCounts <- create_caseCounts(date_vec = sample_dates,
                                location_vec = sample_location,
                                cases_vec = sample_cases)
head(caseCounts)

## ----dpi=100, fig.height=2.75, fig.width=6.75,dev='png'-----------------------
plot(caseCounts)

## ----gethead,dpi=100, fig.height=2.75, fig.width=6.75,dev='png'---------------
caseCounts <- caseCounts[1:80, ]
plot(caseCounts)

## ----serial,dpi=100, fig.height=2.75, fig.width=6.75,dev='png'----------------
sip <- si(14, 4.29, 1.18)
plot(sip, type = 'l')

## ----backcalc1, eval=FALSE----------------------------------------------------
#  out_list_demo <- run_backnow(caseCounts,
#                          MAX_ITER = as.integer(2000),
#                          norm_sigma = 0.2,
#                          sip = sip,
#                          NB_maxdelay = as.integer(20),
#                          NB_size = as.integer(6),
#                          printProgress = 1,
#                          reportF_missP = 0.6)

## ----fakeplot, eval=FALSE-----------------------------------------------------
#  plot(out_list_demo, "est")

## ----plot1, dpi=100, fig.height=2.75, fig.width=6.75,dev='png',echo=FALSE-----
data("out_list_demo")
# plot
oldpar <- par(mfrow = c(1,2))
par(oma = c(0, 0, 0, 0), mar = c(4, 4, 1, 1))
plot(out_list_demo, "est")
par(oldpar)

## ----fakeplot2, eval=FALSE----------------------------------------------------
#  plot(out_list_demo, "rt")

## ----plot2, dpi=100, fig.height=2.75, fig.width=6.75,dev='png',echo=FALSE-----
data("out_list_demo")
# plot
oldpar <- par(mfrow = c(1,2))
par(oma = c(0, 0, 0, 0), mar = c(4, 4, 1, 1))
plot(out_list_demo, "rt")
par(oldpar)

## -----------------------------------------------------------------------------
data("sample_report_dates")
data("sample_onset_dates")

## ----casecounts2--------------------------------------------------------------
caseCounts_line <- create_linelist(report_dates = sample_report_dates,
                                onset_dates = sample_onset_dates)
head(caseCounts_line)

## ----backcalc2, eval=FALSE----------------------------------------------------
#  out_list_demo <- run_backnow(caseCounts_line,
#                          MAX_ITER = as.integer(2000),
#                          norm_sigma = 0.2,
#                          sip = sip,
#                          NB_maxdelay = as.integer(20),
#                          NB_size = as.integer(6),
#                          printProgress = 1)

## -----------------------------------------------------------------------------
my_linelist <- convert_to_linelist(caseCounts, 
                                   reportF = rnbinom, 
                                   reportF_args = list(size = 3, mu = 9),
                                   reportF_missP = 0.6)
head(my_linelist)

