#' Sample datasets
#'
#' These data sets provide tests for use in converting between case counts
#' and line list data.
#'
#' `sample_cases` provides a vector of daily case counts.
#' `sample_dates` are the dates of the sample case counts.
#' `sample_location` are the locations of the sample case counts.
#' `sample_onset_dates` are the same information as in sample_cases, but 
#' with one entry per case indicating the date of symptom onset.
#' `sample_report_dates` are the same information as in sample_cases, but 
#' with one entry per case indicating the date of symptom reporting.
#' `out_list_demo` is a precomputed output from run_backnow, useful for
#' plotting in the vignettes
#'
#' @format Either vectors or a list object (out_list_demo)
#' @examples
#' sample_cases
#' sample_dates
#' sample_location
#' sample_onset_dates
#' sample_report_dates
"sample_cases"

#' @rdname sample_cases
#' @format NULL
"sample_dates"

#' @rdname sample_cases
#' @format NULL
"sample_location"

#' @rdname sample_cases
#' @format NULL
"sample_onset_dates"

#' @rdname sample_cases
#' @format NULL
"sample_report_dates"

#' @rdname sample_cases
#' @format NULL
"out_list_demo"