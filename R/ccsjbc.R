#' @title 
#' Japanese Breast Cancer Data (JBC) - Case-Control Data Example
#' @description 
#' Data on breast cancer cases and cases+control in prefects in Japan across 5 time periods. 
#' 
#' @docType data
#' 
#' 
#' 
#' @format An object of class "data.frame" with 1040 rows and 10 variables:
#' \itemize{
#'     \item id: identifier for each centroid location.
#'     \item period: time period at which each centroid location was measured (there are 5 total time periods, coded).
#'     \item numcases: count of the number of breast cancer cases in each centroid-time period.
#'     \item n: total number of cases+controls (number of trials) in each centroid-time period.
#'     \item covar1: simulated covariate 1 (will be unpenalized in model).
#'     \item covar2: simulated covariate 2 (will be unpenalized in model).
#'     \item covar3: simulated covariate 3 (will be unpenalized in model).
#'     \item covar4: simulated covariate 4 (will be unpenalized in model).
#'     \item covar5: simulated covariate 5 (will be unpenalized in model).
#'     \item prop: proportion of cases/(cases+controls) (used in exploratory analysis).
#' }
#' 
#' @keywords datasets
#' 
#' @examples 
#' data(ccsjbc)
"ccsjbc"