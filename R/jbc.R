#' @title 
#' Japanese Breast Cancer Data (JBC)
#' @description 
#' Data on breast cancer incidence in prefects in Japan across 5 time periods. 
#' 
#' @docType data
#' 
#' 
#' 
#' @format An object of class "data.frame" with 1040 rows and 9 variables:
#' \itemize{
#'     \item id: identifier for each centroid location.
#'     \item period: time period at which each centroid location was measured (there are 5 total time periods, coded).
#'     \item death: count of the number of observed incident breast cancer deaths in each centroid-time period.
#'     \item expdeath: expected number of incident breast cancer deaths in each centroid-time period (age-standardized).
#'     \item covar1: simulated covariate 1 (will be unpenalized in model).
#'     \item covar2: simulated covariate 2 (will be unpenalized in model).
#'     \item covar3: simulated covariate 3 (will be unpenalized in model).
#'     \item covar4: simulated covariate 4 (will be unpenalized in model).
#'     \item covar5: simulated covariate 5 (will be unpenalized in model).
#' }
#' 
#' @keywords datasets
#' 
#' @examples 
#' data(jbc)
"jbc"