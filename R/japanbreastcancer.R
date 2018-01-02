#' @title 
#' Japanese Breast Cancer Data
#' @description 
#' Data on breast cancer incidence in prefects in Japan across 5 time periods. UTM coordinates and polygon data included.
#' 
#' @docType data
#' 
#' @usage data(japanbreastcancer)
#' 
#' 
#' @format An object of class "data.frame". Includes four separate dataframes: 
#' \itemize{
#'  \item 1) dframe.poly2 - coordinates for plotting polygons
#'  \item 2) dframe.prefect2 - coordinates for plotting prefect boundaries
#'  \item 3) japanbreastcancer - data on observed death from breast cancer and expected death by polygon-time period.
#'  \item 4) utmJapan - UTM coordinates
#' }
#' 
#' @keywords datasets
#' 
#' @aliases dframe.poly2 dframe.prefect2 utmJapan
#' 
#' @examples 
#' data(japanbreastcancer)
"japanbreastcancer"