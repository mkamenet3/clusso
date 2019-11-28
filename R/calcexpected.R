#'Estimate expected counts if data not available
#'
#'@title
#'calcexpected
#' @description
#' Expected counts can be estimated if not available, following selected distribution function
#'@param observed Vector of observed counts.
#'@param population Vector of population sizes in each geographic unit.
#'@param periods Vector of time periods matching observed counts. 
#'@param ids Vector of centroid identifiers matching the observed counts.
#'@export
#'@return Returns a dataframe of geographic ids, time periods, population counts, observed counts, and expected counts.
#'@examples
#'\donttest{
#'data(japanbreastcancer)
#'calcexpected(japanbreastcancer$observed, japanbreastcancer$period, 
#'    japanbreastcancer$id, family="poisson") }

calcexpected <- function(observed, population,periods, ids){
    if(is.character(periods)){
        periods <- as.factor(periods)
    }
    if(is.character(ids)){
        ids <- as.factor(ids)
    }
    rate <- sum(observed)/sum(population)
    expected <- population*rate
    dat <- cbind.data.frame(ids, periods, population, observed, expected)
    return(dat)
}
