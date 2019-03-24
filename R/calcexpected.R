#'Estimate expected counts if data not available
#'
#'@title
#'calcexpected
#' @description
#' Expected counts can be estimated if not available, following selected distribution function
#'@param observed Vector of observed counts
#'@param periods Cector of time periods matching observed counts 
#'@param ids Vector of centroid identifiers matching the observed counts.
#'@param family Must be provided. Can either be specified as \code{poisson} or \code{binomial}.
#'@export
#'@return returns vector of expected counts
#'@examples
#'\donttest{
#'data(japanbreastcancer)
#'calcexpected(japanbreastcancer$observed, japanbreastcancer$period, japanbreastcancer$id, family="poisson") }

calcexpected <- function(observed, periods, ids,family=c("poisson","binomial")){
    if(missing(family)){
        stop("You must specify family to be either 'poisson' or 'binomial'")
    }
    if(is.character(periods)){
        periods <- as.factor(periods)
    }
    if(is.character(ids)){
        ids <- as.factor(ids)
    }
    family <- match.arg(family, several.ok = FALSE)
    switch(family, 
           poisson=fitted(glm(observed ~ ids + periods, family="poisson")),
           binomial=fitted(glm(observed ~ ids + periods, family="binomial")))
}
