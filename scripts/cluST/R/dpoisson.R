#' Poisson Distribution Function
#' 
#' This is the main distriution function for our model. This assumes we have a Poisson fixed effect and Gamma random effect. In order to deal with constraints from the Lasso function, we use the Poisson distirbution function here and account for overdispersion in the QIC.
#' @param y observed values
#' @param lambda vector of expected outcomes * exp(each column of each potential path)
#' @param log whether or not the log-likelihood should be returned or the likelihood. Default is to be TRUE
#' @return returns a matrix 

# dpoisson <- function(y, lambda, log = FALSE) {
#     if(log == FALSE) 
#         return(lambda^y * exp(-lambda)/factorial(y))
#     else
#         return(y*ifelse(lambda==0,1,log(lambda))-lambda)
# }



dpoisson <- function(x, lambda, log = FALSE) {
    if(log == FALSE) 
        return(lambda^x * exp(-lambda)/factorial(x))
    else
        return(x*ifelse(lambda==0,1,log(lambda))-lambda-log(factorial(x)))
}