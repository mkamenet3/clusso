#' Poisson Function
#' 
#' Poisson function. 
#' @param y observed values
#' @param lambda vector of exp(each column of each potential path). Assumes you exponentiate prior to using this function
#' @param E0 expected counts (the offset). If in simulation, then this should be the E0 counts for each simulation. These E0 have been standardized such that E0 = E*(sum(y)/sum(E)) with the scale_sim()/scale() function depending on if a simulation is being run or not
#' @return returns a a vector (or a list of vectors) of the Poisson log-likelihood for each proposed path of the Lasso tuning parameter values 
#' @export
#' 
#' @examples 
#' set.seed(1)
#' E <- MASS::rnegbin(n = 20,mu = 15,theta = 1000)
#' theta = 1000
#' y <- MASS::rnegbin(E, theta=theta)
#' E0 <- E*(sum(y)/sum(E))
#' lambda <- exp(c(rep(0,15),rep(0.4,5)))
#' dpoisson(y, lambda, E0)

dpoisson <- function(y, lambda, E0) {
    if (any(lambda == 0)) stop("Element of lambda is zero - log of zero will return -Inf. Make sure you exponentiated already")
    if (any(E0 == 0)) warning("At least one element of expected E0 is zero")
    loglik_i <- y*log(lambda*E0) - (lambda*E0)
    return(sum(loglik_i))
}

