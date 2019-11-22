#' Binomial Function
#' 
#' @title dbinom
#' @param linpred
#' @return Returns vector of Binomial log-likelihood for each proposed path of the LASSO tuning parameter values.
#' 
#' @example 
#' set.seed(1)
#' 
dbin <- function(x, n, phat){
    loglik_i <- x*log(phat) + (n-x)*log(1-phat) 
    return(sum(loglik_i))
}

invlogit <- function(linpred){
    phat<- exp(linpred)/(1+exp(linpred))
    return(phat)
}

# logL <- function(Yx, Ex,p){
#     sum(dbinom(Yx, Ex, p,log=TRUE))
# }
# # 
# # pseq <- seq(0.01, 0.99, 0.01)
# # logL(pseq)
# # 
# logL(0.5)
# logL(0.4)
# # test <- sapply(1:341, function(i) logL(phat[,i]))
