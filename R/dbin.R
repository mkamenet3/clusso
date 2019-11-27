#' Binomial Function
#' @title dbinom
#' @param Yx Number of cases.
#' @param Ex Total number of trials (or cases + controls).
#' @param p Vector of estimated probabilities for each space-time location.
#' @return Returns vector of Binomial log-likelihood for each proposed path of the LASSO tuning parameter values.
#' @example 
#' 
dbin <- function(Yx, Ex,p){
    sum(dbinom(Yx, Ex, p,log=TRUE))
}


#' Convert linear predictor to probability scale
#' @title pihat
#' @param xbeta Sparse matrix of potential space-time clusters times coefficient for each lasso path.
#' @return Returns matrix of estimated probabilities for each space-time location by each lasso path.
pihat <- function(xbeta){
    phat <- 1/(1+exp(-xbeta))
    return(phat)
}
