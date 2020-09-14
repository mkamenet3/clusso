#' Bernoullii Function
#' @title dbern
#' @param Yx Cases (0 or 1).
#' @param p Vector of estimated probabilities for each individual.
#' @return Returns vector of Bernoulli log-likelihood for each proposed path of the LASSO tuning parameter values.
dbern <- function(Yx,p){
    sum(dbinom(Yx, size=1, p,log=TRUE))
}
#' Convert linear predictor to probability scale
#' @title pihat
#' @param xbeta Sparse matrix of potential space-time clusters times coefficient for each lasso path.
#' @return Returns matrix of estimated probabilities for each space-time location by each lasso path.
pihat <- function(xbeta){
    phat <- 1/(1+exp(-xbeta))
    return(phat)
}
