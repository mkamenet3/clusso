#' Color-mapping of risk ratios to red-blue gradient
#' 
#' redblue
#' 
#' This function establishes the spread of reds and blues for the risk ratios to be mapped to. Higher risk ratios will be deeper red colors and lower risk ratios will be deeper blue colors.
#' @param x this will be the risk ratios shrunk to be on the scale of half risk to twice the risk as end points.
#' @return colors
redblue=function(x,...) { 
    y=colorRamp(brewer.pal(11,"RdBu")[11:1])(x); rgb(y[,1],y[,2],y[,3],max=255) 
}


#' colormapping
#' 
#' This function establishes the spread of reds and blues for the risk ratios to be mapped to. Higher risk ratios will be deeper red colors and lower risk ratios will be deeper blue colors.
#' @param riskratios this will be the risk ratios shrunk to be on the scale of half risk to twice the risk as end points.
#' @param Time how many time periods are in the model? If this is only a spatial model, then time is set to 1
#' @return returns vectors ofcolors for each time period, where risk ratios have been constrained to be between half risk and twice the risk
#' @export
#'@examples
#'set.seed(2)
#'riskratios <- list(RRobs = rnorm(20, mean = 1, sd = 0.10)
#'RRbic <- c(rep(1,18), c(2.1,2.2))
#'RRaic <- c(rep(1,15), rnorm(5, 2, 0.05))
#'RRaicc <- c(rep(1,15), rnorm(5, 2, 0.05))
#'Time = 2
#'colormapping(riskratios,Time)

colormapping <- function(riskratios,Time,...) {
    if(max(riskratios$RRbic)>2) {warning("Max riskratios from BIC greater than 2")}
    if(max(riskratios$RRaic)>2) {warning("Max riskratios from AIC greater than 2")}
    if(max(riskratios$RRaicc)>2) {warning("Max riskratios from AICc greater than 2")}
    color.obs <- sapply(1:Time, function(i) redblue(log(2*pmax(1/2,pmin(riskratios$RRobs[,i],2)))/log(4)))
    color.qbic <- sapply(1:Time, function(i) redblue(log(2*pmax(1/2,pmin(riskratios$RRbic[,i],2)))/log(4))) 
    color.qaic <- sapply(1:Time, function(i) redblue(log(2*pmax(1/2,pmin(riskratios$RRaic[,i],2)))/log(4)))
    color.qaicc <- sapply(1:Time, function(i) redblue(log(2*pmax(1/2,pmin(riskratios$RRaicc[,i],2)))/log(4)))
    return(list(colors.obs = color.obs, color.qbic = color.qbic, color.qaic = color.qaic, color.qaicc = color.qaicc)) 
}
