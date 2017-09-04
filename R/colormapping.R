#' Color-mapping of risk ratios to red-blue gradient
#' 
#' redblue
#' 
#' This function establishes the spread of reds and blues for the risk ratios to be mapped to. Higher risk ratios will be deeper red colors and lower risk ratios will be deeper blue colors.
#' @param x this will be the risk ratios shrunk to be on the scale of half risk to twice the risk as end points.
#' @return colors
redblue=function(x) { 
    y=colorRamp(RColorBrewer::brewer.pal(11,"RdBu")[11:1])(x); rgb(y[,1],y[,2],y[,3],maxColorValue=255) 
}


#' colormapping
#' 
#' This function establishes the spread of reds and blues for the risk ratios to be mapped to. Higher risk ratios will be deeper red colors and lower risk ratios will be deeper blue colors.
#' @param rates this will be the risk ratios shrunk to be on the scale of half risk to twice the risk as end points.
#' @param Time how many time periods are in the model? If this is only a spatial model, then time is set to 1
#' @param cv triggered if cross-validation model run on real data
#' @return returns vectors ofcolors for each time period, where risk ratios have been constrained to be between half risk and twice the risk
#' @export

colormapping <- function(rates,Time, cv) {
    if(!is.null(cv)){
        cl <- match.call()
        rate <- cl[[2]]
        #assign
        rcv <- eval(rate)[[1]]
        robs <- eval(rate)[[2]]
        
        color.obs <- sapply(1:Time, function(i) redblue(log(2*pmax(1/2,pmin(robs[,i],2)))/log(4)))
        color.cv <- sapply(1:Time, function(i) redblue(log(2*pmax(1/2,pmin(rcv[,i],2)))/log(4)))
        res <- list(colors.obs = color.obs, color.cv = color.cv)
    }
    else{
        #order: RRbic,RRaic, RRaicc, RRobs
        cl <- match.call()        
        rate <- cl[[2]]
        #assign
        rbic <- matrix(eval(rate)[[1]], ncol=Time)
        raic <- matrix(eval(rate)[[2]], ncol=Time)
        raicc <- matrix(eval(rate)[[3]], ncol=Time)
        robs <- matrix(eval(rate)[[4]], ncol=Time)
        print(str(rbic))
        #print(c(str(robs), str(rbic)))
        
        if(max(rbic)>2) {warning("Max riskratios from BIC greater than 2")}
        if(max(raic)>2) {warning("Max riskratios from AIC greater than 2")}
        if(max(raicc)>2) {warning("Max riskratios from AICc greater than 2")}
        color.obs <- sapply(1:Time, function(i) redblue(log(2*pmax(1/2,pmin(robs[,i],2)))/log(4)))
        color.qbic <- sapply(1:Time, function(i) redblue(log(2*pmax(1/2,pmin(rbic[,i],2)))/log(4)))
        color.qaic <- sapply(1:Time, function(i) redblue(log(2*pmax(1/2,pmin(raic[,i],2)))/log(4)))
        color.qaicc <- sapply(1:Time, function(i) redblue(log(2*pmax(1/2,pmin(raicc[,i],2)))/log(4)))
        res <- list(colors.obs = color.obs, color.qbic = color.qbic, color.qaic = color.qaic, color.qaicc = color.qaicc)
    }
    return(res) 
}

