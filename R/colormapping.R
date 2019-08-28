#' @title
#' Color-mapping of risk ratios to red-blue gradient
#' @description 
#' This function establishes the spread of reds and blues for the risk ratios to be mapped to. Higher risk ratios will be deeper red colors and lower risk ratios will be deeper blue colors.
#' @param x this will be the risk ratios shrunk to be on the scale of half risk to twice the risk as end points.
#' @export
#' @return colors
redblue=function(x) { 
    y=colorRamp(RColorBrewer::brewer.pal(11,"RdBu")[11:1])(x); rgb(y[,1],y[,2],y[,3],maxColorValue=255) 
}


#' Color-mapping of probabilities to grey scale
#' 
#' greys
#' 
#' This function establishes the spread of grey scale for probabilities to be mapped to. Higher probabilities will be darker and lower probabilities will be lighter colors.
#' @param x vectors of probabilities to be mapped
#' @return colors
greys=function(x) { 
    y=colorRamp(RColorBrewer::brewer.pal(9,"Greys")[1:9])(x); rgb(y[,1],y[,2],y[,3],maxColorValue = 255) 
}

#' Color-mapping of red
#' 
#' reds
#' 
#' This function shows scale of red; used to create the legend for rr maps
#' @param x vector to be mapped
#' @return colors
reds=function(x) { 
    y=colorRamp(RColorBrewer::brewer.pal(9,"Reds")[1:9])(x); rgb(y[,1],y[,2],y[,3],maxColorValue = 255)
}

#' Color-mapping of blue
#' 
#' blues
#' 
#' This function shows scale of blue; used to create the legend for rr maps
#' @param x vector to be mapped
#' @return colors
blues=function(x) { 
    y=colorRamp(RColorBrewer::brewer.pal(9,"Blues")[1:9])(x); rgb(y[,1],y[,2],y[,3],maxColorValue = 255) 
}


#' colormapping
#' 
#' This function establishes the spread of reds and blues for the risk ratios to be mapped to. Higher risk ratios will be deeper red colors and lower risk ratios will be deeper blue colors.
#' @param rate this will be the risk ratios shrunk to be on the scale of half risk to twice the risk as end points.
#' @param Time how many time periods are in the model? If this is only a spatial model, then time is set to 1
#' @param cv triggered if cross-validation model run on real data
#' @param prob TRUE or FALSE depending on if you want to map probabilities or not (default is to map RR)
#' @return returns vectors ofcolors for each time period, where risk ratios have been constrained to be between half risk and twice the risk
#' @export

colormapping <- function(rate,Time,cv, prob) {
    if(!is.null(cv)){
        #assign
        rcv <- eval(rate)[[1]]
        robs <- eval(rate)[[2]]
        
        color.obs <- sapply(1:Time, function(i) redblue(log(2*pmax(1/2,pmin(robs[,i],2)))/log(4)))
        color.cv <- sapply(1:Time, function(i) redblue(log(2*pmax(1/2,pmin(rcv[,i],2)))/log(4)))
        res <- list(colors.obs = color.obs, color.cv = color.cv)
    }
    else{
        if(prob==TRUE){
            
            rbic <- matrix(rate[[1]], ncol=Time)
            raic <- matrix(rate[[2]], ncol=Time)
            raicc <- matrix(rate[[3]], ncol=Time)
            robs <- matrix(rate[[4]], ncol=Time)
            
            if(max(rbic)>1) {warning("Max probability from BIC greater than 1")}
            if(max(raic)>1) {warning("Max probability from AIC greater than 1")}
            if(max(raicc)>1) {warning("Max probability from AICc greater than 1")}
            color.obs <- sapply(1:Time, function(i) greys(robs[,i]))
            color.qbic <- sapply(1:Time, function(i) greys(rbic[,i]))
            color.qaic <- sapply(1:Time, function(i) greys(raic[,i]))
            color.qaicc <- sapply(1:Time, function(i) greys(raicc[,i]))
            res <- list(colors.obs = color.obs, color.qbic = color.qbic, color.qaic = color.qaic, color.qaicc = color.qaicc)
            
        }
        else{
            rbic <- matrix(rate[[1]], ncol=Time)
            raic <- matrix(rate[[2]], ncol=Time)
            raicc <- matrix(rate[[3]], ncol=Time)
            robs <- matrix(rate[[4]], ncol=Time)
            
            if(max(rbic)>2) {warning("Max riskratios from BIC greater than 2")}
            if(max(raic)>2) {warning("Max riskratios from AIC greater than 2")}
            if(max(raicc)>2) {warning("Max riskratios from AICc greater than 2")}
            color.obs <- sapply(1:Time, function(i) redblue(log(2*pmax(1/2,pmin(robs[,i],2)))/log(4)))
            color.qbic <- sapply(1:Time, function(i) redblue(log(2*pmax(1/2,pmin(rbic[,i],2)))/log(4)))
            color.qaic <- sapply(1:Time, function(i) redblue(log(2*pmax(1/2,pmin(raic[,i],2)))/log(4)))
            color.qaicc <- sapply(1:Time, function(i) redblue(log(2*pmax(1/2,pmin(raicc[,i],2)))/log(4)))
            # color.obs <- sapply(1:Time, function(i) redblue(log(1.5 * pmax((1/1.5), pmin(robs[, i], 1.5)))/log(2.25)))
            # color.qbic <- sapply(1:Time, function(i) redblue(log(1.5 * pmax((1/1.5), pmin(rbic[, i], 1.5)))/log(2.25)))
            # color.qaic <- sapply(1:Time, function(i) redblue(log(1.5 * pmax((1/1.5), pmin(raic[, i], 1.5)))/log(2.25)))
            # color.qaicc <- sapply(1:Time, function(i) redblue(log(1.5 * pmax((1/1.5), pmin(raicc[, i], 1.5)))/log(2.25)))
            # color.obs <- sapply(1:Time, function(i) redblue(log(1.2 * pmax((1/1.2), pmin(robs[, i], 1.2)))/log(1.44)))
            # color.qbic <- sapply(1:Time, function(i) redblue(log(1.2 * pmax((1/1.2), pmin(rbic[, i], 1.2)))/log(1.44)))
            # color.qaic <- sapply(1:Time, function(i) redblue(log(1.2 * pmax((1/1.2), pmin(raic[, i], 1.2)))/log(1.44)))
            # color.qaicc <- sapply(1:Time, function(i) redblue(log(1.2 * pmax((1/1.2), pmin(raicc[, i], 1.2)))/log(1.44)))
            res <- list(colors.obs = color.obs, color.qbic = color.qbic, color.qaic = color.qaic, color.qaicc = color.qaicc)
            
        }
        
    }
    return(res) 
}

