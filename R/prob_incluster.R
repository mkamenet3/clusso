
#' @title
#' get_prob
#' @description 
#' Finds the probability of being in the cluster for BIC, AIC, and AICc based on the expected risk ratio
#' @param lassoresult result of either space-time lasso simulation or space-only simulation
#' @param ncentroids number of centroids or centers
#' @param Time number of time periods
#' @return returns list of probabilities.
#' 
get_prob <- function(lassoresult, spacetime,init, E1, ncentroids, Time){
    prob.bic <- probs(lassoresult$E.qbic, lassoresult$probs.qbic)
    prob.aic <- probs(lassoresult$E.qaic, lassoresult$probs.qaic)
    prob.aicc <- probs(lassoresult$E.qaicc, lassoresult$probs.qaicc)
    if(spacetime == FALSE){
        obs <-as.vector(matrix(as.vector(E1)/as.vector(init$E0),ncol=Time))
    }
    else{
         obs <- matrix(as.vector(E1)/as.vector(init$E0),ncol=Time)
    }
    prob.obs <- ifelse(obs>1,1,0)
    res <- list(prob.bic = prob.bic, prob.aic = prob.aic, prob.aicc = prob.aicc, prob.obs = prob.obs)
    return(res)
}

#' @title 
#' prob_incluster
#' @description
#'Mapping colors to create a black/white probability map. This function will create a probability map based on simulation data. In each simulation, it identifies where a cluster was selected,
#'compared to the background rate. It then average over the number of simulations, giving us a matrix which ranges from 0 to 1 in probability.
#'To map this probabilities into a color scheme, please see the $colormapping$ function and select probmap=TRUE. TODO integrate all of this
#'into a workflow and extend to observed data, not only simulated data.
#'@param lassoresult List of QBIC, QAIC, QAICc estimates from the mylasso.sim function
#'@param vectors.sim dataframe of initial vectors of the observed and expected counts that went into mylasso.sim function
#'@param rr risk ratio matrix that was used in the simulation
#'@param nsim number of simulations
#'@param Time number of time period
#'@param colormap default is FALSE. Signals whether the probabilities should directly be mapped to the red-blue color scheme. If this is false,
#'only the probability values will be returned. If true, then the probability values and the mapped colors will be returned in a list.
#'@return returns vector which calculated the number of time the cluster was correctly identified out of the simulations
prob_incluster <- function(lassoresult_selectmu, ncentroids, Time, nsim){
    # if(Time == 1){
    #     print("space-only")
    #     
    # }
    vec <- rep(0, ncentroids * Time)
    position <- list(vec)[rep(1, nsim)]
    ix <- lapply(1:nsim, function(x) which(lassoresult_selectmu[[x]] >1))
    
    #quick function to recode
    reval <- function(probs, ix){
        probs[ix] <-1
        return(probs)
    }
    simindicator <- mapply(reval, position, ix)
    probs <- Matrix::rowSums(simindicator)/nsim
    
    return(probs)
}    


#' @title
#' probs
#' @description 
#' Helper function
#' @param lassoresult_E
#' @param lassoresult_probs
#' @return returns vector of probabilities mapped of being included in the cluster based on IC or cv.

probs <- function(lassoresult_E, lassoresult_probs){
    out <- log(as.vector(lassoresult_E))/lassoresult_probs
    ix <- which(is.na(out) | is.infinite(out) | is.nan(out) | (out < 0))
    out[ix] <- 0
    return(out)
}



