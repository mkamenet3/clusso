#' @title
#' get_prob
#' @description 
#' Finds the probability of being in the cluster for BIC, AIC, and AICc based on the expected risk ratio
#' @param lassoresult result of either space-time lasso simulation or space-only simulation
#' @param init list of initial vector values
#' @param E1 standardized expected counts
#' @param ncentroids number of centroids or centers
#' @param Time number of time periods
#' @param nsim number of simulations
#' @param threshold vector of two threshold values TODO: allow flexibility for number of threshold
#' @return returns list of probabilities.
#' 
get_prob <- function(lassoresult,init, E1, ncentroids, Time, nsim, threshold){
    #nsim <- nsim
    prob.bic <- prob_incluster(lassoresult$select_mu.qbic, ncentroids, Time, nsim)
    prob.aic <- prob_incluster(lassoresult$select_mu.qaic, ncentroids, Time, nsim)
    prob.aicc <- prob_incluster(lassoresult$select_mu.qaicc, ncentroids, Time, nsim)
    obs <- matrix(as.vector(E1)/as.vector(init$E0),ncol=Time)
    prob.obs <- ifelse(obs>1,1,0)
    #thresholding
    ##BIC
    bic_thresh1 <- ifelse(prob.bic>threshold[1],1,0)
    attr(bic_thresh1, 'thresh') <- threshold[1]
    bic_thresh2 <- ifelse(prob.bic>threshold[2],1,0)
    attr(bic_thresh2, 'thresh') <- threshold[2]
    
    ##AIC
    aic_thresh1 <- ifelse(prob.aic>threshold[1],1,0)
    attr(aic_thresh1, 'thresh') <- threshold[1]
    aic_thresh2 <- ifelse(prob.aic>threshold[2],1,0)
    attr(aic_thresh2, 'thresh') <- threshold[2]
    
    ##AICc
    aicc_thresh1 <- ifelse(prob.aicc>threshold[1],1,0)
    attr(aicc_thresh1, 'thresh') <- threshold[1]
    aicc_thresh2 <- ifelse(prob.aicc>threshold[2],1,0)
    attr(aicc_thresh2, 'thresh') <- threshold[2]
    
    probs <- list(prob.bic = prob.bic, prob.aic = prob.aic, prob.aicc = prob.aicc, prob.obs = prob.obs)
    probs.thresh1 <- list( prob.bic.thresh1 = bic_thresh1,  prob.aic.thresh1 = aic_thresh1, 
                          prob.aicc.thresh1 = aicc_thresh1, prob.obs = prob.obs) 
    probs.thresh2  <- list(prob.bic.thresh2 = bic_thresh2, prob.aic.thresh2 = aic_thresh2,
                           prob.aicc.thresh2 = aicc_thresh2, prob.obs = prob.obs)                         
                          
    
    res <- list(probs = probs, probs.thresh1 = probs.thresh1, probs.thresh2 = probs.thresh2)
    return(res)
}

#' @title 
#' prob_incluster
#' @description
#'Mapping colors to create a  probability map. This function will create a probability map based on simulation data. In each simulation, it identifies where a cluster was selected,
#'compared to the background rate. It then average over the number of simulations, giving us a matrix which ranges from 0 to 1 in probability.
#'To map this probabilities into a color scheme, please see the $colormapping$ function and select probmap=TRUE. TODO integrate all of this
#'into a workflow and extend to observed data, not only simulated data.
#'@param select_mu List of selected mu vectors selected by the respective information criteria.
#'@param ncentroids number of centroids
#'@param Time number of time period
#'@param nsim number of simulations
#'@return returns vector which calculated the number of time the cluster was correctly identified out of the simulations
prob_incluster <- function(select_mu, ncentroids, Time, nsim){
    vec <- rep(0, ncentroids * Time)
    position <- list(vec)[rep(1, nsim)]
    bgRate = lapply(1:nsim, function(i) as.numeric(names(which.max(table(select_mu[[i]])))))
    ix <- lapply(1:nsim, function(i) which(abs(select_mu[[i]] - bgRate[[i]])>=10^-6))
    #apply(rate, 1, function(x) sum(abs(x-bgRate)>=10^-6)/ncol(rate))
    
    
    #ix <- lapply(1:nsim, function(i) which(as.vector(select_mu[[i]]) > 1))
    #quick function to recode
    reval <- function(probs, ix){
        probs[ix] <-1
        return(probs)
    }
    simindicator <- mapply(reval, position, ix)
    probs <- Matrix::rowSums(simindicator)/nsim
    return(probs)
}    




