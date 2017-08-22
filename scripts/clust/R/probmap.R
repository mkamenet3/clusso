
#' Mapping colors to create a black/white probability map
#'
#'This function will create a probability map based on simulation data. In each simulation, it identifies where a cluster was selected,
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
#'@export
probmap <- function(lassoresult, vectors.sim, rr, nsim, Time, colormap=FALSE){
    prob.simBIC <- lapply(1:nsim, function(x) matrix(0, nrow(rr)*Time))
    prob.simAIC <- lapply(1:nsim, function(x) matrix(0, nrow(rr)*Time))
    prob.simAICc <- lapply(1:nsim, function(x) matrix(0, nrow(rr)*Time))
    indx <- which(rr !=1)
    rr.simBIC <- lapply(1:nsim, function(i) lassoresult$select_mu.qbic[[i]]/vectors.sim$E0_fit)
    rr.simAIC <- lapply(1:nsim, function(i) lassoresult$select_mu.qaic[[i]]/vectors.sim$E0_fit)
    rr.simAICc <- lapply(1:nsim, function(i) lassoresult$select_mu.qaicc[[i]]/vectors.sim$E0_fit)
    alphaBIC <- lapply(1:nsim, function(i) lapply(1:Time, function(k) 
        sort(table(matrix(rr.simBIC[[i]], ncol=Time)[,k]),decreasing=TRUE)[1]))
    alphaAIC <- lapply(1:nsim, function(i) lapply(1:Time, function(k) 
        sort(table(matrix(rr.simAIC[[i]], ncol=Time)[,k]),decreasing=TRUE)[1]))
    alphaAICc <- lapply(1:nsim, function(i) lapply(1:Time, function(k) 
        sort(table(matrix(rr.simAICc[[i]], ncol=Time)[,k]),decreasing=TRUE)[1]))
    
    for(j in 1:length(prob.simBIC)){
        for(i in 1:length(indx)){
            if (rr.simBIC[[j]][indx[i]] >= as.numeric(attributes(alphaBIC[[j]][[1]])[[1]]))  {
                prob.simBIC[[j]][indx[i]] <- 1
            }
            else {
                prob.simBIC[[j]][indx[i]] <- 0
            }
        }
    }
    for(j in 1:length(prob.simAIC)){
        for(i in 1:length(indx)){
            if (rr.simAIC[[j]][indx[i]] >= as.numeric(attributes(alphaAIC[[j]][[1]])[[1]]))  {
                prob.simAIC[[j]][indx[i]] <- 1
            }
            else {
                prob.simAIC[[j]][indx[i]] <- 0
            }
        }
    }
    for(j in 1:length(prob.simAICc)){
        for(i in 1:length(indx)){
            if (rr.simAICc[[j]][indx[i]] >= as.numeric(attributes(alphaAICc[[j]][[1]])[[1]]))  {
                prob.simAICc[[j]][indx[i]] <- 1
            }
            else {
                prob.simAICc[[j]][indx[i]] <- 0
            }
        }
    }
    prob <- as.vector(rr)
    prob[prob==1] <-0
    prob[prob!=0] <-1
    probBIC <- Reduce("+", prob.simBIC)/nsim
    probAIC <- Reduce("+", prob.simAIC)/nsim
    probAICc <- Reduce("+", prob.simAICc)/nsim
    if (colormap==TRUE){
        color.probmap <- sapply(1:Time, function(i) redblue(log(2*pmax(1/2,pmin(matrix(prob,ncol=Time)[,i],2)))/log(4)))    
        color.probmapBIC <- sapply(1:Time, function(i) redblue(log(2*pmax(1/2,pmin(matrix(probBIC,ncol=Time)[,i],2)))/log(4)))    
        color.probmapAIC <- sapply(1:Time, function(i) redblue(log(2*pmax(1/2,pmin(matrix(probAIC,ncol=Time)[,i],2)))/log(4)))    
        color.probmapAICc <- sapply(1:Time, function(i) redblue(log(2*pmax(1/2,pmin(matrix(probAICc,ncol=Time)[,i],2)))/log(4)))    
        return(list(prob = prob, probBIC = probBIC, probAIC = probAIC, probAICc = probAICc, 
                    color.probmap = color.probmap,
                    color.probmapBIC = color.probmapBIC,
                    color.probmapAIC = color.probmapAIC,
                    color.probmapAICc = color.probmapAICc))
    }
    else{
        return(list(prob=prob, probBIC = probBIC,probAIC = probAIC, probAICc = probAICc ))
    }
}    