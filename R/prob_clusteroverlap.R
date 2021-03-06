#' @title 
#' prob_clusteroverlap
#' @description
#' Finds the probability of any overlap with true cluster based on BIC, AIC, and AICc based on the expected risk ratio
#'@param sparseMAT large sparse matrix created in \code{clusso_sim} function.
#'@param lassoresult List of QBIC, QAIC, QAICc estimates from the mylasso.sim function
#'@param rr risk ratio matrix that was used in the simulation
#'@param risk.ratio Risk ratio that was set for cluster in simulation
#'@param x x-coordinates
#'@param y y-coordinates
#'@param rMax Maximum radius for threshold in simulation
#'@param nsim number of simulations
#'@param Time number of time period
#'@param thresh thresholds for simulation detection (TODO)
#'@param ncentroids number of centroids
#'@param nullmod Is model run the the null model with no cluster to detect? Default is that model has a cluster to detect (\code{nullmod=NULL})
#'@return List of detection probabilities (in vs. out of cluster) by criterion. Returns vector which calculated the number of time the cluster was correctly identified out of the simulations.
prob_clusteroverlap <- function(sparseMAT,lassoresult,rr, risk.ratio,x,y,rMax,nsim,Time, thresh,ncentroids, nullmod=NULL){
    #DEFINE TRUTH
    if(risk.ratio==1){
        warning("Risk.ratio was set to 1")
        rrmatvec <- rep(0,length(rr))
    }
    else{
        rrmatvec <- ifelse(as.vector(rr)==risk.ratio,1,0)    
    }
    if(!missing(nullmod)){
        nullmod = TRUE
    }
    else{
        nullmod = FALSE
    }
    if(!is.null(thresh)){
        message("Threshold detection function is still in development. Stay tuned!")
    }
    ##(Q)BIC
    selected <- lassoresult$selections$select.qbic
    #Take out the time vectors - only keep cluster part of matrix
    sparseMAT_clusteronly <- sparseMAT[,-c(ncol(sparseMAT)-Time+1:ncol(sparseMAT))]
    #select out my betas for each sim
    betaselect <- lapply(1:nsim, function(i) lassoresult$coefs.lasso.all[[i]][,selected[[i]]])
    #binarize
    betaselect_bin <- lapply(1:nsim, function(i) ifelse(abs(betaselect[[i]] >= 10e-6),1,0))
    #only take the clusters betas
    betaselect_bin_clusteronly <- lapply(1:nsim, function(i) betaselect_bin[[i]][-c(ncol(sparseMAT)-Time+1:ncol(sparseMAT))])
    ##INCLUSTER
    if (nullmod == TRUE){
        in.qbic <- lapply(1:nsim, function(i) any(betaselect_bin_clusteronly[[i]]!=0))
        inperc.qbic <- sum(unlist(in.qbic))/nsim    
    }
    else{
        clusteroverlap <- t(rrmatvec) %*% sparseMAT_clusteronly
        clusteroverlap_bin <- ifelse(clusteroverlap !=0,1,0)
        incluster_sim <- lapply(1:nsim, function(i) clusteroverlap_bin %*% betaselect_bin_clusteronly[[i]])
        incluster_sim_bin <- lapply(1:nsim, function(i) ifelse(incluster_sim[[i]] !=0,1,0))
        inperc.qbic <- paste0((sum(unlist(incluster_sim_bin))/nsim)*100, "%")
    }
    ##OUTCLUSTER
    if (nullmod == TRUE){
        outperc.qbic <- NA
    }
    else{
        clusteroverlap_bin <- ifelse(clusteroverlap !=0,0,1)
        outcluster_sim <- lapply(1:nsim, function(i) clusteroverlap_bin %*% betaselect_bin_clusteronly[[i]])
        outcluster_sim_bin <- lapply(1:nsim, function(i) ifelse(outcluster_sim[[i]] !=0,1,0))
        outperc.qbic <- paste0((sum(unlist(outcluster_sim_bin))/nsim)*100, "%")
    }
    ##(Q)AIC
    selected <- lassoresult$selections$select.qaic
    #Take out the time vectors - only keep cluster part of matrix
    sparseMAT_clusteronly <- sparseMAT[,-c(ncol(sparseMAT)-Time+1:ncol(sparseMAT))]
    #select out my betas for each sim
    betaselect <- lapply(1:nsim, function(i) lassoresult$coefs.lasso.all[[i]][,selected[[i]]])
    #binarize
    betaselect_bin <- lapply(1:nsim, function(i) ifelse(abs(betaselect[[i]] >= 10e-6),1,0))
    #only take the clusters betas
    betaselect_bin_clusteronly <- lapply(1:nsim, function(i) betaselect_bin[[i]][-c(ncol(sparseMAT)-Time+1:ncol(sparseMAT))])
    ##INCLUSTER
    if (nullmod == TRUE){
        in.qbic <- lapply(1:nsim, function(i) any(betaselect_bin_clusteronly[[i]]!=0))
        inperc.qaic <- sum(unlist(in.qbic))/nsim    
    }
    else{
        clusteroverlap <- t(rrmatvec) %*% sparseMAT_clusteronly
        clusteroverlap_bin <- ifelse(clusteroverlap !=0,1,0)
        incluster_sim <- lapply(1:nsim, function(i) clusteroverlap_bin %*% betaselect_bin_clusteronly[[i]])
        incluster_sim_bin <- lapply(1:nsim, function(i) ifelse(incluster_sim[[i]] !=0,1,0))
        inperc.qaic <- paste0((sum(unlist(incluster_sim_bin))/nsim)*100, "%")
    }
    ##OUTCLUSTER
    if (nullmod == TRUE){
        outperc.qaic <- NA
    }
    else{
        clusteroverlap_bin <- ifelse(clusteroverlap !=0,0,1)
        outcluster_sim <- lapply(1:nsim, function(i) clusteroverlap_bin %*% betaselect_bin_clusteronly[[i]])
        outcluster_sim_bin <- lapply(1:nsim, function(i) ifelse(outcluster_sim[[i]] !=0,1,0))
        outperc.qaic <- paste0((sum(unlist(outcluster_sim_bin))/nsim)*100, "%")
    }
    #(Q)AICc
    selected <- lassoresult$select.qaicc
    #Take out the time vectors - only keep cluster part of matrix
    sparseMAT_clusteronly <- sparseMAT[,-c(ncol(sparseMAT)-Time+1:ncol(sparseMAT))]
    #select out my betas for each sim
    betaselect <- lapply(1:nsim, function(i) lassoresult$coefs.lasso.all[[i]][,selected[[i]]])
    #binarize
    betaselect_bin <- lapply(1:nsim, function(i) ifelse(abs(betaselect[[i]] >= 10e-6),1,0))
    #only take the clusters betas
    betaselect_bin_clusteronly <- lapply(1:nsim, function(i) betaselect_bin[[i]][-c(ncol(sparseMAT)-Time+1:ncol(sparseMAT))])
    ##INCLUSTER
    if (nullmod == TRUE){
        in.qbic <- lapply(1:nsim, function(i) any(betaselect_bin_clusteronly[[i]]!=0))
        inperc.qaicc <- sum(unlist(in.qbic))/nsim    
    }
    else{
        clusteroverlap <- t(rrmatvec) %*% sparseMAT_clusteronly
        clusteroverlap_bin <- ifelse(clusteroverlap !=0,1,0)
        incluster_sim <- lapply(1:nsim, function(i) clusteroverlap_bin %*% betaselect_bin_clusteronly[[i]])
        incluster_sim_bin <- lapply(1:nsim, function(i) ifelse(incluster_sim[[i]] !=0,1,0))
        inperc.qaicc <- paste0((sum(unlist(incluster_sim_bin))/nsim)*100, "%")
    }
    ##OUTCLUSTER
    if (nullmod == TRUE){
        outperc.qaicc <- NA
    }
    else{
        clusteroverlap_bin <- ifelse(clusteroverlap !=0,0,1)
        outcluster_sim <- lapply(1:nsim, function(i) clusteroverlap_bin %*% betaselect_bin_clusteronly[[i]])
        outcluster_sim_bin <- lapply(1:nsim, function(i) ifelse(outcluster_sim[[i]] !=0,1,0))
        outperc.qaicc <- paste0((sum(unlist(outcluster_sim_bin))/nsim)*100, "%")
    }
    
    
    return(list(outperc.qbic = outperc.qbic, outperc.qaic = outperc.qaic, outperc.qaicc = outperc.qaicc,
                inperc.qbic = inperc.qbic, inperc.qaic = inperc.qaic, inperc.qaicc = inperc.qaicc))
    
}
    
    
#' @title
#'get_prob
#'@description 
#'Finds the probability of being in the cluster for BIC, AIC, and AICc based on the expected risk ratio
#'@param lassoresult result of either space-time lasso simulation or space-only simulation
#'@param init list of initial vector values
#'@param E1 standardized expected counts
#'@param ncentroids number of centroids or centers
#'@param Time number of time periods
#'@param nsim number of simulations
#'@param threshold vector of two threshold values TODO: allow flexibility for number of threshold
#'@return returns list of probabilities.
 
get_prob <- function(lassoresult,init, E1, ncentroids, Time, nsim, threshold){
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
#'@param background option to return background in addition to probabilities. Default is NULL (just return probabilities, not background)
#'@return returns vector which calculated the number of time the cluster was correctly identified out of the simulations
prob_incluster <- function(select_mu, ncentroids, Time, nsim, background=NULL){
    #prob.bic <- prob_incluster(lassoresult$select_mu.qbic, ncentroids, Time, nsim)
    vec <- rep(0, ncentroids * Time)
    position <- list(vec)[rep(1, nsim)]
    bgRate_i <- lapply(1:nsim, function(i) sapply(1:Time,
                                                  function(j) as.numeric(names(which.max(table(matrix(as.vector(select_mu[[i]]),
                                                                                                      ncol=Time)[,j]))))))
    bgRate <- lapply(1:nsim, function(i) rep(bgRate_i[[i]], each = ncentroids))
    ix <- lapply(1:nsim, function(i) which(abs(as.vector(select_mu[[i]]) - bgRate[[i]])>=10^-3))
    #quick function to recode
    reval <- function(probs, ix){
        probs[ix] <-1
        return(probs)
    }
    simindicator <- mapply(reval, position, ix)
    probs <- Matrix::rowSums(simindicator)/nsim
    if(!is.null(background)){
        out <- list(bgRate = bgRate_i)
       return(out)
    }
    else{
        return(probs)
    }
}






