#' @title 
#' prob_clusteroverlap
#' @description
#' Finds the probability of any overlap with true cluster based on BIC, AIC, and AICc based on the expected risk ratio
#'@param sparseMAT large sparse matrix created in \code{clust_sim} function.
#'@param lassoresult List of QBIC, QAIC, QAICc estimates from the mylasso.sim function
#'@param rr risk ratio matrix that was used in the simulation
#'@param risk.ratio Risk ratio that was set for cluster in simulation
#'@param x x-coordinates
#'@param y y-coordinates
#'@param rMax Maximum radius for threshold in simulation
#'@param nsim number of simulations
#'@param Time number of time period
#'@param thresh Default is NULL; vector or value as threshold for cluster detection
#'@return returns vector which calculated the number of time the cluster was correctly identified out of the simulations
prob_clusteroverlap <- function(sparseMAT,lassoresult,rr, risk.ratio,x,y,rMax,nsim,Time, thresh){
    #DEFINE TRUTH
    if(risk.ratio==1){
     warning("Risk.ratio was set to 1")
        rrmatvec <- rep(0,length(rr))
    }
    else{
        rrmatvec <- ifelse(as.vector(rr)==risk.ratio,1,0)    
    }
    
    
    ##############################################
    bgRate_i
    
    
    
    ##############################################
    vec <- rep(0, ncentroids * Time) #ncentroids <- 208
    position <- list(vec)[rep(1, nsim)]
    select_mu <- lassoresult$select_mu.qaic
    
    bgRate_i <- lapply(1:nsim, 
                       function(i) sapply(1:Time,
                                                  function(j) as.numeric(names(which.max(table(matrix(as.vector(select_mu[[i]]),ncol=Time)[,j]))))))
    bgRate <- lapply(1:nsim, function(i) rep(bgRate_i[[i]], each = ncentroids))
    ix <- lapply(1:nsim, function(i) which(abs(log(as.vector(select_mu[[i]])) - log(bgRate[[i]]))>=10^-3))
    #quick function to recode
    reval <- function(probs, ix){
        probs[ix] <-1
        return(probs)
    }
    simindicator <- mapply(reval, position, ix)
    probs <- Matrix::rowSums(simindicator)/nsim
    

    ###################################
    #BIC
    ###################################
    select_mu <- as.vector(lassoresult$select.qbic)
    betaMat <- sapply(1:nsim, function(i) lassoresult$lasso[[i]]$beta)
    betaSelect <- sapply(1:nsim, function(i) betaMat[[i]][,select_mu[[i]]])
    #select what's not in the background
    betaSelect_bin <- lapply(1:nsim, function(i) ifelse(is.element(betaSelect[,i],log(bgRate_i[[i]])),0,1))
    clusteroverlap <- rrmatvec %*% sparseMAT #non-zeros are good - those touch the cluster
    clusteroverlap_bin <- ifelse(clusteroverlap !=0,1,0) #1= incluster, 0=not in cluster
    #select incluster
    clustin_sim <- lapply(1:nsim, function(i) clusteroverlap_bin %*% betaSelect_bin[[i]])
    clustin_sim_bin <- lapply(1:nsim, function(i) ifelse(clustin_sim[[i]]>0,1,0))
    inperc.bic <- paste0((sum(unlist(clustin_sim_bin))/nsim)*100,"%")
    #select not in cluster
    notincluster <- ifelse(clusteroverlap_bin==1,0,1)
    notinclust_sim <- lapply(1:nsim, function(i) notincluster %*% betaSelect_bin[[i]])
    notinclust_sim_bin <- ifelse(notinclust_sim!=0,1,0)
    #notinclust
    notinperc.bic<- paste0((sum(unlist(notinclust_sim_bin))/nsim)*100,"%")

    
    # ################################################
    # #print(table(rrmatvec))
    # #GO through what was detected 
    # #Let A = true cluster (clusteroverlap), B = detected cluster (betaSelect_bin)
    # ###################################
    # #BIC
    # ###################################
    # select_mu <- as.vector(lassoresult$select.qbic)
    # betaMat <- sapply(1:nsim, function(i) lassoresult$lasso[[i]]$beta)
    # betaSelect <- sapply(1:nsim, function(i) betaMat[[i]][,select_mu[[i]]])
    # betaSelect_bin <- lapply(1:nsim, function(i) ifelse(betaSelect[,i]>0,1,0))
    # 
    # ##########################
    # clusteroverlap <- rrmatvec %*% sparseMAT #non-zeros are good - those touch the cluster
    # clusteroverlap_bin <- ifelse(clusteroverlap !=0,1,0) #1= incluster, 0=not in cluster
    # notincluster <- ifelse(clusteroverlap_bin==1,0,1)
    # notinclust_sim <- lapply(1:nsim, function(i) notincluster %*% betaSelect_bin[[i]])
    # notinclust_sim_bin <- ifelse(notinclust_sim!=0,1,0)
    # #notinclust
    # notinperc.bic<- paste0((sum(unlist(notinclust_sim_bin))/nsim)*100,"%")
    # #inclust
    # #clustin_sim <- clusteroverlap %*% betaSelect
    # clustin_sim <- lapply(1:nsim, function(i) clusteroverlap %*% betaSelect_bin[[i]])
    # clustin_sim_bin <- lapply(1:nsim, function(i) ifelse(clustin_sim[[i]]>0,1,0))
    # inperc.bic <- paste0((sum(unlist(clustin_sim_bin))/nsim)*100,"%")

    if(!is.null(thresh)){
    ##Diagnostics with thresh
        truth <- which(rrmatvec!=0) #this is true location of cluster
        detected <- sapply(1:nsim, function(i) sparseMAT %*% betaSelect_bin[[i]]) 
        
        truth_and_detected <- sapply(1:nsim, function(i) union(which(as.vector(detected[[i]])!=0), truth))
        intersect_AandB <- sapply(1:nsim, function(i) ifelse(length(clustin_sim)/length(truth_and_detected[[i]])>thresh,1,0))
        percintersect_AandB.bic <- paste0(mean(unlist(intersect_AandB)), "%")
        
        ##|(A and B)|/|B|?
        intersect_AandB <- sapply(1:nsim, function(i) ifelse(length(clustin_sim)/length(which(as.vector(detected[[i]])!=0))>thresh,1,0))
        percintersect_B.bic <- paste0(mean(unlist(intersect_AandB)), "%")
        
        ##|(A and B)|/|A|?
        intersect_AandB <- sapply(1:nsim, function(i) ifelse(length(clustin_sim)/length(truth)>thresh,1,0))
        percintersect_A.bic <- paste0(mean(unlist(intersect_AandB)), "%")
    }
    else{
        message("No threshold diagnostics - BIC")
    }
    
    ###################################
    #AIC
    ###################################
    select_mu <- as.vector(lassoresult$select.qaic)
    betaMat <- sapply(1:nsim, function(i) lassoresult$lasso[[i]]$beta)
    betaSelect <- sapply(1:nsim, function(i) betaMat[[i]][,select_mu[[i]]])
    #select what's not in the background
    betaSelect_bin <- lapply(1:nsim, function(i) ifelse(is.element(betaSelect[,i],log(bgRate_i[[i]])),0,1))
    clusteroverlap <- rrmatvec %*% sparseMAT #non-zeros are good - those touch the cluster
    clusteroverlap_bin <- ifelse(clusteroverlap !=0,1,0) #1= incluster, 0=not in cluster
    #select incluster
    clustin_sim <- lapply(1:nsim, function(i) clusteroverlap_bin %*% betaSelect_bin[[i]])
    clustin_sim_bin <- lapply(1:nsim, function(i) ifelse(clustin_sim[[i]]>0,1,0))
    inperc.bic <- paste0((sum(unlist(clustin_sim_bin))/nsim)*100,"%")
    #select not in cluster
    notincluster <- ifelse(clusteroverlap_bin==1,0,1)
    notinclust_sim <- lapply(1:nsim, function(i) notincluster %*% betaSelect_bin[[i]])
    notinclust_sim_bin <- ifelse(notinclust_sim!=0,1,0)
    #notinclust
    notinperc.bic<- paste0((sum(unlist(notinclust_sim_bin))/nsim)*100,"%")
    
    # select_mu <- lassoresult$select.qaic
    # betaMat <- sapply(1:nsim, function(i) lassoresult$lasso[[i]]$beta)
    # betaSelect <- sapply(1:nsim, function(i) betaMat[[i]][,select_mu[[i]]])
    # betaSelect_bin <- lapply(1:nsim, function(i) ifelse(betaSelect[,i]>0,1,0))
    # 
    # ##########################
    # clusteroverlap <- rrmatvec %*% sparseMAT #non-zeros are good - those touch the cluster
    # clusteroverlap_bin <- ifelse(clusteroverlap !=0,1,0) #1= incluster, 0=not in cluster
    # notincluster <- ifelse(clusteroverlap_bin==1,0,1)
    # notinclust_sim <- lapply(1:nsim, function(i) notincluster %*% betaSelect_bin[[i]])
    # notinclust_sim_bin <- ifelse(notinclust_sim!=0,1,0)
    # #notinclust
    # notinperc.aic<- paste0((sum(unlist(notinclust_sim_bin))/nsim)*100,"%")
    # #inclust
    # clustin_sim <- lapply(1:nsim, function(i) clusteroverlap_bin %*% betaSelect_bin[[i]])
    # clustin_sim_bin <- lapply(1:nsim, function(i) ifelse(clustin_sim[[i]]>0,1,0))
    # inperc.aic <- paste0((sum(unlist(clustin_sim_bin))/nsim)*100,"%")
    
    ###########################
    if(!is.null(thresh)){
        ##Diagnostics with thresh
        ##|(A and B)|/|A U B|? 
        ##AandB = clustin_sim
        ##AUB = truth_and_detected
        truth <- which(rrmatvec!=0) #this is true location of cluster
        detected <- sapply(1:nsim, function(i) sparseMAT %*% betaSelect_bin[[i]]) 
        truth_and_detected <- sapply(1:nsim, function(i) union(which(as.vector(detected[[i]])!=0), truth))
        intersect_AandB <- sapply(1:nsim, function(i) ifelse(length(clustin_sim)/length(truth_and_detected[[i]])>thresh,1,0))
        percintersect_AandB.aic <- paste0(mean(unlist(intersect_AandB)), "%")
        
        ##|(A and B)|/|B|?
        intersect_AandB <- sapply(1:nsim, function(i) ifelse(length(clustin_sim)/length(which(as.vector(detected[[i]])!=0))>thresh,1,0))
        percintersect_B.aic <- paste0(mean(unlist(intersect_AandB)), "%")
        
        ##|(A and B)|/|A|?
        intersect_AandB <- sapply(1:nsim, function(i) ifelse(length(clustin_sim)/length(truth)>thresh,1,0))
        percintersect_A.aic <- paste0(mean(unlist(intersect_AandB)), "%")
    }
    else{
        message("No threshold diagnostics - AIC")
    }
    
    ###################################
    #AICc
    ###################################
    select_mu <- as.vector(lassoresult$select.qaicc)
    betaMat <- sapply(1:nsim, function(i) lassoresult$lasso[[i]]$beta)
    betaSelect <- sapply(1:nsim, function(i) betaMat[[i]][,select_mu[[i]]])
    #select what's not in the background
    betaSelect_bin <- lapply(1:nsim, function(i) ifelse(is.element(betaSelect[,i],log(bgRate_i[[i]])),0,1))
    clusteroverlap <- rrmatvec %*% sparseMAT #non-zeros are good - those touch the cluster
    clusteroverlap_bin <- ifelse(clusteroverlap !=0,1,0) #1= incluster, 0=not in cluster
    #select incluster
    clustin_sim <- lapply(1:nsim, function(i) clusteroverlap_bin %*% betaSelect_bin[[i]])
    clustin_sim_bin <- lapply(1:nsim, function(i) ifelse(clustin_sim[[i]]>0,1,0))
    inperc.bic <- paste0((sum(unlist(clustin_sim_bin))/nsim)*100,"%")
    #select not in cluster
    notincluster <- ifelse(clusteroverlap_bin==1,0,1)
    notinclust_sim <- lapply(1:nsim, function(i) notincluster %*% betaSelect_bin[[i]])
    notinclust_sim_bin <- ifelse(notinclust_sim!=0,1,0)
    #notinclust
    notinperc.bic<- paste0((sum(unlist(notinclust_sim_bin))/nsim)*100,"%")
    # select_mu <- lassoresult$select.qaicc
    # betaMat <- sapply(1:nsim, function(i) lassoresult$lasso[[i]]$beta)
    # betaSelect <- sapply(1:nsim, function(i) betaMat[[i]][,select_mu[[i]]])
    # betaSelect_bin <- lapply(1:nsim, function(i) ifelse(betaSelect[,i]>0,1,0))
    # 
    # ##########################
    # 
    # clusteroverlap <- rrmatvec %*% sparseMAT #non-zeros are good - those touch the cluster
    # clusteroverlap_bin <- ifelse(clusteroverlap !=0,1,0) #1= incluster, 0=not in cluster
    # notincluster <- ifelse(clusteroverlap_bin==1,0,1)
    # notinclust_sim <- lapply(1:nsim, function(i) notincluster %*% betaSelect_bin[[i]])
    # notinclust_sim_bin <- ifelse(notinclust_sim!=0,1,0)
    # #notinclust
    # notinperc.aicc<- paste0((sum(unlist(notinclust_sim_bin))/nsim)*100,"%")
    # #inclust
    # #clustin_sim <- clusteroverlap %*% betaSelect
    # clustin_sim <- lapply(1:nsim, function(i) clusteroverlap %*% betaSelect_bin[[i]])
    # #clustin_sim_bin <- ifelse(clustin_sim>0,1,0)
    # clustin_sim_bin <- lapply(1:nsim, function(i) ifelse(clustin_sim[[i]]>0,1,0))
    # inperc.aicc <- paste0((sum(unlist(clustin_sim_bin))/nsim)*100,"%")
   
    if(!is.null(thresh)){
        ##Diagnostics with thresh
        ##|(A and B)|/|A U B|? 
        ##AandB = clustin_sim
        ##AUB = truth_and_detected
        truth <- which(rrmatvec!=0) #this is true location of cluster
        detected <- sapply(1:nsim, function(i) sparseMAT %*% betaSelect_bin[[i]]) 
        truth_and_detected <- sapply(1:nsim, function(i) union(which(as.vector(detected[[i]])!=0), truth))
        intersect_AandB <- sapply(1:nsim, function(i) ifelse(length(clustin_sim)/length(truth_and_detected[[i]])>thresh,1,0))
        percintersect_AandB.aicc <- paste0(mean(unlist(intersect_AandB)), "%")
        
        ##|(A and B)|/|B|?
        intersect_AandB <- sapply(1:nsim, function(i) ifelse(length(clustin_sim)/length(which(as.vector(detected[[i]])!=0))>thresh,1,0))
        percintersect_B.aicc <- paste0(mean(unlist(intersect_AandB)), "%")
        
        ##|(A and B)|/|A|?
        intersect_AandB <- sapply(1:nsim, function(i) ifelse(length(clustin_sim)/length(truth)>thresh,1,0))
        percintersect_A.aicc <- paste0(mean(unlist(intersect_AandB)), "%")
    }
    else{
        message("No threshold diagnostics - AICc")
    }
    if(is.null(thresh)){
        return(list(notinperc.bic =notinperc.bic, notinperc.aic = notinperc.aic, notinperc.aicc = notinperc.aicc,
                    inperc.bic = inperc.bic, inperc.aic = inperc.aic, inperc.aicc=inperc.aicc))
    }
    else{
        #new return with additional thresh calcs
        return(list(notinperc.bic =notinperc.bic, notinperc.aic = notinperc.aic, notinperc.aicc = notinperc.aicc,
                    inperc.bic = inperc.bic, inperc.aic = inperc.aic, inperc.aicc=inperc.aicc, 
                    thresh = list(c(percintersect_AandB.bic,  percintersect_AandB.aic, percintersect_AandB.aicc,
                    percintersect_A.bic, percintersect_A.aic, percintersect_A.aicc,
                    percintersect_B.bic, percintersect_B.aic, percintersect_B.aicc),
                    c("AandB.bic","AandB.aic","AandB.aicc","A.bic","A.aic","A.aicc","B.bic","B.aic", "B.aicc"))))
        }
}    
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
#'@param option to return background in addition to probabilities. Default is NULL (just return probabilities, not background)
#'@return returns vector which calculated the number of time the cluster was correctly identified out of the simulations
prob_incluster <- function(select_mu, ncentroids, Time, nsim, background=NULL){
    vec <- rep(0, ncentroids * Time)
    position <- list(vec)[rep(1, nsim)]
    bgRate_i <- lapply(1:nsim, function(i) sapply(1:Time,
                                                  function(j) as.numeric(names(which.max(table(matrix(as.vector(select_mu[[i]]),ncol=Time)[,j]))))))
    bgRate <- lapply(1:nsim, function(i) rep(bgRate_i[[i]], each = ncentroids))
    ix <- lapply(1:nsim, function(i) which(abs(log(as.vector(select_mu[[i]])) - log(bgRate[[i]]))>=10^-3))
    #quick function to recode
    reval <- function(probs, ix){
        probs[ix] <-1
        return(probs)
    }
    simindicator <- mapply(reval, position, ix)
    probs <- Matrix::rowSums(simindicator)/nsim
    if(!is.null(background)){
        out <- list(probs = probs, bgRate = bgRate_i)
       return(out) 
    }
    else{
        return(probs)    
    }
}  

