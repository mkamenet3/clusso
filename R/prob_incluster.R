#' @title 
#' prob_inoutcluster
#' @description
#' Finds the probability of being in the cluster for BIC, AIC, and AICc based on the expected risk ratio
#'@param select_mu List of selected mu vectors selected by the respective information criteria.
#'@param ncentroids number of centroids
#'@param Time number of time period
#'@param nsim number of simulations
#'@return returns vector which calculated the number of time the cluster was correctly identified out of the simulations
prob_inoutcluster <- function(lassoresult,rr, risk.ratio,x,y,rMax,nsim){
    #############SPACE-TIME
    #DEFINE TRUTH
    rrmatvec <- ifelse(as.vector(rr)!=risk.ratio,0,1)
    clusters <- clusters2df(x,y,rMax, utm=TRUE, length(x))
    n <- length(unique(clusters$center))
    potClus <- n
    numCenters <- n
    sparseMAT <- spacetimeMat(clusters, numCenters, Time)
    #GO through what was detected 
    #BIC
    select_mu <- as.vector(lassoresult$select.qbic)
    betaMat <- sapply(1:nsim, function(i) lassoresult$lasso[[i]]$beta)
    betaSelect <- sapply(1:nsim, function(i) betaMat[[i]][,select_mu[[i]]])
    betaSelect_bin <- ifelse(betaSelect!=0,1,0)
    clusteroverlap <- rrmatvec %*% sparseMAT #non-zeros are good - those touch the cluster
    clusteroverlap_bin <- ifelse(clusteroverlap !=0,1,0) #1= incluster, 0=not in cluster
    notincluster <- ifelse(clusteroverlap!=0,0,1) 
    notinclust_sim <- notincluster %*% betaSelect_bin
    notinclust_sim_bin <- ifelse(notinclust_sim!=0,1,0)
    #notinclust
    notinperc.bic<- paste0((sum(unlist(notinclust_sim_bin))/nsim)*100,"%")
    #inclust
    clustin_sim <- clusteroverlap %*% betaSelect_bin
    clustin_sim_bin <- ifelse(clustin_sim!=0,1,0)
    inperc.bic <- paste0((sum(unlist(clustin_sim_bin))/nsim)*100,"%")
    
    #AIC
    select_mu <- lassoresult$select.qaic
    betaMat <- sapply(1:nsim, function(i) lassoresult$lasso[[i]]$beta)
    betaSelect <- sapply(1:nsim, function(i) betaMat[[i]][,select_mu[[i]]])
    betaSelect_bin <- ifelse(betaSelect!=0,1,0)
    clusteroverlap <- rrmatvec %*% sparseMAT #non-zeros are good - those touch the cluster
    clusteroverlap_bin <- ifelse(clusteroverlap !=0,1,0) #1= incluster, 0=not in cluster
    notincluster <- ifelse(clusteroverlap!=0,0,1) 
    notinclust_sim <- notincluster %*% betaSelect_bin
    notinclust_sim_bin <- ifelse(notinclust_sim!=0,1,0)
    #notinclust
    notinperc.aic<- paste0((sum(unlist(notinclust_sim_bin))/nsim)*100,"%")
    #inclust
    clustin_sim <- clusteroverlap %*% betaSelect_bin
    clustin_sim_bin <- ifelse(clustin_sim!=0,1,0)
    inperc.aic <- paste0((sum(unlist(clustin_sim_bin))/nsim)*100,"%")
    
    #AICc
    select_mu <- lassoresult$select.qaicc
    betaMat <- sapply(1:nsim, function(i) lassoresult$lasso[[i]]$beta)
    betaSelect <- sapply(1:nsim, function(i) betaMat[[i]][,select_mu[[i]]])
    betaSelect_bin <- ifelse(betaSelect!=0,1,0)
    clusteroverlap <- rrmatvec %*% sparseMAT #non-zeros are good - those touch the cluster
    clusteroverlap_bin <- ifelse(clusteroverlap !=0,1,0) #1= incluster, 0=not in cluster
    notincluster <- ifelse(clusteroverlap!=0,0,1) 
    notinclust_sim <- notincluster %*% betaSelect_bin
    notinclust_sim_bin <- ifelse(notinclust_sim!=0,1,0)
    #notinclust
    notinperc.aicc<- paste0((sum(unlist(notinclust_sim_bin))/nsim)*100,"%")
    #inclust
    clustin_sim <- clusteroverlap %*% betaSelect_bin
    clustin_sim_bin <- ifelse(clustin_sim!=0,1,0)
    inperc.aicc <- paste0((sum(unlist(clustin_sim_bin))/nsim)*100,"%")
    
    return(list(notinperc.bic =notinperc.bic, notinperc.aic = notinperc.aic, notinperc.aicc = notinperc.aicc,
                inperc.bic = inperc.bic, inperc.aic = inperc.aic, inperc.aicc=inperc.aicc))
    
}    


