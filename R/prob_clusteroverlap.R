#11-1-18 testing new detection function
prob_clusteroverlap2 <- function(sparseMAT,lassoresult,rr, risk.ratio,x,y,rMax,nsim,Time, ncentroids, nullmod=NULL){
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
    ##(Q)BIC
    selected <- lassoresult$select.qbic
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
    selected <- lassoresult$select.qaic
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
    
    # return(list(outperc.qbic =outperc.qbic,# outperc.qaic = outperc.qaic, outperc.qaicc = outperc.qaicc,
    #             inperc.qbic = inperc.qbic),#, inperc.qaic = inperc.qaic, inperc.qaicc=inperc.qaicc, 
    #                           c("AandB.bic","AandB.aic","AandB.aicc","A.bic","A.aic","A.aicc","B.bic","B.aic", "B.aicc"))

    
    
}
    


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
#'@param ncentroids number of centroids
#'@return returns vector which calculated the number of time the cluster was correctly identified out of the simulations
#prob_clusteroverlap <- function(sparseMAT,lassoresult,rr, risk.ratio,x,y,rMax,nsim,Time, thresh, ncentroids){
prob_clusteroverlap <- function(sparseMAT,lassoresult,rr, risk.ratio,x,y,rMax,nsim,Time, thresh,ncentroids, nullmod=NULL){
    #DEFINE TRUTH
    if(risk.ratio==1){
     warning("Risk.ratio was set to 1")
        rrmatvec <- rep(0,length(rr))
    }
    else{
        rrmatvec <- ifelse(as.vector(rr)==risk.ratio,1,0)    
    }
    sparseMAT_clusteronly <- sparseMAT[,-c(ncol(sparseMAT)-Time+1:ncol(sparseMAT))]
    
    # ###################################
    # #TEST 10-18-18 - this works thank goodness
    # ###################################
    # #test <- t(sparseMAT)%*%rrmatvec
    # sparseMAT_clusteronly <- sparseMAT[,-c(ncol(sparseMAT)-Time+1:ncol(sparseMAT))]
    # 
    # ###INCLUSTER
    # clusteroverlap <- t(rrmatvec) %*% sparseMAT_clusteronly 
    # clusteroverlap_bin <- ifelse(clusteroverlap !=0,1,0) 
    # selected <- lassoresult.qp.st$select.qbic
    # betaselect <- lapply(1:nsim, function(i) lassoresult.qp.st$coefs.lasso.all[[i]][,selected[[i]]])
    # #betaselect_bin[[1]][-c(ncol(sparseMAT)-Time+1:ncol(sparseMAT))]
    # betaselect_bin <- lapply(1:nsim, function(i) ifelse(betaselect[[i]]!=0,1,0))
    # betaselect_bin_clusteronly <- lapply(1:nsim, function(i) betaselect_bin[[i]][-c(ncol(sparseMAT)-Time+1:ncol(sparseMAT))])
    # #incluster_sim <- lapply(1:nsim, function(i) clusteroverlap_bin %*% betaselect_bin[[i]][-c(ncol(sparseMAT)-Time+1:ncol(sparseMAT))])
    # incluster_sim <- lapply(1:nsim, function(i) clusteroverlap_bin %*% betaselect_bin_clusteronly[[i]])
    # incluster_sim_bin <- lapply(1:nsim, function(i) ifelse(incluster_sim[[i]] !=0,1,0)) 
    # inperc <- paste0((sum(unlist(incluster_sim_bin))/nsim)*100, "%")
    # 
    # ##OUTCLUSTER
    # clusteroverlap_bin <- ifelse(clusteroverlap !=0,0,1) 
    # selected <- lassoresult.qp.st$select.qbic
    # betaselect <- lapply(1:nsim, function(i) lassoresult.qp.st$coefs.lasso.all[[i]][,selected[[i]]])
    # betaselect_bin <- lapply(1:nsim, function(i) ifelse(betaselect[[i]]!=0,1,0))
    # betaselect_bin_clusteronly <- lapply(1:nsim, function(i) betaselect_bin[[i]][-c(ncol(sparseMAT)-Time+1:ncol(sparseMAT))])
    # outcluster_sim <- lapply(1:nsim, function(i) clusteroverlap_bin %*% betaselect_bin_clusteronly[[i]])
    # outcluster_sim_bin <- lapply(1:nsim, function(i) ifelse(outcluster_sim[[i]] !=0,1,0)) 
    # outperc <- paste0((sum(unlist(outcluster_sim_bin))/nsim)*100, "%")
    # 
    #lassoresult <- lassoresult.qp.st
    ###################################
    #(Q)BIC
    ###################################
    #selected <- lassoresult$select.qbic
    #selected <- lassoresultnull$select.qaic
    selected <- lassoresult$select.qaic
    betaselect <- lapply(1:nsim, function(i) lassoresult$coefs.lasso.all[[i]][,selected[[i]]])
    betaselect_bin <- lapply(1:nsim, function(i) ifelse(betaselect[[i]]!=0,1,0))
    
    ###
    betaselect <- lapply(1:nsim, function(i) lassoresult$coefs.lasso.all[[i]][,selected[[i]]])
    betaselect_bin <- lapply(1:nsim, function(i) ifelse(abs(betaselect[[i]] >= 10e-3),1,0))
    ###
    betaselect_bin_clusteronly <- lapply(1:nsim, function(i) betaselect_bin[[i]][-c(ncol(sparseMAT)-Time+1:ncol(sparseMAT))])
    ###INCLUSTER
    clusteroverlap <- t(rrmatvec) %*% sparseMAT_clusteronly
    clusteroverlap_bin <- ifelse(clusteroverlap !=0,1,0)
    incluster_sim <- lapply(1:nsim, function(i) clusteroverlap_bin %*% betaselect_bin_clusteronly[[i]])
    incluster_sim_bin <- lapply(1:nsim, function(i) ifelse(incluster_sim[[i]] !=0,1,0))
    inperc.qbic <- paste0((sum(unlist(incluster_sim_bin))/nsim)*100, "%")
    if(!is.null(nullmod)){
        a <- lapply(1:nsim, function(i) any(betaselect_bin_clusteronly[[i]]!=0))
        inperc.qbic <- sum(unlist(a))/nsim
    }
    ##OUTCLUSTER
    clusteroverlap_bin <- ifelse(clusteroverlap !=0,0,1)
    # selected <- lassoresult.qp.st$select.qbic
    # betaselect <- lapply(1:nsim, function(i) lassoresult.qp.st$coefs.lasso.all[[i]][,selected[[i]]])
    # betaselect_bin <- lapply(1:nsim, function(i) ifelse(betaselect[[i]]!=0,1,0))
    # betaselect_bin_clusteronly <- lapply(1:nsim, function(i) betaselect_bin[[i]][-c(ncol(sparseMAT)-Time+1:ncol(sparseMAT))])
    outcluster_sim <- lapply(1:nsim, function(i) clusteroverlap_bin %*% betaselect_bin_clusteronly[[i]])
    outcluster_sim_bin <- lapply(1:nsim, function(i) ifelse(outcluster_sim[[i]] !=0,1,0))
    outperc.qbic <- paste0((sum(unlist(outcluster_sim_bin))/nsim)*100, "%")
    
    # select_mu <- as.vector(lassoresult$select.qbic) #selected path for each sim
    # 
    # #bgRate_i <- prob_incluster(lassoresult$select_mu.qbic, ncentroids, Time, nsim, background = TRUE)
    # betaMat <-  sapply(1:nsim, function(i) lassoresult$lasso[[i]]$beta) #rows are potential clusters; columns are lambdas
    # #betaMat <- sapply(1:nsim, function(i) lassoresult.qp.st$lasso[[i]]$beta) #rows are potential clusters; columns are lambdas
    # 
    # betaSelect <- sapply(1:nsim, function(i) betaMat[[i]][,select_mu[[i]]]) 
    # #these are columns of betas selected for each n in nsim
    # 
    # #select what's not in the background
    # betaSelect_bin <- lapply(1:nsim, function(i) ifelse(is.element(betaSelect[,i],log(bgRate_i$bgRate[[i]])),0,1))
    # 
    # 
    # 
    # clusteroverlap <- t(rrmatvec) %*% sparseMAT #non-zeros are good - those touch the cluster
    # clusteroverlap_bin <- ifelse(clusteroverlap !=0,1,0) #1= incluster, 0=not in cluster
    # #select incluster
    # clustin_sim <- lapply(1:nsim, function(i) t(clusteroverlap_bin) %*% betaSelect_bin[[i]]) #something weird here
    # clustin_sim_bin <- lapply(1:nsim, function(i) ifelse(clustin_sim[[i]]>0,1,0))
    # inperc.bic <- paste0((sum(unlist(clustin_sim_bin))/nsim)*100,"%")
    # #select not in cluster
    # notincluster <- ifelse(clusteroverlap_bin==1,0,1)
    # notinclust_sim <- lapply(1:nsim, function(i) t(notincluster) %*% betaSelect_bin[[i]])
    # notinclust_sim_bin <- ifelse(notinclust_sim!=0,1,0)
    # #notinclust
    # notinperc.bic<- paste0((sum(unlist(notinclust_sim_bin))/nsim)*100,"%")


    # if(!is.null(thresh)){
    # ##Diagnostics with thresh
    #     
    #     truth <- which(rrmatvec!=0) #this is true location of cluster
    #     detected <- sapply(1:nsim, function(i) sparseMAT %*% betaSelect_bin[[i]]) 
    #     detectedincluster <- lapply(1:nsim, function(i) rrmatvec %*% ifelse(detected[[i]]@x!=0,1,0))
    #     
    #     ##|(A and B)|/|A U B|?
    #     truth_and_detected <- sapply(1:nsim, function(i) union(which(as.vector(detected[[i]]@x)!=0), truth))
    #     intersect_AandB <- sapply(1:nsim, function(i) ifelse(detectedincluster[[i]]/length(truth_and_detected[[i]])>thresh,1,0))
    #     percintersect_AandB.bic <- paste0(mean(unlist(intersect_AandB))*100, "%")
    #     
    #     ##|(A and B)|/|B|?
    #     intersect_AandB <- sapply(1:nsim, function(i) ifelse((detectedincluster[[i]]/length(which(detected[[i]]@x!=0)))>thresh,1,0))
    #     percintersect_B.bic <- paste0(mean(unlist(intersect_AandB))*100, "%")
    #     
    #     ##|(A and B)|/|A|?
    #     intersect_AandB <- sapply(1:nsim, function(i) ifelse((detectedincluster[[i]]/length(truth))>thresh,1,0))
    #     percintersect_A.bic <- paste0(mean(unlist(intersect_AandB))*100, "%")
    # }
    # else{
    #     message("No threshold diagnostics - BIC")
    # }
    
    ###################################
    #AIC
    ###################################
    selected <- lassoresult$select.qaic
    betaselect <- lapply(1:nsim, function(i) lassoresult$coefs.lasso.all[[i]][,selected[[i]]])
    betaselect_bin <- lapply(1:nsim, function(i) ifelse(betaselect[[i]]!=0,1,0))
    betaselect_bin_clusteronly <- lapply(1:nsim, function(i) betaselect_bin[[i]][-c(ncol(sparseMAT)-Time+1:ncol(sparseMAT))])
    if(!is.null(nullmod)){
        a <- lapply(1:nsim, function(i) any(betaselect_bin_clusteronly[[i]]!=0))
        inperc.qaic <- sum(unlist(a))/nsim
    }
    ###INCLUSTER
    clusteroverlap <- t(rrmatvec) %*% sparseMAT_clusteronly
    clusteroverlap_bin <- ifelse(clusteroverlap !=0,1,0)
    incluster_sim <- lapply(1:nsim, function(i) clusteroverlap_bin %*% betaselect_bin_clusteronly[[i]])
    incluster_sim_bin <- lapply(1:nsim, function(i) ifelse(incluster_sim[[i]] !=0,1,0))
    inperc.qaic <- paste0((sum(unlist(incluster_sim_bin))/nsim)*100, "%")
    ##OUTCLUSTER
    clusteroverlap_bin <- ifelse(clusteroverlap !=0,0,1)
    # selected <- lassoresult.qp.st$select.qbic
    # betaselect <- lapply(1:nsim, function(i) lassoresult.qp.st$coefs.lasso.all[[i]][,selected[[i]]])
    # betaselect_bin <- lapply(1:nsim, function(i) ifelse(betaselect[[i]]!=0,1,0))
    # betaselect_bin_clusteronly <- lapply(1:nsim, function(i) betaselect_bin[[i]][-c(ncol(sparseMAT)-Time+1:ncol(sparseMAT))])
    outcluster_sim <- lapply(1:nsim, function(i) clusteroverlap_bin %*% betaselect_bin_clusteronly[[i]])
    outcluster_sim_bin <- lapply(1:nsim, function(i) ifelse(outcluster_sim[[i]] !=0,1,0))
    outperc.qaic <- paste0((sum(unlist(outcluster_sim_bin))/nsim)*100, "%")
    
    # select_mu <- as.vector(lassoresult$select.qaic)
    # bgRate_i <- prob_incluster(lassoresult$select_mu.qaic, ncentroids, Time, nsim, background = TRUE)
    # betaMat <- sapply(1:nsim, function(i) lassoresult$lasso[[i]]$beta)
    # betaSelect <- sapply(1:nsim, function(i) betaMat[[i]][,select_mu[[i]]])
    # #select what's not in the background
    # betaSelect_bin <- lapply(1:nsim, function(i) ifelse(is.element(betaSelect[,i],log(bgRate_i$bgRate[[i]])),0,1))
    # clusteroverlap <- rrmatvec %*% sparseMAT #non-zeros are good - those touch the cluster
    # clusteroverlap_bin <- ifelse(clusteroverlap !=0,1,0) #1= incluster, 0=not in cluster
    # #select incluster
    # clustin_sim <- lapply(1:nsim, function(i) clusteroverlap_bin %*% betaSelect_bin[[i]])
    # clustin_sim_bin <- lapply(1:nsim, function(i) ifelse(clustin_sim[[i]]>0,1,0))
    # inperc.aic <- paste0((sum(unlist(clustin_sim_bin))/nsim)*100,"%")
    # #select not in cluster
    # notincluster <- ifelse(clusteroverlap_bin==1,0,1)
    # notinclust_sim <- lapply(1:nsim, function(i) notincluster %*% betaSelect_bin[[i]])
    # notinclust_sim_bin <- ifelse(notinclust_sim!=0,1,0)
    # #notinclust
    # notinperc.aic<- paste0((sum(unlist(notinclust_sim_bin))/nsim)*100,"%")
    # 
   
    
    ###########################
    # if(!is.null(thresh)){
    #     
    #     truth <- which(rrmatvec!=0) #this is true location of cluster
    #     detected <- sapply(1:nsim, function(i) sparseMAT %*% betaSelect_bin[[i]]) 
    #     detectedincluster <- lapply(1:nsim, function(i) rrmatvec %*% ifelse(detected[[i]]@x!=0,1,0))
    #     
    #     ##|(A and B)|/|A U B|?
    #     truth_and_detected <- sapply(1:nsim, function(i) union(which(as.vector(detected[[i]]@x)!=0), truth))
    #     intersect_AandB <- sapply(1:nsim, function(i) ifelse(detectedincluster[[i]]/length(truth_and_detected[[i]])>thresh,1,0))
    #     percintersect_AandB.aic <- paste0(mean(unlist(intersect_AandB))*100, "%")
    #     
    #     ##|(A and B)|/|B|?
    #     intersect_AandB <- sapply(1:nsim, function(i) ifelse((detectedincluster[[i]]/length(which(detected[[i]]@x!=0)))>thresh,1,0))
    #     percintersect_B.aic <- paste0(mean(unlist(intersect_AandB))*100, "%")
    #     
    #     ##|(A and B)|/|A|?
    #     intersect_AandB <- sapply(1:nsim, function(i) ifelse((detectedincluster[[i]]/length(truth))>thresh,1,0))
    #     percintersect_A.aic <- paste0(mean(unlist(intersect_AandB))*100, "%")
    # }
    # else{
    #     message("No threshold diagnostics - AIC")
    # }
    
    ###################################
    #AICc
    ###################################
    selected <- lassoresult$select.qaicc
    betaselect <- lapply(1:nsim, function(i) lassoresult$coefs.lasso.all[[i]][,selected[[i]]])
    betaselect_bin <- lapply(1:nsim, function(i) ifelse(betaselect[[i]]!=0,1,0))
    betaselect_bin_clusteronly <- lapply(1:nsim, function(i) betaselect_bin[[i]][-c(ncol(sparseMAT)-Time+1:ncol(sparseMAT))])
    ###INCLUSTER
    clusteroverlap <- t(rrmatvec) %*% sparseMAT_clusteronly
    clusteroverlap_bin <- ifelse(clusteroverlap !=0,1,0)
    incluster_sim <- lapply(1:nsim, function(i) clusteroverlap_bin %*% betaselect_bin_clusteronly[[i]])
    incluster_sim_bin <- lapply(1:nsim, function(i) ifelse(incluster_sim[[i]] !=0,1,0))
    inperc.qaicc <- paste0((sum(unlist(incluster_sim_bin))/nsim)*100, "%")
    ##OUTCLUSTER
    clusteroverlap_bin <- ifelse(clusteroverlap !=0,0,1)
    # selected <- lassoresult.qp.st$select.qbic
    # betaselect <- lapply(1:nsim, function(i) lassoresult.qp.st$coefs.lasso.all[[i]][,selected[[i]]])
    # betaselect_bin <- lapply(1:nsim, function(i) ifelse(betaselect[[i]]!=0,1,0))
    # betaselect_bin_clusteronly <- lapply(1:nsim, function(i) betaselect_bin[[i]][-c(ncol(sparseMAT)-Time+1:ncol(sparseMAT))])
    outcluster_sim <- lapply(1:nsim, function(i) clusteroverlap_bin %*% betaselect_bin_clusteronly[[i]])
    outcluster_sim_bin <- lapply(1:nsim, function(i) ifelse(outcluster_sim[[i]] !=0,1,0))
    outperc.qaicc <- paste0((sum(unlist(outcluster_sim_bin))/nsim)*100, "%")
    # select_mu <- as.vector(lassoresult$select.qaicc)
    # bgRate_i <- prob_incluster(lassoresult$select_mu.qaicc, ncentroids, Time, nsim, background = TRUE)
    # betaMat <- sapply(1:nsim, function(i) lassoresult$lasso[[i]]$beta)
    # betaSelect <- sapply(1:nsim, function(i) betaMat[[i]][,select_mu[[i]]])
    # #select what's not in the background
    # betaSelect_bin <- lapply(1:nsim, function(i) ifelse(is.element(betaSelect[,i],log(bgRate_i$bgRate[[i]])),0,1))
    # clusteroverlap <- rrmatvec %*% sparseMAT #non-zeros are good - those touch the cluster
    # clusteroverlap_bin <- ifelse(clusteroverlap !=0,1,0) #1= incluster, 0=not in cluster
    # #select incluster
    # clustin_sim <- lapply(1:nsim, function(i) clusteroverlap_bin %*% betaSelect_bin[[i]])
    # clustin_sim_bin <- lapply(1:nsim, function(i) ifelse(clustin_sim[[i]]>0,1,0))
    # inperc.aicc <- paste0((sum(unlist(clustin_sim_bin))/nsim)*100,"%")
    # #select not in cluster
    # notincluster <- ifelse(clusteroverlap_bin==1,0,1)
    # notinclust_sim <- lapply(1:nsim, function(i) notincluster %*% betaSelect_bin[[i]])
    # notinclust_sim_bin <- ifelse(notinclust_sim!=0,1,0)
    # #notinclust
    # notinperc.aicc<- paste0((sum(unlist(notinclust_sim_bin))/nsim)*100,"%")
    
   
    # if(!is.null(thresh)){
    #     ##Diagnostics with thresh
    #     
    #     truth <- which(rrmatvec!=0) #this is true location of cluster
    #     detected <- sapply(1:nsim, function(i) sparseMAT %*% betaSelect_bin[[i]]) 
    #     detectedincluster <- lapply(1:nsim, function(i) rrmatvec %*% ifelse(detected[[i]]@x!=0,1,0))
    #     
    #     ##|(A and B)|/|A U B|?
    #     truth_and_detected <- sapply(1:nsim, function(i) union(which(as.vector(detected[[i]]@x)!=0), truth))
    #     intersect_AandB <- sapply(1:nsim, function(i) ifelse(detectedincluster[[i]]/length(truth_and_detected[[i]])>thresh,1,0))
    #     percintersect_AandB.aicc <- paste0(mean(unlist(intersect_AandB))*100, "%")
    #     
    #     ##|(A and B)|/|B|?
    #     intersect_AandB <- sapply(1:nsim, function(i) ifelse((detectedincluster[[i]]/length(which(detected[[i]]@x!=0)))>thresh,1,0))
    #     percintersect_B.aicc <- paste0(mean(unlist(intersect_AandB))*100, "%")
    #     
    #     ##|(A and B)|/|A|?
    #     intersect_AandB <- sapply(1:nsim, function(i) ifelse((detectedincluster[[i]]/length(truth))>thresh,1,0))
    #     percintersect_A.aicc <- paste0(mean(unlist(intersect_AandB))*100, "%")
    # }
    # else{
    #     message("No threshold diagnostics - AICc")
    # }
    if(is.null(thresh)){
        return(list(outperc.qbic =outperc.qbic, outperc.qaic = outperc.qaic, outperc.qaicc = outperc.qaicc,
                    inperc.qbic = inperc.qbic, inperc.qaic = inperc.qaic, inperc.qaicc=inperc.qaicc))
    }
    else{
        #new return with additional thresh calcs
        return(list(outperc.bic =outperc.bic, outperc.aic = outperc.aic, outperc.aicc = outperc.aicc,
                    inperc.bic = inperc.bic, inperc.aic = inperc.aic, inperc.aicc=inperc.aicc, 
                    thresh = list(c(percintersect_AandB.bic,  percintersect_AandB.aic, percintersect_AandB.aicc,
                    percintersect_A.bic, percintersect_A.aic, percintersect_A.aicc,
                    percintersect_B.bic, percintersect_B.aic, percintersect_B.aicc),
                    c("AandB.bic","AandB.aic","AandB.aicc","A.bic","A.aic","A.aicc","B.bic","B.aic", "B.aicc"))))
        }
}    
#' #' @title
#' #' get_prob
#' #' @description 
#' #' Finds the probability of being in the cluster for BIC, AIC, and AICc based on the expected risk ratio
#' #' @param lassoresult result of either space-time lasso simulation or space-only simulation
#' #' @param init list of initial vector values
#' #' @param E1 standardized expected counts
#' #' @param ncentroids number of centroids or centers
#' #' @param Time number of time periods
#' #' @param nsim number of simulations
#' #' @param threshold vector of two threshold values TODO: allow flexibility for number of threshold
#' #' @return returns list of probabilities.
#' #' 
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
    #prob.bic <- prob_incluster(lassoresult$select_mu.qbic, ncentroids, Time, nsim)
    vec <- rep(0, ncentroids * Time)
    position <- list(vec)[rep(1, nsim)]
    bgRate_i <- lapply(1:nsim, function(i) sapply(1:Time,
                                                  function(j) as.numeric(names(which.max(table(matrix(as.vector(select_mu[[i]]),
                                                                                                      ncol=Time)[,j]))))))
    bgRate <- lapply(1:nsim, function(i) rep(bgRate_i[[i]], each = ncentroids))
    #ix <- lapply(1:nsim, function(i) which(abs(log(as.vector(select_mu[[i]])) - log(bgRate[[i]]))>=10^-3))
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






