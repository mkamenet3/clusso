#' Suite of function for cluster detection diagnostics - incluster diagnostics
#'
#'clust.diagnostics
#'@param incluster inherited object from detect.incluster.ic
#'@param threshold detection threshold set by user
#'@return list of detection
#'@example 
#'set <- detect.set(lassoresult.qp.s, vectors.sim, rr, Time, x, y, rMax, center, radius)
#'incluster.qp.s <- detect.incluster(lassoresult.qp.s, vectors.sim, rr, set, timeperiod, Time, nsim, x, y, rMax, center, radius, IC = "ic")
#'detect.qp.s <- list(clust.diagnostics(incluster.qp.s, threshold[1]), clust.diagnostics(incluster.qp.s , threshold[2]))
#'@export
#'
#'
clust.diagnostics <- function(incluster, threshold){
    #thresholding of prop.alldetect
    alldetect.aic <- paste0((length(which(unlist(incluster$prop.alldetect.aic)> threshold))/nsim)*100, "%")
    alldetect.aicc <- paste0((length(which(unlist(incluster$prop.alldetect.aicc)> threshold))/nsim)*100, "%")
    alldetect.bic <- paste0((length(which(unlist(incluster$prop.alldetect.bic)> threshold))/nsim)*100, "%")
    
    #thresholding of prop.wasdetect
    potentialclusterdetect.aic <- paste0((length(which(unlist(incluster$prop.wasdetect.aic)> threshold))/nsim)*100, "%")
    potentialclusterdetect.aicc <- paste0((length(which(unlist(incluster$prop.wasdetect.aicc)> threshold))/nsim)*100, "%")
    potentialclusterdetect.bic <- paste0((length(which(unlist(incluster$prop.wasdetect.bic)> threshold))/nsim)*100, "%")
    
    #thresholding of prop.shoulddetect
    trueclusterdetect.aic <- paste0((length(which(unlist(incluster$prop.shoulddetect.aic)> threshold))/nsim)*100, "%")
    trueclusterdetect.aicc <- paste0((length(which(unlist(incluster$prop.shoulddetect.aicc)> threshold))/nsim)*100, "%")
    trueclusterdetect.bic <- paste0((length(which(unlist(incluster$prop.shoulddetect.bic)> threshold))/nsim)*100, "%")
    
    return(list(incluster.any.aic = incluster$incluster.any.aic, incluster.any.aicc = incluster$incluster.any.aicc,
                incluster.any.bic = incluster$incluster.any.bic,
                alldetect.aic = alldetect.aic, alldetect.aicc = alldetect.aicc, alldetect.bic = alldetect.bic,
                potentialclusterdetect.aic = potentialclusterdetect.aic, potentialclusterdetect.aicc = potentialclusterdetect.aicc,
                potentialclusterdetect.bic = potentialclusterdetect.bic,
                trueclusterdetect.aic = trueclusterdetect.aic, trueclusterdetect.aicc = trueclusterdetect.aicc,
                trueclusterdetect.bic = trueclusterdetect.bic
    ))
}




#'detect.incluster.ic
#'This function will calculate detection based on the three information criterion. You can run detection on a null model (no cluster),
#'a model with an under-estimated cluster (artificial cluster <1), and on an elevated relative risk cluster. This is called by the more general
#'function *detect.incluster* which will switch to this function (TODO - Allow for detection options for AIC, AICc, and BIC only). 
#'@param lassoresult List of QBIC, QAIC, QAICc estimates from the mylasso.sim function
#'@param vectors.sim  dataframe of initial vectors of the observed and expected counts that went into simulation function
#'@param rr risk ratio matrix that was used in the simulation
#'@param res result of detect_set function
#'@param timeperiod Time period span of the cluster
#'@param Time number of time periods total in the model
#'@param nsim number of simulations run
#'@param under Default is NULL. If not null, then it will estimate detection based on a relative risk which is < 1. In other words, it will consider 
#'clusters to be identified where the estimated values are less than the background rate, not more than the background rate as is the case in the 
#'elevated relative risk models.
#'@param nullmod Default is NULL. If not null, then it will estimate detection based on the null model where there is no cluster. 
detect.incluster.ic <- function(lassoresult, vectors.sim, rr, set, period, Time, nsim, under=NULL,nullmod = NULL,...){
    prob.simBIC <- lapply(1:nsim, function(x) matrix(0, nrow(rr)*Time))
    prob.simAIC <- lapply(1:nsim, function(x) matrix(0, nrow(rr)*Time))
    prob.simAICc <- lapply(1:nsim, function(x) matrix(0, nrow(rr)*Time))
    #print(period)
    if(tail(period, n=1) == Time | tail(period, n=1)==1){
        maxTime = tail(period, n=1)    
    }
    else {
        maxTime = tail(period, n=1) +1
    }
    if(period[1] == 1){
        minTime = 1
    }
    else{
        minTime = period[1]-1
    }
    #(Q)AIC
    #extract things that are not the background rate
    if(!is.null(under)){
        ix <- lapply(1:nsim, function(i) which(round(matrix(set$rr.simAIC[[i]],ncol=Time),6) < round(set$alphaAIC[[i]],6), arr.ind=TRUE))
    }
    else{
        ix <- lapply(1:nsim, function(i) which(round(matrix(set$rr.simAIC[[i]],ncol=Time),6) > round(set$alphaAIC[[i]],6), arr.ind=TRUE))    
    }
    
    #1) Did it find anything in the cluster?
    ##NULL MODEL
    if(!is.null(nullmod)){
        clust <- lapply(1:nsim, function(i) is.element(ix[[i]][,1],set$neighs$cluster))
        in.clust.any <- lapply(1:nsim, function(i) if(isTRUE(length(clust[[i]])==0)) in.clust.any=0 else in.clust.any = 1)
        incluster.any.aic <- paste0((sum(unlist(in.clust.any))/nsim)*100,"%")
    }
    ##CLUSTER MODEL
    else{
        clust <- lapply(1:nsim, function(i) is.element(ix[[i]][,1],set$neighs$cluster))
        in.clust.any <- lapply(1:nsim, function(i) if(any(clust[[i]]==TRUE)) in.clust.any=1 else in.clust.any = 0)
        incluster.any.aic <- paste0((sum(unlist(in.clust.any))/nsim)*100,"%")    
    }
    
    #2) |(A and B)|/|A U B|?
    ##Calculate
    wasDetected <- lapply(1:nsim, function(i) which(round(matrix(set$rr.simAIC[[i]],ncol=Time),6) > round(set$alphaAIC[[i]],6), arr.ind=FALSE))    
    shouldDetected <- which(rr!=1, arr.ind=FALSE)
    ##Numerator
    AandB <- lapply(1:nsim, function(i) length(intersect(wasDetected[[i]], shouldDetected)))
    ##Denominator
    AuB <- lapply(1:nsim, function(i) length(union(wasDetected[[i]], shouldDetected)))
    ##Divide
    prop.alldetect.aic <- lapply(1:nsim, function(i) AandB[[i]]/AuB[[i]])
    
    
    #3) |(A and B)|/|A|? Proportion of what was detected was in overlap?
    ##Calculate length of what was detected
    A <- lapply(1:nsim, function(i) length(wasDetected[[i]]))
    ##Divide
    prop.wasdetect.aic <- lapply(1:nsim, function(i) AandB[[i]]/A[[i]])
    
    #4) |(A and B)|/|B|? Proportion of what should be detected was in overlap?
    B <- length(shouldDetected)
    prop.shoulddetect.aic <- lapply(1:nsim, function(i) AandB[[i]]/B)
    
    #5) ONLY FOR NULL MODEL - DID IT FIND ANYTHING?
    if(!is.null(nullmod)){
        null.any.aic <- length(unlist(ix))    
    }
    
    
    ########################################################################
    #(Q)AICc
    #extract things that are not the background rate
    if(!is.null(under)){
        ix <- lapply(1:nsim, function(i) which(round(matrix(set$rr.simAICc[[i]],ncol=Time),6) < round(set$alphaAICc[[i]],6), arr.ind=TRUE))
    }
    else{
        ix <- lapply(1:nsim, function(i) which(round(matrix(set$rr.simAICc[[i]],ncol=Time),6) > round(set$alphaAICc[[i]],6), arr.ind=TRUE))    
    }
    
    
    #1) Did it find anything in the cluster?
    ##NULL MODEL
    if(!is.null(nullmod)){
        clust <- lapply(1:nsim, function(i) is.element(ix[[i]][,1],set$neighs$cluster))
        in.clust.any <- lapply(1:nsim, function(i) if(isTRUE(length(clust[[i]])==0)) in.clust.any=0 else in.clust.any = 1)
        incluster.any.aicc <- paste0((sum(unlist(in.clust.any))/nsim)*100,"%")
    }
    else{
        clust <- lapply(1:nsim, function(i) is.element(ix[[i]][,1],set$neighs$cluster))
        in.clust.any <- lapply(1:nsim, function(i) if(any(clust[[i]]==TRUE)) in.clust.any=1 else in.clust.any = 0)
        incluster.any.aicc <- paste0((sum(unlist(in.clust.any))/nsim)*100,"%")    
    }
    #2) |(A and B)|/|A U B|?
    ##Calculate
    wasDetected <- lapply(1:nsim, function(i) which(round(matrix(set$rr.simAICc[[i]],ncol=Time),6) > round(set$alphaAICc[[i]],6), arr.ind=FALSE))    
    shouldDetected <- which(rr!=1, arr.ind=FALSE)
    ##Numerator
    AandB <- lapply(1:nsim, function(i) length(intersect(wasDetected[[i]], shouldDetected)))
    ##Denominator
    AuB <- lapply(1:nsim, function(i) length(union(wasDetected[[i]], shouldDetected)))
    ##Divide
    prop.alldetect.aicc <- lapply(1:nsim, function(i) AandB[[i]]/AuB[[i]])
    
    
    #3) |(A and B)|/|A|? Proportion of what was detected was in overlap?
    ##Calculate length of what was detected
    A <- lapply(1:nsim, function(i) length(wasDetected[[i]]))
    ##Divide
    prop.wasdetect.aicc <- lapply(1:nsim, function(i) AandB[[i]]/A[[i]])
    
    #4) |(A and B)|/|B|? Proportion of what should be detected was in overlap?
    B <- length(shouldDetected)
    prop.shoulddetect.aicc <- lapply(1:nsim, function(i) AandB[[i]]/B)
    
    #5) ONLY FOR NULL MODEL - DID IT FIND ANYTHING?
    if(!is.null(nullmod)){
        null.any.aicc <- length(unlist(ix))    
    }
    
    ################################################################
    #(Q)BIC
    #extract things that are not the background rate
    if(!is.null(under)){
        ix <- lapply(1:nsim, function(i) which(round(matrix(set$rr.simBIC[[i]],ncol=Time),6) < round(set$alphaBIC[[i]],6), arr.ind=TRUE))
    }
    else{
        ix <- lapply(1:nsim, function(i) which(round(matrix(set$rr.simBIC[[i]],ncol=Time),6) > round(set$alphaBIC[[i]],6), arr.ind=TRUE))    
    }
    
    #1) Did it find anything in the cluster?
    ##NULL MODEL
    if(!is.null(nullmod)){
        clust <- lapply(1:nsim, function(i) is.element(ix[[i]][,1],set$neighs$cluster))
        in.clust.any <- lapply(1:nsim, function(i) if(isTRUE(length(clust[[i]])==0)) in.clust.any=0 else in.clust.any = 1)
        incluster.any.bic <- paste0((sum(unlist(in.clust.any))/nsim)*100,"%")
    }
    else{
        clust <- lapply(1:nsim, function(i) is.element(ix[[i]][,1],set$neighs$cluster))
        in.clust.any <- lapply(1:nsim, function(i) if(any(clust[[i]]==TRUE)) in.clust.any=1 else in.clust.any = 0)
        incluster.any.bic <- paste0((sum(unlist(in.clust.any))/nsim)*100,"%")    
    }
    #2) |(A and B)|/|A U B|?
    ##Calculate
    wasDetected <- lapply(1:nsim, function(i) which(round(matrix(set$rr.simBIC[[i]],ncol=Time),6) > round(set$alphaBIC[[i]],6), arr.ind=FALSE))    
    shouldDetected <- which(rr!=1, arr.ind=FALSE)
    ##Numerator
    AandB <- lapply(1:nsim, function(i) length(intersect(wasDetected[[i]], shouldDetected)))
    ##Denominator
    AuB <- lapply(1:nsim, function(i) length(union(wasDetected[[i]], shouldDetected)))
    ##Divide
    prop.alldetect.bic <- lapply(1:nsim, function(i) AandB[[i]]/AuB[[i]])
    
    
    #3) |(A and B)|/|A|? Proportion of what was detected was in overlap?
    ##Calculate length of what was detected
    A <- lapply(1:nsim, function(i) length(wasDetected[[i]]))
    ##Divide
    prop.wasdetect.bic <- lapply(1:nsim, function(i) AandB[[i]]/A[[i]])
    
    #4) |(A and B)|/|B|? Proportion of what should be detected was in overlap?
    B <- length(shouldDetected)
    prop.shoulddetect.bic <- lapply(1:nsim, function(i) AandB[[i]]/B)
    
    #5) ONLY FOR NULL MODEL - DID IT FIND ANYTHING?
    if(!is.null(nullmod)){
        null.any.bic <- length(unlist(ix))    
    }
    if(exists("null.any.aic") & exists("null.any.bic") & exists("null.any.aicc")){
        return(list(
            incluster.any.aic = incluster.any.aic, incluster.any.aicc = incluster.any.aicc,incluster.any.bic = incluster.any.bic,
            prop.alldetect.aic = prop.alldetect.aic, prop.alldetect.aicc = prop.alldetect.aicc, prop.alldetect.bic = prop.alldetect.bic,
            prop.wasdetect.aic = prop.wasdetect.aic, prop.wasdetect.aicc = prop.wasdetect.aicc, prop.wasdetect.bic = prop.wasdetect.bic,
            prop.shoulddetect.aic = prop.shoulddetect.aic, prop.shoulddetect.aicc = prop.shoulddetect.aicc, prop.shoulddetect.bic = prop.shoulddetect.bic,
            null.any.aic = null.any.aic, null.any.aicc = null.any.aicc, null.any.bic = null.any.bic ))
    }
    else{
        return(list(
            incluster.any.aic = incluster.any.aic, incluster.any.aicc = incluster.any.aicc,incluster.any.bic = incluster.any.bic,
            prop.alldetect.aic = prop.alldetect.aic, prop.alldetect.aicc = prop.alldetect.aicc, prop.alldetect.bic = prop.alldetect.bic,
            prop.wasdetect.aic = prop.wasdetect.aic, prop.wasdetect.aicc = prop.wasdetect.aicc, prop.wasdetect.bic = prop.wasdetect.bic,
            prop.shoulddetect.aic = prop.shoulddetect.aic, prop.shoulddetect.aicc = prop.shoulddetect.aicc, prop.shoulddetect.bic = prop.shoulddetect.bic))
    }
    
}




#'detect.incluster
#'
#'This function will calculate the percent of simulations which correctly identify elements in cluster based on (Q)AIC, (Q)AICc, and (Q)BIC. The user can specify
#'if they want to only return one of these criterion or all three for further analysis.
#'@param lassoresult List of QBIC, QAIC, QAICc estimates from the mylasso.sim function
#'@param vectors.sim  dataframe of initial vectors of the observed and expected counts that went into simulation function
#'@param rr risk ratio matrix that was used in the simulation
#'@param res result of detect_set function
#'@param period_start time period where the cluster  starts in the simulation
#'@param period_end time period where cluster ends in the simulation
#'@param multi_period FALSE by default meaning that period_start and period_end are two unique periods. For example, if period_start = 2 and period_end =5, then
#'we will only be looking at periods 2 and 5. If multi_period is TRUE, then we will instead consider period_start through period_end (period_start:period_end). Following the same example,
#'this would mean we look at periods 2, 3, 4, and 5.
#'@param IC the information criteria you would like to be returned. Options are: IC = aic or IC = qaic; IC = aicc or IC = qaicc; IC = bic or IC = qbic; IC = ic or IC = qic 
#'(for aic/aicc/bic and qaic/qaicc/qbic, respectively).
#'TODO add cluster detection in case where risk ratio is less than background rate
#'@return returns
detect.incluster <- function(lassoresult, vectors.sim, rr, set,timeperiod, Time, nsim, x, y, rMax, center, 
                             radius, IC = c("aic","aicc","bic","ic"),nullmod=NULL,...){
    period = timeperiod
    message("Detection Results for:\n"
            , "\t Time Period: ", period,
            "\n \t Num. simulations: ", nsim,
            "\n \t Cluster center: ", center,
            "\n \t Cluster radius: ", radius,
            "\n \t Cluster rel.risk: ",unique(as.vector(rr))[2])
    IC <- match.arg(IC, several.ok= TRUE)
    switch(IC,
           aic = detect.incluster.aic(lassoresult, vectors.sim, rr, set, period, Time, nsim, under=NULL),
           aicc = detect.incluster.aicc(lassoresult, vectors.sim, rr, set, period, Time, nsim, under=NULL),
           bic = detect.incluster.bic(lassoresult, vectors.sim, rr, set, period, Time, nsim, under=NULL),
           ic = detect.incluster.ic(lassoresult, vectors.sim, rr, set, period, Time, nsim, under=NULL, nullmod=NULL))
} 


