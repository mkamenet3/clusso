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
clust.diagnostics <- function(incluster, threshold, nullmod,...){
    if(is.null(nullmod)){
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
        
        #summary stats by IC
        ##alldetect
        ###AIC
        alldetect.summary.mean.aic <- mean(unlist(incluster$prop.alldetect.aic))
        alldetect.summary.median.aic <- median(unlist(incluster$prop.alldetect.aic))
        alldetect.summary.sd.aic <- sd(unlist(incluster$prop.alldetect.aic))
        
        ###AICC
        alldetect.summary.mean.aicc <- mean(unlist(incluster$prop.alldetect.aicc))
        alldetect.summary.median.aicc <- median(unlist(incluster$prop.alldetect.aicc))
        alldetect.summary.sd.aicc <- sd(unlist(incluster$prop.alldetect.aicc))
        
        ###BIC
        alldetect.summary.mean.bic <- mean(unlist(incluster$prop.alldetect.bic))
        alldetect.summary.median.bic <- median(unlist(incluster$prop.alldetect.bic))
        alldetect.summary.sd.bic <- sd(unlist(incluster$prop.alldetect.bic))
        
        ##potentialclusterdetect
        ###AIC
        potentialclusterdetect.summary.mean.aic <- mean(unlist(incluster$prop.wasdetect.aic))
        potentialclusterdetect.summary.median.aic <- median(unlist(incluster$prop.wasdetect.aic))
        potentialclusterdetect.summary.sd.aic <- sd(unlist(incluster$prop.wasdetect.aic))
        
        ###AICC
        potentialclusterdetect.summary.mean.aicc <- mean(unlist(incluster$prop.wasdetect.aicc))
        potentialclusterdetect.summary.median.aicc <- median(unlist(incluster$prop.wasdetect.aicc))
        potentialclusterdetect.summary.sd.aicc <- sd(unlist(incluster$prop.wasdetect.aicc))
        
        ###BIC
        potentialclusterdetect.summary.mean.bic <- mean(unlist(incluster$prop.wasdetect.bic))
        potentialclusterdetect.summary.median.bic <- median(unlist(incluster$prop.wasdetect.bic))
        potentialclusterdetect.summary.sd.bic <- sd(unlist(incluster$prop.wasdetect.bic))
        
        ##truclusterdetect
        ###AIC
        trueclusterdetect.summary.mean.aic <- mean(unlist(incluster$prop.shoulddetect.aic))
        trueclusterdetect.summary.median.aic <- median(unlist(incluster$prop.shoulddetect.aic))
        trueclusterdetect.summary.sd.aic <- sd(unlist(incluster$prop.shoulddetect.aic))
        
        ###AICC
        trueclusterdetect.summary.mean.aicc <- mean(unlist(incluster$prop.shoulddetect.aicc))
        trueclusterdetect.summary.median.aicc <- median(unlist(incluster$prop.shoulddetect.aicc))
        trueclusterdetect.summary.sd.aicc <- sd(unlist(incluster$prop.shoulddetect.aicc))
        
        ###BIC
        trueclusterdetect.summary.mean.bic <- mean(unlist(incluster$prop.shoulddetect.bic))
        trueclusterdetect.summary.median.bic <- median(unlist(incluster$prop.shoulddetect.bic))
        trueclusterdetect.summary.sd.bic <- sd(unlist(incluster$prop.shoulddetect.bic))

        return(list(incluster.any.aic = incluster$incluster.any.aic, incluster.any.aicc = incluster$incluster.any.aicc,
                        incluster.any.bic = incluster$incluster.any.bic,
                    outcluster.any.aic = incluster$outcluster.any.aic, outcluster.any.aicc = incluster$outcluster.any.aicc,
                        outcluster.any.bic= incluster$outcluster.any.bic, 
                    alldetect.aic = alldetect.aic, alldetect.aicc = alldetect.aicc, alldetect.bic = alldetect.bic,
                    potentialclusterdetect.aic = potentialclusterdetect.aic, potentialclusterdetect.aicc = potentialclusterdetect.aicc,
                     potentialclusterdetect.bic = potentialclusterdetect.bic,
                    trueclusterdetect.aic = trueclusterdetect.aic, trueclusterdetect.aicc = trueclusterdetect.aicc,
                        trueclusterdetect.bic = trueclusterdetect.bic,
                    alldetect.summary.mean.aic = alldetect.summary.mean.aic, alldetect.summary.mean.aicc = alldetect.summary.mean.aicc, 
                        alldetect.summary.mean.bic = alldetect.summary.mean.bic,
                    alldetect.summary.median.aic = alldetect.summary.median.aic, alldetect.summary.median.aicc = alldetect.summary.median.aicc,
                        alldetect.summary.median.bic = alldetect.summary.median.bic,
                    alldetect.summary.sd.aic = alldetect.summary.sd.aic, alldetect.summary.sd.aicc = alldetect.summary.sd.aicc, alldetect.summary.sd.bic = alldetect.summary.sd.bic,
                    
                    potentialclusterdetect.summary.mean.aic = potentialclusterdetect.summary.mean.aic, potentialclusterdetect.summary.mean.aicc = potentialclusterdetect.summary.mean.aicc,
                        potentialclusterdetect.summary.mean.bic = potentialclusterdetect.summary.mean.bic,
                    potentialclusterdetect.summary.median.aic = potentialclusterdetect.summary.median.aic, potentialclusterdetect.summary.median.aicc = potentialclusterdetect.summary.median.aicc,
                        potentialclusterdetect.summary.median.bic = potentialclusterdetect.summary.median.bic,
                    potentialclusterdetect.summary.sd.aic = potentialclusterdetect.summary.sd.aic, potentialclusterdetect.summary.sd.aicc = potentialclusterdetect.summary.sd.aicc,
                        potentialclusterdetect.summary.sd.bic = potentialclusterdetect.summary.sd.bic,
                    
                    trueclusterdetect.summary.mean.aic = trueclusterdetect.summary.mean.aic, trueclusterdetect.summary.mean.aicc = trueclusterdetect.summary.mean.aicc, 
                        trueclusterdetect.summary.mean.bic = trueclusterdetect.summary.mean.bic,
                    trueclusterdetect.summary.median.aic = trueclusterdetect.summary.median.aic, trueclusterdetect.summary.median.aicc = trueclusterdetect.summary.median.aicc, 
                        trueclusterdetect.summary.median.bic = trueclusterdetect.summary.median.bic,
                    trueclusterdetect.summary.sd.aic = trueclusterdetect.summary.sd.aic, trueclusterdetect.summary.sd.aicc = trueclusterdetect.summary.sd.aicc,
                    trueclusterdetect.summary.sd.bic = trueclusterdetect.summary.sd.bic))
    }
    else{
        message("Returning Diagnostics for Null Model")
        #prop.null.aic
        null.prop.aic <- lapply(1:nsim, function(i) if(any(incluster$null.aic[[i]] != 0)) null.prop.aic=1 else null.prop.aic = 0)
        prop.null.aic <- paste0((sum(unlist(null.prop.aic))/nsim)*100 , "%")
        
        #prop.null.aicc
        null.prop.aicc <- lapply(1:nsim, function(i) if(any(incluster$null.aicc[[i]] != 0)) null.prop.aicc=1 else null.prop.aicc = 0)
        prop.null.aicc <- paste0((sum(unlist(null.prop.aicc))/nsim)*100 , "%")
        
        #prop.null.bic
        null.prop.bic <- lapply(1:nsim, function(i) if(any(incluster$null.bic[[i]] != 0)) null.prop.bic=1 else null.prop.bic = 0)
        prop.null.bic <- paste0((sum(unlist(null.prop.bic))/nsim)*100 , "%")
        
        # #summary on numbser of cells detected
        # ##null.summary.aic
        # null.summary.mean.aic <- mean(unlist(incluster$null.aic))
        # null.summary.median.aic <- median(unlist(incluster$null.aic))
        # null.summary.sd.aic <- sd(unlist(incluster$null.aic))
        # 
        # #null.summary.aicc
        # null.summary.mean.aicc <- mean(unlist(incluster$null.aicc))
        # null.summary.median.aicc <- median(unlist(incluster$null.aicc))
        # null.summary.sd.aicc <- sd(unlist(incluster$null.aicc))
        # 
        # #null.summary.bic
        # null.summary.mean.bic <- mean(unlist(incluster$null.bic))
        # null.summary.median.bic <- median(unlist(incluster$null.bic))
        # null.summary.sd.bic <- sd(unlist(incluster$null.bic))
        
        
        return(list(#null.any.aic = incluster$null.any.aic, null.any.aicc = incluster$null.any.aicc, null.any.bic = incluster$null.any.bic,
                    null.summary.mean.aic = incluster$null.summary.mean.aic, null.summary.mean.aicc = incluster$null.summary.mean.aicc, null.summary.mean.bic = incluster$null.summary.mean.bic,
                    null.summary.median.aic = incluster$null.summary.median.aic, null.summary.median.aicc = incluster$null.summary.median.aicc, null.summary.median.bic = incluster$null.summary.median.bic,
                    null.summary.sd.aic = incluster$null.summary.sd.aic, null.summary.sd.aicc = incluster$null.summary.sd.aicc, null.summary.sd.bic = incluster$null.summary.sd.bic,
                    prop.null.aic = prop.null.aic, prop.null.aicc = prop.null.aicc, prop.null.bic = prop.null.bic))
    }
    
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
#'@param under Default is FALSE. If TRUE, then it will estimate detection based on a relative risk which is < 1. In other words, it will consider 
#'clusters to be identified where the estimated values are less than the background rate, not more than the background rate as is the case in the 
#'elevated relative risk models.
#'@param nullmod Default is FALSE. If TRUE, then it will estimate detection based on the null model where there is no cluster. 
detect.incluster.ic <- function(lassoresult, vectors.sim, rr, set, period, Time, nsim, under=FALSE,nullmod){
    prob.simBIC <- lapply(1:nsim, function(x) matrix(0, nrow(rr)*Time))
    prob.simAIC <- lapply(1:nsim, function(x) matrix(0, nrow(rr)*Time))
    prob.simAICc <- lapply(1:nsim, function(x) matrix(0, nrow(rr)*Time))
    print(nullmod)
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
    if(under==TRUE){
        #ix <- lapply(1:nsim, function(i) which(round(matrix(set$rr.simAIC[[i]],ncol=Time),6) < round(set$alphaAIC[[i]],6)))
        ix <- lapply(1:nsim, function(i) which(round(set$rr.simAIC[[i]],6) < round(as.vector(set$alphaAIC[[i]]),6)))
    }
    else{
        print("ok1")
        #print(str(set))
        ix <- lapply(1:nsim, function(i) which(round(set$rr.simAIC[[i]],6) > round(as.vector(set$alphaAIC[[i]]),6)))
        #ix <- lapply(1:nsim, function(i) which(round(as.vector(set$rr.simAIC[[i]]),6) > round(as.vector(set$alphaAIC[[i]]),6)))    
        print("ok2")
    }
    
    #1)a) Did it find anything in the cluster?
    
    clust <- lapply(1:nsim, function(i) is.element(ix[[i]],set$indx.clust.truth))
    in.clust.any <- lapply(1:nsim, function(i) if(any(clust[[i]]==TRUE)) in.clust.any=1 else in.clust.any = 0)
    incluster.any.aic <- paste0((sum(unlist(in.clust.any))/nsim)*100,"%")    
    
    
    #1)b) Did it find something OUTSIDE of the cluster?
    outclust <- lapply(1:nsim, function(i) setdiff(ix[[i]],set$indx.clust.truth))
    out.clust.any <- lapply(1:nsim, function(i) if(isTRUE(length(setdiff(ix[[i]],set$indx.clust.truth))==0) == FALSE) out.clust.any=1 else out.clust.any = 0)
    outcluster.any.aic <- paste0((sum(unlist(out.clust.any))/nsim)*100,"%")    
    
    
    #2) |(A and B)|/|A U B|?
    ##Calculate
    #wasDetected <- lapply(1:nsim, function(i) which(round(matrix(set$rr.simAIC[[i]],ncol=Time),6) > round(set$alphaAIC[[i]],6), arr.ind=FALSE))    
    wasDetected <- lapply(1:nsim, function(i) which(round(set$rr.simAIC[[i]],6) > round(as.vector(set$alphaAIC[[i]]),6)))
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
        #of sims that find something, how many cells on average are being detected
        index.nonzero <- unlist(lapply(1:nsim, function(i) if(isTRUE(length(ix[[i]])!=0)) index.nonzero = i))
        nonzero <- NULL
        for(i in index.nonzero){
            nonzero <- append(nonzero, length(ix[[i]]))
        }
        null.summary.mean.aic <- mean(nonzero)
        null.summary.median.aic <- median(nonzero)
        null.summary.sd.aic <- sd(nonzero)
        
        #propr detected
        null.aic <- lapply(1:nsim, function(i) length(unlist(ix[[i]])))    
        #null.any <- lapply(1:nsim, function(i) if(isTRUE(null.aic[[i]] == 0)) null.any = 0 else null.any = 1)
        #null.any.aic <- paste0((sum(unlist(null.any))/nsim)*100,"%")   
    }
    
    
    ########################################################################
    #(Q)AICc
    #extract things that are not the background rate
    if(under==TRUE){
        #ix <- lapply(1:nsim, function(i) which(round(matrix(set$rr.simAICc[[i]],ncol=Time),6) < round(set$alphaAICc[[i]],6)))
        ix <- lapply(1:nsim, function(i) which(round(set$rr.simAIC[[i]],6) < round(as.vector(set$alphaAIC[[i]]),6)))
    }
    else{
        #ix <- lapply(1:nsim, function(i) which(round(matrix(set$rr.simAICc[[i]],ncol=Time),6) > round(set$alphaAICc[[i]],6)))    
        ix <- lapply(1:nsim, function(i) which(round(set$rr.simAICc[[i]],6) > round(as.vector(set$alphaAICc[[i]]),6)))
    }
    
    #1)a) Did it find anything in the cluster?
    
    clust <- lapply(1:nsim, function(i) is.element(ix[[i]],set$indx.clust.truth))
    in.clust.any <- lapply(1:nsim, function(i) if(any(clust[[i]]==TRUE)) in.clust.any=1 else in.clust.any = 0)
    incluster.any.aicc <- paste0((sum(unlist(in.clust.any))/nsim)*100,"%")    
    
    
    #1)b) Did it find something OUTSIDE of the cluster?
    outclust <- lapply(1:nsim, function(i) setdiff(ix[[i]],set$indx.clust.truth))
    out.clust.any <- lapply(1:nsim, function(i) if(isTRUE(length(setdiff(ix[[i]],set$indx.clust.truth))==0) == FALSE) out.clust.any=1 else out.clust.any = 0)
    outcluster.any.aicc <- paste0((sum(unlist(out.clust.any))/nsim)*100,"%")    
    
    
    #2) |(A and B)|/|A U B|?
    ##Calculate
    #wasDetected <- lapply(1:nsim, function(i) which(round(matrix(set$rr.simAICc[[i]],ncol=Time),6) > round(set$alphaAICc[[i]],6), arr.ind=FALSE))    
    wasDetected <- lapply(1:nsim, function(i) which(round(set$rr.simAICc[[i]],6) > round(as.vector(set$alphaAICc[[i]]),6)))
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
        #of sims that find something, how many cells on average are being detected
        index.nonzero <- unlist(lapply(1:nsim, function(i) if(isTRUE(length(ix[[i]])!=0)) index.nonzero = i))
        nonzero <- NULL
        for(i in index.nonzero){
            nonzero <- append(nonzero, length(ix[[i]]))
        }
        null.summary.mean.aicc <- mean(nonzero)
        null.summary.median.aicc <- median(nonzero)
        null.summary.sd.aicc <- sd(nonzero)
        
        #propr detected
        null.aicc <- lapply(1:nsim, function(i) length(unlist(ix[[i]])))    
        #null.any <- lapply(1:nsim, function(i) if(isTRUE(null.aic[[i]] == 0)) null.any = 0 else null.any = 1)
    }
    ################################################################
    #(Q)BIC
    #extract things that are not the background rate
    if(under==TRUE){
        #ix <- lapply(1:nsim, function(i) which(round(matrix(set$rr.simBIC[[i]],ncol=Time),6) < round(set$alphaBIC[[i]],6)))
        ix <- lapply(1:nsim, function(i) which(round(set$rr.simAIC[[i]],6) < round(as.vector(set$alphaAIC[[i]]),6)))
    }
    else{
        #ix <- lapply(1:nsim, function(i) which(round(matrix(set$rr.simBIC[[i]],ncol=Time),6) > round(set$alphaBIC[[i]],6)))    
        ix <- lapply(1:nsim, function(i) which(round(set$rr.simBIC[[i]],6) > round(as.vector(set$alphaBIC[[i]]),6)))
    }
    
    #1)a) Did it find anything in the cluster?
    
    clust <- lapply(1:nsim, function(i) is.element(ix[[i]],set$indx.clust.truth))
    in.clust.any <- lapply(1:nsim, function(i) if(any(clust[[i]]==TRUE)) in.clust.any=1 else in.clust.any = 0)
    incluster.any.bic <- paste0((sum(unlist(in.clust.any))/nsim)*100,"%")    
    
    
    #1)b) Did it find something OUTSIDE of the cluster?
    outclust <- lapply(1:nsim, function(i) setdiff(ix[[i]],set$indx.clust.truth))
    out.clust.any <- lapply(1:nsim, function(i) if(isTRUE(length(setdiff(ix[[i]],set$indx.clust.truth))==0) == FALSE) out.clust.any=1 else out.clust.any = 0)
    outcluster.any.bic <- paste0((sum(unlist(out.clust.any))/nsim)*100,"%")    
    
    
    #2) |(A and B)|/|A U B|?
    ##Calculate
    #wasDetected <- lapply(1:nsim, function(i) which(round(matrix(set$rr.simBIC[[i]],ncol=Time),6) > round(set$alphaBIC[[i]],6), arr.ind=FALSE))    
    wasDetected <- lapply(1:nsim, function(i) which(round(set$rr.simBIC[[i]],6) > round(as.vector(set$alphaBIC[[i]]),6)))
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
        #of sims that find something, how many cells on average are being detected
        index.nonzero <- unlist(lapply(1:nsim, function(i) if(isTRUE(length(ix[[i]])!=0)) index.nonzero = i))
        nonzero <- NULL
        for(i in index.nonzero){
            nonzero <- append(nonzero, length(ix[[i]]))
        }
        null.summary.mean.bic <- mean(nonzero)
        null.summary.median.bic <- median(nonzero)
        null.summary.sd.bic <- sd(nonzero)
        
        #propr detected
        null.bic <- lapply(1:nsim, function(i) length(unlist(ix[[i]])))    
        #null.any <- lapply(1:nsim, function(i) if(isTRUE(null.aic[[i]] == 0)) null.any = 0 else null.any = 1)
        #null.any.aic <- paste0((sum(unlist(null.any))/nsim)*100,"%")   
    }
    
    
    #Returns
    if(exists("null.aic") & exists("null.bic") & exists("null.aicc")){
        return(list(
            null.aic = null.aic, null.aicc = null.aicc, null.bic = null.bic,
            null.summary.mean.aic = null.summary.mean.aic, null.summary.mean.aicc = null.summary.mean.aicc, null.summary.mean.bic = null.summary.mean.bic,
            null.summary.median.aic = null.summary.median.aic, null.summary.median.aicc = null.summary.median.aicc, null.summary.median.bic = null.summary.median.bic,
            null.summary.sd.aic = null.summary.sd.aic, null.summary.sd.aicc = null.summary.sd.aicc, null.summary.sd.bic = null.summary.sd.bic))
    }
    else{
        return(list(
            incluster.any.aic = incluster.any.aic, incluster.any.aicc = incluster.any.aicc,incluster.any.bic = incluster.any.bic,
            outcluster.any.aic = outcluster.any.aic, outcluster.any.aicc = outcluster.any.aicc, outcluster.any.bic = outcluster.any.bic,
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
detect.incluster <- function(lassoresult, vectors.sim, rr, set, timeperiod, Time, nsim, x, y, rMax, center, 
                             radius, IC = c("aic","aicc","bic","ic"),under=FALSE, nullmod){
    period = timeperiod
    message("Detection Results for:\n"
            , "\t Time Period: ", period,
            "\n \t Num. simulations: ", nsim,
            "\n \t Cluster center: ", center,
            "\n \t Cluster radius: ", radius,
            "\n \t Cluster rel.risk: ",unique(as.vector(rr))[2])
    IC <- match.arg(IC, several.ok= TRUE)
    switch(IC,
           aic = detect.incluster.aic(lassoresult, vectors.sim, rr, set, period, Time, nsim, under=FALSE),
           aicc = detect.incluster.aicc(lassoresult, vectors.sim, rr, set, period, Time, nsim, under=FALSE),
           bic = detect.incluster.bic(lassoresult, vectors.sim, rr, set, period, Time, nsim, under=FALSE),
           ic = detect.incluster.ic(lassoresult, vectors.sim, rr, set, period, Time, nsim, under=FALSE, nullmod))
} 


