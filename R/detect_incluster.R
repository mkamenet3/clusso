
#'@title
#'detect_incluster_ic
#'@description
#'This function will calculate detection based on the three information criterion. You can run detection on a null model (no cluster),
#'a model with an under-estimated cluster (artificial cluster <1), and on an elevated relative risk cluster. This is called by the more general
#'function *detect_incluster* which will switch to this function (TODO - Allow for detection options for AIC, AICc, and BIC only). 
#'@param lassoresult List of QBIC, QAIC, QAICc estimates from the mylasso.sim function
#'@param vectors.sim  dataframe of initial vectors of the observed and expected counts that went into simulation function
#'@param rr risk ratio matrix that was used in the simulation
#'@param set output from detet_set function
#'@param timeperiod Time period span of the cluster
#'@param Time number of time periods total in the model
#'@param nsim number of simulations run
#'@param under Default is FALSE. If TRUE, then it will estimate detection based on a relative risk which is < 1. In other words, it will consider 
#'clusters to be identified where the estimated values are less than the background rate, not more than the background rate as is the case in the 
#'elevated relative risk models.
#'@param nullmod Default is NULL. If TRUE, then it will estimate detection based on the null model where there is no cluster. 
detect_incluster_ic <- function(lassoresult, vectors.sim, rr, set, timeperiod, Time, nsim, under=FALSE,nullmod, risk.ratio,x,y,rMax,thresh){
    if(is.null(nullmod)){
        message("Returning results for simulation model")
        print(thresh)
        ############################
        ##Prob in/out cluster function as function of all potential clusters
        if(is.null(thresh)){
            clusterdetectionrates <- prob_inoutcluster(lassoresult,rr,risk.ratio,x,y,rMax,nsim,Time,thresh)
        }
        else{
            if(length(thresh)>1){
                print("yah1")
                clusterdetectionrates <- lapply(1:length(thresh), function(i) prob_inoutcluster(lassoresult,rr,risk.ratio,x,y,rMax,nsim,Time,thresh[[i]]))
            }
            else{
                print("yah2")
                clusterdetectionrates <- prob_inoutcluster(lassoresult,rr,risk.ratio,x,y,rMax,nsim,Time,thresh)    
            }
            
        }
        
        
        #Prob in/out for individual cells
        ############################
        ##AIC
        ############################
        #extract things that are not the background rate
        if(under==TRUE){
            if(!is.null(under)) stop("Specify if you want under or null estimates")
            ix <- lapply(1:nsim, function(i) which(round(set$rr.simAIC[[i]],6) < round(as.vector(set$alphaAIC[[i]]),6)))
            message("Running < 1 risk model")
        }
        else{
            ix <- lapply(1:nsim, function(i) which(round(set$rr.simAIC[[i]],6) > round(as.vector(set$alphaAIC[[i]]),6)))
            stopifnot(all.equal(ix, lapply(1:nsim, function(i) which(round(as.vector(set$rr.simAIC[[i]]),6) > round(as.vector(set$alphaAIC[[i]]),6)))))
        }
        #1)a) Did it find anything in the cluster?
        clust <- lapply(1:nsim, function(i) is.element(ix[[i]],set$indx.clust.truth))
        in.clust.any <- lapply(1:nsim, function(i) if(any(clust[[i]]==TRUE)) in.clust.any=1 else in.clust.any = 0)
        incluster.any.aic <- paste0((sum(unlist(in.clust.any))/nsim)*100,"%")    
        
        
        #1)b) Did it find something OUTSIDE of the cluster?
        outclust <- lapply(1:nsim, function(i) setdiff(ix[[i]],set$indx.clust.truth))
        out.clust.any <- lapply(1:nsim, function(i) if(isTRUE(length(setdiff(ix[[i]],set$indx.clust.truth))==0) == FALSE) out.clust.any=1 else out.clust.any = 0)
        outcluster.any.aic <- paste0((sum(unlist(out.clust.any))/nsim)*100,"%")    
        
        ############################
        ##AICc
        ############################
        if(under==TRUE){
            ix <- lapply(1:nsim, function(i) which(round(set$rr.simAICc[[i]],6) < round(as.vector(set$alphaAICc[[i]]),6)))
        }
        else{
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
        
        ############################
        ##BIC
        ############################
        if(under==TRUE){
            ix <- lapply(1:nsim, function(i) which(round(set$rr.simBIC[[i]],6) < round(as.vector(set$alphaBIC[[i]]),6)))
        }
        else{
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
          
        # if(is.null(thresh)){
        #     print("yapppe")
        #     return(list(
        #     
        #         incluster.any.aic = incluster.any.aic, incluster.any.aicc = incluster.any.aicc,incluster.any.bic = incluster.any.bic,
        #         outcluster.any.aic = outcluster.any.aic, outcluster.any.aicc = outcluster.any.aicc, outcluster.any.bic = outcluster.any.bic,
        #         notinperc.aic = clusterdetectionrates$notinperc.aic, notinperc.aicc = clusterdetectionrates$notinperc.aicc, notinperc.bic = clusterdetectionrates$notinperc.bic,
        #         inperc.aic = clusterdetectionrates$inperc.aic, inperc.aicc = clusterdetectionrates$inperc.aicc, inperc.bic = clusterdetectionrates$inperc.bic))    
        # }
        print(length(thresh))
        if (length(thresh) == 1){
            print("yeck")
            return(list(
                incluster.any.aic = incluster.any.aic, incluster.any.aicc = incluster.any.aicc,incluster.any.bic = incluster.any.bic,
                outcluster.any.aic = outcluster.any.aic, outcluster.any.aicc = outcluster.any.aicc, outcluster.any.bic = outcluster.any.bic,
                notinperc.aic = clusterdetectionrates$notinperc.aic, notinperc.aicc = clusterdetectionrates$notinperc.aicc, notinperc.bic = clusterdetectionrates$notinperc.bic,
                inperc.aic = clusterdetectionrates$inperc.aic, inperc.aicc = clusterdetectionrates$inperc.aicc, inperc.bic = clusterdetectionrates$inperc.bic,
                threshresults = clusterdetectionrates$thresh)) 
        }
        else if (length(thresh)>1){
            print("yeppers")
            return(list(
                incluster.any.aic = incluster.any.aic, incluster.any.aicc = incluster.any.aicc,incluster.any.bic = incluster.any.bic,
                outcluster.any.aic = outcluster.any.aic, outcluster.any.aicc = outcluster.any.aicc, outcluster.any.bic = outcluster.any.bic,
                notinperc.aic = clusterdetectionrates[[1]]$notinperc.aic, notinperc.aicc = clusterdetectionrates[[1]]$notinperc.aicc, notinperc.bic = clusterdetectionrates[[1]]$notinperc.bic,
                inperc.aic = clusterdetectionrates[[1]]$inperc.aic, inperc.aicc = clusterdetectionrates[[1]]$inperc.aicc, inperc.bic = clusterdetectionrates[[1]]$inperc.bic,
                threshresults = matrix(unlist(lapply(1:length(thresh), function(i) clusterdetectionrates[[i]]$thresh)),nrow=length(thresh),byrow=TRUE)))
        }
        else{
            print("yapppe")
            return(list(
                
                incluster.any.aic = incluster.any.aic, incluster.any.aicc = incluster.any.aicc,incluster.any.bic = incluster.any.bic,
                outcluster.any.aic = outcluster.any.aic, outcluster.any.aicc = outcluster.any.aicc, outcluster.any.bic = outcluster.any.bic,
                notinperc.aic = clusterdetectionrates$notinperc.aic, notinperc.aicc = clusterdetectionrates$notinperc.aicc, notinperc.bic = clusterdetectionrates$notinperc.bic,
                inperc.aic = clusterdetectionrates$inperc.aic, inperc.aicc = clusterdetectionrates$inperc.aicc, inperc.bic = clusterdetectionrates$inperc.bic))
        }
    }
    ###NULL Diagnostics
    else{
        message("Returning Diagnostics for Null Model")
        ########
        #AIC
        ########
        if(under==TRUE){
            if(!is.null(under)) stop("Specify if you want under or null estimates")
            ix <- lapply(1:nsim, function(i) which(round(set$rr.simAIC[[i]],6) < round(as.vector(set$alphaAIC[[i]]),6)))
            message("Running < 1 risk model")
        }
        else{
            ix <- lapply(1:nsim, function(i) which(round(set$rr.simAIC[[i]],6) > round(as.vector(set$alphaAIC[[i]]),6)))
        }
        null.aic <- lapply(1:nsim, function(i) length(unlist(ix[[i]]))) 
        #prop.null.aic
        null.prop.aic <- lapply(1:nsim, function(i) if(any(null.aic[[i]] != 0)) null.prop.aic=1 else null.prop.aic = 0)
        prop.null.aic <- paste0((sum(unlist(null.prop.aic))/nsim)*100 , "%")
        
        
        ########
        #AICc
        ########
        if(under==TRUE){
            if(!is.null(under)) stop("Specify if you want under or null estimates")
            ix <- lapply(1:nsim, function(i) which(round(set$rr.simAICc[[i]],6) < round(as.vector(set$alphaAICc[[i]]),6)))
            message("Running < 1 risk model")
        }
        else{
            ix <- lapply(1:nsim, function(i) which(round(set$rr.simAICc[[i]],6) > round(as.vector(set$alphaAICc[[i]]),6)))
        }
        null.aicc <- lapply(1:nsim, function(i) length(unlist(ix[[i]]))) 
        #prop.null.aicc
        null.prop.aicc <- lapply(1:nsim, function(i) if(any(null.aicc[[i]] != 0)) null.prop.aicc=1 else null.prop.aicc = 0)
        prop.null.aicc <- paste0((sum(unlist(null.prop.aicc))/nsim)*100 , "%")
        
        ########
        #BIC
        ########
        if(under==TRUE){
            if(!is.null(under)) stop("Specify if you want under or null estimates")
            ix <- lapply(1:nsim, function(i) which(round(set$rr.simBIC[[i]],6) < round(as.vector(set$alphaBIC[[i]]),6)))
            message("Running < 1 risk model")
        }
        else{
            ix <- lapply(1:nsim, function(i) which(round(set$rr.simBIC[[i]],6) > round(as.vector(set$alphaBIC[[i]]),6)))
        }
        null.bic <- lapply(1:nsim, function(i) length(unlist(ix[[i]]))) 
        #prop.null.bic
        null.prop.bic <- lapply(1:nsim, function(i) if(any(null.bic[[i]] != 0)) null.prop.bic=1 else null.prop.bic = 0)
        prop.null.bic <- paste0((sum(unlist(null.prop.bic))/nsim)*100 , "%")
        
        return(list(
            prop.null.aic = prop.null.aic, prop.null.aicc = prop.null.aicc, prop.null.bic = prop.null.bic))
    }
    
}
    
    
    
#'detect_incluster
#'
#'This function will calculate the percent of simulations which correctly identify elements in cluster based on (Q)AIC, (Q)AICc, and (Q)BIC. The user can specify
#'if they want to only return one of these criterion or all three for further analysis.
#'@param lassoresult List of QBIC, QAIC, QAICc estimates from the mylasso.sim function
#'@param vectors.sim  dataframe of initial vectors of the observed and expected counts that went into simulation function
#'@param rr risk ratio matrix that was used in the simulation
#'@param set result of detect_set function
#'@param timeperiod span of timeperiods in simulation
#'@param Time total number of time periods
#'@param nsim number of simulations
#'@param x x coordinates
#'@param y y coordinates
#'@param rMax max radius
#'@param center center of centroid for cluster
#'@param radius radius of the cluster
#'@param IC the information criteria you would like to be returned. Options are: IC = aic or IC = qaic; IC = aicc or IC = qaicc; IC = bic or IC = qbic; IC = ic or IC = qic 
#'(for aic/aicc/bic and qaic/qaicc/qbic, respectively).
#'@param under default is FALSE. If risk.ratio is less than one (under-risk)
#'@param nullmod default is NULL. If not null, then null model results will be estimated and returned.
#'@return returns
detect_incluster <- function(lassoresult, vectors.sim, rr, set, timeperiod, Time, nsim, x, y, rMax, center, 
                             radius, under=FALSE, nullmod, risk.ratio,thresh){
    #period = timeperiod
    message("Detection Results for:\n"
            , "\t Time Period: ", timeperiod,
            "\n \t Num. simulations: ", nsim,
            "\n \t Cluster center: ", center,
            "\n \t Cluster radius: ", radius,
            "\n \t Cluster rel.risk: ",ifelse(length(unique(as.vector(rr)))==1,unique(as.vector(rr))[1],unique(as.vector(rr))[2]))
    ic <- detect_incluster_ic(lassoresult, vectors.sim, rr, set, timeperiod, Time, nsim, under=FALSE, nullmod, risk.ratio,x,y,rMax,thresh)
    return(ic)

} 


