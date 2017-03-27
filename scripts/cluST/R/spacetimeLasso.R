#' Spatial and Spatio-Temporal Lasso
#' 
#' This function runs the Lasso regularization technique on our large sparse matric of potential space or space-time clusters.
#' @param clusters clusters dataframe from (cluster.df function) that includes the center, x,y, r (radius), n (counter), and last (last observation in potential cluster)
#' @param vectors takes in the list of expected and observed counts from setVectors function
#' @param Time number of time periods in the dataset
#' @param spacetime indicator of whether the cluster detection method should be run on all space-time clusters(default) or on only the potential space clusters.
#' @param pois whether or not the Quasi-Poisson or Poisson model should be run. Default is pois=FALSE (default is Quasi-Poisson model is to be run)
#' @return This function will return a list with the expected counts as selected by QBIC, QAIC, QAICc, a list of original expected counts (Ex),
#' a list of observed counts (Yx), the lasso object, a list of K values (number of unique values in each decision path), and n (length of unique centers in the clusters dataframe)
#' @export
#' @example 
#' potentialclusters <- clusters2df(lat, long, utm=FALSE, length(lat))
#' myvectors <- setVectors(period, expected, observed, Time, byrow=TRUE)
#' myresults <- spacetimeLasso(potentialclusters, myvectors, spacetime=TRUE, pois=FALSE)
#' 
spacetimeLasso<- function(clusters, vectors, Time, spacetime=TRUE,pois=FALSE,...){
    n <- length(unique(clusters$center))
    potClus <- n
    numCenters <- n
    message("Creating space-time matrix")
    if(spacetime==TRUE){
        sparseMAT <- spacetimeMat(clusters, numCenters, Time)
    }
    else{
        sparseMAT <- spaceMat(clusters, numCenters)
    }
    message("Space-time matrix created")
    Ex <- vectors$E0
    Yx <- vectors$Y.vec
    Period <- vectors$Period
    message("Running Lasso - stay tuned")
    lasso <- glmnet(sparseMAT, Yx, family=("poisson"), alpha=1, offset=log(Ex), nlambda = 2000, standardize = FALSE, dfmax = 10)
    message("Lasso complete - extracting estimates and paths")
    coefs.lasso.all <- coef(lasso)
    intercept <- rep(1, dim(sparseMAT)[1])
    sparseMAT <- cBind(intercept, sparseMAT)
    xbetaPath<- sparseMAT%*%coefs.lasso.all
    mu <- sapply(1:length(lasso$lambda), function(i) exp(xbetaPath[,i]))    
    loglike <- sapply(1:length(lasso$lambda), function(i) sum(dpoisson(Yx, mu[,i],log=TRUE)))
    K <- lasso$df + 1
    if(spacetime==TRUE & pois == FALSE){
        message("returning results for space-time Quasi-Poisson model")
        offset_reg <- glm(Yx ~ 1 + as.factor(vectors$Period) + offset(log(Ex)),family=poisson)
        overdisp.est <- overdisp(offset_reg)
        message("Selecting best paths")
        
        #QBIC
        PLL.qbic  <- (loglike/overdisp.est)-log(n*Time)/2*K
        select.qbic <- which.max(PLL.qbic)
        E.qbic <- mu[,select.qbic]
        
        #QAIC
        PLL.qaic = (loglike/overdisp.est) - K
        select.qaic <- which.max(PLL.qaic)
        E.qaic <- mu[,select.qaic]
        
        #QAICc
        PLL.qaicc=(loglike/overdisp.est)- ((K*n*Time)/(n*Time-K-1))
        select.qaicc <- which.max(PLL.qaicc)
        E.qaicc <- mu[,select.qaic]
        message("Returning results")
    }
    #########################################################
    #Space-Time, Poisson only
    #########################################################                
    else if(spacetime==TRUE & pois==TRUE){
        message("returning results for space-time Poisson model")
        #QBIC
        PLL.qbic  <- (loglike)-log(n*Time)/2*K
        select.qbic <- which.max(PLL.qbic)
        E.qbic <- mu[,select.qbic]
        
        #QAIC
        PLL.qaic = (loglike) - K
        select.qaic <- which.max(PLL.qaic)
        E.qaic <- mu[,select.qaic]
        
        #QAICc
        PLL.qaicc=(loglike)- ((K*n*Time)/(n*Time-K-1))
        select.qaicc <- which.max(PLL.qaicc)
        E.qaicc <- mu[,select.qaic]
        message("Returning results")
    }
    #########################################################
    #Space-Only, Quasi-Poisson
    #########################################################
    else if(spacetime==FALSE & pois==FALSE){
        message("Returning results for space-only  Quasi-Poisson model")
        offset_reg <- glm(Yx ~ 1 + offset(log(Ex)),family=poisson)
        overdisp.est <- overdisp(offset_reg)
        message("Selecting best paths")
        
        #QBIC
        PLL.qbic  <- (loglike/overdisp.est)-log(n*Time)/2*K
        select.qbic <- which.max(PLL.qbic)
        E.qbic <- mu[,select.qbic]
        
        #QAIC
        PLL.qaic = (loglike/overdisp.est) - K
        select.qaic <- which.max(PLL.qaic)
        E.qaic <- mu[,select.qaic]
        
        #QAICc
        PLL.qaicc=(loglike/overdisp.est)- ((K*n*Time)/(n*Time-K-1))
        select.qaicc <- which.max(PLL.qaicc)
        E.qaicc <- mu[,select.qaic]
        message("Returning results")
    }
    #########################################################
    #Space-only, Poisson only
    #########################################################
    else if(spacetime==FALSE & pois == TRUE){
        message("Returning results for space-only  Poisson model")
        #QBIC
        PLL.qbic  <- (loglike)-log(n*Time)/2*K
        select.qbic <- which.max(PLL.qbic)
        E.qbic <- mu[,select.qbic]
        
        #QAIC
        PLL.qaic = (loglike) - K
        select.qaic <- which.max(PLL.qaic)
        E.qaic <- mu[,select.qaic]
        
        #QAICc
        PLL.qaicc=(loglike)- ((K*n*Time)/(n*Time-K-1))
        select.qaicc <- which.max(PLL.qaicc)
        E.qaicc <- mu[,select.qaic]
        message("Returning results")
    }
    
    return(list(E.qbic = E.qbic, E.qaic = E.qaic, E.qaicc = E.qaicc, Ex = Ex, Yx = Yx, lasso = lasso, K = K, n = n))  
}

