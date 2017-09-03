#' Simulations forSpatial and Spatio-Temporal Lasso
#' 
#' This function runs the Lasso regularization technique on our large sparse matric of potential space-time clusters.It is specifically created to use with simulations.
#' @param clusters clusters dataframe from (cluster.df function) that includes the center, x,y, r (radius), n (counter), and last (last observation in potential cluster)
#' @param vectors.sim takes in the list of expected and observed counts from setVectors function
#' @param Time number of time periods in the dataset
#' @param spacetime indicator of whether the cluster detection method should be run on all space-time clusters(default) or on only the potential space clusters.
#' @param pois whether or not the Quasi-Poisson or Poisson model should be run. Default is pois=FALSE (default is Quasi-Poisson model is to be run)
#' @param nsim number of simulations
#' @param YSIM vector of simulated observed counts
#' @param overdispfloor default is TRUE. When TRUE, it limits phi (overdispersion parameter) to be greater or equal to 1. If FALSE, will allow for under dispersion.
#' @return This function will return a list with the expected counts as selected by QBIC, QAIC, QAICc, a list of original expected counts (Ex),
#' a list of observed counts (Yx), the lasso object, a list of K values (number of unique values in each decision path), and n (length of unique centers in the clusters dataframe)
#' @export

spacetimeLasso_sim <- function(clusters, vectors.sim, Time, spacetime,pois, nsim,YSIM,overdispfloor){
    n <- length(unique(clusters$center))
    potClus <- n
    numCenters <- n
    covars <- vectors.sim$covars
    #message("Creating space-time matrix")
    if(spacetime==TRUE){
        sparseMAT <- spacetimeMat(clusters, numCenters, Time)
        message("Space-time matrix created")
    }
    else{
        sparseMAT <- spaceMat(clusters, numCenters)
        message("Spatial matrix created")
    }
    if(!is.null(covars)){
        covarMAT <- Matrix::Matrix(data.matrix(covars), sparse=TRUE)
        sparseMAT <- Matrix::cBind(sparseMAT, covarMAT)
        message("Running with covariates")
    }
    else{
        message("No covariates found")
    }
    
    message(paste("Number of potential clusters to scan through: ", dim(sparseMAT)[2]))
    Ex <- vectors.sim$Ex
    Yx <- YSIM
    Period <- vectors.sim$Period
    message("Running Lasso - stay tuned")
    lasso <- lapply(1:nsim, function(i) glmnet::glmnet(sparseMAT, Yx[[i]], family=("poisson"), alpha=1, offset=log(Ex[[i]]), 
                                               nlambda = 2000, standardize = FALSE, intercept = FALSE, dfmax = 10))
    message("Lasso complete - extracting estimates and paths")
    coefs.lasso.all <- lapply(1:nsim, function(i) coef(lasso[[i]])[-1,]) #do not take intercept
    xbetaPath<- lapply(1:nsim, function(i) sparseMAT%*%coefs.lasso.all[[i]])
    mu <- lapply(1:nsim, function(j) sapply(1:length(lasso[[j]]$lambda), 
                                            function(i) exp(xbetaPath[[j]][,i])))    
    loglike <- lapply(1:nsim, function(k) sapply(1:length(lasso[[k]]$lambda), 
                                                 function(i) sum(dpoisson(Yx[[k]], mu[[k]][,i],Ex[[k]]))))
    K <- lapply(1:nsim, function(i) lasso[[i]]$df)   
    message("Selecting best paths")
    
    #########################################################
    #Space-Time, Quasi-Poisson only (yes overdispersion)
    #########################################################
    if(spacetime==TRUE & pois == FALSE){
        message("returning results for space-time Quasi-Poisson model")
        #calculate max overdispersion
        if(!is.null(covars)){
            offset_reg <- lapply(1:nsim, 
                                 function(i) glm(Yx[[i]] ~ . + as.factor(vectors.sim$Period) +offset(log(Ex[[i]])),
                                                         data = covars, family=quasipoisson))

        }
        else{
            offset_reg <- lapply(1:nsim, function(i) glm(Yx[[i]] ~ 1 + as.factor(vectors.sim$Period) +offset(log(Ex[[i]])),
                                                         family=quasipoisson))
        }
        overdisp.est <- overdisp(offset_reg, sim=TRUE, overdispfloor = overdispfloor)
        message(paste("Overdispersion estimate:", overdisp.est))
        if(pois == FALSE & is.null(overdisp.est)) warning("No overdispersion for quasi-Poisson model. Please check.")

        #QBIC
        PLL.qbic <- lapply(1:nsim, function(i) -2*(loglike[[i]]/overdisp.est) + ((K[[i]])*log(n*Time)))
        select.qbic <- lapply(1:nsim, function(i) which.min(unlist(PLL.qbic[[i]])))
        select_mu.qbic <- lapply(1:nsim, function(i) sapply(select.qbic[[i]], function(j) mu[[i]][,j]))
        select_muRR.qbic <- Reduce("+", select_mu.qbic)/nsim
        E.qbic <- select_muRR.qbic
        numclust.qbic <- lapply(1:nsim, function(i) length(unique(coefs.lasso.all[[i]][,select.qbic[[i]]]))-1)
        ##find probability of being detected = number of times RR was elevated/nsim
        probs.qbic <- prob_incluster(select_mu.qbic, numCenters, Time, nsim)
        
         #QAIC
        PLL.qaic <- lapply(1:nsim, function(i) 2*(K[[i]]) - 2*(loglike[[i]]/overdisp.est))
        select.qaic <- lapply(1:nsim, function(i) which.min(unlist(PLL.qaic[[i]])))
        select_mu.qaic <- lapply(1:nsim, function(i) sapply(select.qaic[[i]], function(j) mu[[i]][,j]))
        select_muRR.qaic <- Reduce("+", select_mu.qaic)/nsim
        E.qaic <- select_muRR.qaic
        numclust.qaic <- lapply(1:nsim, function(i) length(unique(coefs.lasso.all[[i]][,select.qaic[[i]]]))-1)
        ##find probability of being detected = number of times RR was elevated/nsim
        probs.qaic <- prob_incluster(select_mu.qaic, numCenters, Time, nsim)
        


         #QAICc
        PLL.qaicc <- lapply(1:nsim, function(i) 2*(K[[i]]) - 2*(loglike[[i]]/overdisp.est) + 
                                ((2*K[[i]]*(K[[i]] + 1))/(n*Time - K[[i]] - 1)))
        select.qaicc <- lapply(1:nsim, function(i) which.min(unlist(PLL.qaicc[[i]])))
        select_mu.qaicc <- lapply(1:nsim, function(i) sapply(select.qaicc[[i]], function(j) mu[[i]][,j]))
        select_muRR.qaicc <- Reduce("+", select_mu.qaicc)/nsim
        E.qaicc <- select_muRR.qaicc
        numclust.qaicc <- lapply(1:nsim, function(i) length(unique(coefs.lasso.all[[i]][,select.qaicc[[i]]]))-1)
        ##find probability of being detected = number of times RR was elevated/nsim
        probs.qaicc <- prob_incluster(select_mu.qaicc, numCenters, Time, nsim)
    }

    #########################################################
    #Space-Time, Poisson only (no overdispersion)
    #########################################################
    if(spacetime==TRUE & pois == TRUE){
        message("returning results for space-time Poisson model")
        #if(pois == FALSE & !is.null(overdisp.est)) stop("Overdispersion parameter estimated - model is no longer Poisson")
        #QBIC
        PLL.qbic <- lapply(1:nsim, function(i) -2*(loglike[[i]]) + ((K[[i]])*log(n*Time)))
        select.qbic <- lapply(1:nsim, function(i) which.min(unlist(PLL.qbic[[i]])))
        select_mu.qbic <- lapply(1:nsim, function(i) sapply(select.qbic[[i]], function(j) mu[[i]][,j]))
        select_muRR.qbic <- Reduce("+", select_mu.qbic)/nsim
        E.qbic <- select_muRR.qbic
        numclust.qbic <- lapply(1:nsim, function(i) length(unique(coefs.lasso.all[[i]][,select.qbic[[i]]]))-1)
        ##find probability of being detected = number of times RR was elevated/nsim
        probs.qbic <- prob_incluster(select_mu.qbic, numCenters, Time, nsim)
        
      
        #QAIC
        PLL.qaic <- lapply(1:nsim, function(i) 2*(K[[i]]) - 2*(loglike[[i]]))
        select.qaic <- lapply(1:nsim, function(i) which.min(unlist(PLL.qaic[[i]])))
        select_mu.qaic <- lapply(1:nsim, function(i) sapply(select.qaic[[i]], function(j) mu[[i]][,j]))
        select_muRR.qaic <- Reduce("+", select_mu.qaic)/nsim
        E.qaic <- select_muRR.qaic
        numclust.qaic <- lapply(1:nsim, function(i) length(unique(coefs.lasso.all[[i]][,select.qaic[[i]]]))-1)
        ##find probability of being detected = number of times RR was elevated/nsim
        probs.qaic <- prob_incluster(select_mu.qaic, numCenters, Time, nsim)
        
        #QAICc
        PLL.qaicc <- lapply(1:nsim, function(i) 2*(K[[i]]) - 2*(loglike[[i]]) + 
                                ((2*K[[i]]*(K[[i]] + 1))/(n*Time - K[[i]] - 1)))
        select.qaicc <- lapply(1:nsim, function(i) which.min(unlist(PLL.qaicc[[i]])))
        select_mu.qaicc <- lapply(1:nsim, function(i) sapply(select.qaicc[[i]], function(j) mu[[i]][,j]))
        select_muRR.qaicc <- Reduce("+", select_mu.qaicc)/nsim
        E.qaicc <- select_muRR.qaicc
        numclust.qaicc <- lapply(1:nsim, function(i) length(unique(coefs.lasso.all[[i]][,select.qaicc[[i]]]))-1)
        ##find probability of being detected = number of times RR was elevated/nsim
        probs.qaicc <- prob_incluster(select_mu.qaicc, numCenters, Time, nsim)
        
    }
    
    
    #########################################################
    #Space-Only, Quasi-Poisson
    #########################################################
    else if(spacetime==FALSE & pois == FALSE){
        message("Returning results for space-only  Quasi-Poisson model")
        if(!is.null(covars)){
            offset_reg <- lapply(1:nsim, function(i) glm(as.vector(Yx[[i]]) ~ . +offset(log(as.vector(Ex[[i]]))),
                                                         data = covars, family=quasipoisson))
        }
        else{
            offset_reg <- lapply(1:nsim, function(i) glm(Yx[[i]] ~ 1  +offset(log(Ex[[i]])),
                                                         family=quasipoisson))
        }
        overdisp.est <- overdisp(offset_reg, sim=TRUE, overdispfloor = overdispfloor)
        message(paste("Overdispersion estimate:", overdisp.est))
        if(pois == FALSE & is.null(overdisp.est)) warning("No overdispersion for quasi-Poisson model. Please check.")

        #QBIC
        PLL.qbic <- lapply(1:nsim, function(i) -2*(loglike[[i]]/overdisp.est) + ((K[[i]])*log(n)))
        select.qbic <- lapply(1:nsim, function(i) which.min(unlist(PLL.qbic[[i]])))
        select_mu.qbic <- lapply(1:nsim, function(i) sapply(select.qbic[[i]], function(j) mu[[i]][,j]))
        select_muRR.qbic <- Reduce("+", select_mu.qbic)/nsim
        E.qbic <- select_muRR.qbic
        numclust.qbic <- lapply(1:nsim, function(i) length(unique(coefs.lasso.all[[i]][,select.qbic[[i]]]))-1)
        ##find probability of being detected = number of times RR was elevated/nsim
        probs.qbic <- prob_incluster(select_mu.qbic, numCenters, Time, nsim)
        
        

        #QAIC
        PLL.qaic <- lapply(1:nsim, function(i) 2*(K[[i]]) - 2*(loglike[[i]]/overdisp.est))
        select.qaic <- lapply(1:nsim, function(i) which.min(unlist(PLL.qaic[[i]])))
        select_mu.qaic <- lapply(1:nsim, function(i) sapply(select.qaic[[i]], function(j) mu[[i]][,j]))
        select_muRR.qaic <- Reduce("+", select_mu.qaic)/nsim
        E.qaic <- select_muRR.qaic
        numclust.qaic <- lapply(1:nsim, function(i) length(unique(coefs.lasso.all[[i]][,select.qaic[[i]]]))-1)
        ##find probability of being detected = number of times RR was elevated/nsim
        probs.qaic <- prob_incluster(select_mu.qaic, numCenters, Time, nsim)
    


        #QAICc
        PLL.qaicc <- lapply(1:nsim, function(i) 2*(K[[i]]) - 2*(loglike[[i]]/overdisp.est) + 
                                ((2*K[[i]]*(K[[i]] + 1))/(n - K[[i]] - 1)))
        select.qaicc <- lapply(1:nsim, function(i) which.min(unlist(PLL.qaicc[[i]])))
        select_mu.qaicc <- lapply(1:nsim, function(i) sapply(select.qaicc[[i]], function(j) mu[[i]][,j]))
        select_muRR.qaicc <- Reduce("+", select_mu.qaicc)/nsim
        E.qaicc <- select_muRR.qaicc
        numclust.qaicc <- lapply(1:nsim, function(i) length(unique(coefs.lasso.all[[i]][,select.qaicc[[i]]]))-1)
        ##find probability of being detected = number of times RR was elevated/nsim
        probs.qaicc <- prob_incluster(select_mu.qaicc, numCenters, Time, nsim)
     }



    #########################################################
    #Space-only, Poisson only
    #########################################################
    else if(spacetime==FALSE & pois == TRUE){
        message("Returning results for space-only  Poisson model")
        #if(pois == FALSE & !is.null(overdisp.est)) stop("Overdispersion parameter estimated - model is no longer Poisson")
        #QBIC
        PLL.qbic <- lapply(1:nsim, function(i) -2*(loglike[[i]]) + ((K[[i]])*log(n)))
        select.qbic <- lapply(1:nsim, function(i) which.min(unlist(PLL.qbic[[i]])))
        select_mu.qbic <- lapply(1:nsim, function(i) sapply(select.qbic[[i]], function(j) mu[[i]][,j]))
        select_muRR.qbic <- Reduce("+", select_mu.qbic)/nsim
        E.qbic <- select_muRR.qbic
        numclust.qbic <- lapply(1:nsim, function(i) length(unique(coefs.lasso.all[[i]][,select.qbic[[i]]]))-1)
        ##find probability of being detected = number of times RR was elevated/nsim
        probs.qbic <- prob_incluster(select_mu.qbic, numCenters, Time, nsim)
        
    
        #QAIC
        PLL.qaic <- lapply(1:nsim, function(i) 2*(K[[i]]) - 2*(loglike[[i]]))
        select.qaic <- lapply(1:nsim, function(i) which.min(unlist(PLL.qaic[[i]])))
        select_mu.qaic <- lapply(1:nsim, function(i) sapply(select.qaic[[i]], function(j) mu[[i]][,j]))
        select_muRR.qaic <- Reduce("+", select_mu.qaic)/nsim
        E.qaic <- select_muRR.qaic
        numclust.qaic <- lapply(1:nsim, function(i) length(unique(coefs.lasso.all[[i]][,select.qaic[[i]]]))-1)
        ##find probability of being detected = number of times RR was elevated/nsim
        probs.qaic <- prob_incluster(select_mu.qaic, numCenters, Time, nsim)
        
        
         
        #QAICc
        PLL.qaicc <- lapply(1:nsim, function(i) 2*(K[[i]]) - 2*(loglike[[i]]) + 
                                ((2*K[[i]]*(K[[i]] + 1))/(n - K[[i]] - 1)))
        select.qaicc <- lapply(1:nsim, function(i) which.min(unlist(PLL.qaicc[[i]])))
        select_mu.qaicc <- lapply(1:nsim, function(i) sapply(select.qaicc[[i]], function(j) mu[[i]][,j]))
        select_muRR.qaicc <- Reduce("+", select_mu.qaicc)/nsim
        E.qaicc <- select_muRR.qaicc
        numclust.qaicc <- lapply(1:nsim, function(i) length(unique(coefs.lasso.all[[i]][,select.qaicc[[i]]]))-1)
        ##find probability of being detected = number of times RR was elevated/nsim
        probs.qaicc <- prob_incluster(select_mu.qaicc, numCenters, Time, nsim)
     }
    
    return(list(nsim = nsim, E.qbic = E.qbic, E.qaic = E.qaic, E.qaicc = E.qaicc,Ex = Ex,mu = mu, Yx = Yx, PLL.qbic = PLL.qbic, 
                PLL.qaic = PLL.qaic, PLL.qaicc = PLL.qaicc, select.qbic = select.qbic, select.qaic = select.qaic, 
                select.qaicc = select.qaicc, select_mu.qbic = select_mu.qbic, select_mu.qaic = select_mu.qaic, 
                select_mu.qaicc = select_mu.qaicc, xbetaPath = xbetaPath, coefs.lasso.all = coefs.lasso.all,
                numclust.qaic = numclust.qaic, numclust.qaicc = numclust.qaicc, numclust.qbic = numclust.qbic,
                probs.qbic = probs.qbic, probs.qaic = probs.qaic, probs.qaicc = probs.qaicc))    
}
