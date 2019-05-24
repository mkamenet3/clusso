#' Simulations for Spatial and Spatio-Temporal Lasso
#' 
#' This function runs the Lasso regularization technique on our large sparse matric of potential space-time clusters.It is specifically created to use with simulations.
#' @param sparseMAT large sparse matrix created in \code{clust_sim} function.
#' @param n_uniq number of unique polygons (ex: counties, zip code, etc). Inherited from \code{clust_sim}.
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
spacetimeLasso_sim <- function(sparseMAT, n_uniq ,vectors.sim, Time, spacetime,pois, nsim,YSIM,maxclust,overdispfloor){
    #check for covariates
    covars <- vectors.sim$covars
    if(!is.null(covars)){
        message("Running with covariates")
        covarMAT <- Matrix::Matrix(data.matrix(covars), sparse=TRUE)
        sparseMAT <- cbind(sparseMAT, covarMAT)
    }
    else{
        message("No covariates found")
    }
    #set initial
    message(paste("Number of potential clusters to scan through: ", dim(sparseMAT)[2]))
    Ex <- vectors.sim$Ex
    Yx <- YSIM
    Period <- vectors.sim$Period
    penalty <- c(rep(1,(ncol(sparseMAT)-Time)), rep(0,Time))
    message("Running Lasso - stay tuned")
    lasso <- lapply(1:nsim, function(i) glmnet::glmnet(sparseMAT, Yx[[i]], family=("poisson"), alpha=1, offset=log(Ex[[i]]), 
                                               nlambda = 2000, standardize = FALSE, intercept = FALSE, dfmax = maxclust, 
                                               penalty.factor = penalty))
    message("Lasso complete - extracting estimates and paths")
    coefs.lasso.all <- lapply(1:nsim, function(i) coef(lasso[[i]])[-1,]) #do not take intercept
    xbetaPath<- lapply(1:nsim, function(i) sparseMAT%*%coefs.lasso.all[[i]])
    mu <- lapply(1:nsim, function(j) sapply(1:length(lasso[[j]]$lambda),
                                            function(i) exp(xbetaPath[[j]][,i])))
    loglike <- lapply(1:nsim, function(k) sapply(1:length(lasso[[k]]$lambda),
                                                 function(i) sum(dpoisson(Yx[[k]], mu[[k]][,i],Ex[[k]]))))
    K <- lapply(1:nsim, function(i) lasso[[i]]$df)   
    message("Selecting best paths")
    
    ############################################
    
    # ############################################
    
    
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
        message(paste("Overdispersion estimate:", round(overdisp.est,4)))
        if(pois == FALSE & is.null(overdisp.est)) warning("No overdispersion for quasi-Poisson model. Please check.")

        #QBIC
        PLL.qbic <- lapply(1:nsim, function(i) -2*(loglike[[i]]/overdisp.est) + ((K[[i]])*log(n_uniq*Time)))
        select.qbic <- lapply(1:nsim, function(i) which.min(unlist(PLL.qbic[[i]])))
        #select_mu.qbic <- lapply(1:nsim, function(i) sapply(select.qbic[[i]], function(j) mu[[i]][,j]))
        select_mu.qbic <- lapply(1:nsim, function(i) sapply(select.qbic[[i]], function(j) log(mu[[i]][,j])))
        select_muRR.qbic <- Reduce("+", select_mu.qbic)/nsim
        #E.qbic <- select_muRR.qbic
        E.qbic <- exp(select_muRR.qbic)
        # exp_coefs_qbic <- c(exp(unique(lasso$beta[,select.qbic])),
        #                     exp(lasso$beta[(nrow(lasso$beta)-Time+1):nrow(lasso$beta),select.qbic]))
        #numclust.qbic <- length(unique(exp_coefs_qbic))-Time
        numclust.qbic <- lapply(1:nsim, function(i) length(unique(coefs.lasso.all[[i]][,select.qbic[[i]]]))-(Time+1))

         #QAIC
        PLL.qaic <- lapply(1:nsim, function(i) 2*(K[[i]]) - 2*(loglike[[i]]/overdisp.est))
        select.qaic <- lapply(1:nsim, function(i) which.min(unlist(PLL.qaic[[i]])))
        #select_mu.qaic <- lapply(1:nsim, function(i) sapply(select.qaic[[i]], function(j) mu[[i]][,j]))
        select_mu.qaic <- lapply(1:nsim, function(i) sapply(select.qaic[[i]], function(j) log(mu[[i]][,j])))
        select_muRR.qaic <- Reduce("+", select_mu.qaic)/nsim
        #E.qaic <- select_muRR.qaic
        E.qaic <- exp(select_muRR.qaic)
        #numclust.qaic <- lapply(1:nsim, function(i) length(unique(coefs.lasso.all[[i]][,select.qaic[[i]]]))-1)
        numclust.qaic <- lapply(1:nsim, function(i) length(unique(coefs.lasso.all[[i]][,select.qaic[[i]]]))-(Time+1))

         #QAICc
        PLL.qaicc <- lapply(1:nsim, function(i) 2*(K[[i]]) - 2*(loglike[[i]]/overdisp.est) + 
                                ((2*K[[i]]*(K[[i]] + 1))/(n_uniq*Time - K[[i]] - 1)))
        select.qaicc <- lapply(1:nsim, function(i) which.min(unlist(PLL.qaicc[[i]])))
        #select_mu.qaicc <- lapply(1:nsim, function(i) sapply(select.qaicc[[i]], function(j) mu[[i]][,j]))
        select_mu.qaicc <- lapply(1:nsim, function(i) sapply(select.qaicc[[i]], function(j) log(mu[[i]][,j])))
        select_muRR.qaicc <- Reduce("+", select_mu.qaicc)/nsim
        #E.qaicc <- select_muRR.qaicc
        E.qaicc <- exp(select_muRR.qaicc)
        #numclust.qaicc <- lapply(1:nsim, function(i) length(unique(coefs.lasso.all[[i]][,select.qaicc[[i]]]))-1)
        numclust.qaicc <- lapply(1:nsim, function(i) length(unique(coefs.lasso.all[[i]][,select.qaicc[[i]]]))-(Time+1))

    }

    #########################################################
    #Space-Time, Poisson only (no overdispersion)
    #########################################################
    if(spacetime==TRUE & pois == TRUE){
        message("returning results for space-time Poisson model")
        #QBIC
        PLL.qbic <- lapply(1:nsim, function(i) -2*(loglike[[i]]) + ((K[[i]])*log(n_uniq*Time)))
        select.qbic <- lapply(1:nsim, function(i) which.min(unlist(PLL.qbic[[i]])))
        #select_mu.qbic <- lapply(1:nsim, function(i) sapply(select.qbic[[i]], function(j) mu[[i]][,j]))
        select_mu.qbic <- lapply(1:nsim, function(i) sapply(select.qbic[[i]], function(j) log(mu[[i]][,j])))
        select_muRR.qbic <- Reduce("+", select_mu.qbic)/nsim
        #E.qbic <- select_muRR.qbic
        E.qbic <- exp(select_muRR.qbic)
        #numclust.qbic <- lapply(1:nsim, function(i) length(unique(coefs.lasso.all[[i]][,select.qbic[[i]]]))-1)
        numclust.qbic <- lapply(1:nsim, function(i) length(unique(coefs.lasso.all[[i]][,select.qbic[[i]]]))-(Time+1))
        
        #QAIC
        PLL.qaic <- lapply(1:nsim, function(i) 2*(K[[i]]) - 2*(loglike[[i]]))
        select.qaic <- lapply(1:nsim, function(i) which.min(unlist(PLL.qaic[[i]])))
        select_mu.qaic <- lapply(1:nsim, function(i) sapply(select.qaic[[i]], function(j) log(mu[[i]][,j])))
        select_muRR.qaic <- Reduce("+", select_mu.qaic)/nsim
        E.qaic <- exp(select_muRR.qaic)
        #numclust.qaic <- lapply(1:nsim, function(i) length(unique(coefs.lasso.all[[i]][,select.qaic[[i]]]))-1)
        numclust.qaic <- lapply(1:nsim, function(i) length(unique(coefs.lasso.all[[i]][,select.qaic[[i]]]))-(Time+1))
       
        #QAICc
        PLL.qaicc <- lapply(1:nsim, function(i) 2*(K[[i]]) - 2*(loglike[[i]]) + 
                                ((2*K[[i]]*(K[[i]] + 1))/(n_uniq*Time - K[[i]] - 1)))
        select.qaicc <- lapply(1:nsim, function(i) which.min(unlist(PLL.qaicc[[i]])))
        #select_mu.qaicc <- lapply(1:nsim, function(i) sapply(select.qaicc[[i]], function(j) mu[[i]][,j]))
        select_mu.qaicc <- lapply(1:nsim, function(i) sapply(select.qaicc[[i]], function(j) log(mu[[i]][,j])))
        select_muRR.qaicc <- Reduce("+", select_mu.qaicc)/nsim
        #E.qaicc <- select_muRR.qaicc
        E.qaicc <- exp(select_muRR.qaicc)
        #numclust.qaicc <- lapply(1:nsim, function(i) length(unique(coefs.lasso.all[[i]][,select.qaicc[[i]]]))-1)
        numclust.qaicc <- lapply(1:nsim, function(i) length(unique(coefs.lasso.all[[i]][,select.qaicc[[i]]]))-(Time+1))
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
        message(paste("Overdispersion estimate:", round(overdisp.est,4)))
        if(pois == FALSE & is.null(overdisp.est)) warning("No overdispersion for quasi-Poisson model. Please check.")

        #QBIC
        PLL.qbic <- lapply(1:nsim, function(i) -2*(loglike[[i]]/overdisp.est) + ((K[[i]])*log(n_uniq)))
        select.qbic <- lapply(1:nsim, function(i) which.min(unlist(PLL.qbic[[i]])))
        #select_mu.qbic <- lapply(1:nsim, function(i) sapply(select.qbic[[i]], function(j) mu[[i]][,j]))
        select_mu.qbic <- lapply(1:nsim, function(i) sapply(select.qbic[[i]], function(j) log(mu[[i]][,j])))
        select_muRR.qbic <- Reduce("+", select_mu.qbic)/nsim
        #E.qbic <- select_muRR.qbic
        E.qbic <- exp(select_muRR.qbic)
        #numclust.qbic <- lapply(1:nsim, function(i) length(unique(coefs.lasso.all[[i]][,select.qbic[[i]]]))-1)
        numclust.qbic <- lapply(1:nsim, function(i) length(unique(coefs.lasso.all[[i]][,select.qbic[[i]]]))-(Time+1))
       
        #QAIC
        PLL.qaic <- lapply(1:nsim, function(i) 2*(K[[i]]) - 2*(loglike[[i]]/overdisp.est))
        select.qaic <- lapply(1:nsim, function(i) which.min(unlist(PLL.qaic[[i]])))
        #select_mu.qaic <- lapply(1:nsim, function(i) sapply(select.qaic[[i]], function(j) mu[[i]][,j]))
        select_mu.qaic <- lapply(1:nsim, function(i) sapply(select.qaic[[i]], function(j) log(mu[[i]][,j])))
        select_muRR.qaic <- Reduce("+", select_mu.qaic)/nsim
        #E.qaic <- select_muRR.qaic
        E.qaic <- exp(select_muRR.qaic)
        #numclust.qaic <- lapply(1:nsim, function(i) length(unique(coefs.lasso.all[[i]][,select.qaic[[i]]]))-1)
        numclust.qaic <- lapply(1:nsim, function(i) length(unique(coefs.lasso.all[[i]][,select.qaic[[i]]]))-(Time+1))
       
        #QAICc
        PLL.qaicc <- lapply(1:nsim, function(i) 2*(K[[i]]) - 2*(loglike[[i]]/overdisp.est) + 
                                ((2*K[[i]]*(K[[i]] + 1))/(n_uniq - K[[i]] - 1)))
        select.qaicc <- lapply(1:nsim, function(i) which.min(unlist(PLL.qaicc[[i]])))
        #select_mu.qaicc <- lapply(1:nsim, function(i) sapply(select.qaicc[[i]], function(j) mu[[i]][,j]))
        select_mu.qaicc <- lapply(1:nsim, function(i) sapply(select.qaicc[[i]], function(j) log(mu[[i]][,j])))
        select_muRR.qaicc <- Reduce("+", select_mu.qaicc)/nsim
        #E.qaicc <- select_muRR.qaicc
        E.qaicc <- exp(select_muRR.qaicc)
        #numclust.qaicc <- lapply(1:nsim, function(i) length(unique(coefs.lasso.all[[i]][,select.qaicc[[i]]]))-1)
        numclust.qaicc <- lapply(1:nsim, function(i) length(unique(coefs.lasso.all[[i]][,select.qaicc[[i]]]))-(Time+1))
     }

    #########################################################
    #Space-only, Poisson only
    #########################################################
    else if(spacetime==FALSE & pois == TRUE){
        message("Returning results for space-only  Poisson model")
        #QBIC
        PLL.qbic <- lapply(1:nsim, function(i) -2*(loglike[[i]]) + ((K[[i]])*log(n_uniq)))
        select.qbic <- lapply(1:nsim, function(i) which.min(unlist(PLL.qbic[[i]])))
        #select_mu.qbic <- lapply(1:nsim, function(i) sapply(select.qbic[[i]], function(j) mu[[i]][,j]))
        select_mu.qbic <- lapply(1:nsim, function(i) sapply(select.qbic[[i]], function(j) log(mu[[i]][,j])))
        select_muRR.qbic <- Reduce("+", select_mu.qbic)/nsim
        #E.qbic <- select_muRR.qbic
        E.qbic <- exp(select_muRR.qbic)
        #numclust.qbic <- lapply(1:nsim, function(i) length(unique(coefs.lasso.all[[i]][,select.qbic[[i]]]))-1)
        numclust.qbic <- lapply(1:nsim, function(i) length(unique(coefs.lasso.all[[i]][,select.qbic[[i]]]))-(Time+1))
       
        #QAIC
        PLL.qaic <- lapply(1:nsim, function(i) 2*(K[[i]]) - 2*(loglike[[i]]))
        select.qaic <- lapply(1:nsim, function(i) which.min(unlist(PLL.qaic[[i]])))
        #select_mu.qaic <- lapply(1:nsim, function(i) sapply(select.qaic[[i]], function(j) mu[[i]][,j]))
        select_mu.qaic <- lapply(1:nsim, function(i) sapply(select.qaic[[i]], function(j) log(mu[[i]][,j])))
        select_muRR.qaic <- Reduce("+", select_mu.qaic)/nsim
        #E.qaic <- select_muRR.qaic
        E.qaic <- exp(select_muRR.qaic)
        #numclust.qaic <- lapply(1:nsim, function(i) length(unique(coefs.lasso.all[[i]][,select.qaic[[i]]]))-1)
        numclust.qaic <- lapply(1:nsim, function(i) length(unique(coefs.lasso.all[[i]][,select.qaic[[i]]]))-(Time+1))
         
        #QAICc
        PLL.qaicc <- lapply(1:nsim, function(i) 2*(K[[i]]) - 2*(loglike[[i]]) + 
                                ((2*K[[i]]*(K[[i]] + 1))/(n_uniq - K[[i]] - 1)))
        select.qaicc <- lapply(1:nsim, function(i) which.min(unlist(PLL.qaicc[[i]])))
        #select_mu.qaicc <- lapply(1:nsim, function(i) sapply(select.qaicc[[i]], function(j) mu[[i]][,j]))
        select_mu.qaicc <- lapply(1:nsim, function(i) sapply(select.qaicc[[i]], function(j) log(mu[[i]][,j])))
        select_muRR.qaicc <- Reduce("+", select_mu.qaicc)/nsim
        #E.qaicc <- select_muRR.qaicc
        E.qaicc <- exp(select_muRR.qaicc)
        #numclust.qaicc <- lapply(1:nsim, function(i) length(unique(coefs.lasso.all[[i]][,select.qaicc[[i]]]))-1)
        numclust.qaicc <- lapply(1:nsim, function(i) length(unique(coefs.lasso.all[[i]][,select.qaicc[[i]]]))-(Time+1))
    }
    
    #Return only changepoints from lasso
    changepoints_ix <- lapply(1:nsim, function(k) which(diff(K[[k]])!=0)) #Find lambda where new coef introduced
    lambda_changepoint <- lapply(1:nsim, function(k) lasso[[k]]$lambda[changepoints_ix[[k]]])
    #QIC
    coefs_qbic <- lapply(1:nsim, 
                         function(k) coefs.lasso.all[[k]][which(coefs.lasso.all[[k]][,select.qbic[[k]]]!=0), changepoints_ix[[k]]])
    coefs_qaic <- lapply(1:nsim, 
                         function(k) coefs.lasso.all[[k]][which(coefs.lasso.all[[k]][,select.qaic[[k]]]!=0), changepoints_ix[[k]]])
    coefs_qaicc <- lapply(1:nsim, 
                         function(k) coefs.lasso.all[[k]][which(coefs.lasso.all[[k]][,select.qaicc[[k]]]!=0), changepoints_ix[[k]]])

    lasso_out <- list(
        lambdas = lambda_changepoint,
        coefs_bic = coefs_qbic,
        coefs_aic = coefs_qaic,
        coefs_aicc = coefs_qaicc
    )
    
    
    return(list(lasso = lasso, lasso_out=lasso_out, nsim = nsim, E.qbic = E.qbic, E.qaic = E.qaic, E.qaicc = E.qaicc,Ex = Ex,mu = mu, Yx = Yx, PLL.qbic = PLL.qbic, 
                PLL.qaic = PLL.qaic, PLL.qaicc = PLL.qaicc, select.qbic = select.qbic, select.qaic = select.qaic, 
                select.qaicc = select.qaicc, select_mu.qbic = select_mu.qbic, select_mu.qaic = select_mu.qaic, 
                select_mu.qaicc = select_mu.qaicc, xbetaPath = xbetaPath, coefs.lasso.all = coefs.lasso.all,
                numclust.qaic = numclust.qaic, numclust.qaicc = numclust.qaicc, numclust.qbic = numclust.qbic))    
}
