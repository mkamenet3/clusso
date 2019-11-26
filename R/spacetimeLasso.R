#' Spatial and Spatio-Temporal Cluster Detection Using the LASSO
#' @title 
#' spacetimeLasso
#' 
#' @description 
#' This function runs the LASSO regularization technique on the large sparse matrix of potential space or space-time clusters.
#'@param model A string specifying which model to use, Poisson or binomial. For Poisson, specify \code{"poisson"} and both the Poisson and quasi-Poisson model results are returned. For binomial, specify \code{"binomial"} and both the binomial and quasi-binomial model results are returned.
#' @param sparseMAT Large sparse matrix of indicators of space-time potential clusters.
#' @param n_uniq Number of unique geographic subregions/centroids (ex: counties, zip code, ,etc). Calculated as the number of unique geographic identifiers.
#' @param vectors List of expected and observed counts inherited from \code{setVectors}.
#' @param Time Number of timeperiods in the dataset. This is calculated based on the number of unique timeperiods (factor levels) supplied to \code{clusso}.
#' @param quasi Whether or not quasi- models for the Poisson and Binomial models should be executed. Default is quasi=FALSE.
#' @param maxclust Upper limit on the maximum number of clusters you expect to find in the region. This equivalent to setting \code{dfmax} in \code{glmnet}. If none supplied, default is \code{11}.
#' @param overdispfloor Default is \code{TRUE}. When \code{TRUE}, it limits \eqn{\phi} (overdispersion parameter) to be greater or equal to 1. If \code{FALSE}, it will allow for under-dispersion in the model.
#' @param cv Numeric argument for the number of folds to use if using k-fold cross-validation. Default is \code{NULL}, indicating that cross-validation should not be performed in favor of \code{clusso}.
#' @return Returns a list of lists with the 1) estimated relative risks for each geographic subregion-time period based on selection by (Q)BIC, (Q)AIC, (Q)AICc, 2) number of identified clusters by (Q)BIC, (Q)AIC, (Q)AICc, 3) original expected counts (Ex), 4) original observed counts (Yx),5) the LASSO object, 6)  list of K values (number of unique values in each LASSO path), 7) and the estimated coefficients from all estimated lambdas.

spacetimeLasso<- function(model, sparseMAT, n_uniq, vectors,Time, quasi,maxclust, overdispfloor, cv){
    #check for covariates
    covars <- vectors$covars
    if(!is.null(covars)){
        message("Running with covariates")
        covarMAT <- Matrix::Matrix(data.matrix(covars), sparse=TRUE)
        sparseMAT <- cbind(sparseMAT, covarMAT)
        message(paste("Number of potential clusters to scan through: ", (dim(sparseMAT)[2]-(Time+ncol(covarMAT)))))
    }
    else{
        message("No covariates found")
        message(paste("Number of potential clusters to scan through: ", (dim(sparseMAT)[2]-Time)))
    }
    #Set initial
    Ex <- vectors$Ex
    Yx <- vectors$Y.vec
    Period <- vectors$Period
    
    if(!is.null(cv)){
        message("Path selection: cross-validation")
        if(!is.null(covars)){
            penalty <- c(rep(1,(ncol(sparseMAT)-Time)), rep(0,(Time+ncol(covarMAT))))   
            # #print(table(penalty))
        }
        else{
            penalty <- c(rep(1,(ncol(sparseMAT)-Time)), rep(0,Time))    
        }
        if(model=="poisson"){
            lasso <- glmnet::cv.glmnet(sparseMAT, Yx, family=("poisson"), alpha=1, offset=log(Ex), nlambda = 2000,
                                       standardize = FALSE, intercept=FALSE,dfmax = maxclust,
                                       nfolds = cv,
                                       penalty.factor = penalty)
        }
        else if(model=="binomial"){
            lasso <- glmnet::cv.glmnet(sparseMAT, cbind((Ex-Yx),Yx), family=("binomial"), alpha=1, nlambda = 2000,
                                       standardize = FALSE, intercept=FALSE,dfmax = maxclust,
                                       nfolds = cv,
                                       penalty.factor = penalty)
        }
    }
    else{
        message("Path selection: information criteria")
        if(!is.null(covars)){
            penalty <- c(rep(1,(ncol(sparseMAT)-Time)), rep(0,(Time+ncol(covarMAT))))   
        }
        else{
            penalty <- c(rep(1,(ncol(sparseMAT)-Time)), rep(0,Time))    
        }
        if(model=="poisson"){
            #print("Model poisson selected in spacetimeLasso")
            lasso <- glmnet::glmnet(sparseMAT, Yx, family=("poisson"), alpha=1, offset=log(Ex), 
                                    nlambda = 2000,
                                    standardize = FALSE, intercept=FALSE,dfmax = maxclust,
                                    penalty.factor = penalty)    
        }
        else if(model=="binomial"){
            #print("Model binomial selected in spacetimeLasso")
            lasso <- glmnet::glmnet(sparseMAT, cbind((Ex-Yx),Yx), family=("binomial"), alpha=1, 
                                    nlambda = 2000,
                                    standardize = FALSE, intercept=FALSE,dfmax = maxclust,
                                    penalty.factor = penalty)    
        }
        else{
            stop("No model specified.")
        }
        
        
    }
    message("Lasso complete - extracting estimates and paths")
    coefs.lasso.all <- coef(lasso)[-1,]
    xbetaPath<- sparseMAT%*%coefs.lasso.all
    K <- lasso$df
    #if running cross-validation version:
    if(!is.null(cv)){
        res <- stLasso_cv(lasso, sparseMAT, Ex, Yx, Time)
    }
    #information criteria selection version:
    else{
        #print("Pois post-process")
        #print(paste0("quasi: ",quasi))
        if(model=="poisson"){
            mu <- sapply(1:length(lasso$lambda), function(i) exp(xbetaPath[,i]))
            loglike <- sapply(1:length(lasso$lambda), function(i) sum(dpoisson(Yx, mu[,i],Ex)))
            res <- spacetimeLassoPois(lasso, coefs.lasso.all,loglike,mu,K, quasi, 
                                      covars,Yx, Ex, Period, Time, n_uniq, overdispfloor, maxclust)
        }
        else if(model=="binomial"){
            #print("Binom post-process")
            phat <-  sapply(1:length(lasso$lambda), function(i) pihat(xbetaPath[,i]))
            loglike <- sapply(1:length(lasso$lambda), function(i) sum(dbin(Yx, Ex, phat[,i])))
            res <- spacetimeLassoBinom(lasso,coefs.lasso.all, loglike, K, quasi,covars,
                                       Yx, Ex, Period, Time, n_uniq, overdispfloor,maxclust)
        }
        return(res)
    }
}


#'@title spacetimeLassoPois
#'@description
#'@param loglike Loglikelihood for Poisson model
#'@param coefs.lasso.all Matrix of coefficient estimates for every lambda in lasso path.
#'@param mu Estimated relative risk estimates (\code{\rho_i\E_i})
#'@param K Vector of the number of K parameters estimated for every lambda in lasso path.
#'@param quasi Boolean. \code{TRUE} indicates a quasi-Poisson model that accounts for overdispersion. \code{FALSE} indicates a Poisson model without adjustment for overdispersion.
#'@param covars Dataframe of additional covariates to be included in the model that are un-penalized by the LASSO.
#'@param Yx Number of observed cases for each space-time location.
#'@param Ex Expected number of cases for each space-time location.
#'@param Period Vector of time periods for each space-time location.
#'@param Time Integer. Number of time periods in the study.
#'@param n_uniq Number of unique centroids in the study area.
#'@param overdispfloor Default is \code{TRUE}. When TRUE, it limits \eqn{\phi} (overdispersion parameter) to be greater or equal to 1. If FALSE, will allow for under-dispersion in the model.  
#'@param maxclust  Upper limit on the maximum number of clusters you expect to find in the region. This equivalent to setting \code{dfmax} in the lasso. If none supplied, default is \code{11}.
spacetimeLassoPois <- function(lasso, coefs.lasso.all, loglike,mu, K, quasi, covars,Yx, Ex, Period, Time, n_uniq, overdispfloor, maxclust){    
    #########################################################
    #Quasi-Poisson only (yes overdispersion)
    #########################################################
    if(quasi == TRUE){
        #message("Returning results for space-time Quasi-Poisson model")
        if(!is.null(covars)){
            offset_reg <- glm(Yx ~ . + as.factor(Period) + offset(log(Ex)),
                              data = covars,family=quasipoisson)
        }
        else{
            offset_reg <- glm(Yx ~ 1 + as.factor(Period) + offset(log(Ex)),
                              family=quasipoisson)
        }
        overdisp.est <- overdisp(offset_reg, sim = FALSE, overdispfloor = overdispfloor)
        message(paste("Overdispersion estimate:", round(overdisp.est,4)))
        if(quasi == TRUE & is.null(overdisp.est)) warning("No overdispersion for quasi-Poisson model. Please check.")
        
        #QBIC
        PLL.qbic  <- -2*(loglike/overdisp.est) + ((K)*log(sum(Yx)))
        select.qbic <- which.min(PLL.qbic)
        E.qbic <- mu[,select.qbic]
        numclust.qbic <- K[select.qbic]-Time 
        
        #QAIC
        PLL.qaic <-  2*(K) - 2*(loglike/overdisp.est)
        select.qaic <- which.min(PLL.qaic)
        E.qaic <- mu[,select.qaic]
        numclust.qaic <- K[select.qaic]-Time 
        
        #QAICc
        PLL.qaicc <- 2*(K) - 2*(loglike/overdisp.est) +
            ((2*K*(K + 1))/(sum(Yx) - K - 1))
        select.qaicc <- which.min(PLL.qaicc)
        E.qaicc <- mu[,select.qaicc]
        numclust.qaicc <- K[select.qaicc]-Time 

    }
    #########################################################
    #Poisson only (no overdispersion)
    #########################################################
    else if(quasi==FALSE){
        #QBIC
        PLL.qbic  <- -2*(loglike) + ((K)*log(sum(Yx)))
        select.qbic <- which.min(PLL.qbic)
        E.qbic <- mu[,select.qbic]
        numclust.qbic <- K[select.qbic]-Time 
        
        #QAIC
        PLL.qaic <-  2*(K) - 2*(loglike)
        select.qaic <- which.min(PLL.qaic)
        E.qaic <- mu[,select.qaic]
        numclust.qaic <- K[select.qaic]-Time 
        
        #QAICc
        PLL.qaicc <- 2*(K) - 2*(loglike) +
            ((2*K*(K + 1))/(sum(Yx) - K - 1))
        select.qaicc <- which.min(PLL.qaicc)
        E.qaicc <- mu[,select.qaicc]
        numclust.qaicc <- K[select.qaicc]-Time 
    }
    #warning if numclust similar to maxclust
    if (numclust.qaic==(maxclust-Time)){
        message("The number of clusters selected by at least one criterion is equal to maxclust. You may want to increase maxclust.")
    }
    res <- list(E.qbic = E.qbic, E.qaic = E.qaic, E.qaicc = E.qaicc, numclust.qaic = numclust.qaic,
                numclust.qaicc = numclust.qaicc, numclust.qbic= numclust.qbic, Ex = Ex, Yx = Yx, 
                lasso = lasso, K = K, coefs.lasso.all = coefs.lasso.all)
}




#'@title spacetimeLassoBinom
#'@description
#'@param loglike Loglikelihood for binomial model.
#'@param coefs.lasso.all Matrix of coefficient estimates for every lambda in lasso path.
#'@param K Vector of the number of K parameters estimated for every lambda in lasso path.
#'@param quasi Boolean. \code{TRUE} indicates a quasi-Poisson model that accounts for overdispersion. \code{FALSE} indicates a Poisson model without adjustment for overdispersion.
#'@param covars Dataframe of additional covariates to be included in the model that are un-penalized by the LASSO.
#'@param Yx Number of observed cases for each space-time location.
#'@param Ex Total number of trials (cases + controls) for each space-time location.
#'@param Period Vector of time periods for each space-time location.
#'@param Time Integer. Number of time periods in the study.
#'@param n_uniq Number of unique centroids in the study area.
#'@param overdispfloor Default is \code{TRUE}. When TRUE, it limits \eqn{\phi} (overdispersion parameter) to be greater or equal to 1. If FALSE, will allow for under-dispersion in the model.  
#'@param maxclust  Upper limit on the maximum number of clusters you expect to find in the region. This equivalent to setting \code{dfmax} in the lasso. If none supplied, default is \code{11}.
spacetimeLassoBinom <- function(lasso, coefs.lasso.all, loglike, K, quasi,covars, Yx, Ex, Period, Time, n_uniq, overdispfloor,maxclust){   
    if(quasi==TRUE){
        #message("Returning results for space-time Quasi-Poisson model")
        if(!is.null(covars)){
            offset_reg <- glm(cbind(Yx, (Ex-Yx))~ . + as.factor(Period),
                              data = covars,family=quasibinomial)
        }
        else{
            offset_reg <- glm(cbind(Yx, (Ex-Yx)) ~ 1 + as.factor(Period),
                              family=quasibinomial)
        }
        overdisp.est <- overdisp(offset_reg, sim = FALSE, overdispfloor = overdispfloor)
        message(paste("Overdispersion estimate:", round(overdisp.est,4)))
        if(quasi == TRUE & is.null(overdisp.est)) warning("No overdispersion for quasi-Poisson model. Please check.")
        #########################################################
        #Space-Time, Binomial
        #########################################################
        #QBIC
        PLL.qbic  <- -2*(loglike/overdisp.est) + ((K)*log(min(sum(Yx),sum(Ex-Yx))))
        select.qbic <- which.min(PLL.qbic)
        E.qbic <- mu[,select.qbic]
        numclust.qbic <- K[select.qbic]-Time 
        
        #QAIC
        PLL.qaic <-  2*(K) - 2*(loglike/overdisp.est)
        select.qaic <- which.min(PLL.qaic)
        E.qaic <- mu[,select.qaic]
        numclust.qaic <- K[select.qaic]-Time 
        
        #QAICc
        PLL.qaicc <- 2*(K) - 2*(loglike/overdisp.est) +
            ((2*K*(K + 1))/(min(sum(Yx),sum(Ex-Yx)) - K - 1))
        select.qaicc <- which.min(PLL.qaicc)
        E.qaicc <- mu[,select.qaicc]
        numclust.qaicc <- K[select.qaicc]-Time 

    }
    else if (quasi==FALSE){
        #########################################################
        #Space-Time, Binomial
        #########################################################
        #QBIC
        PLL.qbic  <- -2*(loglike) + ((K)*log(min(sum(Yx),sum(Ex-Yx))))
        select.qbic <- which.min(PLL.qbic)
        E.qbic <- mu[,select.qbic]
        numclust.qbic <- K[select.qbic]-Time 
        
        #QAIC
        PLL.qaic <-  2*(K) - 2*(loglike)
        select.qaic <- which.min(PLL.qaic)
        E.qaic <- mu[,select.qaic]
        numclust.qaic <- K[select.qaic]-Time 
        
        #QAICc
        PLL.qaicc <- 2*(K) - 2*(loglike) +
            ((2*K*(K + 1))/(min(sum(Yx),sum(Ex-Yx)) - K - 1))
        select.qaicc <- which.min(PLL.qaicc)
        E.qaicc <- mu[,select.qaicc]

        numclust.qaicc <- K[select.qaicc]-Time 

    }
    
    
    #warning if numclust similar to maxclust
    if (numclust.qaic==(maxclust-Time)){
        message("The number of clusters selected by at least one criterion is equal to maxclust. You may want to increase maxclust.")
    }
    #Add if quasi=TRUE for quasi-binomial
    #warning if numclust similar to maxclust
    if (numclust.qaic==(maxclust-Time)){
        message("The number of clusters selected by at least one criterion is equal to maxclust. You may want to increase maxclust.")
    }

    res <- list(E.qbic = E.qbic, E.qaic = E.qaic, E.qaicc = E.qaicc, numclust.qaic = numclust.qaic,
                numclust.qaicc = numclust.qaicc, numclust.qbic= numclust.qbic, Ex = Ex, Yx = Yx, 
                lasso = lasso,K = K, coefs.lasso.all = coefs.lasso.all)
}




#' Selection via cross-validation
#' @title 
#' stLasso_cv
#' 
#' @description 
#' This function will output results from cross-validation results.
#' @param lasso results of cv.glmnet
#' @param sparseMAT sparse spatial or spatio-temporal design matrix
#' @param Ex expected counts (inherits from setVectors)
#' @param Yx observed counts (inherits from setVectors)
#' @param Time number of time periods
#' @return list of expected values, number of clusters, Ex, Yx, and lasso object
#' 
stLasso_cv <- function(lasso, sparseMAT, Ex, Yx, Time){
    #cv version
    ix <- which(lasso$lambda == lasso$lambda.min)
    E.cv <- lasso$glmnet.fit$beta[,ix]
    #mu.cv <- exp(E.cv)
    xbetaPath<- sparseMAT%*%lasso$glmnet.fit$beta[,ix]
    mu <- exp(xbetaPath)
    numclust.cv <- length(unique(E.cv))-Time
    return(list(E.cv = mu[,1], numclust.cv = numclust.cv,
                Ex = Ex, Yx = Yx, lasso = lasso))
}


