#' Spatial and Spatio-Temporal Lasso
#' @title 
#' spacetimeLasso
#' 
#' @description 
#' This function runs the LASSOregularization technique on our large sparse matric of potential space or space-time clusters.
#'@param model A string specifying which model to use, Poisson or binomial. For Poisson, specify \code{"poisson"} and both the Poisson and quasi-Poisson model results are returned. 
#' @param sparseMAT large sparse matrix created in \code{clusso} function.
#' @param n_uniq number of unique polygons (ex: counties, zip code, etc). Inherited from \code{clusso}.
#' @param vectors takes in the list of expected and observed counts from setVectors function
#' @param Time number of time periods in the dataset
#' @param model A string specifying which model to use, Poisson or binomial. For Poisson, specify \code{"poisson"} and both the Poisson and quasi-Poisson model results are returned. For binomial, specify \code{"binomial"}.
#' @param quasi Whether or not the Quasi-Poisson or Poisson model should be run. Default is quasi=FALSE (default is Poisson model is to be run)
#' @param maxclust Upper limit on the maximum number of clusters you expect to find in the region. This equivalent to setting \code{dfmax} in the lasso. If none supplied, default is \code{11}.
#' @param overdispfloor default is TRUE. If TRUE, does not allow for underdispersion. If FALSE, allows for underdispersion (phi < 1)
#' @param cv option for cross-validation instead of AIC/BIC. Default is set to FALSE
#' @return This function will return a list with the expected counts as selected by QBIC, QAIC, QAICc, a list of original expected counts (Ex),
#' a list of observed counts (Yx), the lasso object, a list of K values (number of unique values in each decision path), and n (length of unique centers in the clusters dataframe)
#' @export

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
            # print(table(penalty))
        }
        else{
            penalty <- c(rep(1,(ncol(sparseMAT)-Time)), rep(0,Time))    
        }
        print("TODO")
        #lasso <- glmnet::cv.glmnet(sparseMAT, Yx, family=("poisson"), alpha=1, offset=log(Ex), nlambda = 2000, 

                   #                standardize = FALSE, intercept=FALSE,dfmax = maxclust, 
                #                   nfolds = cv, 
                 #                  penalty.factor = penalty) 

        #                standardize = FALSE, intercept=FALSE,dfmax = maxclust, 
        #                   nfolds = cv, 
        #                  penalty.factor = penalty) 
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
            print("Model poisson selected in spacetimeLasso")
            lasso <- glmnet::glmnet(sparseMAT, Yx, family=("poisson"), alpha=1, offset=log(Ex), 
                                    nlambda = 2000,
                                    standardize = FALSE, intercept=FALSE,dfmax = maxclust,
                                    penalty.factor = penalty)    
        }
        else if(model=="binomial"){
            print("Model binomial selected in spacetimeLasso")
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
        print("Pois post-process")
        print(paste0("quasi: ",quasi))
        if(model=="poisson"){
            mu <- sapply(1:length(lasso$lambda), function(i) exp(xbetaPath[,i]))
            loglike <- sapply(1:length(lasso$lambda), function(i) sum(dpoisson(Yx, mu[,i],Ex)))
            res <- spacetimeLassoPois(lasso, coefs.lasso.all,loglike,mu,K, quasi, 
                                      covars,Yx, Ex, Period, Time, n_uniq, overdispfloor, maxclust)
        }
        else if(model=="binomial"){
            print("Binom post-process")
            phat <-  sapply(1:length(lasso$lambda), function(i) pihat(xbetaPath[,i]))
            loglike <- sapply(1:length(lasso$lambda), function(i) sum(dbin(Yx, Ex, phat[,i])))
            res <- spacetimeLassoBinom(lasso,coefs.lasso.all, loglike, mu, K, covars,
                                       Yx, Ex, Period, Time, n_uniq, maxclust)
        }
        return(res)
    }
}


#'@title spacetimeLassoPois
#'@description
#'@param loglike
#'@param coefs.lasso.all
#'@param mu
#'@param K
#'@param quasi
#'@param covars
#'@param Yx
#'@param Ex
#'@param Period
#'@param Time
#'@param n_uniq
#'@param overdispfloor    
#'@param maxclust
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
        exp_coefs_qbic <- c(exp(unique(lasso$beta[,select.qbic])),
                            exp(lasso$beta[(nrow(lasso$beta)-Time+1):nrow(lasso$beta),select.qbic]))
        numclust.qbic <- K[select.qbic]-Time 
        
        #QAIC
        PLL.qaic <-  2*(K) - 2*(loglike/overdisp.est)
        select.qaic <- which.min(PLL.qaic)
        E.qaic <- mu[,select.qaic]
        exp_coefs_qaic <- c(exp(unique(lasso$beta[,select.qaic])),
                            exp(lasso$beta[(nrow(lasso$beta)-Time+1):nrow(lasso$beta),select.qaic]))
        numclust.qaic <- K[select.qaic]-Time 
        
        #QAICc
        PLL.qaicc <- 2*(K) - 2*(loglike/overdisp.est) +
            ((2*K*(K + 1))/(sum(Yx) - K - 1))
        select.qaicc <- which.min(PLL.qaicc)
        E.qaicc <- mu[,select.qaicc]
        exp_coefs_qaicc <- c(exp(unique(lasso$beta[,select.qaicc])),
                             exp(lasso$beta[(nrow(lasso$beta)-Time+1):nrow(lasso$beta),select.qaicc]))
        numclust.qaicc <- K[select.qaicc]-Time 
        
        selections <- list(select.qbic=select.qbic,
                           select.qaic = select.qaic,
                           select.qaicc = select.qaic)
    }
    #########################################################
    #Poisson only (no overdispersion)
    #########################################################
    else if(quasi==FALSE){
        #QBIC
        PLL.qbic  <- -2*(loglike) + ((K)*log(sum(Yx)))
        select.qbic <- which.min(PLL.qbic)
        E.qbic <- mu[,select.qbic]
        exp_coefs_qbic <- c(exp(unique(lasso$beta[,select.qbic])),
                            exp(lasso$beta[(nrow(lasso$beta)-Time+1):nrow(lasso$beta),select.qbic]))
        numclust.qbic <- K[select.qbic]-Time 
        
        #QAIC
        PLL.qaic <-  2*(K) - 2*(loglike)
        select.qaic <- which.min(PLL.qaic)
        E.qaic <- mu[,select.qaic]
        exp_coefs_qaic <- c(exp(unique(lasso$beta[,select.qaic])),
                            exp(lasso$beta[(nrow(lasso$beta)-Time+1):nrow(lasso$beta),select.qaic]))
        numclust.qaic <- K[select.qaic]-Time 
        
        #QAICc
        PLL.qaicc <- 2*(K) - 2*(loglike) +
            ((2*K*(K + 1))/(sum(Yx) - K - 1))
        select.qaicc <- which.min(PLL.qaicc)
        E.qaicc <- mu[,select.qaicc]
        exp_coefs_qaicc <- c(exp(unique(lasso$beta[,select.qaicc])),
                             exp(lasso$beta[(nrow(lasso$beta)-Time+1):nrow(lasso$beta),select.qaicc]))
        numclust.qaicc <- K[select.qaicc]-Time 
        
        selections <- list(select.qbic=select.qbic,
                           select.qaic = select.qaic,
                           select.qaicc = select.qaic)
    }
    #warning if numclust similar to maxclust
    if (numclust.qaic==(maxclust-Time)){
        message("The number of clusters selected by at least one criterion is equal to maxclust. You may want to increase maxclust.")
    }
    #Return only changepoints from lasso
    changepoints_ix <- which(diff(K)!=0) #Find lambda where new coef introduced
    lambda_changepoint <- lasso$lambda[changepoints_ix]
    #QIC
    coefs_qbic <- coefs.lasso.all[which(coefs.lasso.all[,select.qbic]!=0), changepoints_ix]
    coefs_qaic <-  coefs.lasso.all[which(coefs.lasso.all[,select.qaic]!=0), changepoints_ix]
    coefs_qaicc <- coefs.lasso.all[which(coefs.lasso.all[,select.qaicc]!=0), changepoints_ix]
    
    lasso_out <- list(
        lambdas = lambda_changepoint,
        coefs_bic = coefs_qbic,
        coefs_aic = coefs_qaic,
        coefs_aicc = coefs_qaicc
    )
    
    res <- list(E.qbic = E.qbic, E.qaic = E.qaic, E.qaicc = E.qaicc, numclust.qaic = numclust.qaic,
                numclust.qaicc = numclust.qaicc, numclust.qbic= numclust.qbic, Ex = Ex, Yx = Yx, 
                lasso = lasso, lasso_out=lasso_out,K = K, coefs.lasso.all = coefs.lasso.all,
                exp_coefs_qbic = exp_coefs_qbic, exp_coefs_qaic = exp_coefs_qaic, exp_coefs_qaicc = exp_coefs_qaicc,
                selections = selections)
}




#'@title spacetimeLassoBinom
#'@description
#'@param loglike
#'@param coefs.lasso.all
#'@param mu
#'@param K
#'@param covars
#'@param Yx
#'@param Ex
#'@param Period
#'@param Time
#'@param n_uniq
#'@param maxclust
spacetimeLassoBinom <- function(lasso, coefs.lasso.all, loglike, mu, K, covars, Yx, Ex, Period, Time, n_uniq, maxclust){    
    #########################################################
    #Space-Time, Binomial
    #########################################################
    #QBIC
    PLL.qbic  <- -2*(loglike) + ((K)*log(min(sum(Yx),sum(Ex-Yx))))
    select.qbic <- which.min(PLL.qbic)
    E.qbic <- mu[,select.qbic]
    exp_coefs_qbic <- c(exp(unique(lasso$beta[,select.qbic])),
                        exp(lasso$beta[(nrow(lasso$beta)-Time+1):nrow(lasso$beta),select.qbic]))
    numclust.qbic <- K[select.qbic]-Time 
    
    #QAIC
    PLL.qaic <-  2*(K) - 2*(loglike)
    select.qaic <- which.min(PLL.qaic)
    E.qaic <- mu[,select.qaic]
    exp_coefs_qaic <- c(exp(unique(lasso$beta[,select.qaic])),
                        exp(lasso$beta[(nrow(lasso$beta)-Time+1):nrow(lasso$beta),select.qaic]))
    numclust.qaic <- K[select.qaic]-Time 
    
    #QAICc
    PLL.qaicc <- 2*(K) - 2*(loglike) +
        ((2*K*(K + 1))/(min(sum(Yx),sum(Ex-Yx)) - K - 1))
    select.qaicc <- which.min(PLL.qaicc)
    E.qaicc <- mu[,select.qaicc]
    exp_coefs_qaicc <- c(exp(unique(lasso$beta[,select.qaicc])),
                         exp(lasso$beta[(nrow(lasso$beta)-Time+1):nrow(lasso$beta),select.qaicc]))
    numclust.qaicc <- K[select.qaicc]-Time 
    
    selections <- list(select.qbic=select.qbic,
                       select.qaic = select.qaic,
                       select.qaicc = select.qaic)
    
    #warning if numclust similar to maxclust
    if (numclust.qaic==(maxclust-Time)){
        message("The number of clusters selected by at least one criterion is equal to maxclust. You may want to increase maxclust.")
    }
    #Add if quasi=TRUE for quasi-binomial
    #warning if numclust similar to maxclust
    if (numclust.qaic==(maxclust-Time)){
        message("The number of clusters selected by at least one criterion is equal to maxclust. You may want to increase maxclust.")
    }
    #Return only changepoints from lasso
    changepoints_ix <- which(diff(K)!=0) #Find lambda where new coef introduced
    lambda_changepoint <- lasso$lambda[changepoints_ix]
    #QIC
    coefs_qbic <- coefs.lasso.all[which(coefs.lasso.all[,select.qbic]!=0), changepoints_ix]
    coefs_qaic <-  coefs.lasso.all[which(coefs.lasso.all[,select.qaic]!=0), changepoints_ix]
    coefs_qaicc <- coefs.lasso.all[which(coefs.lasso.all[,select.qaicc]!=0), changepoints_ix]
    
    lasso_out <- list(
        lambdas = lambda_changepoint,
        coefs_bic = coefs_qbic,
        coefs_aic = coefs_qaic,
        coefs_aicc = coefs_qaicc
    )
    
    res <- list(E.qbic = E.qbic, E.qaic = E.qaic, E.qaicc = E.qaicc, numclust.qaic = numclust.qaic,
                numclust.qaicc = numclust.qaicc, numclust.qbic= numclust.qbic, Ex = Ex, Yx = Yx, 
                lasso = lasso, lasso_out=lasso_out,K = K, coefs.lasso.all = coefs.lasso.all,
                exp_coefs_qbic = exp_coefs_qbic, exp_coefs_qaic = exp_coefs_qaic, exp_coefs_qaicc = exp_coefs_qaicc,
                selections = selections)
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


