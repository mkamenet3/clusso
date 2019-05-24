#' Spatial and Spatio-Temporal Lasso
#' @title 
#' spacetimeLasso
#' 
#' @description 
#' This function runs the Lasso regularization technique on our large sparse matric of potential space or space-time clusters.
#' @param sparseMAT large sparse matrix created in \code{clust} function.
#' @param n_uniq number of unique polygons (ex: counties, zip code, etc). Inherited from \code{clust}.
#' @param vectors takes in the list of expected and observed counts from setVectors function
#' @param Time number of time periods in the dataset
#' @param spacetime indicator of whether the cluster detection method should be run on all space-time clusters(default) or on only the potential space clusters.
#' @param pois whether or not the Quasi-Poisson or Poisson model should be run. Default is pois=FALSE (default is Quasi-Poisson model is to be run)
#' @param maxclust Upper limit on the maximum number of clusters you expect to find in the region. This equivalent to setting \code{dfmax} in the lasso. If none supplied, default is \code{10}.
#' @param overdispfloor default is TRUE. If TRUE, does not allow for underdispersion. If FALSE, allows for underdispersion (phi < 1)
#' @param cv option for cross-validation instead of AIC/BIC. Default is set to FALSE
#' @return This function will return a list with the expected counts as selected by QBIC, QAIC, QAICc, a list of original expected counts (Ex),
#' a list of observed counts (Yx), the lasso object, a list of K values (number of unique values in each decision path), and n (length of unique centers in the clusters dataframe)
#' @export
 
spacetimeLasso<- function(sparseMAT, n_uniq, vectors,Time, spacetime=TRUE,pois=FALSE,maxclust, overdispfloor, cv){
    #check for covariates
    covars <- vectors$covars
#    message(str(covars))
    if(!is.null(covars)){
        message("Running with covariates")
        covarMAT <- Matrix::Matrix(data.matrix(covars), sparse=TRUE)
        sparseMat <- cbind(sparseMAT, covarMAT)
    }
    else{
        message("No covariates found")
    }
    message(paste("Number of potential clusters to scan through: ", (dim(sparseMAT)[2]-Time)))
    #Set initial
    Ex <- vectors$Ex
    Yx <- vectors$Y.vec
    Period <- vectors$Period
    
    ############################################
    ####Create time matrix - not lasso'd
    # time_period <- factor(rep(1:Time, each=n_uniq))
    # timeMat <- Matrix(model.matrix(~ time_period - 1), sparse=TRUE)
    # #add this to sparsemat
    # sparseMAT <- cbind(sparseMAT, timeMat)
    ############################################
    #Run
    message("Running Lasso - stay tuned")
    if(!is.null(cv)){
        message("Path selection: cross-validation")
        penalty <- c(rep(1,(ncol(sparseMAT)-Time)), rep(0,Time))
        lasso <- glmnet::cv.glmnet(sparseMAT, Yx, family=("poisson"), alpha=1, offset=log(Ex), nlambda = 2000, 
                                   standardize = FALSE, intercept=FALSE,dfmax = maxclust, 
                                   #standardize = FALSE, intercept=FALSE,dfmax = 10, 
                                   nfolds = cv, 
                                   penalty.factor = penalty) 
    }
    else{
        message("Path selection: information criteria")
        #message(ncol(sparseMAT))
        penalty <- c(rep(1,(ncol(sparseMAT)-Time)), rep(0,Time))
       # message(str(penalty))
        #message(tail(penalty))
        lasso <- glmnet::glmnet(sparseMAT, Yx, family=("poisson"), alpha=1, offset=log(Ex), nlambda = 2000,
                                #standardize = FALSE, intercept=FALSE,dfmax = 10,
                                standardize = FALSE, intercept=FALSE,dfmax = maxclust,
                                penalty.factor = penalty)
        # lasso <- glmnet::glmnet(sparseMAT, Yx, family=("poisson"), alpha=1, offset=log(Ex), nlambda = 2000,
        #                         standardize = FALSE, intercept=FALSE,dfmax = 10,
        #                         exclude=(ncol(sparseMAT)-(Time-1)):ncol(sparseMAT))
    }
    message("Lasso complete - extracting estimates and paths")
    coefs.lasso.all <- coef(lasso)[-1,]
    xbetaPath<- sparseMAT%*%coefs.lasso.all
    
    #if running cross-validation version:
    if(!is.null(cv)){
        res <- stLasso_cv(lasso, sparseMAT, Ex, Yx)
    }
    #information criteria selection version:
    else{
        mu <- sapply(1:length(lasso$lambda), function(i) exp(xbetaPath[,i]))
        loglike <- sapply(1:length(lasso$lambda), function(i) sum(dpoisson(Yx, mu[,i],Ex)))
        K <- lasso$df
        message("Selecting best paths")
        
        #########################################################
        #Space-Time, Quasi-Poisson only (yes overdispersion)
        #########################################################
        if(spacetime==TRUE & pois == FALSE){
            message("Returning results for space-time Quasi-Poisson model")
            if(!is.null(covars)){
                offset_reg <- glm(Yx ~ . + as.factor(vectors$Period) + offset(log(Ex)),
                                  data = covars,family=quasipoisson)
            }
            else{
                offset_reg <- glm(Yx ~ 1 + as.factor(vectors$Period) + offset(log(Ex)),
                                  family=quasipoisson)
            }
            overdisp.est <- overdisp(offset_reg, sim = FALSE, overdispfloor = overdispfloor)
            message(paste("Overdispersion estimate:", round(overdisp.est,4)))
            if(pois == FALSE & is.null(overdisp.est)) warning("No overdispersion for quasi-Poisson model. Please check.")
            
            #QBIC
            PLL.qbic  <- -2*(loglike/overdisp.est) + ((K)*log(n_uniq*Time))
            select.qbic <- which.min(PLL.qbic)
            E.qbic <- mu[,select.qbic]
            #numclust.qbic <- length(unique(coefs.lasso.all[,select.qbic]))-1
            #numclust.qbic <- sort(c(exp(unique(lasso$beta[,select.qbic])),exp(lasso$beta[66871:66875,select.qbic])))
            #message(print((nrow(lasso$beta)-Time+1):nrow(lasso$beta)))
            exp_coefs_qbic <- c(exp(unique(lasso$beta[,select.qbic])),
                                    exp(lasso$beta[(nrow(lasso$beta)-Time+1):nrow(lasso$beta),select.qbic]))
            numclust.qbic <- K[select.qbic]-Time #length(unique(exp_coefs_qbic))-Time
            #test <- t(sparseMAT)%*%xbetaPath[,select.qbic]
            
            
            #QAIC
            PLL.qaic <-  2*(K) - 2*(loglike/overdisp.est)
            select.qaic <- which.min(PLL.qaic)
            E.qaic <- mu[,select.qaic]
            #numclust.qaic <- length(unique(exp_coefs_lasso.all[,select.qaic]))-1
            # numclust.qaic <- sort(c(exp(unique(lasso$beta[,select.qaic])),
            #                         exp(lasso$beta[66871:66875,select.qaic])))
            exp_coefs_qaic <- c(exp(unique(lasso$beta[,select.qaic])),
                                    exp(lasso$beta[(nrow(lasso$beta)-Time+1):nrow(lasso$beta),select.qaic]))
            numclust.qaic <- K[select.qaic]-Time #length(unique(exp_coefs_qaic)) - Time
            
            
            #QAICc
            PLL.qaicc <- 2*(K) - 2*(loglike/overdisp.est) +
                ((2*K*(K + 1))/(n_uniq*Time - K - 1))
            select.qaicc <- which.min(PLL.qaicc)
            E.qaicc <- mu[,select.qaicc]
            #numclust.qaicc <- length(unique(coefs.lasso.all[,select.qaicc]))-1
            #numclust.qaicc <- length(unique(coefs.lasso.all[,select.qaicc]))-1
            exp_coefs_qaicc <- c(exp(unique(lasso$beta[,select.qaicc])),
                                    exp(lasso$beta[(nrow(lasso$beta)-Time+1):nrow(lasso$beta),select.qaicc]))
            numclust.qaicc <- K[select.qaicc]-Time #length(unique(exp_coefs_qaicc))-Time
            
            selections <- list(select.qbic=select.qbic,
                               select.qaic = select.qaic,
                               select.qaicc = select.qaic)
            
        }
        #########################################################
        #Space-Time, Poisson only (no overdispersion)
        #########################################################
        else if(spacetime==TRUE & pois==TRUE){
            message("Returning results for space-time Poisson model")
            #QBIC
            PLL.qbic  <- -2*(loglike) + ((K)*log(n_uniq*Time))
            select.qbic <- which.min(PLL.qbic)
            E.qbic <- mu[,select.qbic]
            exp_coefs_qbic <- c(exp(unique(lasso$beta[,select.qbic])),
                           exp(lasso$beta[(nrow(lasso$beta)-Time+1):nrow(lasso$beta),select.qbic]))
            numclust.qbic <- K[select.qbic]-Time #length(unique(exp_coefs_qbic)) - Time
            #numclust.qbic <- length(unique(coefs.lasso.all[,select.qbic]))-1
            
            #QAIC
            PLL.qaic <-  2*(K) - 2*(loglike)
            select.qaic <- which.min(PLL.qaic)
            E.qaic <- mu[,select.qaic]
            #numclust.qaic <- length(unique(coefs.lasso.all[,select.qaic]))-1
            exp_coefs_qaic <- c(exp(unique(lasso$beta[,select.qaic])),
                            exp(lasso$beta[(nrow(lasso$beta)-Time+1):nrow(lasso$beta),select.qaic]))
            numclust.qaic <- K[select.qaic]-Time #length(unique(exp_coefs_qaic)) - Time
            
            
            #QAICc
            PLL.qaicc <- 2*(K) - 2*(loglike) +
                ((2*K*(K + 1))/(n_uniq*Time - K - 1))
            select.qaicc <- which.min(PLL.qaicc)
            E.qaicc <- mu[,select.qaicc]
            #numclust.qaicc <- length(unique(coefs.lasso.all[,select.qaicc]))-1
            exp_coefs_qaicc <- c(exp(unique(lasso$beta[,select.qaicc])),
                            exp(lasso$beta[(nrow(lasso$beta)-Time+1):nrow(lasso$beta),select.qaicc]))
            numclust.qaicc <- K[select.qaicc]-Time #length(unique(exp_coefs_qaicc)) - Time
            
            selections <- list(select.qbic=select.qbic,
                               select.qaic = select.qaic,
                               select.qaicc = select.qaic)
            
        }
        #########################################################
        #Space-Only, Quasi-Poisson
        #########################################################
        else if(spacetime==FALSE & pois==FALSE){
            message("Returning results for space-only Quasi-Poisson model")
            if(!is.null(covars)){
                offset_reg <- glm(Yx ~ . + offset(log(Ex)),data = covars,family=quasipoisson)
            }
            else{
                offset_reg <- glm(Yx ~ 1 + offset(log(Ex)),family=quasipoisson)
            }
            offset_reg <- glm(Yx ~ 1 + offset(log(Ex)),family=poisson)
            overdisp.est <- overdisp(offset_reg, sim = FALSE, overdispfloor = overdispfloor)
            message(paste("Overdispersion estimate:", round(overdisp.est,4)))
            if(pois==FALSE & is.null(overdisp.est)) warning("No overdispersion for quasi-Poisson model. Please check.")
            
            #QBIC
            PLL.qbic  <- -2*(loglike/overdisp.est) + ((K)*log(n_uniq*Time))
            select.qbic <- which.min(PLL.qbic)
            E.qbic <- mu[,select.qbic]
            #numclust.qbic <- length(unique(coefs.lasso.all[,select.qbic]))-1
            exp_coefs_qbic <- c(exp(unique(lasso$beta[,select.qbic])),
                            exp(lasso$beta[(nrow(lasso$beta)-Time+1):nrow(lasso$beta),select.qbic]))
            numclust.qbic <- K[select.qbic]-Time #length(unique(exp_coefs_qbic)) - Time
            
            #QAIC
            PLL.qaic <-  2*(K) - 2*(loglike/overdisp.est)
            select.qaic <- which.min(PLL.qaic)
            E.qaic <- mu[,select.qaic]
            #numclust.qaic <- length(unique(coefs.lasso.all[,select.qaic]))-1
            exp_coefs_qaic <- c(exp(unique(lasso$beta[,select.qaic])),
                            exp(lasso$beta[(nrow(lasso$beta)-Time+1):nrow(lasso$beta),select.qaic]))
            numclust.qaic <- K[select.qaic]-Time #length(unique(exp_coefs_qaic)) - Time
            
            
            #QAICc
            PLL.qaicc <- 2*(K) - 2*(loglike/overdisp.est) +
                ((2*K*(K + 1))/(n_uniq*Time - K - 1))
            select.qaicc <- which.min(PLL.qaicc)
            E.qaicc <- mu[,select.qaicc]
            #numclust.qaicc <- length(unique(coefs.lasso.all[,select.qaicc]))-1
            exp_coefs_qaicc <- c(exp(unique(lasso$beta[,select.qaicc])),
                            exp(lasso$beta[(nrow(lasso$beta)-Time+1):nrow(lasso$beta),select.qaicc]))
            numclust.qaicc <- K[select.qaicc]-Time  #length(unique(exp_coefs_qaicc)) - Time
            
            selections <- list(select.qbic=select.qbic,
                               select.qaic = select.qaic,
                               select.qaicc = select.qaic)
            
        }
        #########################################################
        #Space-only, Poisson only
        #########################################################
        else if(spacetime==FALSE & pois == TRUE){
            message("Returning results for space-only  Poisson model")
            #QBIC
            PLL.qbic  <- -2*(loglike) + ((K)*log(n_uniq*Time))
            select.qbic <- which.min(PLL.qbic)
            E.qbic <- mu[,select.qbic]
            #numclust.qbic <- length(unique(coefs.lasso.all[,select.qbic]))-1
            exp_coefs_qbic <- c(exp(unique(lasso$beta[,select.qbic])),
                            exp(lasso$beta[(nrow(lasso$beta)-Time+1):nrow(lasso$beta),select.qbic]))
            numclust.qbic <- K[select.qbic]-Time #length(unique(exp_coefs_qbic)) - Time
            
            #QAIC
            PLL.qaic <-  2*(K) - 2*(loglike)
            select.qaic <- which.min(PLL.qaic)
            E.qaic <- mu[,select.qaic]
            #numclust.qaic <- length(unique(coefs.lasso.all[,select.qaic]))-1
            exp_coefs_qaic <- c(exp(unique(lasso$beta[,select.qaic])),
                            exp(lasso$beta[(nrow(lasso$beta)-Time+1):nrow(lasso$beta),select.qaic]))
            numclust.qaic <- K[select.qaic]-Time #length(unique(exp_coefs_qaic)) - Time
            
            #QAICc
            PLL.qaicc <- 2*(K) - 2*(loglike) +
                ((2*K*(K + 1))/(n_uniq*Time - K - 1))
            select.qaicc <- which.min(PLL.qaicc)
            E.qaicc <- mu[,select.qaicc]
            #numclust.qaicc <- length(unique(coefs.lasso.all[,select.qaicc]))-1
            exp_coefs_qaicc <- c(exp(unique(lasso$beta[,select.qaicc])),
                            exp(lasso$beta[(nrow(lasso$beta)-Time+1):nrow(lasso$beta),select.qaicc]))
            numclust.qaicc <- K[select.qaicc]-Time #length(unique(exp_coefs_qaicc)) - Time
            
            selections <- list(select.qbic=select.qbic,
                               select.qaic = select.qaic,
                               select.qaicc = select.qaic)
            
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
                    penalty = penalty,
                    selections = selections)
    }
    return(res)
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


