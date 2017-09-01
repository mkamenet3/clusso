#' Spatial and Spatio-Temporal Lasso
#' @title 
#' spacetimeLasso
#' 
#' @description 
#' This function runs the Lasso regularization technique on our large sparse matric of potential space or space-time clusters.
#' @param clusters clusters dataframe from (cluster.df function) that includes the center, x,y, r (radius), n (counter), and last (last observation in potential cluster)
#' @param vectors takes in the list of expected and observed counts from setVectors function
#' @param Time number of time periods in the dataset
#' @param spacetime indicator of whether the cluster detection method should be run on all space-time clusters(default) or on only the potential space clusters.
#' @param pois whether or not the Quasi-Poisson or Poisson model should be run. Default is pois=FALSE (default is Quasi-Poisson model is to be run)
#' @param overdispfloor default is TRUE. If TRUE, does not allow for underdispersion. If FALSE, allows for underdispersion (phi < 1)
#' @param cv option for cross-validation instead of AIC/BIC. Default is set to FALSE
#' @return This function will return a list with the expected counts as selected by QBIC, QAIC, QAICc, a list of original expected counts (Ex),
#' a list of observed counts (Yx), the lasso object, a list of K values (number of unique values in each decision path), and n (length of unique centers in the clusters dataframe)
#' @export
#' @examples 
#' potentialclusters <- clusters2df(lat, long, utm=FALSE, length(lat))
#' myvectors <- setVectors(period, expected, observed, Time, byrow=TRUE)
#' myresults <- spacetimeLasso(potentialclusters, myvectors, spacetime=TRUE, pois=FALSE)
#' 
spacetimeLasso<- function(clusters, vectors,Time, spacetime=TRUE,pois=FALSE,overdispfloor, cv){
    n <- length(unique(clusters$center))
    potClus <- n
    numCenters <- n
    if(spacetime == FALSE && nrow(vectors$covars) ==0){
        covars<- NULL
    }
    else{
        covars <- vectors$covars    
    }
    if(spacetime==TRUE){
        message("Creating space-time matrix")
        sparseMAT <- spacetimeMat(clusters, numCenters, Time)
        
        #set initial
        Ex <- vectors$Ex
        Yx <- vectors$Y.vec
        Period <- vectors$Period
    }
    else{
        message("Creating space-only matrix")
        sparseMAT <- spaceMat(clusters, numCenters)
        #set initial
        Ex <- as.vector(vectors$Ex)
        Yx <- as.vector(vectors$Y.vec)
        Period <- as.factor(as.vector(vectors$Period))
    }
    if(!is.null(covars)){
        str(covars)
        print(nrow(covars))
        print(isTRUE(nrow(covars)!=0))
        message("Running with covariates")
        covarMAT <- Matrix::Matrix(data.matrix(covars), sparse=TRUE)
        dim(sparseMAT)
        sparseMat <- Matrix::cBind(sparseMAT, covarMAT)
    }
    else{
        message("No covariates found")
    }
    print(paste("Number of potential clusters to scan through: ", dim(sparseMAT)[2]))
    message("Running Lasso - stay tuned")
    if(!is.null(cv)){
        message("Path selection: cross-validation")
        lasso <- glmnet::cv.glmnet(sparseMAT, Yx, family=("poisson"), alpha=1, offset=log(Ex), nlambda = 2000, 
                                   standardize = FALSE, intercept=FALSE,dfmax = 10, nfolds = cv) 
    }
    else{
        message("Path selection: information criteria")
        lasso <- glmnet::glmnet(sparseMAT, Yx, family=("poisson"), alpha=1, offset=log(Ex), nlambda = 2000, 
                                standardize = FALSE, intercept=FALSE,dfmax = 10) 
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
        
        #########################################################
        #Space-Time, Quasi-Poisson only (yes overdispersion)
        #########################################################
        if(spacetime==TRUE & pois == FALSE){
            message("returning results for space-time Quasi-Poisson model")
            if(!is.null(covars)){
                offset_reg <- glm(Yx ~ . + as.factor(vectors$Period) + offset(log(Ex)),
                                  data = covars,family=quasipoisson)
            }
            else{
                offset_reg <- glm(Yx ~ 1 + as.factor(vectors$Period) + offset(log(Ex)),
                                  family=quasipoisson)
            }
            overdisp.est <- overdisp(offset_reg, sim = FALSE, overdispfloor = overdispfloor)
            message(paste("Overdispersion estimate:", overdisp.est))
            #overdisp.est <- max(unlist(deviance(offset_reg)/df.residual(offset_reg)))
            message("Selecting best paths")
            
            #QBIC
            PLL.qbic  <- -2*(loglike/overdisp.est) + ((K)*log(n*Time))
            select.qbic <- which.min(PLL.qbic)
            E.qbic <- mu[,select.qbic]
            numclust.qbic <- length(unique(coefs.lasso.all[,select.qbic]))-1
            
            #QAIC
            PLL.qaic <-  2*(K) - 2*(loglike/overdisp.est)
            select.qaic <- which.min(PLL.qaic)
            E.qaic <- mu[,select.qaic]
            numclust.qaic <- length(unique(coefs.lasso.all[,select.qaic]))-1
            
            
            #QAICc
            PLL.qaicc <- 2*(K) - 2*(loglike/overdisp.est) +
                ((2*K*(K + 1))/(n*Time - K - 1))
            select.qaicc <- which.min(PLL.qaicc)
            E.qaicc <- mu[,select.qaic]
            numclust.qaicc <- length(unique(coefs.lasso.all[,select.qaicc]))-1
            
        }
        #########################################################
        #Space-Time, Poisson only (no overdispersion)
        #########################################################
        else if(spacetime==TRUE & pois==TRUE){
            message("returning results for space-time Poisson model")
            #QBIC
            PLL.qbic  <- -2*(loglike) + ((K)*log(n*Time))
            select.qbic <- which.min(PLL.qbic)
            E.qbic <- mu[,select.qbic]
            numclust.qbic <- length(unique(coefs.lasso.all[,select.qbic]))-1
            
            #QAIC
            PLL.qaic <-  2*(K) - 2*(loglike)
            select.qaic <- which.min(PLL.qaic)
            E.qaic <- mu[,select.qaic]
            numclust.qaic <- length(unique(coefs.lasso.all[,select.qaic]))-1
            
            
            #QAICc
            PLL.qaicc <- 2*(K) - 2*(loglike) +
                ((2*K*(K + 1))/(n*Time - K - 1))
            select.qaicc <- which.min(PLL.qaicc)
            E.qaicc <- mu[,select.qaic]
            numclust.qaicc <- length(unique(coefs.lasso.all[,select.qaicc]))-1
        }
        #########################################################
        #Space-Only, Quasi-Poisson
        #########################################################
        else if(spacetime==FALSE & pois==FALSE){
            message("Returning results for space-only  Quasi-Poisson model")
            if(!is.null(covars)){
                offset_reg <- glm(Yx ~ . + as.factor(vectors$Period) + offset(log(Ex)),data = covars,family=quasipoisson)
            }
            else{
                offset_reg <- glm(Yx ~ 1 + as.factor(vectors$Period) + offset(log(Ex)),family=quasipoisson)
            }
            offset_reg <- glm(Yx ~ 1 + offset(log(Ex)),family=poisson)
            #overdisp.est <- max(unlist(deviance(offset_reg)/df.residual(offset_reg)))
            overdisp.est <- overdisp(offset_reg, sim = FALSE, overdispfloor = overdispfloor)
            message(paste("Overdispersion estimate:", overdisp.est))
            message("Selecting best paths")
            if(pois==FALSE & is.null(overdisp.est)) warning("No overdispersion for quasi-Poisson model. Please check.")
            
            #QBIC
            PLL.qbic  <- -2*(loglike/overdisp.est) + ((K)*log(n*Time))
            select.qbic <- which.min(PLL.qbic)
            E.qbic <- mu[,select.qbic]
            numclust.qbic <- length(unique(coefs.lasso.all[,select.qbic]))-1
            
            #QAIC
            PLL.qaic <-  2*(K) - 2*(loglike/overdisp.est)
            select.qaic <- which.min(PLL.qaic)
            E.qaic <- mu[,select.qaic]
            numclust.qaic <- length(unique(coefs.lasso.all[,select.qaic]))-1
            
            
            #QAICc
            PLL.qaicc <- 2*(K) - 2*(loglike/overdisp.est) +
                ((2*K*(K + 1))/(n*Time - K - 1))
            select.qaicc <- which.min(PLL.qaicc)
            E.qaicc <- mu[,select.qaic]
            numclust.qaicc <- length(unique(coefs.lasso.all[,select.qaicc]))-1
        }
        #########################################################
        #Space-only, Poisson only
        #########################################################
        else if(spacetime==FALSE & pois == TRUE){
            message("Returning results for space-only  Poisson model")
            #QBIC
            PLL.qbic  <- -2*(loglike) + ((K)*log(n*Time))
            select.qbic <- which.min(PLL.qbic)
            E.qbic <- mu[,select.qbic]
            numclust.qbic <- length(unique(coefs.lasso.all[,select.qbic]))-1
            
            #QAIC
            PLL.qaic <-  2*(K) - 2*(loglike)
            select.qaic <- which.min(PLL.qaic)
            E.qaic <- mu[,select.qaic]
            numclust.qaic <- length(unique(coefs.lasso.all[,select.qaic]))-1
            
            
            #QAICc
            PLL.qaicc <- 2*(K) - 2*(loglike) +
                ((2*K*(K + 1))/(n*Time - K - 1))
            select.qaicc <- which.min(PLL.qaicc)
            E.qaicc <- mu[,select.qaic]
            numclust.qaicc <- length(unique(coefs.lasso.all[,select.qaicc]))-1
        }
        res <- list(E.qbic = E.qbic, E.qaic = E.qaic, E.qaicc = E.qaicc, numclust.qaic = numclust.qaic,
                    numclust.qaicc = numclust.qaicc, numclust.qbic= numclust.qbic, Ex = Ex, Yx = Yx, lasso = lasso, K = K)
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
#' #' @return list of expected values, number of clusters, Ex, Yx, and lasso object
#' 
stLasso_cv <- function(lasso, sparseMAT, Ex, Yx){
    #cv version
    ix <- which(lasso$lambda == lasso$lambda.min)
    print(dim(sparseMAT), dim(lasso$glmnet.fit$beta[,ix]))
    xbetaPath<- sparseMAT%*%lasso$glmnet.fit$beta[,ix]
    mu <- exp(xbetaPath)
    numclust.cv <- length(unique(mu@x))-1
    return(list(E.cv = mu[,1], numclust.cv = numclust.cv,
                Ex = Ex, Yx = Yx, lasso = lasso))
}


