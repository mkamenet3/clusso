#' Create pretty table from \code{clusso} output    
#' 
#'@title 
#'clussopretty
#'@description 
#'This function takes the cluster detection results from the \code{clusso} function and creates a nice table of detected clusters.  
#'@param outclusso Object with output from \code{clusso}.
#'@param analysis A string specifying if the spatial (\code{"space"}), spatio-temporal (\code{"spacetime"}), or both spatial and spatio-temporal (\code{"both"}) analysis should be executed.  
#'@param model A string specifying which model to use, Poisson or binomial. For Poisson, specify \code{"poisson"} and both the Poisson and quasi-Poisson model results are returned. For binomial, specify \code{"binomial"} and both the binomial and quasi-binomial model results are returned.
#'@param clusteridentify Whether specific clusters should be identified in output; default is \code{FALSE}.
#'@param clusterRR Risk ratio cut off for a cluster to be identified. Alternatively, a risk ratio can be provided (ex: 1.0, 1.25) that would serve as a cut off value for identifiying a cluster. Current default is 1.0. We define a cluster as being different from the background rate. 
#'@param cv Boolean, default is \code{FALSE}. If \code{TRUE}, then this will assume that the cross-validation (and not information criteria) model selection was performed.
#'@param covars Default is FALSE (time main effects and covariate estimates suppressed). If \code{TRUE}, then results reported for time main effects and covariate estimates. Results reported as both raw and exponentiated estimates.
#'@return Data frame of output that can be sent to a comma-separated values file (csv) or list of data frames (if \code{clusteridentify} set to \code{TRUE}).
#'@export
#'@examples
#'\donttest{
#'#load data
#'data("jbc")
#'data("utmJapan")
#'data("japan.poly2")
#'data("japan.prefect2")
#'#Initial inputs
#'x <- utmJapan$utmx/1000
#'y <- utmJapan$utmy/1000
#'rMax <- 20 
#'system.time(resreal <- clusso(df=jbc, expected = expdeath, observed=death,
#'    timeperiod = factor(period), covars=FALSE,x= x,y = y, rMax =  rMax, 
#'    utm=TRUE, analysis="both", model="poisson",maxclust=11))
#'clussopretty(resreal, analysis="both", model="poisson",clusteridentify=FALSE)
#'
#'#example with including covariates
#'system.time(resreal <- clusso(df=jbc, expected = expdeath, observed=death,
#'    timeperiod = factor(period), covars=TRUE, id=id,x= x,y = y, rMax =  rMax, 
#'    utm=TRUE, analysis="both", model="poisson",maxclust=11))
#'}

clussopretty <- function(outclusso, analysis="both", model = c("poisson", "binomial"), clusteridentify=FALSE, clusterRR, cv=FALSE, covars=FALSE){
    err <- 1e-4
    if (length(model)>1){
        stop("You must select either `poisson` or `binomial`")
    }
    else if(model=="poisson"){
        model <- c("Poisson", "Quasi-Poisson")
    }
    else if(model=="binomial"){
        model <- c("Binomial", "Quasi-Binomial")
    }
    else {stop("Unknown model type. If you think this was by error, please submit an issue.")}
    if(cv==TRUE){
    #cv version
        if(analysis=="space"){
            #model <- c("Poisson", "Quasi-Poisson")
            analysistype <- rep("Space",2)
            numclust.cv <- c(outclusso$lassoresult.p.s$numclust.cv,outclusso$lassoresult.qp.s$numclust.cv)
        }
        else if(analysis=="spacetime"){
            #model <- c("Poisson", "Quasi-Poisson")
            analysistype <- rep("Space-Time",2)
            numclust.cv <- c(outclusso$lassoresult.p.st$numclust.cv, outclusso$lassoresult.qp.st$numclust.cv)
        }
        else{
            #both
            #model <- c(rep("Poisson",2), rep("Quasi-Poisson",2))
            model <- rep(model, each=2)
            analysistype <- rep(c("Space", "Space-Time"),2)
            numclust.cv <- c(outclusso$lassoresult.p.s$numclust.cv, outclusso$lassoresult.p.st$numclust.cv, 
                              outclusso$lassoresult.qp.s$numclust.cv, outclusso$lassoresult.qp.st$numclust.cv)
        }
        table.clusters <- cbind.data.frame(model, analysistype, numclust.cv)
    }
    else{
        if(analysis=="space"){
            #model <- c("Poisson", "Quasi-Poisson")
            analysistype <- rep("Space",2)
            numclust.AIC <- c(outclusso$lassoresult.p.s$numclust.qaic,outclusso$lassoresult.qp.s$numclust.qaic)
            numclust.AICc <- c(outclusso$lassoresult.p.s$numclust.qaicc, outclusso$lassoresult.qp.s$numclust.qaicc)
            numclust.BIC <- c(outclusso$lassoresult.p.s$numclust.qbic, outclusso$lassoresult.qp.s$numclust.qbic)
        }
        else if(analysis=="spacetime"){
            #model <- c("Poisson", "Quasi-Poisson")
            analysistype <- rep("Space-Time",2)
            numclust.AIC <- c(outclusso$lassoresult.p.st$numclust.qaic, outclusso$lassoresult.qp.st$numclust.qaic)
            numclust.AICc <- c(outclusso$lassoresult.p.st$numclust.qaicc, outclusso$lassoresult.qp.st$numclust.qaicc)
            numclust.BIC <- c(outclusso$lassoresult.p.st$numclust.qbic, outclusso$lassoresult.qp.st$numclust.qbic)
        }
        else{
            #both
            model <- rep(model, each=2)
            analysistype <- rep(c("Space", "Space-Time"),2)
            numclust.AIC <- c(outclusso$lassoresult.p.s$numclust.qaic, outclusso$lassoresult.p.st$numclust.qaic, 
                              outclusso$lassoresult.qp.s$numclust.qaic, outclusso$lassoresult.qp.st$numclust.qaic)
            numclust.AICc <- c(outclusso$lassoresult.p.s$numclust.qaicc, outclusso$lassoresult.p.st$numclust.qaicc, 
                               outclusso$lassoresult.qp.s$numclust.qaicc, outclusso$lassoresult.qp.st$numclust.qaicc)
            numclust.BIC <- c(outclusso$lassoresult.p.s$numclust.qbic, outclusso$lassoresult.p.st$numclust.qbic, 
                              outclusso$lassoresult.qp.s$numclust.qbic, outclusso$lassoresult.qp.st$numclust.qbic)
        }
        table.clusters <- cbind.data.frame(model, analysistype, numclust.AIC, numclust.AICc, numclust.BIC)
    }
    ####
    if(covars==TRUE){
        if(analysis=="space"){
            analysistype <- rep("Space",2)
            #AIC
            ix.AIC <- which(names(outclusso$lassoresult.p.s$coefs.lasso.all[,outclusso$lassoresult.p.s$selections$select.qaic]) != "")
            ix.AIC.q <- which(names(outclusso$lassoresult.qp.s$coefs.lasso.all[,outclusso$lassoresult.qp.s$selections$select.qaic]) != "")
            coefs.AIC <- rbind(outclusso$lassoresult.p.s$coefs.lasso.all[ix.AIC, outclusso$lassoresult.p.s$selections$select.qaic],
                               outclusso$lassoresult.qp.s$coefs.lasso.all[ix.AIC.q, outclusso$lassoresult.qp.s$selections$select.qaic])
            #AICc
            ix.AICc <- which(names(outclusso$lassoresult.p.s$coefs.lasso.all[,outclusso$lassoresult.p.s$selections$select.qaicc]) != "")
            ix.AICc.q <- which(names(outclusso$lassoresult.qp.s$coefs.lasso.all[,outclusso$lassoresult.qp.s$selections$select.qaicc]) != "")
            coefs.AICc <- rbind(outclusso$lassoresult.p.s$coefs.lasso.all[ix.AICc, outclusso$lassoresult.p.s$selections$select.qaicc],
                               outclusso$lassoresult.qp.s$coefs.lasso.all[ix.AICc.q, outclusso$lassoresult.qp.s$selections$select.qaicc])
            #BIC
            ix.BIC <- which(names(outclusso$lassoresult.p.s$coefs.lasso.all[,outclusso$lassoresult.p.s$selections$select.qbic]) != "")
            ix.BIC.q <- which(names(outclusso$lassoresult.qp.s$coefs.lasso.all[,outclusso$lassoresult.qp.s$selections$select.qbic]) != "")
            coefs.BIC <- rbind(outclusso$lassoresult.p.s$coefs.lasso.all[ix.BIC, outclusso$lassoresult.p.s$selections$select.qbic],
                               outclusso$lassoresult.qp.s$coefs.lasso.all[ix.BIC.q, outclusso$lassoresult.qp.s$selections$select.qbic])
        }
        #table.coefs <- cbind(model, analysistype,rbind.data.frame(coefs.AIC, coefs.AICc, coefs.BIC))
        #outclusso
        #ix <- which(names(resrealcovars$lassoresult.p.st$coefs.lasso.all[,resrealcovars$lassoresult.p.st$selections$select.qaic]) != "")
        else if(analysis=="spacetime"){
            analysistype <- rep("Spacetime",2)
            #AIC
            ix.AIC <- which(names(outclusso$lassoresult.p.st$coefs.lasso.all[,outclusso$lassoresult.p.st$selections$select.qaic]) != "")
            ix.AIC.q <- which(names(outclusso$lassoresult.qp.st$coefs.lasso.all[,outclusso$lassoresult.qp.st$selections$select.qaic]) != "")
            coefs.AIC <- rbind(outclusso$lassoresult.p.st$coefs.lasso.all[ix.AIC, outclusso$lassoresult.p.st$selections$select.qaic],
                               outclusso$lassoresult.qp.st$coefs.lasso.all[ix.AIC.q, outclusso$lassoresult.qp.st$selections$select.qaic])
            #AICc
            ix.AICc <- which(names(outclusso$lassoresult.p.st$coefs.lasso.all[,outclusso$lassoresult.p.st$selections$select.qaicc]) != "")
            ix.AICc.q <- which(names(outclusso$lassoresult.qp.st$coefs.lasso.all[,outclusso$lassoresult.qp.st$selections$select.qaicc]) != "")
            coefs.AICc <- rbind(outclusso$lassoresult.p.st$coefs.lasso.all[ix.AICc, outclusso$lassoresult.p.st$selections$select.qaicc],
                                outclusso$lassoresult.qp.st$coefs.lasso.all[ix.AICc.q, outclusso$lassoresult.qp.st$selections$select.qaicc])
            #BIC
            ix.BIC <- which(names(outclusso$lassoresult.p.st$coefs.lasso.all[,outclusso$lassoresult.p.st$selections$select.qbic]) != "")
            ix.BIC.q <- which(names(outclusso$lassoresult.qp.st$coefs.lasso.all[,outclusso$lassoresult.qp.st$selections$select.qbic]) != "")
            coefs.BIC <- rbind(outclusso$lassoresult.p.st$coefs.lasso.all[ix.BIC, outclusso$lassoresult.p.st$selections$select.qbic],
                               outclusso$lassoresult.qp.st$coefs.lasso.all[ix.BIC.q, outclusso$lassoresult.qp.st$selections$select.qbic])
        }
        else{
            model <- rep(rep(model, each=2),times=3)
            analysistype <- rep(c("Space", "Space-Time"),6)
            #AIC
            ix.AIC <- which(names(outclusso$lassoresult.p.s$coefs.lasso.all[,outclusso$lassoresult.p.s$selections$select.qaic]) != "")
            ix.AIC.q <- which(names(outclusso$lassoresult.qp.s$coefs.lasso.all[,outclusso$lassoresult.qp.s$selections$select.qaic]) != "")
            ixt.AIC <- which(names(outclusso$lassoresult.p.st$coefs.lasso.all[,outclusso$lassoresult.p.st$selections$select.qaic]) != "")
            ixt.AIC.q <- which(names(outclusso$lassoresult.qp.st$coefs.lasso.all[,outclusso$lassoresult.qp.st$selections$select.qaic]) != "")
            
            coefs.AIC <- rbind(outclusso$lassoresult.p.s$coefs.lasso.all[ix.AIC, outclusso$lassoresult.p.s$selections$select.qaic],
                               outclusso$lassoresult.qp.s$coefs.lasso.all[ix.AIC.q, outclusso$lassoresult.qp.s$selections$select.qaic],
                               outclusso$lassoresult.p.st$coefs.lasso.all[ixt.AIC, outclusso$lassoresult.p.st$selections$select.qaic],
                               outclusso$lassoresult.qp.st$coefs.lasso.all[ixt.AIC.q, outclusso$lassoresult.qp.st$selections$select.qaic])
            #AICc
            ix.AICc <- which(names(outclusso$lassoresult.p.s$coefs.lasso.all[,outclusso$lassoresult.p.s$selections$select.qaicc]) != "")
            ix.AICc.q <- which(names(outclusso$lassoresult.qp.s$coefs.lasso.all[,outclusso$lassoresult.qp.s$selections$select.qaicc]) != "")
            ixt.AICc <- which(names(outclusso$lassoresult.p.st$coefs.lasso.all[,outclusso$lassoresult.p.st$selections$select.qaicc]) != "")
            ixt.AICc.q <- which(names(outclusso$lassoresult.qp.st$coefs.lasso.all[,outclusso$lassoresult.qp.st$selections$select.qaicc]) != "")
            
            coefs.AICc <- rbind(outclusso$lassoresult.p.s$coefs.lasso.all[ix.AICc, outclusso$lassoresult.p.s$selections$select.qaicc],
                                outclusso$lassoresult.qp.s$coefs.lasso.all[ix.AICc.q, outclusso$lassoresult.qp.s$selections$select.qaicc],
                                outclusso$lassoresult.p.st$coefs.lasso.all[ixt.AICc, outclusso$lassoresult.p.st$selections$select.qaicc],
                                outclusso$lassoresult.qp.st$coefs.lasso.all[ixt.AICc.q, outclusso$lassoresult.qp.st$selections$select.qaicc])
            #BIC
            ix.BIC <- which(names(outclusso$lassoresult.p.s$coefs.lasso.all[,outclusso$lassoresult.p.s$selections$select.qbic]) != "")
            ix.BIC.q <- which(names(outclusso$lassoresult.qp.s$coefs.lasso.all[,outclusso$lassoresult.qp.s$selections$select.qbic]) != "")
            ixt.BIC <- which(names(outclusso$lassoresult.p.st$coefs.lasso.all[,outclusso$lassoresult.p.st$selections$select.qbic]) != "")
            ixt.BIC.q <- which(names(outclusso$lassoresult.qp.st$coefs.lasso.all[,outclusso$lassoresult.qp.st$selections$select.qbic]) != "")
            
            coefs.BIC <- rbind(outclusso$lassoresult.p.s$coefs.lasso.all[ix.BIC, outclusso$lassoresult.p.s$selections$select.qbic],
                               outclusso$lassoresult.qp.s$coefs.lasso.all[ix.BIC.q, outclusso$lassoresult.qp.s$selections$select.qbic],
                               outclusso$lassoresult.p.st$coefs.lasso.all[ixt.BIC, outclusso$lassoresult.p.st$selections$select.qbic],
                               outclusso$lassoresult.qp.st$coefs.lasso.all[ixt.BIC.q, outclusso$lassoresult.qp.st$selections$select.qbic])
  
            
        }
        table.coefs <- cbind(IC=rep(c("(Q)AIC", "(Q)AICc", "(Q)BIC"),each=4),model, analysistype,rbind.data.frame(coefs.AIC, coefs.AICc, coefs.BIC))
        table.expcoefs <- cbind(IC=rep(c("(Q)AIC", "(Q)AICc", "(Q)BIC"),each=4),model, analysistype,exp(rbind.data.frame(coefs.AIC, coefs.AICc, coefs.BIC)))
        return(list(table.clusters = table.clusters,
                    table.coefs= table.coefs,
                    table.expcoefs = table.expcoefs))
    }
    
    ####
    if(clusteridentify==FALSE){
        return(table.clusters)
    }
    if(clusteridentify==TRUE){
        if(missing(clusterRR)){
            message("Using default risk ratio threshold for cluster of > 1.0. To specify a different risk ratio threshold, specify clusterRR. To specify different from background, specify clusterRR='background'")
            clusterRR <- 1.0
        }
        else if (clusterRR=="background"){
            clusterRR <- "background" 
            #as.numeric(apply(clustout$riskratios$riskratios.qp.s$RRbic,2,function(x) names(sort(table(x),decreasing=TRUE))[[1]]))
        }
        else{
            clusterRR <- clusterRR
        }
        #todo set default to background
        if(analysis=="space"){
            varnames <- paste0("Period",1:length(unique(outclusso$init.vec.s$Period))) 
            if(clusterRR=="background"){
                #background rates
                clusterRRbic.qp <- as.numeric(apply(outclusso$riskratios$riskratios.qp.s$RRbic,2,
                                 function(x) names(sort(table(x),decreasing=TRUE))[[1]]))
                clusterRRaic.qp <- as.numeric(apply(outclusso$riskratios$riskratios.qp.s$RRaic,2,
                                                 function(x) names(sort(table(x),decreasing=TRUE))[[1]]))
                clusterRRaicc.qp <- as.numeric(apply(outclusso$riskratios$riskratios.qp.s$RRaicc,2,
                                                 function(x) names(sort(table(x),decreasing=TRUE))[[1]]))
                clusterRRbic.p <- as.numeric(apply(outclusso$riskratios$riskratios.p.s$RRbic,2,
                                                    function(x) names(sort(table(x),decreasing=TRUE))[[1]]))
                clusterRRaic.p <- as.numeric(apply(outclusso$riskratios$riskratios.p.s$RRaic,2,
                                                    function(x) names(sort(table(x),decreasing=TRUE))[[1]]))
                clusterRRaicc.p <- as.numeric(apply(outclusso$riskratios$riskratios.p.s$RRaicc,2,
                                                     function(x) names(sort(table(x),decreasing=TRUE))[[1]]))
                #calc identified
                #identBIC.qp <- lapply(1:ncol(outclusso$riskratios$riskratios.qp.s$RRbic),
                 #                     function(i) which(round(outclusso$riskratios$riskratios.qp.s$RRbic[,i],4)>clusterRRbic.qp[[i]])) 
                identBIC.qp <- lapply(1:ncol(outclusso$riskratios$riskratios.qp.s$RRbic),
                                      function(i) which(abs(outclusso$riskratios$riskratios.qp.s$RRbic[,i] - clusterRRbic.qp[[i]]) > err)) 
                
                
                identBIC.p <- lapply(1:ncol(outclusso$riskratios$riskratios.p.s$RRbic),
                                     function(i) which(abs(outclusso$riskratios$riskratios.p.s$RRbic[,i]-clusterRRbic.p[[i]])> err)) 
                identAIC.qp <- lapply(1:ncol(outclusso$riskratios$riskratios.qp.s$RRaic),
                                      function(i) which(abs(outclusso$riskratios$riskratios.qp.s$RRaic[,i]-clusterRRaic.qp[[i]])> err)) 
                identAIC.p <- lapply(1:ncol(outclusso$riskratios$riskratios.p.s$RRaic),
                                     function(i) which(abs(outclusso$riskratios$riskratios.p.s$RRaic[,i]-clusterRRaic.p[[i]])> err)) 
                identAICc.qp <- lapply(1:ncol(outclusso$riskratios$riskratios.qp.s$RRaicc),
                                       function(i) which(abs(outclusso$riskratios$riskratios.qp.s$RRaicc[,i]-clusterRRaicc.qp[[i]])> err)) 
                identAICc.p <- lapply(1:ncol(outclusso$riskratios$riskratios.p.s$RRaicc),
                                      function(i) which(abs(outclusso$riskratios$riskratios.p.s$RRaicc[,i]>clusterRRaicc.p[[i]])> err)) 
            }
            else{
                identBIC.qp <- lapply(1:ncol(outclusso$riskratios$riskratios.qp.s$RRbic),
                                      function(i) which(outclusso$riskratios$riskratios.qp.s$RRbic[,i]>clusterRR)) 
                identBIC.p <- lapply(1:ncol(outclusso$riskratios$riskratios.p.s$RRbic),
                                     function(i) which(outclusso$riskratios$riskratios.p.s$RRbic[,i]>clusterRR)) 
                identAIC.qp <- lapply(1:ncol(outclusso$riskratios$riskratios.qp.s$RRaic),
                                      function(i) which(outclusso$riskratios$riskratios.qp.s$RRaic[,i]>clusterRR)) 
                identAIC.p <- lapply(1:ncol(outclusso$riskratios$riskratios.p.s$RRaic),
                                     function(i) which(outclusso$riskratios$riskratios.p.s$RRaic[,i]>clusterRR)) 
                identAICc.qp <- lapply(1:ncol(outclusso$riskratios$riskratios.qp.s$RRaicc),
                                       function(i) which(outclusso$riskratios$riskratios.qp.s$RRaicc[,i]>clusterRR)) 
                identAICc.p <- lapply(1:ncol(outclusso$riskratios$riskratios.p.s$RRaicc),
                                      function(i) which(outclusso$riskratios$riskratios.p.s$RRaicc[,i]>clusterRR)) 
            }
            QBIC = as.data.frame(stringi::stri_list2matrix(identBIC.qp)); colnames(QBIC) <- varnames
            BIC =  as.data.frame(stringi::stri_list2matrix(identBIC.p)); colnames(BIC) <- varnames
            QAIC = as.data.frame(stringi::stri_list2matrix(identAIC.qp)); colnames(QAIC) <- varnames
            AIC = as.data.frame(stringi::stri_list2matrix(identAIC.p)); colnames(AIC) <- varnames
            QAICc = as.data.frame(stringi::stri_list2matrix(identAICc.qp)); colnames(QAICc) <- varnames
            AICc = as.data.frame(stringi::stri_list2matrix(identAICc.p)); colnames(AICc) <- varnames
            identQP <- list(QBIC = QBIC,
                            QAIC = QAIC, 
                            QAICc = QAICc)
            identP <- list(BIC = BIC, 
                           AIC = AIC, 
                           AICc = AICc)
            table.identified <- list(identQuasiPoisson = identQP, identPoisson=identP)
            return(list(table.clusters = table.clusters, table.identified = table.identified))
        }
        else if(analysis=="spacetime"){
            varnames <- paste0("Period",1:length(unique(outclusso$init.vec$Period))) 
            if(clusterRR=="background"){
                #background rates
                clusterRRbic.qp <- as.numeric(apply(outclusso$riskratios$riskratios.qp.st$RRbic,2,
                                                 function(x) names(sort(table(x),decreasing=TRUE))[[1]]))
                clusterRRaic.qp <- as.numeric(apply(outclusso$riskratios$riskratios.qp.st$RRaic,2,
                                                 function(x) names(sort(table(x),decreasing=TRUE))[[1]]))
                clusterRRaicc.qp <- as.numeric(apply(outclusso$riskratios$riskratios.qp.st$RRaicc,2,
                                                  function(x) names(sort(table(x),decreasing=TRUE))[[1]]))
                clusterRRbic.p <- as.numeric(apply(outclusso$riskratios$riskratios.p.st$RRbic,2,
                                                    function(x) names(sort(table(x),decreasing=TRUE))[[1]]))
                clusterRRaic.p <- as.numeric(apply(outclusso$riskratios$riskratios.p.st$RRaic,2,
                                                    function(x) names(sort(table(x),decreasing=TRUE))[[1]]))
                clusterRRaicc.p <- as.numeric(apply(outclusso$riskratios$riskratios.p.st$RRaicc,2,
                                                     function(x) names(sort(table(x),decreasing=TRUE))[[1]]))
                #calc identified
                identBIC.qp <- lapply(1:ncol(outclusso$riskratios$riskratios.qp.st$RRbic),
                                      function(i) which(abs(outclusso$riskratios$riskratios.qp.st$RRbic[,i]-clusterRRbic.qp[[i]])> err)) 
                identBIC.p <- lapply(1:ncol(outclusso$riskratios$riskratios.p.st$RRbic),
                                     function(i) which(abs(outclusso$riskratios$riskratios.p.st$RRbic[,i]-clusterRRbic.p[[i]])> err)) 
                identAIC.qp <- lapply(1:ncol(outclusso$riskratios$riskratios.qp.st$RRaic),
                                      function(i) which(abs(outclusso$riskratios$riskratios.qp.st$RRaic[,i]-clusterRRaic.qp[[i]])> err)) 
                identAIC.p <- lapply(1:ncol(outclusso$riskratios$riskratios.p.st$RRaic),
                                     function(i) which(abs(outclusso$riskratios$riskratios.p.st$RRaic[,i]-clusterRRaic.p[[i]])> err)) 
                identAICc.qp <- lapply(1:ncol(outclusso$riskratios$riskratios.qp.st$RRaicc),
                                       function(i) which(abs(outclusso$riskratios$riskratios.qp.st$RRaicc[,i]-clusterRRaicc.qp[[i]])> err)) 
                identAICc.p <- lapply(1:ncol(outclusso$riskratios$riskratios.p.st$RRaicc),
                                      function(i) which(abs(outclusso$riskratios$riskratios.p.st$RRaicc[,i]-clusterRRaicc.p[[i]])> err)) 
            }
            else{
                identBIC.qp <- lapply(1:ncol(outclusso$riskratios$riskratios.qp.st$RRbic),
                                      function(i) which(outclusso$riskratios$riskratios.qp.st$RRbic[,i]>clusterRR)) 
                identBIC.p <- lapply(1:ncol(outclusso$riskratios$riskratios.p.st$RRbic),
                                     function(i) which(outclusso$riskratios$riskratios.p.st$RRbic[,i]>clusterRR)) 
                identAIC.qp <- lapply(1:ncol(outclusso$riskratios$riskratios.qp.st$RRaic),
                                      function(i) which(outclusso$riskratios$riskratios.qp.st$RRaic[,i]>clusterRR)) 
                identAIC.p <- lapply(1:ncol(outclusso$riskratios$riskratios.p.st$RRaic),
                                     function(i) which(outclusso$riskratios$riskratios.p.st$RRaic[,i]>clusterRR)) 
                identAICc.qp <- lapply(1:ncol(outclusso$riskratios$riskratios.qp.st$RRaicc),
                                       function(i) which(outclusso$riskratios$riskratios.qp.st$RRaicc[,i]>clusterRR)) 
                identAICc.p <- lapply(1:ncol(outclusso$riskratios$riskratios.p.st$RRaicc),
                                      function(i) which(outclusso$riskratios$riskratios.p.st$RRaicc[,i]>clusterRR)) 
            }
            QBIC = as.data.frame(stringi::stri_list2matrix(identBIC.qp)); colnames(QBIC) <- varnames
            BIC =  as.data.frame(stringi::stri_list2matrix(identBIC.p)); colnames(BIC) <- varnames
            QAIC = as.data.frame(stringi::stri_list2matrix(identAIC.qp)); colnames(QAIC) <- varnames
            AIC = as.data.frame(stringi::stri_list2matrix(identAIC.p)); colnames(AIC) <- varnames
            QAICc = as.data.frame(stringi::stri_list2matrix(identAICc.qp)); colnames(QAICc) <- varnames
            AICc = as.data.frame(stringi::stri_list2matrix(identAICc.p)); colnames(AICc) <- varnames
            identQP <- list(QBIC = QBIC,
                            QAIC = QAIC, 
                            QAICc = QAICc)
            identP <- list(BIC = BIC, 
                           AIC = AIC, 
                           AICc = AICc)
            table.identified <- list(identQuasiPoisson = identQP, identPoisson=identP)
            return(list(table.clusters = table.clusters, table.identified = table.identified))
        }
        else{
            #both
            varnames <- paste0("Period",1:length(unique(outclusso$init.vec.s$Period))) 
            if(clusterRR=="background"){
                #background rates
                #space
                clusterRRbic.qp.s <- as.numeric(apply(outclusso$riskratios$riskratios.qp.s$RRbic,2,
                                                 function(x) names(sort(table(x),decreasing=TRUE))[[1]]))
                clusterRRaic.qp.s <- as.numeric(apply(outclusso$riskratios$riskratios.qp.s$RRaic,2,
                                                 function(x) names(sort(table(x),decreasing=TRUE))[[1]]))
                clusterRRaicc.qp.s <- as.numeric(apply(outclusso$riskratios$riskratios.qp.s$RRaicc,2,
                                                  function(x) names(sort(table(x),decreasing=TRUE))[[1]]))
                clusterRRbic.p.s <- as.numeric(apply(outclusso$riskratios$riskratios.p.s$RRbic,2,
                                                      function(x) names(sort(table(x),decreasing=TRUE))[[1]]))
                clusterRRaic.p.s <- as.numeric(apply(outclusso$riskratios$riskratios.p.s$RRaic,2,
                                                      function(x) names(sort(table(x),decreasing=TRUE))[[1]]))
                clusterRRaicc.p.s <- as.numeric(apply(outclusso$riskratios$riskratios.p.s$RRaicc,2,
                                                       function(x) names(sort(table(x),decreasing=TRUE))[[1]]))
                #spacetime
                clusterRRbic.qp.st <- as.numeric(apply(outclusso$riskratios$riskratios.qp.st$RRbic,2,
                                                    function(x) names(sort(table(x),decreasing=TRUE))[[1]]))
                clusterRRaic.qp.st <- as.numeric(apply(outclusso$riskratios$riskratios.qp.st$RRaic,2,
                                                    function(x) names(sort(table(x),decreasing=TRUE))[[1]]))
                clusterRRaicc.qp.st <- as.numeric(apply(outclusso$riskratios$riskratios.qp.st$RRaicc,2,
                                                     function(x) names(sort(table(x),decreasing=TRUE))[[1]]))
                clusterRRbic.p.st <- as.numeric(apply(outclusso$riskratios$riskratios.p.st$RRbic,2,
                                                       function(x) names(sort(table(x),decreasing=TRUE))[[1]]))
                clusterRRaic.p.st <- as.numeric(apply(outclusso$riskratios$riskratios.p.st$RRaic,2,
                                                       function(x) names(sort(table(x),decreasing=TRUE))[[1]]))
                clusterRRaicc.p.st <- as.numeric(apply(outclusso$riskratios$riskratios.p.st$RRaicc,2,
                                                        function(x) names(sort(table(x),decreasing=TRUE))[[1]]))
                #calc identified
                #space
                identBIC.qp.s <- lapply(1:ncol(outclusso$riskratios$riskratios.qp.s$RRbic),
                                      function(i) which(abs(outclusso$riskratios$riskratios.qp.s$RRbic[,i]-clusterRRbic.qp.s[[i]])> err)) 
                identBIC.p.s <- lapply(1:ncol(outclusso$riskratios$riskratios.p.s$RRbic),
                                     function(i) which(abs(outclusso$riskratios$riskratios.p.s$RRbic[,i]-clusterRRbic.p.s[[i]])> err)) 
                identAIC.qp.s <- lapply(1:ncol(outclusso$riskratios$riskratios.qp.s$RRaic),
                                      function(i) which(abs(outclusso$riskratios$riskratios.qp.s$RRaic[,i]-clusterRRaic.qp.s[[i]])> err)) 
                identAIC.p.s <- lapply(1:ncol(outclusso$riskratios$riskratios.p.s$RRaic),
                                     function(i) which(abs(outclusso$riskratios$riskratios.p.s$RRaic[,i]-clusterRRaic.p.s[[i]])> err)) 
                identAICc.qp.s <- lapply(1:ncol(outclusso$riskratios$riskratios.qp.s$RRaicc),
                                       function(i) which(abs(outclusso$riskratios$riskratios.qp.s$RRaicc[,i]-clusterRRaicc.qp.s[[i]])> err)) 
                identAICc.p.s <- lapply(1:ncol(outclusso$riskratios$riskratios.p.s$RRaicc),
                                      function(i) which(abs(outclusso$riskratios$riskratios.p.s$RRaicc[,i]-clusterRRaicc.p.s[[i]])> err))
                #spacetime
                identBIC.qp.st <- lapply(1:ncol(outclusso$riskratios$riskratios.qp.st$RRbic),
                                      function(i) which(abs(outclusso$riskratios$riskratios.qp.st$RRbic[,i]-clusterRRbic.qp.st[[i]])> err)) 
                identBIC.p.st <- lapply(1:ncol(outclusso$riskratios$riskratios.p.st$RRbic),
                                     function(i) which(abs(outclusso$riskratios$riskratios.p.st$RRbic[,i]-clusterRRbic.p.st[[i]])> err)) 
                identAIC.qp.st <- lapply(1:ncol(outclusso$riskratios$riskratios.qp.st$RRaic),
                                      function(i) which(abs(outclusso$riskratios$riskratios.qp.st$RRaic[,i]-clusterRRaic.qp.st[[i]])> err)) 
                identAIC.p.st <- lapply(1:ncol(outclusso$riskratios$riskratios.p.st$RRaic),
                                     function(i) which(abs(outclusso$riskratios$riskratios.p.st$RRaic[,i]-clusterRRaic.p.st[[i]])> err)) 
                identAICc.qp.st <- lapply(1:ncol(outclusso$riskratios$riskratios.qp.st$RRaicc),
                                       function(i) which(abs(outclusso$riskratios$riskratios.qp.st$RRaicc[,i]-clusterRRaicc.qp.st[[i]])> err)) 
                identAICc.p.st <- lapply(1:ncol(outclusso$riskratios$riskratios.p.st$RRaicc),
                                      function(i) which(abs(outclusso$riskratios$riskratios.p.st$RRaicc[,i]-clusterRRaicc.p.st[[i]])> err)) 
                
            }
            else{
                #notbackground but set clusterRR
                #Space
                identBIC.qp.s <- lapply(1:ncol(outclusso$riskratios$riskratios.qp.s$RRbic),
                                        function(i) which(outclusso$riskratios$riskratios.qp.s$RRbic[,i]>clusterRR)) 
                identBIC.p.s <- lapply(1:ncol(outclusso$riskratios$riskratios.p.s$RRbic),
                                       function(i) which(outclusso$riskratios$riskratios.p.s$RRbic[,i]>clusterRR)) 
                identAIC.qp.s <- lapply(1:ncol(outclusso$riskratios$riskratios.qp.s$RRaic),
                                        function(i) which(outclusso$riskratios$riskratios.qp.s$RRaic[,i]>clusterRR)) 
                identAIC.p.s <- lapply(1:ncol(outclusso$riskratios$riskratios.p.s$RRaic),
                                       function(i) which(outclusso$riskratios$riskratios.p.s$RRaic[,i]>clusterRR)) 
                identAICc.qp.s <- lapply(1:ncol(outclusso$riskratios$riskratios.qp.s$RRaicc),
                                         function(i) which(outclusso$riskratios$riskratios.qp.s$RRaicc[,i]>clusterRR)) 
                identAICc.p.s <- lapply(1:ncol(outclusso$riskratios$riskratios.p.s$RRaicc),
                                        function(i) which(outclusso$riskratios$riskratios.p.s$RRaicc[,i]>clusterRR)) 
                #spacetime
                #spacetime
                identBIC.qp.st <- lapply(1:ncol(outclusso$riskratios$riskratios.qp.st$RRbic),
                                         function(i) which(outclusso$riskratios$riskratios.qp.st$RRbic[,i]>clusterRR)) 
                identBIC.p.st <- lapply(1:ncol(outclusso$riskratios$riskratios.p.st$RRbic),
                                        function(i) which(outclusso$riskratios$riskratios.p.st$RRbic[,i]>clusterRR)) 
                identAIC.qp.st <- lapply(1:ncol(outclusso$riskratios$riskratios.qp.st$RRaic),
                                         function(i) which(outclusso$riskratios$riskratios.qp.st$RRaic[,i]>clusterRR)) 
                identAIC.p.st <- lapply(1:ncol(outclusso$riskratios$riskratios.p.st$RRaic),
                                        function(i) which(outclusso$riskratios$riskratios.p.st$RRaic[,i]>clusterRR)) 
                identAICc.qp.st <- lapply(1:ncol(outclusso$riskratios$riskratios.qp.s$RRaicc),
                                          function(i) which(outclusso$riskratios$riskratios.qp.st$RRaicc[,i]>clusterRR)) 
                identAICc.p.st <- lapply(1:ncol(outclusso$riskratios$riskratios.p.st$RRaicc),
                                         function(i) which(outclusso$riskratios$riskratios.p.s$RRaicc[,i]>clusterRR))
            }
            #space
            QBIC.s = as.data.frame(stringi::stri_list2matrix(identBIC.qp.s)); colnames(QBIC.s) <- varnames
            BIC.s =  as.data.frame(stringi::stri_list2matrix(identBIC.p.s)); colnames(BIC.s) <- varnames
            QAIC.s = as.data.frame(stringi::stri_list2matrix(identAIC.qp.s)); colnames(QAIC.s) <- varnames
            AIC.s = as.data.frame(stringi::stri_list2matrix(identAIC.p.s)); colnames(AIC.s) <- varnames
            QAICc.s = as.data.frame(stringi::stri_list2matrix(identAICc.qp.s)); colnames(QAICc.s) <- varnames
            AICc.s = as.data.frame(stringi::stri_list2matrix(identAICc.p.s)); colnames(AICc.s) <- varnames
            #spacetime 
            QBIC.st = as.data.frame(stringi::stri_list2matrix(identBIC.qp.st)); colnames(QBIC.st) <- varnames
            BIC.st =  as.data.frame(stringi::stri_list2matrix(identBIC.p.st)); colnames(BIC.st) <- varnames
            QAIC.st = as.data.frame(stringi::stri_list2matrix(identAIC.qp.st)); colnames(QAIC.st) <- varnames
            AIC.st = as.data.frame(stringi::stri_list2matrix(identAIC.p.st)); colnames(AIC.st) <- varnames
            QAICc.st = as.data.frame(stringi::stri_list2matrix(identAICc.qp.st)); colnames(QAICc.st) <- varnames
            AICc.st = as.data.frame(stringi::stri_list2matrix(identAICc.p.st)); colnames(AICc.st) <- varnames
            #Output
            identQP <- list(QBIC.s = QBIC.s, QBIC.st = QBIC.st,
                            QAIC.s = QAIC.s, QAIC.st = QAIC.st, 
                            QAICc.s = QAICc.s, QAICc.st = QAICc.st)
            identP <- list(BIC.s = BIC.s, BIC.st = BIC.st,
                           AIC.s = AIC.s, AIC.st = AIC.st,
                           AICc.s = AICc.s, AICc.st = AICc.st)
            table.identified <- list(identQuasiPoisson = identQP, identPoisson=identP)
            return(list(table.clusters = table.clusters, table.identified = table.identified))
        }
    }
}

