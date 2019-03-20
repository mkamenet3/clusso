#' Create pretty table from \code{clust} output    
#' @title
#'clustpretty
#' @description 
#' This function takes the cluster detection results from the \code{clust} function and creates a nice table of detected clusters.  
#'@param clustout
#'@param analysis  A string specifying if the spatial (\code{"space")), spatio-temporal (\code{"spacetime"}), or both spatial and spatio-temporal (\code{"both"}) analysis should be executed. Default is \code{"both"}. 
#'@param clusteridentify Whether specific clusters should be identified in output; default is FALSE
#'@param clusterRR Risk ratio cut off for a cluster to be identified. Alternatively, a risk ratio can be provided (ex: 1.0, 1.25) that would serve as a cut off value for identifiying a cluster. Current default is 1.0. We define a cluster as being different from the background rate (TODO). 
#'@return Data frame of output that can be sent to a csv or list of data frames (if clusteridentify set to TRUE).
#'@examples
#'\donttest{
#'data(japanbreastcancer)
#'#Set Some Initial Conditions
#'x1=utmJapan$utmx/1000
#'y1=utmJapan$utmy/1000
#'rMax <- 20 
#'Time=5
#'japanbreastcancer <- japanbreastcancer[,-1] #get rid of indicator column
#'clst <- toclust(japanbreastcancer, expected = expdeath, observed=death,timeperiod = period)
#'system.time(resreal <- clust(clst, x,y, rMax, Time, utm=TRUE, analysis="both", maxclust=10))
#'clustpretty(resreal, analysis="both")
#'clustpretty(resreal, analysis="both", clusteridentify=TRUE)
#'clustpretty(resreal, analysis="both", clusteridentify=TRUE, clusterRR = 1.05)
#'  }


clustpretty <- function(clustout, analysis="both", clusteridentify=FALSE, clusterRR){
    if(analysis=="space"){
        model <- c("Poisson", "Quasi-Poisson")
        analysistype <- rep("Space",2)
        numclust.AIC <- c(resreal$lassoresult.p.s$numclust.qaic,resreal$lassoresult.qp.s$numclust.qaic)
        numclust.AICc <- c(resreal$lassoresult.p.s$numclust.qaicc, resreal$lassoresult.qp.s$numclust.qaicc)
        numclust.BIC <- c(resreal$lassoresult.p.s$numclust.qbic, resreal$lassoresult.qp.s$numclust.qbic)
    }
    else if(analysis=="spacetime"){
        model <- c("Poisson", "Quasi-Poisson")
        analysistype <- rep("Space-Time",2)
        numclust.AIC <- c(resreal$lassoresult.p.st$numclust.qaic, resreal$lassoresult.qp.st$numclust.qaic)
        numclust.AICc <- c(resreal$lassoresult.p.st$numclust.qaicc, resreal$lassoresult.qp.st$numclust.qaicc)
        numclust.BIC <- c(resreal$lassoresult.p.st$numclust.qbic, resreal$lassoresult.qp.st$numclust.qbic)
    }
    else{
        #both
        model <- c(rep("Poisson",2), rep("Quasi-Poisson",2))
        analysistype <- rep(c("Space", "Space-Time"),2)
        numclust.AIC <- c(resreal$lassoresult.p.s$numclust.qaic, resreal$lassoresult.p.st$numclust.qaic, 
                          resreal$lassoresult.qp.s$numclust.qaic, resreal$lassoresult.qp.st$numclust.qaic)
        numclust.AICc <- c(resreal$lassoresult.p.s$numclust.qaicc, resreal$lassoresult.p.st$numclust.qaicc, 
                           resreal$lassoresult.qp.s$numclust.qaicc, resreal$lassoresult.qp.st$numclust.qaicc)
        numclust.BIC <- c(resreal$lassoresult.p.s$numclust.qbic, resreal$lassoresult.p.st$numclust.qbic, 
                          resreal$lassoresult.qp.s$numclust.qbic, resreal$lassoresult.qp.st$numclust.qbic)
    }
    table.clusters <- cbind.data.frame(model, analysistype, numclust.AIC, numclust.AICc, numclust.BIC)
    
    if(clusteridentify==FALSE){
        return(table.clusters)
    }
    if(clusteridentify==TRUE){
        if(missing(clusterRR)){
            message("Using default risk ratio threshold for cluster of > 1.0. To specify a different risk ratio threshold, specify clusterRR. To specify different from background, specify clusterRR='background'")
            clusterRR <- 1.0
        }
        else if (clusterRR=="background"){
            clusterRR <- NULL 
            #as.numeric(apply(resreal$riskratios$riskratios.qp.s$RRbic,2,function(x) names(sort(table(x),decreasing=TRUE))[[1]]))
        }
        else{
            clusterRR <- clusterRR
        }
        #todo set default to background
        if(analysis=="space"){
            varnames <- paste0("Period",1:length(unique(resreal$init.vec.s$Period))) 
            identBIC.qp <- lapply(1:ncol(resreal$riskratios$riskratios.qp.s$RRbic),
                               function(i) which(resreal$riskratios$riskratios.qp.s$RRbic[,i]>clusterRR)) 
            identBIC.p <- lapply(1:ncol(resreal$riskratios$riskratios.p.s$RRbic),
                                  function(i) which(resreal$riskratios$riskratios.p.s$RRbic[,i]>clusterRR)) 
            identAIC.qp <- lapply(1:ncol(resreal$riskratios$riskratios.qp.s$RRaic),
                                  function(i) which(resreal$riskratios$riskratios.qp.s$RRaic[,i]>clusterRR)) 
            identAIC.p <- lapply(1:ncol(resreal$riskratios$riskratios.p.s$RRaic),
                                 function(i) which(resreal$riskratios$riskratios.p.s$RRaic[,i]>clusterRR)) 
            identAICc.qp <- lapply(1:ncol(resreal$riskratios$riskratios.qp.s$RRaicc),
                                  function(i) which(resreal$riskratios$riskratios.qp.s$RRaicc[,i]>clusterRR)) 
            identAICc.p <- lapply(1:ncol(resreal$riskratios$riskratios.p.s$RRaicc),
                                 function(i) which(resreal$riskratios$riskratios.p.s$RRaicc[,i]>clusterRR)) 
            QBIC = as.data.frame(stri_list2matrix(identBIC.qp)); colnames(QBIC) <- varnames
            BIC =  as.data.frame(stri_list2matrix(identBIC.p)); colnames(BIC) <- varnames
            QAIC = as.data.frame(stri_list2matrix(identAIC.qp)); colnames(QAIC) <- varnames
            AIC = as.data.frame(stri_list2matrix(identAIC.p)); colnames(AIC) <- varnames
            QAICc = as.data.frame(stri_list2matrix(identAICc.qp)); colnames(QAICc) <- varnames
            AICc = as.data.frame(stri_list2matrix(identAICc.p)); colnames(AICc) <- varnames
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
            varnames <- paste0("Period",1:length(unique(resreal$init.vec$Period))) 
            identBIC.qp <- lapply(1:ncol(resreal$riskratios$riskratios.qp.st$RRbic),
                                  function(i) which(resreal$riskratios$riskratios.qp.st$RRbic[,i]>clusterRR)) 
            identBIC.p <- lapply(1:ncol(resreal$riskratios$riskratios.p.st$RRbic),
                                 function(i) which(resreal$riskratios$riskratios.p.st$RRbic[,i]>clusterRR)) 
            identAIC.qp <- lapply(1:ncol(resreal$riskratios$riskratios.qp.st$RRaic),
                                  function(i) which(resreal$riskratios$riskratios.qp.st$RRaic[,i]>clusterRR)) 
            identAIC.p <- lapply(1:ncol(resreal$riskratios$riskratios.p.st$RRaic),
                                 function(i) which(resreal$riskratios$riskratios.p.st$RRaic[,i]>clusterRR)) 
            identAICc.qp <- lapply(1:ncol(resreal$riskratios$riskratios.qp.st$RRaicc),
                                   function(i) which(resreal$riskratios$riskratios.qp.st$RRaicc[,i]>clusterRR)) 
            identAICc.p <- lapply(1:ncol(resreal$riskratios$riskratios.p.st$RRaicc),
                                  function(i) which(resreal$riskratios$riskratios.p.st$RRaicc[,i]>clusterRR)) 
            QBIC = as.data.frame(stri_list2matrix(identBIC.qp)); colnames(QBIC) <- varnames
            BIC =  as.data.frame(stri_list2matrix(identBIC.p)); colnames(BIC) <- varnames
            QAIC = as.data.frame(stri_list2matrix(identAIC.qp)); colnames(QAIC) <- varnames
            AIC = as.data.frame(stri_list2matrix(identAIC.p)); colnames(AIC) <- varnames
            QAICc = as.data.frame(stri_list2matrix(identAICc.qp)); colnames(QAICc) <- varnames
            AICc = as.data.frame(stri_list2matrix(identAICc.p)); colnames(AICc) <- varnames
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
            varnames <- paste0("Period",1:length(unique(resreal$init.vec.s$Period))) 
            #Space
            identBIC.qp.s <- lapply(1:ncol(resreal$riskratios$riskratios.qp.s$RRbic),
                                  function(i) which(resreal$riskratios$riskratios.qp.s$RRbic[,i]>clusterRR)) 
            identBIC.p.s <- lapply(1:ncol(resreal$riskratios$riskratios.p.s$RRbic),
                                 function(i) which(resreal$riskratios$riskratios.p.s$RRbic[,i]>clusterRR)) 
            identAIC.qp.s <- lapply(1:ncol(resreal$riskratios$riskratios.qp.s$RRaic),
                                  function(i) which(resreal$riskratios$riskratios.qp.s$RRaic[,i]>clusterRR)) 
            identAIC.p.s <- lapply(1:ncol(resreal$riskratios$riskratios.p.s$RRaic),
                                 function(i) which(resreal$riskratios$riskratios.p.s$RRaic[,i]>clusterRR)) 
            identAICc.qp.s <- lapply(1:ncol(resreal$riskratios$riskratios.qp.s$RRaicc),
                                   function(i) which(resreal$riskratios$riskratios.qp.s$RRaicc[,i]>clusterRR)) 
            identAICc.p.s <- lapply(1:ncol(resreal$riskratios$riskratios.p.s$RRaicc),
                                  function(i) which(resreal$riskratios$riskratios.p.s$RRaicc[,i]>clusterRR)) 
            QBIC.s = as.data.frame(stri_list2matrix(identBIC.qp.s)); colnames(QBIC.s) <- varnames
            BIC.s =  as.data.frame(stri_list2matrix(identBIC.p.s)); colnames(BIC.s) <- varnames
            QAIC.s = as.data.frame(stri_list2matrix(identAIC.qp.s)); colnames(QAIC.s) <- varnames
            AIC.s = as.data.frame(stri_list2matrix(identAIC.p.s)); colnames(AIC.s) <- varnames
            QAICc.s = as.data.frame(stri_list2matrix(identAICc.qp.s)); colnames(QAICc.s) <- varnames
            AICc.s = as.data.frame(stri_list2matrix(identAICc.p.s)); colnames(AICc.s) <- varnames
            #spacetime
            identBIC.qp.st <- lapply(1:ncol(resreal$riskratios$riskratios.qp.st$RRbic),
                                    function(i) which(resreal$riskratios$riskratios.qp.st$RRbic[,i]>clusterRR)) 
            identBIC.p.st <- lapply(1:ncol(resreal$riskratios$riskratios.p.st$RRbic),
                                   function(i) which(resreal$riskratios$riskratios.p.st$RRbic[,i]>clusterRR)) 
            identAIC.qp.st <- lapply(1:ncol(resreal$riskratios$riskratios.qp.st$RRaic),
                                    function(i) which(resreal$riskratios$riskratios.qp.st$RRaic[,i]>clusterRR)) 
            identAIC.p.st <- lapply(1:ncol(resreal$riskratios$riskratios.p.st$RRaic),
                                   function(i) which(resreal$riskratios$riskratios.p.st$RRaic[,i]>clusterRR)) 
            identAICc.qp.st <- lapply(1:ncol(resreal$riskratios$riskratios.qp.s$RRaicc),
                                     function(i) which(resreal$riskratios$riskratios.qp.st$RRaicc[,i]>clusterRR)) 
            identAICc.p.st <- lapply(1:ncol(resreal$riskratios$riskratios.p.st$RRaicc),
                                    function(i) which(resreal$riskratios$riskratios.p.s$RRaicc[,i]>clusterRR)) 
            QBIC.st = as.data.frame(stri_list2matrix(identBIC.qp.st)); colnames(QBIC.st) <- varnames
            BIC.st =  as.data.frame(stri_list2matrix(identBIC.p.st)); colnames(BIC.st) <- varnames
            QAIC.st = as.data.frame(stri_list2matrix(identAIC.qp.st)); colnames(QAIC.st) <- varnames
            AIC.st = as.data.frame(stri_list2matrix(identAIC.p.st)); colnames(AIC.st) <- varnames
            QAICc.st = as.data.frame(stri_list2matrix(identAICc.qp.st)); colnames(QAICc.st) <- varnames
            AICc.st = as.data.frame(stri_list2matrix(identAICc.p.st)); colnames(AICc.st) <- varnames
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

