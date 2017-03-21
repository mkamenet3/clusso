#' cluST.R
#' @name cluST
#' @title Spatial and Spatial-Temporal Clustering using Lasso
#' @description This package will host the functions and code for a spatial and spatial-temporal clustering model using Lasso regression. We assume a Poisson distribution for counts. The model includes a spatially-structured random effect to account for heterogeneity in a Poisson model. Model selection is based on QAIC/QAICc/QBIC.
#' @details This package contains C++ code which is compiled at run-time. These functions are located in the "src" file.
#' Dependencies and inputs are documented here.
#' Data used for this analysis comes from MidWest Poverty data.
#' @references Xu, Jiale, Gangnon, Ronald: "Stagewise and Stepwise Methods for Space and Space-Time Cluster Detection"
#' @references Kamenetsky, M., Gangnon, R., Zhu, J., Lee, J. "Space and Space-Time Cluster Detection Using the Lasso"
#' 
#' 
#' 
#' 
#' Unlist2
#' 
#' This function is from https://github.com/Bioconductor-mirror/AnnotationDbi/blob/master/R/unlist2.R
#' It unlists the names without issues
#' 
make.name.tree <- function(x, recursive, what.names)
{
    if (!is.character(what.names) || length(what.names) != 1)
        stop("'what.names' must be a single string")
    what.names <- match.arg(what.names, c("inherited" , "full"))
    .make.name.tree.rec <- function(x, parent_name, depth)
    {
        if (length(x) == 0)
            return(character(0))
        x_names <- names(x)
        if (is.null(x_names))
            x_names <- rep.int(parent_name, length(x))
        else if (what.names == "full")
            x_names <- paste0(parent_name, x_names)
        else
            x_names[x_names == ""] <- parent_name
        if (!is.list(x) || (!recursive && depth >= 1L))
            return(x_names)
        if (what.names == "full")
            x_names <- paste0(x_names, ".")
        lapply(seq_len(length(x)),
               function(i) .make.name.tree.rec(x[[i]], x_names[i], depth + 1L))
    }
    .make.name.tree.rec(x, "", 0L)
}

unlist2 <- function(x, recursive=TRUE, use.names=TRUE, what.names="inherited")
{
    ans <- unlist(x, recursive, FALSE)
    if (!use.names)
        return(ans)
    if (!is.character(what.names) || length(what.names) != 1)
        stop("'what.names' must be a single string")
    what.names <- match.arg(what.names, c("inherited" , "full"))
    names(ans) <- unlist(make.name.tree(x, recursive, what.names), recursive, FALSE)
    ans
}

#' Create the clusters dataframe
#' 
#' @param vector xP x coordinates (easting/latitude); if utm coordinates, scale to km.
#' @param vector yP y coordinates (northing/longitude); if utm coordinates, scale to km.
#' @param r.max set max radius (in km)
#' @param utm TRUE/FALSE as to whether or not the x and y coordinates are in UTM (TRUE) or LAT/LONG(FALSE)
#' @param n Number of coordinate pairs/number of centers
#' @return This function returns a dataframe that contains 
#' @export
#' @examples
#' cluster.df(x1,y1,rMax, utm=TRUE, length(x1))
#' cluster.df(lat, long, utm=FALSE, length(lat))

clusters.df <- function(xP,yP, r.max, utm=FALSE,n){
    indR = (1:n)[!duplicated(cbind(xP,yP))] 
    if(utm==FALSE){
        tmpR <- (as.matrix(distm(cbind(xP, yP), fun=distHaversine))[indR,])/1000    
    } 
    else{
        tmpR = as.matrix(dist(cbind(xP,yP)))[indR,]
    }
    lastR = apply(tmpR, 1, function(x,r) order(x)[1:sum(x<=r)],r=r.max)
    ncR = unlist2(lapply(lastR, length))
    lastR = unlist2(lastR)
    rR=unlist2(apply(tmpR,1, function(x,r) { sort(x[x<=r]) },r=r.max))
    
    clustersR=data.frame(center=rep(indR,ncR),
                         x=xP[rep(indR,ncR)],y=yP[rep(indR,ncR)],
                         r=rR, 
                         n=unlist(lapply(ncR,seq)),
                         last=lastR)    
    return(clustersR)
}

#' Poisson Distribution Function
#' This is the main distriution function for our model. This assumes we have a Poisson fixed effect and Gamma random effect. In order to deal with constraints from the Lasso function, we use the Poisson distirbution function here and account for overdispersion in the QIC.
#' @param y observed values
#' @param lambda vector of expected outcomes * exp(each column of each potential path)
#' @param log whether or not the log-likelihood should be returned or the likelihood. Default is to be TRUE
#' @return returns a matrix 
#' @export

dpoisson <- function(y, lambda, log = FALSE) {
    if(log == FALSE) 
        return(lambda^y * exp(-lambda)/factorial(y))
    else
        return(y*ifelse(lambda==0,1,log(lambda))-lambda)
}


#' Creates a List Arranged by Time Period with Expected and Observed Counts and Time Period
#' 
#' @param period vector of periods or years in dataset. Should be imported as a factor.
#' @param expect vector of expected counts. Expected counts must match up with the year and observed vectors.
#' @param observed vector of observed counts. Observed counts must match up with the year and expected vectors.
#' @param Time Number of time periods or years in your dataset. Must be declared as numeric.
#' @param byrow default is set to TRUE. Data from the dataset should be imported by row. This is most often the case
#' when you have a dataframe ordered by an identifier and then the period/time frame within that id listed chronologically (in panel format by identifier).
#' If you are simulating data and have each observed/expected vector separate and create the period vector with repetitions of each time
#' period by group, this should be set to false.
#' @return This function returns a list of expected and observed counts along with the period. 
#' @export
#' @examples
set.vectors <- function(period, expect, observed,Time, byrow=TRUE,...) {
    if (byrow==TRUE){
        E0=as.vector(matrix(expect, byrow=T, ncol=Time))
        Y.vec <- as.vector(matrix(observed,byrow=T, ncol=Time))
        Year <- as.vector(matrix(period, byrow=T, ncol=Time)) 
    }
    else {
        E0=as.vector(matrix(expect, ncol=Time))
        Y.vec <- as.vector(matrix(observed, ncol=Time))
        Year <- as.vector(matrix(period, ncol=Time))
    }
    return(list(
        E0 = E0,
        Y.vec = Y.vec,
        Year = Year))
}


#' space.mat
#' 
#' This function creates a sparse matrix of 1's of all potential clusters for the Lasso algorithm to cycle over; this incorporates space
#' @param clusters clusters dataframe from (cluster.df function) that includes the center, x,y, r (radius), n (counter), and last (last observation in potential cluster)
#' @param numCenters the number of centers
#' @return returns sparse matrix of 1's
#' @export
space.mat <- function(clusters, numCenters,...){
    potClus <- numCenters
    mymat <- NULL
    for(i in 1:nrow(clusters)){
        myvec <- list(as(max_colCpp(numCenters, i, clusters$n, clusters$last), "sparseVector")) 
        mymat <- c(mymat,myvec) 
    }
    xx <- NULL
    jj <- NULL
    ii <- NULL
    for(k in 1:length(mymat)){
        xx<- c(xx, mymat[[k]]@x)
        jj <- c(jj, mymat[[k]]@i)
        ii <- c(ii, rep(k,length(mymat[[k]]@x)))
    }
    return(t(sparseMatrix(i = ii, j = jj, x =xx, dims = c(length(mymat), numCenters))))
}

#' time.mat
#' 
#' This function creates a sparse matrix of 1's of all of the potential time periods for the cluster to be in. The number of 
#' potential time periods is determined by [(Time*(Time-1)]/2.
#' @param Time Number of time periods in the data
#' @return returns sparse matrix of 1's as indicators of membership in the time cluster
#' @export
time.mat <-function(Time,..){
    block <- Matrix(diag(1,Time),sparse=TRUE)
    master <- block
    for(i in 1:(Time-2)){
        diag(block[(i+1):Time,])<-1
        master <- cBind(master, block[,1:(Time-i)])        
    }
    master <- cBind(master, Matrix(rep(1,Time)))
    return(master)
}
    
    
#' spacetime.mat
#' 
#' This function takes the Kronecker product of the space and time matrices to create the space-time matrix
#' @param time.mat Time matrix
#' @param numCenters number of centroids in space matrix
#' @param Time number of time periods
#' @return Returns sparse space time matrix     
#' @export     
spacetime.mat <- function(clusters, numCenters, Time,...){
    space <- space.mat(clusters, numCenters)
    time <- time.mat(Time)
    spacetime.mat <- kronecker(time, space)
    return(spacetime.mat)
}



#' spacetime.mat.group
#' 
#' This function creates a space-time matrix for the group Lasso procedure, where group is the number of time periods. 
#' This function is used in conjunection with "spacetime.lasso.group".
#' @param time.mat Time matrix
#' @param numCenters number of centroids in space matrix
#' @param Time number of time periods
#' @return Returns sparse space time matrix for group Lasso.         
spacetime.mat.group <- function(clusters, numCenters, Time,...){
    time <- Matrix(diag(1,Time),sparse=TRUE)
    space <- space.mat(clusters, numCenters)
    spacetimematgroup <- kronecker(time, space)
    return(spacetimematgroup)
}

#' get.groups
#' 
#' This function creates a vector of the groups used for the group Lasso.
#' This function is used in conjunection with "spacetime.lasso.group".
#' @param space space-only vector from the spacetime.matGroup function
#' @param Time number of time periods
#' @return Returns sparse space time matrix for group Lasso.         
get.groups <- function(space, Time,...){
    vec <- NULL
    for (i in 1:Time){
        vec <-c(vec,rep(i,(ncol(space)/Time)))
    }
    group <- as.vector(c(vec,Time+1))
    return(group)
}





#' overdisp
#' 
#' This function calculates the overdispersion parameter for the QIC 'c' overdispersion parameter.
#' @param 
#' @return returns sparse matrix of 1's
overdisp <- function(object) {
    with(object,sum((weights * residuals^2)[weights > 0])/df.residual)
}

#' spacetime.lasso
#' 
#' This function runs the Lasso regularization technique on our large sparse matric of potential space-time clusters.
#' @param potClus number of potential clusters. This will usually be the same as 'numCenters'
#' @param clusters clusters dataframe from (cluster.df function) that includes the center, x,y, r (radius), n (counter), and last (last observation in potential cluster)
#' @param numCenters the number of centers
#' @param vectors takes in the list of expected and observed counts from set.vectors function
#' @param Time number of time periods in the dataset
#' @param spacetime indicator of whether the cluster detection method should be run on all space-time clusters(default) or on only the potential space clusters.
#' @return This function will return a list with the expected counts as selected by QBIC, QAIC, QAICc, a list of original expected counts (Ex),
#' a list of observed counts (Yx), the lasso object, a list of K values (number of unique values in each decision path), and n (length of unique centers in the clusters dataframe)
#' @export
spacetime.lasso<- function(clusters, vectors, Time, spacetime=TRUE,pois=FALSE,...){
    n <- length(unique(clusters$center))
    potClus <- n
    numCenters <- n
    message("Creating space-time matrix")
    if(spacetime==TRUE){
        sparseMAT <- spacetime.mat(clusters, numCenters, Time)
    }
    else{
        sparseMAT <- space.mat(clusters, numCenters)
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




#' spacetime.lasso.sim
#' 
#' This function runs the Lasso regularization technique on our large sparse matric of potential space-time clusters.It is specifically created to use with simulations.
#' @param potClus number of potential clusters. This will usually be the same as 'numCenters'
#' @param clusters clusters dataframe from (cluster.df function) that includes the center, x,y, r (radius), n (counter), and last (last observation in potential cluster)
#' @param numCenters the number of centers
#' @param vectors takes in the list of expected and observed counts from set.vectors function
#' @param Time number of time periods in the dataset
#' @param spacetime indicator of whether the cluster detection method should be run on all space-time clusters(default) or on only the potential space clusters.
#' @return This function will return a list with the expected counts as selected by QBIC, QAIC, QAICc, a list of original expected counts (Ex),
#' a list of observed counts (Yx), the lasso object, a list of K values (number of unique values in each decision path), and n (length of unique centers in the clusters dataframe)
#' @param nsim number of simulations
#' @param YSIM vector of simulated observed counts
#' @param getRR default is TRUE. Determines how relative risk is calculated. If getRR=TRUE, then relative risk is first calculated in each simulation and averaged over simulations.
#' If getRR is FALSE, then relative risk will come from averaging the expected counts from each of the Lasso results by the number of simulations. This should then
#' by used in conjunction with the set.rrfunction to get the relative risk. TODO integrate this into one flow.
#' @export
spacetime.lasso.sim <- function(clusters, vectors.sim, Time, spacetime=TRUE,pois=FALSE, nsim,YSIM,...){
    n <- length(unique(clusters$center))
    potClus <- n
    numCenters <- n
    message("Creating space-time matrix")
    if(spacetime==TRUE){
        sparseMAT <- spacetime.mat(clusters, numCenters, Time)
    }
    else{
        sparseMAT <- space.mat(clusters, numCenters)
    }
    message("Space-time matrix created")
    Ex <- vectors.sim$Ex
    Yx <- YSIM
    Period <- vectors.sim$Period
    message("Running Lasso - stay tuned")
    lasso <- lapply(1:nsim, function(i) glmnet(sparseMAT, Yx[[i]], family=("poisson"), alpha=1, offset=log(Ex[[i]]), 
                                               nlambda = 2000, standardize = FALSE, dfmax = 10))
    message("Lasso complete - extracting estimates and paths")
    coefs.lasso.all <- lapply(1:nsim, function(i) coef(lasso[[i]]))
    intercept <- rep(1, dim(sparseMAT)[1])
    sparseMAT <- cBind(intercept, sparseMAT)
    xbetaPath<- lapply(1:nsim, function(i) sparseMAT%*%coefs.lasso.all[[i]])
    mu <- lapply(1:nsim, function(j) sapply(1:length(lasso[[j]]$lambda), 
                                            function(i) exp(xbetaPath[[j]][,i])))    
    loglike <- lapply(1:nsim, function(k) sapply(1:length(lasso[[k]]$lambda), 
                                                 function(i) sum(dpoisson(Yx[[k]], mu[[k]][,i],log=TRUE))))
    #K <- lapply(1:nsim, function(j) sapply(1:length(lasso[[j]]$lambda), function(i) length(unique(xbetaPath[[j]][,i]))))
    
    K <- lapply(1:nsim, function(i) lasso[[i]]$df + 1)
    message("Selecting best paths")
    if(spacetime==TRUE & pois == FALSE){
        message("returning results for space-time Quasi-Poisson model")
        offset_reg <- lapply(1:nsim, function(i) glm(Yx[[i]] ~ 1 + as.factor(vectors.sim$Period) +offset(log(Ex[[i]])),family=poisson))
        overdisp.est <- lapply(1:nsim, function(i) overdisp(offset_reg[[i]]))
        
        #QBIC
        PLL.qbic  <- lapply(1:nsim, function(i) (loglike[[i]]/overdisp.est[[i]])-log(n*Time)/2*K[[i]])
        select.qbic <- lapply(1:nsim, function(i) which.max(unlist(PLL.qbic[[i]])))
        select_mu.qbic <- lapply(1:nsim, function(i) sapply(select.qbic[[i]], function(j) mu[[i]][,j]))
        select_muRR.qbic <- Reduce("+", select_mu.qbic)/nsim
        E.qbic <- select_muRR.qbic
            
        #QAIC
        PLL.qaic = lapply(1:nsim, function(i) (loglike[[i]]/overdisp.est[[i]]) - K[[i]])
        select.qaic <- lapply(1:nsim, function(i) which.max(unlist(PLL.qaic[[i]])))
            select_mu.qaic <- lapply(1:nsim, function(i) sapply(select.qaic[[i]], function(j) mu[[i]][,j]))
            select_muRR.qaic <- Reduce("+", select_mu.qaic)/nsim
            E.qaic <- select_muRR.qaic
            
        
        #QAICc
        PLL.qaicc=lapply(1:nsim, function(i) (loglike[[i]]/overdisp.est[[i]])- ((K[[i]]*n*Time)/(n*Time-K[[i]]-1)))
        select.qaicc <- lapply(1:nsim, function(i) which.max(unlist(PLL.qaicc[[i]])))
        
            select_mu.qaicc <- lapply(1:nsim, function(i) sapply(select.qaicc[[i]], function(j) mu[[i]][,j]))
            select_muRR.qaicc <- Reduce("+", select_mu.qaicc)/nsim
            E.qaicc <- select_muRR.qaicc
        
        
    }
    
    #########################################################
    #Space-Time, Poisson only
    #########################################################
    else if(spacetime==TRUE & pois == TRUE){
        message("returning results for space-time Poisson model")
        #QBIC
        PLL.qbic  <- lapply(1:nsim, function(i) (loglike[[i]]-log(n*Time)/2*K[[i]]))
        select.qbic <- lapply(1:nsim, function(i) which.max(unlist(PLL.qbic[[i]])))
        select_mu.qbic <- lapply(1:nsim, function(i) sapply(select.qbic[[i]], function(j) mu[[i]][,j]))
        select_muRR.qbic <- Reduce("+", select_mu.qbic)/nsim
        E.qbic <- select_muRR.qbic
        
        #QAIC
        PLL.qaic = lapply(1:nsim, function(i) (loglike[[i]] - K[[i]]))
        select.qaic <- lapply(1:nsim, function(i) which.max(unlist(PLL.qaic[[i]])))
        select_mu.qaic <- lapply(1:nsim, function(i) sapply(select.qaic[[i]], function(j) mu[[i]][,j]))
        select_muRR.qaic <- Reduce("+", select_mu.qaic)/nsim
        E.qaic <- select_muRR.qaic
        
        
        #QAICc
        PLL.qaicc=lapply(1:nsim, function(i) (loglike[[i]] - ((K[[i]]*n*Time)/(n*Time-K[[i]]-1))))
        select.qaicc <- lapply(1:nsim, function(i) which.max(unlist(PLL.qaicc[[i]])))
        
        select_mu.qaicc <- lapply(1:nsim, function(i) sapply(select.qaicc[[i]], function(j) mu[[i]][,j]))
        select_muRR.qaicc <- Reduce("+", select_mu.qaicc)/nsim
        E.qaicc <- select_muRR.qaicc
    }
        
    
    #########################################################
    #Space-Only, Quasi-Poisson
    #########################################################
    else if(spacetime==FALSE & pois == FALSE){
        message("Returning results for space-only  Quasi-Poisson model")
        offset_reg <- lapply(1:nsim, function(i) glm(Yx[[i]] ~ 1  +offset(log(Ex[[i]])),family=poisson))
        overdisp.est <- lapply(1:nsim, function(i) overdisp(offset_reg[[i]]))
        
        #QBIC
        PLL.qbic  <- lapply(1:nsim, function(i) (loglike[[i]]/overdisp.est[[i]])-log(n*Time)/2*K[[i]])
        select.qbic <- lapply(1:nsim, function(i) which.max(unlist(PLL.qbic[[i]])))
        select_mu.qbic <- lapply(1:nsim, function(i) sapply(select.qbic[[i]], function(j) mu[[i]][,j]))
        select_muRR.qbic <- Reduce("+", select_mu.qbic)/nsim
        E.qbic <- select_muRR.qbic
        
        #QAIC
        PLL.qaic = lapply(1:nsim, function(i) (loglike[[i]]/overdisp.est[[i]]) - K[[i]])
        select.qaic <- lapply(1:nsim, function(i) which.max(unlist(PLL.qaic[[i]])))
        select_mu.qaic <- lapply(1:nsim, function(i) sapply(select.qaic[[i]], function(j) mu[[i]][,j]))
        select_muRR.qaic <- Reduce("+", select_mu.qaic)/nsim
        E.qaic <- select_muRR.qaic
        
        
        #QAICc
        PLL.qaicc=lapply(1:nsim, function(i) (loglike[[i]]/overdisp.est[[i]])- ((K[[i]]*n*Time)/(n*Time-K[[i]]-1)))
        select.qaicc <- lapply(1:nsim, function(i) which.max(unlist(PLL.qaicc[[i]])))
        select_mu.qaicc <- lapply(1:nsim, function(i) sapply(select.qaicc[[i]], function(j) mu[[i]][,j]))
        select_muRR.qaicc <- Reduce("+", select_mu.qaicc)/nsim
        E.qaicc <- select_muRR.qaicc
    }
    
    
    
    #########################################################
    #Space-only, Poisson only
    #########################################################
    else if(spacetime==FALSE & pois == TRUE){
        message("Returning results for space-only  Poisson model")
        #QBIC
        PLL.qbic  <- lapply(1:nsim, function(i) (loglike[[i]]-log(n*Time)/2*K[[i]]))
        select.qbic <- lapply(1:nsim, function(i) which.max(unlist(PLL.qbic[[i]])))
        select_mu.qbic <- lapply(1:nsim, function(i) sapply(select.qbic[[i]], function(j) mu[[i]][,j]))
        select_muRR.qbic <- Reduce("+", select_mu.qbic)/nsim
        E.qbic <- select_muRR.qbic
        
        #QAIC
        PLL.qaic = lapply(1:nsim, function(i) (loglike[[i]] - K[[i]]))
        select.qaic <- lapply(1:nsim, function(i) which.max(unlist(PLL.qaic[[i]])))
        select_mu.qaic <- lapply(1:nsim, function(i) sapply(select.qaic[[i]], function(j) mu[[i]][,j]))
        select_muRR.qaic <- Reduce("+", select_mu.qaic)/nsim
        E.qaic <- select_muRR.qaic
        
        
        #QAICc
        PLL.qaicc=lapply(1:nsim, function(i) (loglike[[i]] - ((K[[i]]*n*Time)/(n*Time-K[[i]]-1))))
        select.qaicc <- lapply(1:nsim, function(i) which.max(unlist(PLL.qaicc[[i]])))
        
        select_mu.qaicc <- lapply(1:nsim, function(i) sapply(select.qaicc[[i]], function(j) mu[[i]][,j]))
        select_muRR.qaicc <- Reduce("+", select_mu.qaicc)/nsim
        E.qaicc <- select_muRR.qaicc
    }
    
    return(list(nsim = nsim, E.qbic = E.qbic, E.qaic = E.qaic, E.qaicc = E.qaicc,Ex = Ex,mu = mu, Yx = Yx, PLL.qbic = PLL.qbic, 
                PLL.qaic = PLL.qaic, PLL.qaicc = PLL.qaicc, select.qbic = select.qbic, select.qaic = select.qaic, 
                select.qaicc = select.qaicc, select_mu.qbic = select_mu.qbic, select_mu.qaic = select_mu.qaic, 
                select_mu.qaicc = select_mu.qaicc, xbetaPath = xbetaPath, coefs.lasso.all = coefs.lasso.all))    
}


#' spacetime.lasso.group
#' 
#' This function runs the *group* Lasso regularization technique on our large sparse matric of potential space-time clusters. Groups are each of the
#' time periods in the data. For example, if there are 5 time periods, then there are 5 groups.
#' @param potClus number of potential clusters. This will usually be the same as 'numCenters'
#' @param clusters clusters dataframe from (cluster.df function) that includes the center, x,y, r (radius), n (counter), and last (last observation in potential cluster)
#' @param numCenters the number of centers
#' @param vectors takes in the list of expected and observed counts from set.vectors function
#' @param Time number of time periods in the dataset
#' @param spacetime indicator of whether the cluster detection method should be run on all space-time clusters(default) or on only the potential space clusters.
#' @param group option which allows user to specify the model should be run using group Lasso with time periods as the group.
#' @return This function will return a list with the expected counts as selected by QBIC, QAIC, QAICc, a list of original expected counts (Ex),
#' a list of observed counts (Yx), the lasso object, a list of K values (number of unique values in each decision path), and n (length of unique centers in the clusters dataframe)
#' @export
spacetime.lasso.group <- function(potClus, clusters, numCenters, vectors, Time, ...){
    n <- length(unique(clusters$center))
    potClus <- n
    numCenters <- n
    message("Creating space-time matrix grouped on time")
    sparseMAT <- spacetime.matGroup(clusters, numCenters, Time)
    groups <- get.groups(sparseMAT, Time)
    message("Space-time matrix created")
    Ex <- vectors$E0
    Yx <- vectors$Y.vec
    #Add in offset as unpenalized group
    offset <- log(Ex)
    sparseMAT <- cBind(sparseMAT, offset)
    message("Running group Lasso - stay tuned")
    lasso <- grpreg(sparseMAT, Yx, groups, penalty = "grLasso", family=("poisson"))
    message("Group Lasso complete - extracting estimates and paths")
    coefs.lasso.all <- coef(lasso)
    intercept <- rep(1, dim(sparseMAT)[1])
    sparseMAT <- cBind(intercept, sparseMAT)
    xbetaPath<- sparseMAT%*%coefs.lasso.all
    mu <- sapply(1:length(lasso$lambda), function(i) Ex * exp(xbetaPath[,i]))    
    loglike <- sapply(1:length(lasso$lambda), function(i) sum(dpoisson(Yx, mu[,i],log=TRUE)))
    K <- sapply(1:length(lasso$lambda), function(i) length(unique(xbetaPath[,i])))
    
    offset_reg <- glm(Yx ~ offset(log(Ex)),family=poisson)
    overdisp <- overdisp(offset_reg)
    message("Selecting best paths")
    
    #QBIC
    PLL.qbic  <- (loglike/overdisp)-log(n*Time)/2*K
    #PLL.qbic  <- (loglike)-log(n*Time)/2*K
    qbicMax <- which.max(PLL.qbic)
    E.qbic <- mu[,qbicMax]
    
    #QAIC
    #PLL.qaic = (loglike/overdisp) - K
    PLL.qaic = (loglike) - K
    qaicMax <- which.max(PLL.qaic)
    E.qaic <- mu[,qaicMax]
    
    #QAICc
    #PLL.qaicc=(loglike/overdisp)- ((K*n*Time)/(n*Time-K-1))
    PLL.qaicc=(loglike)- ((K*n*Time)/(n*Time-K-1))
    qaiccMax <- which.max(PLL.qaicc)
    E.qaicc <- mu[,qaiccMax]
    
    message("Returning results")
    return(list(E.qbic = E.qbic, E.qaic = E.qaic, E.qaicc = E.qaicc, Ex = Ex, Yx = Yx, lasso = lasso, K = K, n = n))    
}

#' spacetime.lasso.group.sim
#' 
#' This function runs the Lasso regularization technique on our large sparse matric of potential space-time clusters.It is specifically created to use with simulations.
#' @param potClus number of potential clusters. This will usually be the same as 'numCenters'
#' @param clusters clusters dataframe from (cluster.df function) that includes the center, x,y, r (radius), n (counter), and last (last observation in potential cluster)
#' @param numCenters the number of centers
#' @param vectors takes in the list of expected and observed counts from set.vectors function
#' @param Time number of time periods in the dataset
#' @param spacetime indicator of whether the cluster detection method should be run on all space-time clusters(default) or on only the potential space clusters.
#' @return This function will return a list with the expected counts as selected by QBIC, QAIC, QAICc, a list of original expected counts (Ex),
#' a list of observed counts (Yx), the lasso object, a list of K values (number of unique values in each decision path), and n (length of unique centers in the clusters dataframe)
#' @param nsim number of simulations
#' @param YSIM vector of simulated observed counts
#' @param getRR default is TRUE. Determines how relative risk is calculated. If getRR=TRUE, then relative risk is first calculated in each simulation and averaged over simulations.
#' If getRR is FALSE, then relative risk will come from averaging the expected counts from each of the Lasso results by the number of simulations. This should then
#' by used in conjunction with the set.rr function to get the relative risk. TODO integrate this into one flow.
#' @export
spacetime.lasso.group.sim <- function(potClus, clusters, numCenters, vectors, Time, spacetime=TRUE,nsim,YSIM, getRR=TRUE,...){
    n <- length(unique(clusters$center))
    potClus <- n
    numCenters <- n
    message("Creating space-time matrix")
    sparseMAT <- spacetime.matGroup(clusters, numCenters, Time)
    groups <- get.groups(sparseMAT, Time)
    message("Space-time matrix created")
    Ex <- vectors[[2]]
    Yx <- YSIM
    Period <- vectors[[1]]
    
    #add in offset as un-penalized group 0 in design matrix
    offset <- log(Ex[[1]])
    sparseMAT <- cBind(sparseMAT,offset)
    
    message("Running Group Lasso - stay tuned")
    lasso <- lapply(1:nsim, function(i) grpreg(sparseMAT, Yx[,i], penalty="grLasso",family=("poisson"), dfmax=1))
    
    
    message("Group Lasso complete - extracting estimates and paths")
    #coefs.lasso.all <- lapply(1:nsim, function(i) coef(lasso[[i]]))
    #lob off offset estimate
    coefs.lasso.all <- lapply(1:nsim, function(i) lasso[[i]]$beta[-c(44802),])
    
    
    intercept <- rep(1, dim(sparseMAT)[1])
    sparseMAT <- cBind(intercept, sparseMAT)
    #lob off offset columns
    sparseMAT <- sparseMAT[,-c(44802)]
    
    xbetaPath<- lapply(1:nsim, function(i) sparseMAT%*%coefs.lasso.all[[i]])
    mu <- lapply(1:nsim, function(j) sapply(1:length(lasso[[j]]$lambda), 
                                            function(i) exp(xbetaPath[[j]][,i])))    
    if(getRR==FALSE){
        mu <- lapply(1:nsim, function(j) sapply(1:length(lasso[[j]]$lambda), 
                                                function(i) Ex[[i]]*exp(xbetaPath[[j]][,i])))    
    }
    loglike <- lapply(1:nsim, function(k) sapply(1:length(lasso[[k]]$lambda), 
                                                 function(i) sum(dpoisson(Yx[,k], mu[[k]][,i],log=TRUE))))
    K <- lapply(1:nsim, function(j) sapply(1:length(lasso[[j]]$lambda), function(i) length(unique(xbetaPath[[j]][,i]))))
    offset_reg <- lapply(1:nsim, function(i) glm(Yx[,i] ~ 1 + as.factor(vectors$Period) +offset(log(Ex[[i]])),family=poisson))
    overdisp <- lapply(1:nsim, function(i) overdisp(offset_reg[[i]]))
    message("Selecting best paths")
    if(spacetime==TRUE){
        #QBIC
        PLL.qbic  <- lapply(1:nsim, function(i) (loglike[[i]]/overdisp[[i]])-log(n*Time)/2*K[[i]])
        select.qbic <- lapply(1:nsim, function(i) which.max(unlist(PLL.qbic[[i]])))
        if(getRR==FALSE){
            select_mu.qbic <- lapply(1:nsim, function(i) sapply(select.qbic[[i]], function(j) mu[[i]][,j]))    
            E.qbic <- Reduce("+", select_mu.qbic)/nsim
        }
        else{
            select_mu.qbic <- lapply(1:nsim, function(i) sapply(select.qbic[[i]], function(j) mu[[i]][,j]*Ex[[i]]))    
            
            
            #test
            #test <- scale(Yx[,1],select_mu.qbic[[1]],1,5)
            test <- sapply(1:Time, function(j)
                ((matrix(select_mu.qbic[[1]],ncol=5)[,j])*(sum(matrix(Ex[[1]],ncol=5)[,j])))/(sum(matrix(select_mu.qbic[[1]],ncol=5))))
            
            
            
            select_muRR.qbic <- lapply(1:nsim, function(i) select_mu.qbic[[i]]/vectors[[3]])
            E.qbic <- (Reduce("+", select_muRR.qbic)/nsim)/Time
        }
        
        #QAIC
        PLL.qaic = lapply(1:nsim, function(i) (loglike[[i]]/overdisp[[i]]) - K[[i]])
        select.qaic <- lapply(1:nsim, function(i) which.max(unlist(PLL.qaic[[i]])))
        if(getRR==FALSE){
            select_mu.qaic <- lapply(1:nsim, function(i) sapply(select.qaic[[i]], function(j) mu[[i]][,j]))
            E.qaic <- (Reduce("+", select_mu.qaic)/nsim)/risk
        }
        else{
            select_mu.qaic <- lapply(1:nsim, function(i) sapply(select.qaic[[i]], function(j) mu[[i]][,j]*Ex[[i]]))
            select_muRR.qaic <- lapply(1:nsim, function(i) select_mu.qaic[[i]]/vectors[[3]])
            E.qaic <- (Reduce("+", select_muRR.qaic)/nsim)/Time
        }
        
        #QAICc
        PLL.qaicc=lapply(1:nsim, function(i) (loglike[[i]]/overdisp[[i]])- ((K[[i]]*n*Time)/(n*Time-K[[i]]-1)))
        select.qaicc <- lapply(1:nsim, function(i) which.max(unlist(PLL.qaicc[[i]])))
        if(getRR==FALSE){
            select_mu.qaicc <- lapply(1:nsim, function(i) sapply(select.qaicc[[i]], function(j) mu[[i]][,j]))
            E.qaicc <- Reduce("+", select_mu.qaicc)/nsim    
        }
        else{
            select_mu.qaicc <- lapply(1:nsim, function(i) sapply(select.qaicc[[i]], function(j) mu[[i]][,j]*Ex[[i]]))
            select_muRR.qaicc <- lapply(1:nsim, function(i) select_mu.qaicc[[i]]/vectors[[3]])
            E.qaicc <- (Reduce("+", select_muRR.qaicc)/nsim)/Time
        }
        
    }
    else{
        #BIC
        PLL.qbic  <- lapply(1:nsim, function(i) (loglike[[i]])-log(n*Time)/2*K[[i]])
        select.qbic <- lapply(1:nsim, function(i) which.max(unlist(PLL.qbic[[i]])))
        if(getRR==FALSE){
            select_mu.qbic <- lapply(1:nsim, function(i) sapply(select.qbic[[i]], function(j) mu[[i]][,j]))    
            E.qbic <- Reduce("+", select_mu.qbic)/nsim    
        }
        else{
            select_mu.qbic <- lapply(1:nsim, function(i) sapply(select.qbic[[i]], function(j) mu[[i]][,j]*Ex[[i]]))    
            select_muRR.qbic <- lapply(1:nsim, function(i) select_mu.qbic[[i]]/vectors[[3]])
            E.qbic <- Reduce("+", select_mu.qbic)/nsim    
        }
        
        #AIC
        PLL.qaic = lapply(1:nsim, function(i) (loglike[[i]]) - K[[i]])
        select.qaic <- lapply(1:nsim, function(i) which.max(unlist(PLL.qaic[[i]])))
        if(getRR==FALSE){
            select_mu.qaic <- lapply(1:nsim, function(i) sapply(select.qaic[[i]], function(j) mu[[i]][,j]))
            E.qaic <- Reduce("+", select_mu.qaic)/nsim
        }
        else{
            select_mu.qaic <- lapply(1:nsim, function(i) sapply(select.qaic[[i]], function(j) mu[[i]][,j]*Ex[[i]]))
            select_muRR.qaic <- lapply(1:nsim, function(i) select_mu.qaic[[i]]/vectors[[3]])
            E.qaic <- Reduce("+", select_mu.qaic)/nsim    
        }        
        
        #AICc
        PLL.qaicc=lapply(1:nsim, function(i) (loglike[[i]])- ((K[[i]]*n*Time)/(n*Time-K[[i]]-1)))
        select.qaicc <- lapply(1:nsim, function(i) which.max(unlist(PLL.qaicc[[i]])))
        if(getRR==FALSE){
            select_mu.qaicc <- lapply(1:nsim, function(i) sapply(select.qaicc[[i]], function(j) mu[[i]][,j]))
            E.qaicc <- Reduce("+", select_mu.qaicc)/nsim
        }
        else{
            select_mu.qaicc <- lapply(1:nsim, function(i) sapply(select.qaicc[[i]], function(j) mu[[i]][,j]*Ex[[i]]))
            select_muRR.qaicc <- lapply(1:nsim, function(i) select_mu.qaicc[[i]]/vectors[[3]])
            E.qaicc <- Reduce("+", select_mu.qaicc)/nsim
        }        
    }
    message("Returning results")
    return(list(nsim = nsim, E.qbic = E.qbic, E.qaic = E.qaic, E.qaicc = E.qaicc,Ex = Ex,mu = mu, Yx = Yx, PLL.qbic = PLL.qbic, 
                PLL.qaic = PLL.qaic, PLL.qaicc = PLL.qaicc, select.qbic = select.qbic, select.qaic = select.qaic, 
                select.qaicc = select.qaicc, select_mu.qbic = select_mu.qbic, select_mu.qaic = select_mu.qaic, 
                select_mu.qaicc = select_mu.qaicc, xbetaPath = xbetaPath, coefs.lasso.all = coefs.lasso.all))    
}



#' set.rr
#' 
#' This function will create vectors of the risk ratios as determined by observed counts, QBIC, QAIC, and QAICc, respectively.
#' @param lassoresult List of QBIC, QAIC, QAICc estimates from the mylasso function
#' @param vectors dataframe of initial vectors of the observed and expected counts
#' @param Time number of time period
#' @return This returns a list of the risk ratios (observed over expected) as determined by 1) pure observed/expected counts,
#' 2) observed based on QBIC path/expected; 3) observed based on QAIC path/expected; 4) observed based on QAICc path/expected.
#' @export
set.rr<- function(lassoresult, vectors, Time, sim=FALSE,...){
    if(sim==FALSE){
        RRobs <- matrix(as.vector(vectors$Y.vec)/as.vector(vectors$E0),ncol=Time)
        RRbic <- matrix(lassoresult$E.qbic/as.vector(vectors$E0),ncol=Time)
        RRaic <- matrix(lassoresult$E.qaic/as.vector(vectors$E0),ncol=Time)
        RRaicc <- matrix(lassoresult$E.qaicc/as.vector(vectors$E0),ncol=Time)
        message("Relative risks from observed data")
    }
    else{
        E0_avg <- Reduce("+", vectors$E0)/length(vectors$E0)
        RRobs <- matrix(as.vector(E0_avg)/as.vector(vectors$E0_fit),ncol=Time)
        RRbic <- matrix(lassoresult$E.qbic/as.vector(vectors$E0_fit),ncol=Time)
        RRaic <- matrix(lassoresult$E.qaic/as.vector(vectors$E0_fit),ncol=Time)
        RRaicc <- matrix(lassoresult$E.qaicc/as.vector(vectors$E0_fit),ncol=Time) 
        message("Relative risks from simulated data")
    }
    return(list(RRobs=RRobs, RRbic=RRbic, RRaic=RRaic, RRaicc=RRaicc))  
}


#' get.rr
#' 
#' This function will extract the relative risk ratios as determined by observed counts, QBIC, QAIC, and QAICc, respectively. This method determines the 
#' relative risk by averageing the relative risk in each simulation over the number of simulations. This is in contrast to $setRR$ where the expected counts
#' are first average over the number of simulations and then we determine the relative risk to the initial expected counts. Both methods return the same results,
#' however this method returns a smoother background ratio that does not vary across time periods as much
#' @param lassoresult List of QBIC, QAIC, QAICc estimates from the mylasso function
#' @param vectors dataframe of initial vectors of the observed and expected counts
#' @param Time number of time period
#' @param sim default is TRUE. Signals whether these values are to be extracted from a simulated model or real data.
#' @return This returns a list of the risk ratios (observed over expected) as determined by 1) pure observed/expected counts,
#' 2) observed based on QBIC path/expected; 3) observed based on QAIC path/expected; 4) observed based on QAICc path/expected.
#' @export
get.rr <- function(lassoresult,vectors, Time, sim=TRUE,...){
    if(sim==TRUE){
        E0_avg <- Reduce("+", vectors$E0)/length(vectors$E0)
        RRobs <- matrix(as.vector(E0_avg)/as.vector(vectors$E0_fit),ncol=Time)
        message("Relative risk ratios from simulated data - average RR over nsim")   
    }
    if(sim==FALSE){
        RRobs <- matrix(as.vector(vectors$Y.vec)/as.vector(vectors$E0),ncol=Time)
    }
    return(list(RRbic=matrix(lassoresult$E.qbic,ncol=Time),
                RRaic=matrix(lassoresult$E.qaic,ncol=Time),
                RRaicc=matrix(lassoresult$E.qaicc,ncol=Time),
                RRobs= RRobs))
}

#' redblue
#' 
#' This function establishes the spread of reds and blues for the risk ratios to be mapped to. Higher risk ratios will be deeper red colors and lower risk ratios will be deeper blue colors.
#' @param x this will be the risk ratios shrunk to be on the scale of half risk to twice the risk as end points.
#' @return colors
redblue=function(x,...) { 
    y=colorRamp(brewer.pal(11,"RdBu")[11:1])(x); rgb(y[,1],y[,2],y[,3],max=255) 
}


#' colormapping
#' 
#' This function establishes the spread of reds and blues for the risk ratios to be mapped to. Higher risk ratios will be deeper red colors and lower risk ratios will be deeper blue colors.
#' @param x this will be the risk ratios shrunk to be on the scale of half risk to twice the risk as end points.
#' @return returns vectors ofcolors for each time period, where risk ratios have been constrained to be between half risk and twice the risk
#' @export
colormapping <- function(riskratios,Time,...) {
    color.obs <- sapply(1:Time, function(i) redblue(log(2*pmax(1/2,pmin(riskratios$RRobs[,i],2)))/log(4)))
    color.qbic <- sapply(1:Time, function(i) redblue(log(2*pmax(1/2,pmin(riskratios$RRbic[,i],2)))/log(4))) 
    color.qaic <- sapply(1:Time, function(i) redblue(log(2*pmax(1/2,pmin(riskratios$RRaic[,i],2)))/log(4)))
    color.qaicc <- sapply(1:Time, function(i) redblue(log(2*pmax(1/2,pmin(riskratios$RRaicc[,i],2)))/log(4)))
    return(list(colors.obs = color.obs, color.qbic = color.qbic, color.qaic = color.qaic, color.qaicc = color.qaicc)) 
}

#'scale
#'
#'This function standardizes the fitted values by the sum of the observed/the sum of the expected within each time period. Inputs should all be in vector form.
#'The function handles standardizing within time period as long.
#'@param Y.vec vector of observed (in vector format)
#'@param out.sim simulated data based on the simulate(out) step above #TODO integrate this into this step above and streamline into a single function
#'@param nsim number of simulations performed in out.sim
#'@param Time number of time periods
#'@export
scale <- function(Y.vec, out.sim, nsim,Time,...){
    std <- lapply(1:nsim, function(i) sapply(1:Time, function(j) 
        (matrix(out.sim[[i]]$fitted.values,ncol=Time)[,j])*(sum(matrix(Y.vec,ncol=Time)[,j])/sum(matrix(out.sim[[i]]$fitted.values,ncol=Time)[,j]))))
    E0 <- lapply(1:nsim, function(i) as.vector(std[[i]])) 
    return(E0)
}


#'probmap
#'
#'This function will create a probability map based on simulation data. In each simulation, it identifies where a cluster was selected,
#'compared to the background rate. It then average over the number of simulations, giving us a matrix which ranges from 0 to 1 in probability.
#'To map this probabilities into a color scheme, please see the $colormapping$ function and select probmap=TRUE. TODO integrate all of this
#'into a workflow and extend to observed data, not only simulated data.
#'@param lassoresult List of QBIC, QAIC, QAICc estimates from the mylasso.sim function
#'@param vectors dataframe of initial vectors of the observed and expected counts that went into mylasso.sim function
#'@param rr risk ratio matrix that was used in the simulation
#'@param nsim number of simulations
#'@param Time number of time period
#'@param colormap default is FALSE. Signals whether the probabilities should directly be mapped to the red-blue color scheme. If this is false,
#'only the probability values will be returned. If true, then the probability values and the mapped colors will be returned in a list.
#'@return returns vector which calculated the number of time the cluster was correctly identified out of the simulations
#'@export
probmap <- function(lassoresult, vectors.sim, rr, nsim, Time, colormap=FALSE,...){
    prob.simBIC <- lapply(1:nsim, function(x) matrix(0, nrow(rr)*Time))
    prob.simAIC <- lapply(1:nsim, function(x) matrix(0, nrow(rr)*Time))
    prob.simAICc <- lapply(1:nsim, function(x) matrix(0, nrow(rr)*Time))
    indx <- which(rr !=1)
    rr.simBIC <- lapply(1:nsim, function(i) lassoresult$select_mu.qbic[[i]]/vectors.sim$E0_fit)
    rr.simAIC <- lapply(1:nsim, function(i) lassoresult$select_mu.qaic[[i]]/vectors.sim$E0_fit)
    rr.simAICc <- lapply(1:nsim, function(i) lassoresult$select_mu.qaicc[[i]]/vectors.sim$E0_fit)
    alphaBIC <- lapply(1:nsim, function(i) lapply(1:Time, function(k) 
        sort(table(matrix(rr.simBIC[[i]], ncol=Time)[,k]),decreasing=TRUE)[1]))
    alphaAIC <- lapply(1:nsim, function(i) lapply(1:Time, function(k) 
        sort(table(matrix(rr.simAIC[[i]], ncol=Time)[,k]),decreasing=TRUE)[1]))
    alphaAICc <- lapply(1:nsim, function(i) lapply(1:Time, function(k) 
        sort(table(matrix(rr.simAICc[[i]], ncol=Time)[,k]),decreasing=TRUE)[1]))
    
    for(j in 1:length(prob.simBIC)){
        for(i in 1:length(indx)){
            if (rr.simBIC[[j]][indx[i]] >= as.numeric(attributes(alphaBIC[[j]][[1]])[[1]]))  {
                prob.simBIC[[j]][indx[i]] <- 1
            }
            else {
                prob.simBIC[[j]][indx[i]] <- 0
            }
        }
    }
    for(j in 1:length(prob.simAIC)){
        for(i in 1:length(indx)){
            if (rr.simAIC[[j]][indx[i]] >= as.numeric(attributes(alphaAIC[[j]][[1]])[[1]]))  {
                prob.simAIC[[j]][indx[i]] <- 1
            }
            else {
                prob.simAIC[[j]][indx[i]] <- 0
            }
        }
    }
    for(j in 1:length(prob.simAICc)){
        for(i in 1:length(indx)){
            if (rr.simAICc[[j]][indx[i]] >= as.numeric(attributes(alphaAICc[[j]][[1]])[[1]]))  {
                prob.simAICc[[j]][indx[i]] <- 1
            }
            else {
                prob.simAICc[[j]][indx[i]] <- 0
            }
        }
    }
    prob <- as.vector(rr)
    prob[prob==1] <-0
    prob[prob!=0] <-1
    probBIC <- Reduce("+", prob.simBIC)/nsim
    probAIC <- Reduce("+", prob.simAIC)/nsim
    probAICc <- Reduce("+", prob.simAICc)/nsim
    if (colormap==TRUE){
        color.probmap <- sapply(1:Time, function(i) redblue(log(2*pmax(1/2,pmin(matrix(prob,ncol=Time)[,i],2)))/log(4)))    
        color.probmapBIC <- sapply(1:Time, function(i) redblue(log(2*pmax(1/2,pmin(matrix(probBIC,ncol=Time)[,i],2)))/log(4)))    
        color.probmapAIC <- sapply(1:Time, function(i) redblue(log(2*pmax(1/2,pmin(matrix(probAIC,ncol=Time)[,i],2)))/log(4)))    
        color.probmapAICc <- sapply(1:Time, function(i) redblue(log(2*pmax(1/2,pmin(matrix(probAICc,ncol=Time)[,i],2)))/log(4)))    
        return(list(prob = prob, probBIC = probBIC, probAIC = probAIC, probAICc = probAICc, 
                    color.probmap = color.probmap,
                    color.probmapBIC = color.probmapBIC,
                    color.probmapAIC = color.probmapAIC,
                    color.probmapAICc = color.probmapAICc))
    }
    else{
        return(list(prob=prob, probBIC = probBIC,probAIC = probAIC, probAICc = probAICc ))
    }
}    


#'detect.set
#'
#'This function will extract and set the necessary vectors in order to detect what percent of clusters
#'in the simulation were detected (true positives), what percent of clusters were incorrectly in the background (false negatives),
#'and what percent were incorrectly identified as clusters (false positives). 
#'@param lassoresult List of QBIC, QAIC, QAICc estimates from the mylasso.sim function
#'@param vectors.sim  dataframe of initial vectors of the observed and expected counts that went into simulation function
#'@param rr risk ratio matrix that was used in the simulation
#'@return returns a list of 1) indx_truth, which is a vector of indices where our cluster
#'should be, 2) indx, which contains the index of cluster as determined where the risk ratio in the rr matrix is not 1,
#'3) rr.simBIC, rr.simAIC, rr.simAICc - the risk ratios as determined by BIC, AIC, and AICc (respectively), and 4) alphaBIC, alphaAIC, alphaAICc - the calculated
#'background risk rates as determined by BIC, AIC, and AICc or QBIC, QAIC, QAICc. 
detect.set <- function(lassoresult, vectors.sim, rr, Time, x, y, rMax, center, radius,nullmod=NULL,...){
    #indx_truth <- sapply(1:Time, function(i) which(rr[,i] !=1))
    #create neighbors for cluster
    ##Indices of neighbors only in list nbs - space only
    #neighs <- findneighbors(nb, x, y, rMax, center, radius)
    
    #Indices of True cluster only
    if(!is.null(nullmod)){
        indx.clust.truth <- NULL
        message("Returning results for NULL model")
    }
    else{
        indx.clust.truth <- which(as.vector(rr) !=1, arr.ind=TRUE)    
    }
    
    rr.simAIC <- lapply(1:nsim, function(i) lassoresult$select_mu.qaic[[i]])
    rr.simAICc <- lapply(1:nsim, function(i) lassoresult$select_mu.qaicc[[i]])
    rr.simBIC <- lapply(1:nsim, function(i) lassoresult$select_mu.qbic[[i]])
    #Extract background rates as lists
    alpha.list.AIC <- lapply(1:nsim, function(i) lapply(1:Time, function(k)
        sort(table(matrix(rr.simAIC[[i]], ncol=Time)[,k]),decreasing=TRUE)[1]))
    alpha.list.AICc <- lapply(1:nsim, function(i) lapply(1:Time, function(k) 
        sort(table(matrix(rr.simAICc[[i]], ncol=Time)[,k]),decreasing=TRUE)[1]))
    alpha.list.BIC <- lapply(1:nsim, function(i) lapply(1:Time, function(k) 
        sort(table(matrix(rr.simBIC[[i]], ncol=Time)[,k]),decreasing=TRUE)[1]))
    
    #reformat alpha.list.AIC to be a matrix for each simulation
    alphaAIC <- lapply(1:nsim, function(i) matrix(as.numeric(attributes(unlist(alpha.list.AIC[[i]]))$names),
                                                  ncol=Time, byrow=TRUE, nrow=length(x)))
    alphaAICc <- lapply(1:nsim, function(i) matrix(as.numeric(attributes(unlist(alpha.list.AICc[[i]]))$names),
                                                  ncol=Time, byrow=TRUE, nrow=length(x)))
    alphaBIC <- lapply(1:nsim, function(i) matrix(as.numeric(attributes(unlist(alpha.list.BIC[[i]]))$names),
                                                   ncol=Time, byrow=TRUE, nrow=length(x)))
    return(list(
        indx.clust.truth = indx.clust.truth,
        #neighs = neighs,
        rr.simBIC = rr.simBIC,
        rr.simAIC = rr.simAIC,
        rr.simAICc = rr.simAICc,
        alphaBIC = alphaBIC,
        alphaAIC = alphaAIC,
        alphaAICc = alphaAICc))
}


#'detect.incluster.ic
#'This function will calculate detection based on the three information criterion. You can run detection on a null model (no cluster),
#'a model with an under-estimated cluster (artificial cluster <1), and on an elevated relative risk cluster. This is called by the more general
#'function *detect.incluster* which will switch to this function (TODO - Allow for detection options for AIC, AICc, and BIC only). 
#'@param lassoresult List of QBIC, QAIC, QAICc estimates from the mylasso.sim function
#'@param vectors.sim  dataframe of initial vectors of the observed and expected counts that went into simulation function
#'@param rr risk ratio matrix that was used in the simulation
#'@param res result of detect.set function
#'@param timeperiod Time period span of the cluster
#'@param Time number of time periods total in the model
#'@param nsim number of simulations run
#'@param under Default is NULL. If not null, then it will estimate detection based on a relative risk which is < 1. In other words, it will consider 
#'clusters to be identified where the estimated values are less than the background rate, not more than the background rate as is the case in the 
#'elevated relative risk models.
#'@param nullmod Default is NULL. If not null, then it will estimate detection based on the null model where there is no cluster. 
detect.incluster.ic <- function(lassoresult, vectors.sim, rr, set, period, Time, nsim, under=NULL,nullmod = NULL,...){
    prob.simBIC <- lapply(1:nsim, function(x) matrix(0, nrow(rr)*Time))
    prob.simAIC <- lapply(1:nsim, function(x) matrix(0, nrow(rr)*Time))
    prob.simAICc <- lapply(1:nsim, function(x) matrix(0, nrow(rr)*Time))
    #print(period)
    if(tail(period, n=1) == Time | tail(period, n=1)==1){
        maxTime = tail(period, n=1)    
    }
    else {
        maxTime = tail(period, n=1) +1
    }
    if(period[1] == 1){
        minTime = 1
    }
    else{
        minTime = period[1]-1
    }
    #(Q)AIC
    #extract things that are not the background rate
    if(!is.null(under)){
        ix <- lapply(1:nsim, function(i) which(round(matrix(set$rr.simAIC[[i]],ncol=Time),6) < round(set$alphaAIC[[i]],6), arr.ind=TRUE))
    }
    else{
        ix <- lapply(1:nsim, function(i) which(round(matrix(set$rr.simAIC[[i]],ncol=Time),6) > round(set$alphaAIC[[i]],6), arr.ind=TRUE))    
    }
    
    #1) Did it find anything in the cluster?
    ##NULL MODEL
    if(!is.null(nullmod)){
        clust <- lapply(1:nsim, function(i) is.element(ix[[i]][,1],set$neighs$cluster))
        in.clust.any <- lapply(1:nsim, function(i) if(isTRUE(length(clust[[i]])==0)) in.clust.any=0 else in.clust.any = 1)
        incluster.any.aic <- paste0((sum(unlist(in.clust.any))/nsim)*100,"%")
    }
    ##CLUSTER MODEL
    else{
        clust <- lapply(1:nsim, function(i) is.element(ix[[i]][,1],set$neighs$cluster))
        in.clust.any <- lapply(1:nsim, function(i) if(any(clust[[i]]==TRUE)) in.clust.any=1 else in.clust.any = 0)
        incluster.any.aic <- paste0((sum(unlist(in.clust.any))/nsim)*100,"%")    
    }

    #2) |(A and B)|/|A U B|?
    ##Calculate
    wasDetected <- lapply(1:nsim, function(i) which(round(matrix(set$rr.simAIC[[i]],ncol=Time),6) > round(set$alphaAIC[[i]],6), arr.ind=FALSE))    
    shouldDetected <- which(rr!=1, arr.ind=FALSE)
    ##Numerator
    AandB <- lapply(1:nsim, function(i) length(intersect(wasDetected[[i]], shouldDetected)))
    ##Denominator
    AuB <- lapply(1:nsim, function(i) length(union(wasDetected[[i]], shouldDetected)))
    ##Divide
    prop.alldetect.aic <- lapply(1:nsim, function(i) AandB[[i]]/AuB[[i]])
    
    
    #3) |(A and B)|/|A|? Proportion of what was detected was in overlap?
    ##Calculate length of what was detected
    A <- lapply(1:nsim, function(i) length(wasDetected[[i]]))
    ##Divide
    prop.wasdetect.aic <- lapply(1:nsim, function(i) AandB[[i]]/A[[i]])
    
    #4) |(A and B)|/|B|? Proportion of what should be detected was in overlap?
    B <- length(shouldDetected)
    prop.shoulddetect.aic <- lapply(1:nsim, function(i) AandB[[i]]/B)
    
    #5) ONLY FOR NULL MODEL - DID IT FIND ANYTHING?
    if(!is.null(nullmod)){
        null.any.aic <- length(unlist(ix))    
    }
    
    
    ########################################################################
    #(Q)AICc
    #extract things that are not the background rate
    if(!is.null(under)){
        ix <- lapply(1:nsim, function(i) which(round(matrix(set$rr.simAICc[[i]],ncol=Time),6) < round(set$alphaAICc[[i]],6), arr.ind=TRUE))
    }
    else{
        ix <- lapply(1:nsim, function(i) which(round(matrix(set$rr.simAICc[[i]],ncol=Time),6) > round(set$alphaAICc[[i]],6), arr.ind=TRUE))    
    }
    
    
    #1) Did it find anything in the cluster?
    ##NULL MODEL
    if(!is.null(nullmod)){
        clust <- lapply(1:nsim, function(i) is.element(ix[[i]][,1],set$neighs$cluster))
        in.clust.any <- lapply(1:nsim, function(i) if(isTRUE(length(clust[[i]])==0)) in.clust.any=0 else in.clust.any = 1)
        incluster.any.aicc <- paste0((sum(unlist(in.clust.any))/nsim)*100,"%")
    }
    else{
        clust <- lapply(1:nsim, function(i) is.element(ix[[i]][,1],set$neighs$cluster))
        in.clust.any <- lapply(1:nsim, function(i) if(any(clust[[i]]==TRUE)) in.clust.any=1 else in.clust.any = 0)
        incluster.any.aicc <- paste0((sum(unlist(in.clust.any))/nsim)*100,"%")    
    }
    #2) |(A and B)|/|A U B|?
    ##Calculate
    wasDetected <- lapply(1:nsim, function(i) which(round(matrix(set$rr.simAICc[[i]],ncol=Time),6) > round(set$alphaAICc[[i]],6), arr.ind=FALSE))    
    shouldDetected <- which(rr!=1, arr.ind=FALSE)
    ##Numerator
    AandB <- lapply(1:nsim, function(i) length(intersect(wasDetected[[i]], shouldDetected)))
    ##Denominator
    AuB <- lapply(1:nsim, function(i) length(union(wasDetected[[i]], shouldDetected)))
    ##Divide
    prop.alldetect.aicc <- lapply(1:nsim, function(i) AandB[[i]]/AuB[[i]])
    
    
    #3) |(A and B)|/|A|? Proportion of what was detected was in overlap?
    ##Calculate length of what was detected
    A <- lapply(1:nsim, function(i) length(wasDetected[[i]]))
    ##Divide
    prop.wasdetect.aicc <- lapply(1:nsim, function(i) AandB[[i]]/A[[i]])
    
    #4) |(A and B)|/|B|? Proportion of what should be detected was in overlap?
    B <- length(shouldDetected)
    prop.shoulddetect.aicc <- lapply(1:nsim, function(i) AandB[[i]]/B)
    
    #5) ONLY FOR NULL MODEL - DID IT FIND ANYTHING?
    if(!is.null(nullmod)){
        null.any.aicc <- length(unlist(ix))    
    }

    ################################################################
    #(Q)BIC
    #extract things that are not the background rate
    if(!is.null(under)){
        ix <- lapply(1:nsim, function(i) which(round(matrix(set$rr.simBIC[[i]],ncol=Time),6) < round(set$alphaBIC[[i]],6), arr.ind=TRUE))
    }
    else{
        ix <- lapply(1:nsim, function(i) which(round(matrix(set$rr.simBIC[[i]],ncol=Time),6) > round(set$alphaBIC[[i]],6), arr.ind=TRUE))    
    }
    
    #1) Did it find anything in the cluster?
    ##NULL MODEL
    if(!is.null(nullmod)){
        clust <- lapply(1:nsim, function(i) is.element(ix[[i]][,1],set$neighs$cluster))
        in.clust.any <- lapply(1:nsim, function(i) if(isTRUE(length(clust[[i]])==0)) in.clust.any=0 else in.clust.any = 1)
        incluster.any.bic <- paste0((sum(unlist(in.clust.any))/nsim)*100,"%")
    }
    else{
        clust <- lapply(1:nsim, function(i) is.element(ix[[i]][,1],set$neighs$cluster))
        in.clust.any <- lapply(1:nsim, function(i) if(any(clust[[i]]==TRUE)) in.clust.any=1 else in.clust.any = 0)
        incluster.any.bic <- paste0((sum(unlist(in.clust.any))/nsim)*100,"%")    
    }
    #2) |(A and B)|/|A U B|?
    ##Calculate
    wasDetected <- lapply(1:nsim, function(i) which(round(matrix(set$rr.simBIC[[i]],ncol=Time),6) > round(set$alphaBIC[[i]],6), arr.ind=FALSE))    
    shouldDetected <- which(rr!=1, arr.ind=FALSE)
    ##Numerator
    AandB <- lapply(1:nsim, function(i) length(intersect(wasDetected[[i]], shouldDetected)))
    ##Denominator
    AuB <- lapply(1:nsim, function(i) length(union(wasDetected[[i]], shouldDetected)))
    ##Divide
    prop.alldetect.bic <- lapply(1:nsim, function(i) AandB[[i]]/AuB[[i]])
    
    
    #3) |(A and B)|/|A|? Proportion of what was detected was in overlap?
    ##Calculate length of what was detected
    A <- lapply(1:nsim, function(i) length(wasDetected[[i]]))
    ##Divide
    prop.wasdetect.bic <- lapply(1:nsim, function(i) AandB[[i]]/A[[i]])
    
    #4) |(A and B)|/|B|? Proportion of what should be detected was in overlap?
    B <- length(shouldDetected)
    prop.shoulddetect.bic <- lapply(1:nsim, function(i) AandB[[i]]/B)
    
    #5) ONLY FOR NULL MODEL - DID IT FIND ANYTHING?
    if(!is.null(nullmod)){
        null.any.bic <- length(unlist(ix))    
    }
    if(exists("null.any.aic") & exists("null.any.bic") & exists("null.any.aicc")){
        return(list(
            incluster.any.aic = incluster.any.aic, incluster.any.aicc = incluster.any.aicc,incluster.any.bic = incluster.any.bic,
            prop.alldetect.aic = prop.alldetect.aic, prop.alldetect.aicc = prop.alldetect.aicc, prop.alldetect.bic = prop.alldetect.bic,
            prop.wasdetect.aic = prop.wasdetect.aic, prop.wasdetect.aicc = prop.wasdetect.aicc, prop.wasdetect.bic = prop.wasdetect.bic,
            prop.shoulddetect.aic = prop.shoulddetect.aic, prop.shoulddetect.aicc = prop.shoulddetect.aicc, prop.shoulddetect.bic = prop.shoulddetect.bic,
            null.any.aic = null.any.aic, null.any.aicc = null.any.aicc, null.any.bic = null.any.bic ))
    }
    else{
        return(list(
            incluster.any.aic = incluster.any.aic, incluster.any.aicc = incluster.any.aicc,incluster.any.bic = incluster.any.bic,
            prop.alldetect.aic = prop.alldetect.aic, prop.alldetect.aicc = prop.alldetect.aicc, prop.alldetect.bic = prop.alldetect.bic,
            prop.wasdetect.aic = prop.wasdetect.aic, prop.wasdetect.aicc = prop.wasdetect.aicc, prop.wasdetect.bic = prop.wasdetect.bic,
            prop.shoulddetect.aic = prop.shoulddetect.aic, prop.shoulddetect.aicc = prop.shoulddetect.aicc, prop.shoulddetect.bic = prop.shoulddetect.bic))
    }
    
}

  
 
 


#'detect.incluster
#'
#'This function will calculate the percent of simulations which correctly identify elements in cluster based on (Q)AIC, (Q)AICc, and (Q)BIC. The user can specify
#'if they want to only return one of these criterion or all three for further analysis.
#'@param lassoresult List of QBIC, QAIC, QAICc estimates from the mylasso.sim function
#'@param vectors.sim  dataframe of initial vectors of the observed and expected counts that went into simulation function
#'@param rr risk ratio matrix that was used in the simulation
#'@param res result of detect.set function
#'@param period_start time period where the cluster  starts in the simulation
#'@param period_end time period where cluster ends in the simulation
#'@param multi_period FALSE by default meaning that period_start and period_end are two unique periods. For example, if period_start = 2 and period_end =5, then
#'we will only be looking at periods 2 and 5. If multi_period is TRUE, then we will instead consider period_start through period_end (period_start:period_end). Following the same example,
#'this would mean we look at periods 2, 3, 4, and 5.
#'@param IC the information criteria you would like to be returned. Options are: IC = aic or IC = qaic; IC = aicc or IC = qaicc; IC = bic or IC = qbic; IC = ic or IC = qic 
#'(for aic/aicc/bic and qaic/qaicc/qbic, respectively).
#'TODO add cluster detection in case where risk ratio is less than background rate
#'@return returns
detect.incluster <- function(lassoresult, vectors.sim, rr, set,timeperiod, Time, nsim, x, y, rMax, center, radius, IC = c("aic","aicc","bic","ic"),nullmod=NULL,...){
    # if(multi_period==TRUE){
    #     period = c(period_start:period_end)
    # }
    # # else if(period_end - period_start == 0 ){
    # #    period = 1
    # # }
    # else{
    #     period = c(period_start, period_end)
    # }
    period = timeperiod
    message("Detection Results for:\n"
            , "\t Time Period: ", period,
            "\n \t Num. simulations: ", nsim,
            "\n \t Cluster center: ", center,
            "\n \t Cluster radius: ", radius,
            "\n \t Cluster rel.risk: ",unique(as.vector(rr))[2])
    IC <- match.arg(IC, several.ok= TRUE)
    switch(IC,
          aic = detect.incluster.aic(lassoresult, vectors.sim, rr, set, period, Time, nsim, under=NULL),
          aicc = detect.incluster.aicc(lassoresult, vectors.sim, rr, set, period, Time, nsim, under=NULL),
          bic = detect.incluster.bic(lassoresult, vectors.sim, rr, set, period, Time, nsim, under=NULL),
          ic = detect.incluster.ic(lassoresult, vectors.sim, rr, set, period, Time, nsim, under=NULL, nullmod=NULL))
} 


#'detect.falsecluster
#'
#'This function will calculate the percent of simulations which correctly identify elements in cluster based on (Q)AIC, (Q)AICc, and (Q)BIC. The user can specify
#'if they want to only return one of these criterion or all three for further analysis.
#'@param lassoresult List of QBIC, QAIC, QAICc estimates from the mylasso.sim function
#'@param vectors.sim  dataframe of initial vectors of the observed and expected counts that went into simulation function
#'@param rr risk ratio matrix that was used in the simulation
#'@param res result of detect.set function
#'@param period_start time period where the cluster  starts in the simulation
#'@param period_end time period where cluster ends in the simulation
#'@param multi_period FALSE by default meaning that period_start and period_end are two unique periods. For example, if period_start = 2 and period_end =5, then
#'we will only be looking at periods 2 and 5. If multi_period is TRUE, then we will instead consider period_start through period_end (period_start:period_end). Following the same example,
#'this would mean we look at periods 2, 3, 4, and 5.
#'@param IC the information criteria you would like to be returned. Options are: IC = aic or IC = qaic; IC = aicc or IC = qaicc; IC = bic or IC = qbic; IC = ic or IC = qic 
#'(for aic/aicc/bic and qaic/qaicc/qbic, respectively).
#'@return returns list of average false clusters across all simulations and the lowest number of false clusters detected across clusters

detect.falsecluster <- function(lassoresult, vectors.sim, rr, res, period_start, period_end, multi_period = FALSE, IC = c("aic","aicc","bic","ic"), Time,...){
    if(multi_period==TRUE){
        period = period_start:period_end
    }
    else{
        period = c(period_start, period_end)
    }
    IC <- match.arg(IC, several.ok= TRUE)
    switch(IC,
           aic = detect.falsecluster.aic(lassoresult, vectors.sim, rr, res, period, Time),
           aicc = detect.falsecluster.aicc(lassoresult, vectors.sim, rr, res, period, Time),
           bic = detect.falsecluster.bic(lassoresult, vectors.sim, rr, res, period, Time),
           ic = detect.falsecluster.ic(lassoresult, vectors.sim, rr, res, period, Time))
} 
    
#'
#'detect.falsecluster.aic
#'
#'        
#'                        
detect.falsecluster.aic <- function(lassoresult, vectors.sim, rr, res, period, Time){
    prob.simAIC <- lapply(1:nsim, function(x) matrix(0, nrow(rr)*Time))
    ix <- lapply(1: length(prob.simAIC), 
                 function(j) sapply(1:Time, 
                                    function(k) 
                                        which(round(matrix(res$rr.simAIC[[j]], ncol=Time)[,k],6) > round(as.numeric(attributes(res$alphaAIC[[j]][[k]])),6))))
    
    in_fake <- lapply(1:length(prob.simAIC), 
                      function(j) sapply(1:Time, 
                                         function(k) length(which(ix[[j]][[k]] %in% res$indx_truth[[k]] == FALSE))/length(res$indx_truth[[k]])))
    idx_avg <- unlist(lapply(1:length(prob.simAIC), 
                             function(j) mean(in_fake[[j]][period])))
    mean_fp = mean(idx_avg)
    min_fp = min(idx_avg)
    max_fp = max(idx_avg)
    sd_fp = sd(idx_avg)
    message("Results for (Q)AIC")
    return(list(
        mean_fp = mean_fp, min_fp = min_fp, max_fp = max_fp, sd_fp = sd_fp))								
}


#'
#'detect.falsecluster.aicc
#'
#'        
#'                        
detect.falsecluster.aicc <- function(lassoresult, vectors.sim, rr, res, period, Time){
    prob.simAICc <- lapply(1:nsim, function(x) matrix(0, nrow(rr)*Time))
    ix <- lapply(1: length(prob.simAICc), 
                 function(j) sapply(1:Time, 
                                    function(k) 
                                        which(round(matrix(res$rr.simAICc[[j]], ncol=Time)[,k],6) > round(as.numeric(attributes(res$alphaAICc[[j]][[k]])),6))))
    
    in_fake <- lapply(1:length(prob.simAICc), 
                      function(j) sapply(1:Time, 
                                         function(k) length(which(ix[[j]][[k]] %in% res$indx_truth[[k]] == FALSE))/length(res$indx_truth[[k]])))
    idx_avg <- unlist(lapply(1:length(prob.simAICc), 
                             function(j) mean(in_fake[[j]][period])))
    mean_fp = mean(idx_avg)
    min_fp = min(idx_avg)
    max_fp = max(idx_avg)
    sd_fp = sd(idx_avg)
    message("Results for (Q)AICc")
    return(list(
        mean_fp = mean_fp, min_fp = min_fp, max_fp = max_fp, sd_fp = sd_fp))
}



#'
#'detect.falsecluster.bic
#'
#'        
#'                        
detect.falsecluster.bic <- function(lassoresult, vectors.sim, rr, res, period, Time){
    prob.simBIC <- lapply(1:nsim, function(x) matrix(0, nrow(rr)*Time))
    ix <- lapply(1: length(prob.simBIC), 
                 function(j) sapply(1:Time, 
                                    function(k) 
                                        which(round(matrix(res$rr.simBIC[[j]], ncol=Time)[,k],6) > round(as.numeric(attributes(res$alphaBIC[[j]][[k]])),6))))
    
    in_fake <- lapply(1:length(prob.simBIC), 
                      function(j) sapply(1:Time, 
                                         function(k) length(which(ix[[j]][[k]] %in% res$indx_truth[[k]] == FALSE))/length(res$indx_truth[[k]])))
    idx_avg <- unlist(lapply(1:length(prob.simBIC), 
                             function(j) mean(in_fake[[j]][period])))
    mean_fp = mean(idx_avg)
    min_fp = min(idx_avg)
    max_fp = max(idx_avg)
    sd_fp = sd(idx_avg)
    message("Results for (Q)BIC")
    return(list(
        mean_fp = mean_fp, min_fp = min_fp, max_fp = max_fp, sd_fp = sd_fp))
}




#'
#'detect.falsecluster.ic
#'
#'        
#'                        
detect.falsecluster.ic <- function(lassoresult, vectors.sim, rr, res, period, Time){
    prob.simBIC <- lapply(1:nsim, function(x) matrix(0, nrow(rr)*Time))
    prob.simAIC <- lapply(1:nsim, function(x) matrix(0, nrow(rr)*Time))
    prob.simAICc <- lapply(1:nsim, function(x) matrix(0, nrow(rr)*Time))
    #(Q)AIC
    ix <- lapply(1: length(prob.simAIC), 
                 function(j) sapply(1:Time, 
                                    function(k) 
                                        which(round(matrix(res$rr.simAIC[[j]], ncol=Time)[,k],6) > round(as.numeric(attributes(res$alphaAIC[[j]][[k]])),6))))
    
    in_fake <- lapply(1:length(prob.simAICc), 
                      function(j) sapply(1:Time, 
                                         function(k) length(which(ix[[j]][[k]] %in% res$indx_truth[[k]] == FALSE))/length(res$indx_truth[[k]])))
    idx_avg <- unlist(lapply(1:length(prob.simAIC), 
                             function(j) mean(in_fake[[j]][period])))
    
    mean_fp.aic = mean(idx_avg)
    min_fp.aic = min(idx_avg)
    max_fp.aic = max(idx_avg)
    sd_fp.aic = sd(idx_avg)
    
    #(Q)AICc
    ix <- lapply(1: length(prob.simAICc), 
                 function(j) sapply(1:Time, 
                                    function(k) 
                                        which(round(matrix(res$rr.simAICc[[j]], ncol=Time)[,k],6) > round(as.numeric(attributes(res$alphaAICc[[j]][[k]])),6))))
    
    in_fake <- lapply(1:length(prob.simAICc), 
                      function(j) sapply(1:Time, 
                                         function(k) length(which(ix[[j]][[k]] %in% res$indx_truth[[k]] == FALSE))/length(res$indx_truth[[k]])))
    idx_avg <- unlist(lapply(1:length(prob.simAICc), 
                             function(j) mean(in_fake[[j]][period])))
    
    mean_fp.aicc = mean(idx_avg)
    min_fp.aicc = min(idx_avg)
    max_fp.aicc = max(idx_avg)
    sd_fp.aicc = sd(idx_avg)
    
    #(Q)BIC
    ix <- lapply(1: length(prob.simBIC), 
                 function(j) sapply(1:Time, 
                                    function(k) 
                                        which(round(matrix(res$rr.simBIC[[j]], ncol=Time)[,k],6) > round(as.numeric(attributes(res$alphaBIC[[j]][[k]])),6))))
    
    in_fake <- lapply(1:length(prob.simBIC), 
                      function(j) sapply(1:Time, 
                                         function(k) length(which(ix[[j]][[k]] %in% res$indx_truth[[k]] == FALSE))/length(res$indx_truth[[k]])))
    idx_avg <- unlist(lapply(1:length(prob.simBIC), 
                             function(j) mean(in_fake[[j]][period])))
    
    mean_fp.bic = mean(idx_avg)
    min_fp.bic = min(idx_avg)
    max_fp.bic = max(idx_avg)
    sd_fp.bic = sd(idx_avg)
    
    return(list(mean_fp.aic = mean_fp.aic, min_fp.bic = min_fp.aic, max_fp.aic = max_fp.aic, sd_fp.aic = sd_fp.aic,
                mean_fp.aicc = mean_fp.aicc, min_fp.aicc = min_fp.aicc, max_fp.aicc = max_fp.aicc, sd_fp.aicc = sd_fp.aicc,
                mean_fp.bic = mean_fp.bic, min_fp.bic = min_fp.bic, max_fp.bic = max_fp.bic, sd_fp.bic = sd_fp.bic))
}


#'detect.inbackground
#'
#'This function will calculate the percent of simulations which correctly identify elements in cluster based on (Q)AIC, (Q)AICc, and (Q)BIC. The user can specify
#'if they want to only return one of these criterion or all three for further analysis.
#'@param lassoresult List of QBIC, QAIC, QAICc estimates from the mylasso.sim function
#'@param vectors.sim  dataframe of initial vectors of the observed and expected counts that went into simulation function
#'@param rr risk ratio matrix that was used in the simulation
#'@param res result of detect.set function
#'@param period_start time period where the cluster  starts in the simulation
#'@param period_end time period where cluster ends in the simulation
#'@param multi_period FALSE by default meaning that period_start and period_end are two unique periods. For example, if period_start = 2 and period_end =5, then
#'we will only be looking at periods 2 and 5. If multi_period is TRUE, then we will instead consider period_start through period_end (period_start:period_end). Following the same example,
#'this would mean we look at periods 2, 3, 4, and 5.
#'@param IC the information criteria you would like to be returned. Options are: IC = aic or IC = qaic; IC = aicc or IC = qaicc; IC = bic or IC = qbic; IC = ic or IC = qic 
#'(for aic/aicc/bic and qaic/qaicc/qbic, respectively).
#'@return returns list of average false clusters across all simulations and the lowest number of false clusters detected across clusters

detect.inbackground <-function(lassoresult, vectors.sim, rr, res, period_start, period_end, multi_period=FALSE, IC = c("aic","aicc","bic","ic"), Time,...){
    if(multi_period==TRUE){
        period = period_start:period_end
    }
    else{
        period = c(period_start, period_end)
    }
    IC <- match.arg(IC, several.ok= TRUE)
    switch(IC,
           aic = detect.inbackground.aic(lassoresult, vectors.sim, rr, res, period, Time),
           aicc = detect.inbackground.aicc(lassoresult, vectors.sim, rr, res, period, Time),
           bic = detect.inbackground.bic(lassoresult, vectors.sim, rr, res, period, Time),
           ic = detect.inbackground.ic(lassoresult, vectors.sim, rr, res, period, Time))
}

#'detect.inbackground.aic
#'
#'
#'TODO
#'
detect.inbackground.aic <- function(lassoresult, vectors.sim, rr, res, period, Time){
    prob.simAIC <- lapply(1:nsim, function(x) matrix(0, nrow(rr)*Time))
    ix <- lapply(1: length(prob.simAIC), 
                 function(j) sapply(1:Time, 
                                    function(k) 
                                        which(round(matrix(res$rr.simAIC[[j]], ncol=Time)[,k],6) > round(as.numeric(attributes(res$alphaAIC[[j]][[k]])),6))))
    in_background <- lapply(1:length(prob.simAIC), 
                            function(j) sapply(1:Time, 
                                               function(k) 
                                                   ((sum((ix[[j]][[k]] %in% res$indx_truth[[k]])*1) - 
                                                         length(res$indx_truth[[k]]))/length(res$indx_truth[[k]]))*-1))
    in_bsum <- lapply(1:length(prob.simAIC), 
                      function(j) mean(in_background[[j]][period]))
    mean_bkgrd = mean(unlist(in_bsum))
    max_bkgrd = max(unlist(in_bsum))
    min_bkgrd = min(unlist(in_bsum))
    sd_bkgrd = sd(unlist(in_bsum))
    return(list(
        mean_bkgrd = mean_bkgrd,
        max_bkgrd = max_bkgrd,
        min_bkgrd = min_bkgrd,
        sd_bkgrd = sd_bkgrd
    ))
}



#'detect.inbackground.aicc
#'
#'
#'TODO
#'
detect.inbackground.aicc <- function(lassoresult, vectors.sim, rr, res, period, Time){
    prob.simAICc <- lapply(1:nsim, function(x) matrix(0, nrow(rr)*Time))
    ix <- lapply(1: length(prob.simAICc), 
                 function(j) sapply(1:Time, 
                                    function(k) 
                                        which(round(matrix(res$rr.simAICc[[j]], ncol=Time)[,k],6) < round(as.numeric(attributes(res$alphaAICc[[j]][[k]])),6))))
    in_background <- lapply(1:length(prob.simAICc), 
                            function(j) sapply(1:Time, 
                                               function(k) 
                                                   ((sum((ix[[j]][[k]] %in% res$indx_truth[[k]])*1) - 
                                                         length(res$indx_truth[[k]]))/length(res$indx_truth[[k]]))*-1))
    in_bsum <- lapply(1:length(prob.simAICc), 
                      function(j) mean(in_background[[j]][period]))
    mean_bkgrd = mean(unlist(in_bsum))
    max_bkgrd = max(unlist(in_bsum))
    min_bkgrd = min(unlist(in_bsum))
    sd_bkgrd = sd(unlist(in_bsum))
    return(list(
        mean_bkgrd = mean_bkgrd,
        max_bkgrd = max_bkgrd,
        min_bkgrd = min_bkgrd,
        sd_bkgrd = sd_bkgrd
    ))
}	



#'detect.inbackground.bic
#'
#'
#'TODO
#'
detect.inbackground.bic <- function(lassoresult, vectors.sim, rr, res, period, Time){
    prob.simBIC <- lapply(1:nsim, function(x) matrix(0, nrow(rr)*Time))
    ix <- lapply(1: length(prob.simBIC), 
                 function(j) sapply(1:Time, 
                                    function(k) 
                                        which(round(matrix(res$rr.simBIC[[j]], ncol=Time)[,k],6) < round(as.numeric(attributes(res$alphaBIC[[j]][[k]])),6))))
    in_background <- lapply(1:length(prob.simBIC), 
                            function(j) sapply(1:Time, 
                                               function(k) 
                                                   ((sum((ix[[j]][[k]] %in% res$indx_truth[[k]])*1) - 
                                                         length(res$indx_truth[[k]]))/length(res$indx_truth[[k]]))*-1))
    in_bsum <- lapply(1:length(prob.simBIC), 
                      function(j) mean(in_background[[j]][period]))
    mean_bkgrd = mean(unlist(in_bsum))
    max_bkgrd = max(unlist(in_bsum))
    min_bkgrd = min(unlist(in_bsum))
    sd_bkgrd = sd(unlist(in_bsum))
    return(list(
        mean_bkgrd = mean_bkgrd,
        max_bkgrd = max_bkgrd,
        min_bkgrd = min_bkgrd,
        sd_bkgrd = sd_bkgrd
    ))
}	



#'detect.inbackground.ic
#'
#'
#'TODO
#'
detect.inbackground.ic <- function(lassoresult, vectors.sim, rr, res, period, Time){
    prob.simBIC <- lapply(1:nsim, function(x) matrix(0, nrow(rr)*Time))
    prob.simAIC <- lapply(1:nsim, function(x) matrix(0, nrow(rr)*Time))
    prob.simAICc <- lapply(1:nsim, function(x) matrix(0, nrow(rr)*Time))
    #(Q)AIC
    ix <- lapply(1: length(prob.simAIC), 
                 function(j) sapply(1:Time, 
                                    function(k) 
                                        which(round(matrix(res$rr.simAIC[[j]], ncol=Time)[,k],6) < round(as.numeric(attributes(res$alphaAIC[[j]][[k]])),6))))
    in_background <- lapply(1:length(prob.simAIC), 
                            function(j) sapply(1:Time, 
                                               function(k) 
                                                   ((sum((ix[[j]][[k]] %in% res$indx_truth[[k]])*1) - 
                                                         length(res$indx_truth[[k]]))/length(res$indx_truth[[k]]))*-1))
    in_bsum <- lapply(1:length(prob.simAIC), 
                      function(j) mean(in_background[[j]][period]))
    mean_bkgrd.aic = mean(unlist(in_bsum))
    max_bkgrd.aic = max(unlist(in_bsum))
    #(Q)AICc
    ix <- lapply(1: length(prob.simAICc), 
                 function(j) sapply(1:Time, 
                                    function(k) 
                                        which(round(matrix(res$rr.simAICc[[j]], ncol=Time)[,k],6) < round(as.numeric(attributes(res$alphaAICc[[j]][[k]])),6))))
    in_background <- lapply(1:length(prob.simAICc), 
                            function(j) sapply(1:Time, 
                                               function(k) 
                                                   ((sum((ix[[j]][[k]] %in% res$indx_truth[[k]])*1) - 
                                                         length(res$indx_truth[[k]]))/length(res$indx_truth[[k]]))*-1))
    in_bsum <- lapply(1:length(prob.simAICc), 
                      function(j) mean(in_background[[j]][period]))
    mean_bkgrd.aicc = mean(unlist(in_bsum))
    max_bkgrd.aicc = max(unlist(in_bsum))
    #(Q)BIC
    ix <- lapply(1: length(prob.simBIC), 
                 function(j) sapply(1:Time, 
                                    function(k) 
                                        which(round(matrix(res$rr.simBIC[[j]], ncol=Time)[,k],6) < round(as.numeric(attributes(res$alphaBIC[[j]][[k]])),6))))
    in_background <- lapply(1:length(prob.simBIC), 
                            function(j) sapply(1:Time, 
                                               function(k) 
                                                   ((sum((ix[[j]][[k]] %in% res$indx_truth[[k]])*1) - 
                                                         length(res$indx_truth[[k]]))/length(res$indx_truth[[k]]))*-1))
    in_bsum <- lapply(1:length(prob.simBIC), 
                      function(j) mean(in_background[[j]][period]))
    mean_bkgrd.bic = mean(unlist(in_bsum))
    max_bkgrd.bic = max(unlist(in_bsum))
    
    return(list(mean_bkgrd.aic = mean_bkgrd.aic, max_bkgrd.aic = max_bkgrd.aic,
                mean_bkgrd.aicc = mean_bkgrd.aicc, max_bkgrd.aicc = max_bkgrd.aicc,
                mean_bkgrd.bic = mean_bkgrd.bic, max_bkgrd.bic = max_bkgrd.bic
    ))
    
}


#'detect
#'
#'This function will create a probability map based on simulation data. In each simulation, it identifies where a cluster was selected,
#'compared to the background rate. It then average over the number of simulations, giving us a matrix which ranges from 0 to 1 in probability.
#'To map this probabilities into a color scheme, please see the $colormapping$ function and select probmap=TRUE. TODO integrate all of this
#'into a workflow and extend to observed data, not only simulated data.
#'@param lassoresult List of QBIC, QAIC, QAICc estimates from the mylasso.sim function
#'@param vectors.sim dataframe of initial vectors of the observed and expected counts that went into mylasso.sim function
#'@param rr risk ratio matrix that was used in the simulation
#'@param period_start start period of the simulation
#'@param period_end end period of the simulation
#'@param multi_period FALSE by default meaning that period_start and period_end are two unique periods. For example, if period_start = 2 and period_end =5, then
#'we will only be looking at periods 2 and 5. If multi_period is TRUE, then we will instead consider period_start through period_end (period_start:period_end). Following the same example,
#'this would mean we look at periods 2, 3, 4, and 5.
#'@param IC the information criteria you would like to be returned. Options are: IC = aic or IC = qaic; IC = aicc or IC = qaicc; IC = bic or IC = qbic; IC = ic or IC = qic 
#'(for aic/aicc/bic and qaic/qaicc/qbic, respectively).
#'@return returns a list of lists. Lists include proportion of simulations which correctly identify the cluster[in_cluster], proportion of simulations which incorrectly identify
#'a cluster that is not one (false positive) [false_cluster], and proportion of simulations which incorrectly identify a true cluster in the background rate (false negatives)[in_background]. 
#'All of these are reported as determined by (Q)AIC, (Q)AICc, (Q)BIC, or all three IC.
#'@export
#'
detect <- function(lassoresult, vectors.sim, rr, period_start, period_end, multi_period = FALSE, IC = NULL, Time,...){
    #determined time period span
    if(multi_period==TRUE){
        period = period_start:period_end
    }
    else{
        period = c(period_start, period_end)
    }
    #run set-up
    res <- detect.set(lassoresult, vectors.sim, rr, Time)
    IC = IC
    #run detection
    in_cluster <- detect.incluster(lassoresult, vectors.sim, rr, res, period_start, period_end, multi_period, IC, Time)
    false_cluster <- detect.falsecluster(lassoresult, vectors.sim, rr, period_start, period_end, multi_period, IC, Time)
    in_background <- detect.inbackground(lassoresult, vectors.sim, rr, period_start, period_end, multi_period, IC, Time)
    #return as lists
    return(list(in_cluster = in_cluster,
                false_cluster = false_cluster,
                in_background = in_background
                ))
}


#'clust
#'
#'This function runs both the space and space-time Lasso model. This function is to be run on observed data. A separate function (clust.sim) can be used for simulating data and running diagnostics on simulations.
#'@param x x coordinates (easting/latitude); if utm coordinates, scale to km.
#'@param y y coordinates (northing/longitude); if utm coordinates, scale to km.
#'@param rMax set max radius (in km)
#'@param period vector of periods or years in dataset. Should be imported as a factor.
#'@param expected vector of expected counts. Expected counts must match up with the year and observed vectors.
#'@param observed vector of observed counts. Observed counts must match up with the year and expected vectors.
#'@param Time Number of time periods or years in your dataset. Must be declared as numeric.
#'@param spacetime default is TRUE. To run the space-only model, specify `spacetime=FALSE'
#'@return
#'@details Optional functions include:
#'- 1) utm - default is FALSE. If you have utm coordinates, you want to change this to TRUE.
#'@export
clust <- function(x, y, rMax, period, expected, observed, Time, spacetime=TRUE, pois = FALSE,colors=NULL,utm=TRUE, byrow=TRUE,...){
    if(utm==FALSE){
        message("Coordinates are assumed to be in lat/long coordinates. For utm coordinates, please specify 'utm=TRUE'")
        utm=FALSE
    }
    else{
        utm=TRUE
    }
    if(byrow==FALSE){
        row=FALSE
    }
    else{
        row=TRUE
        message("Data assumed to be in panel data. To use vector data instead, please specify 'byrow=FALSE'")
    }
    if(pois==TRUE){
        pois=TRUE
    }
    else{
        pois=FALSE
        message("Running quasi-Poisson model. For Poisson model, please specify 'pois=TRUE'")
    }
    #set up clusters and fitted values
    clusters <- clusters.df(x,y,rMax, utm=utm, length(x))
    n <- length(x)
    #init <- set.vectors(dframe$period, dframe$expdeath, dframe$death, Time=Time, byrow=row)
    init <- set.vectors(period, expected, observed, Time, row)
    E0 <- scale(init, Time, scaler = 1)
    vectors <- list(Period = init$Year, E0=E0, E0_0=init$E0, Y.vec=init$Y.vec)
    lassoresult <- spacetime.lasso(clusters, vectors, Time, spacetime=spacetime,pois=pois)
    riskratios <- get.rr(lassoresult, vectors, Time, sim=FALSE)
    rrcolors <- colormapping(riskratios,Time)
    return(list(lassoresult = lassoresult,
                riskratios = riskratios,
                rrcolors = rrcolors,
                init.vec = init))
}




#'clust.sim
#'
#'This function runs both the space and space-time Lasso model simulations. This function is to be run on simulated data. A separate function (clust) can be used for observed data.
#'@param x x coordinates (easting/latitude); if utm coordinates, scale to km.
#'@param y y coordinates (northing/longitude); if utm coordinates, scale to km.
#'@param rMax set max radius (in km)
#'@param period vector of periods or years in dataset. Should be imported as a factor.
#'@param expected vector of expected counts. Expected counts must match up with the year and observed vectors.
#'@param observed vector of observed counts. Observed counts must match up with the year and expected vectors.
#'@param Time Number of time periods or years in your dataset. Must be declared as numeric.
#'@param spacetime default is TRUE. To run the space-only model, specify `spacetime=FALSE'
#'@param nsim Number of simulations you would like to run
#'@param center can be a single center or for multiple clusters, concatenate them. Max three TODO extend this
#'@param radius radius for the cluster you want in the simulation
#'@param risk.ratio setting for what the risk ratio should be in the cluster to be detected by the simulation
#'@param timeperiod time period where the cluster  starts in the simulation
#'we will only be looking at periods 2 and 5. If multi_period is TRUE, then we will instead consider timeperiod through period_end (timeperiod:period_end). Following the same example,
#'this would mean we look at periods 2, 3, 4, and 5.
#'@return
#'@details Optional functions include:
#'- 1) utm - default is FALSE. If you have utm coordinates, you want to change this to TRUE.
#'@export
#'TODO allow user to change theta parameter in simulation
clust.sim <- function(x, y, rMax, period, expected, observed, Time, spacetime=TRUE, pois = FALSE, nsim, center, radius, risk.ratio, timeperiod,colors=NULL,utm=TRUE, byrow=TRUE,...){
    #initial user setting
    if(utm==FALSE){
        message("Coordinates are assumed to be in lat/long coordinates. For utm coordinates, please specify 'utm=TRUE'")
        utm=FALSE
    }
    else{
        utm=TRUE
    }
    if(byrow==FALSE){
        row=FALSE
    }
    else{
        row=TRUE
        message("Data assumed to be in panel data. To use vector data instead, please specify 'byrow=FALSE'")
    }
    #set up clusters and fitted values
    clusters <- clusters.df(x,y,rMax, utm=TRUE, length(x))
    n <- length(x)
    init <- set.vectors(period, expected, observed, Time=Time, byrow=row)
    
    #TODO change this to be a function and not hard-coded
    if(length(center) == 2){
        tmp <- clusters[clusters$center==center[1] | clusters$center==center[2],]
    }
    else if(length(center) == 3){
        tmp <- clusters[clusters$center==center[1] | clusters$center==center[2] | clusters$center==center[3],]
    }
    else{
        tmp <- clusters[clusters$center==center,]
    }
    cluster <- tmp[(tmp$r <= radius),]
    rr = matrix(1, nrow=n, ncol=Time)
    #TODO change this to be a function and not hard-coded
    if(length(timeperiod) == 3){
        rr[cluster$last, timeperiod[1]] = risk.ratio
        rr[cluster$last, timeperiod[2]] = risk.ratio
        rr[cluster$last, timeperiod[3]] = risk.ratio
        message(paste("Running model for periods",timeperiod[1],"through", timeperiod[3]))
    }
    if(length(timeperiod) == 4){
        rr[cluster$last, timeperiod[1]] = risk.ratio
        rr[cluster$last, timeperiod[2]] = risk.ratio
        rr[cluster$last, timeperiod[3]] = risk.ratio
        rr[cluster$last, timeperiod[4]] = risk.ratio
        message(paste("Running model for periods",timeperiod[1],"through", timeperiod[4]))
    }
    if(length(timeperiod) == 5){
        rr[cluster$last, timeperiod[1]] = risk.ratio
        rr[cluster$last, timeperiod[2]] = risk.ratio
        rr[cluster$last, timeperiod[3]] = risk.ratio
        rr[cluster$last, timeperiod[4]] = risk.ratio
        rr[cluster$last, timeperiod[5]] = risk.ratio
        message(paste("Running model for periods",timeperiod[1],"through", timeperiod[5]))
    }
    if(length(timeperiod) == 2){
        rr[cluster$last, timeperiod[1]] = risk.ratio
        rr[cluster$last, timeperiod[2]] = risk.ratio
        message(paste("Running model for periods",timeperiod[1],"and", timeperiod[2]))
    }
    else if(length(timeperiod)==1){
        rr[cluster$last, timeperiod:Time] = risk.ratio    
    }
    E1 = as.vector(rr)*init$E0
    Period <- init$Year
    #Simulate observed as NB(Eit, theta)
    theta = 1000
    YSIM <- lapply(1:nsim, function(i) rnegbin(E1, theta = theta))
    #Scale YSIM[i] to init$E0
    Ex <- scale.sim(YSIM, init, nsim, Time)
    vectors.sim <- list(Period = Period, Ex = Ex , E0_0 = init$E0, Y.vec=init$Y.vec)
    #set up and run simulation model
    lassoresult <- spacetime.lasso.sim(clusters, vectors.sim, Time, spacetime=spacetime, pois=pois, nsim, YSIM)
    riskratios <- get.rr2(lassoresult, vectors.sim,init, E1,Time, sim=TRUE)
    rrcolors <- colormapping(riskratios,Time)
    # if(!is.null(colors)){
    #     probcolors <- probmap(lassoresult, vectors.sim, rr, nsim,Time, colormap=TRUE)    
    # }
    # else{
    #     probcolors <- probmap(lassoresult, vectors.sim, rr, nsim,Time, colormap=FALSE)    
    # }
    return(list(lassoresult = lassoresult,
                riskratios = riskratios,
                rrcolors = rrcolors,
                rr.mat = rr,
                init.vec = vectors.sim))
}
    


#'clust.sim.all
#'
#'This function runs both the space and space-time Lasso model simulations for all 4 models simulataneously: Quasi-Poisson vs. Poisson in both space and space-time.
#' This function is to be run on simulated data and all four models are run on the same simulated set. 
#'A separate function (clust.sim) can be used for running simulations on individual models and (clust) can be used for observed data.
#'@param x x coordinates (easting/latitude); if utm coordinates, scale to km.
#'@param y y coordinates (northing/longitude); if utm coordinates, scale to km.
#'@param rMax set max radius (in km)
#'@param period vector of periods or years in dataset. Should be imported as a factor.
#'@param expected vector of expected counts. Expected counts must match up with the year and observed vectors.
#'@param observed vector of observed counts. Observed counts must match up with the year and expected vectors.
#'@param Time Number of time periods or years in your dataset. Must be declared as numeric.
#'@param nsim Number of simulations you would like to run
#'@param center can be a single center or for multiple clusters, concatenate them. Max three TODO extend this
#'@param radius radius for the cluster you want in the simulation
#'@param risk.ratio setting for what the risk ratio should be in the cluster to be detected by the simulation
#'@param timeperiod time period where the cluster  starts in the simulation
#'we will only be looking at periods 2 and 5. If multi_period is TRUE, then we will instead consider timeperiod through period_end (timeperiod:period_end). Following the same example,
#'this would mean we look at periods 2, 3, 4, and 5.
#'@return
#'@details Optional functions include:
#'- 1) utm - default is FALSE. If you have utm coordinates, you want to change this to TRUE.
#'@export
#'TODO allow user to change theta parameter in simulation
clust.sim.all <- function(x, y, rMax, period, expected, observed, Time, nsim, center, radius, risk.ratio, 
                          timeperiod,colors=NULL,utm=TRUE, byrow=TRUE,threshold, space=FALSE,...){
    #initial user setting
    if(utm==FALSE){
        message("Coordinates are assumed to be in lat/long coordinates. For utm coordinates, please specify 'utm=TRUE'")
        utm=FALSE
    }
    else{
        utm=TRUE
    }
    if(byrow==FALSE){
        row=FALSE
    }
    else{
        row=TRUE
        message("Data assumed to be in panel data. To use vector data instead, please specify 'byrow=FALSE'")
    }
    if(space==TRUE){
        timeperiod = 1
        Time = 1
        message("Running Space-Only Model")
    }
    #set up clusters and fitted values
    clusters <- clusters.df(x,y,rMax, utm=TRUE, length(x))
    n <- length(x)
    init <- set.vectors(period, expected, observed, Time=Time, byrow=row)
    
    #TODO change this to be a function and not hard-coded
    if(length(center) == 2){
        tmp <- clusters[clusters$center==center[1] | clusters$center==center[2],]
    }
    else if(length(center) == 3){
        tmp <- clusters[clusters$center==center[1] | clusters$center==center[2] | clusters$center==center[3],]
    }
    else{
        tmp <- clusters[clusters$center==center,]
    }
    cluster <- tmp[(tmp$r <= radius),]
    rr = matrix(1, nrow=n, ncol=Time)
    #TODO change this to be a function and not hard-coded
    if(length(timeperiod) == 3){
        rr[cluster$last, timeperiod[1]] = risk.ratio
        rr[cluster$last, timeperiod[2]] = risk.ratio
        rr[cluster$last, timeperiod[3]] = risk.ratio
        message(paste("Running model for periods",timeperiod[1],"through", timeperiod[3]))
    }
    if(length(timeperiod) == 4){
        rr[cluster$last, timeperiod[1]] = risk.ratio
        rr[cluster$last, timeperiod[2]] = risk.ratio
        rr[cluster$last, timeperiod[3]] = risk.ratio
        rr[cluster$last, timeperiod[4]] = risk.ratio
        message(paste("Running model for periods",timeperiod[1],"through", timeperiod[4]))
    }
    if(length(timeperiod) == 5){
        rr[cluster$last, timeperiod[1]] = risk.ratio
        rr[cluster$last, timeperiod[2]] = risk.ratio
        rr[cluster$last, timeperiod[3]] = risk.ratio
        rr[cluster$last, timeperiod[4]] = risk.ratio
        rr[cluster$last, timeperiod[5]] = risk.ratio
        message(paste("Running model for periods",timeperiod[1],"through", timeperiod[5]))
    }
    if(length(timeperiod) == 2){
        rr[cluster$last, timeperiod[1]] = risk.ratio
        rr[cluster$last, timeperiod[2]] = risk.ratio
        message(paste("Running model for periods",timeperiod[1],"and", timeperiod[2]))
    }
    else if(length(timeperiod)==1){
        rr[cluster$last, timeperiod:Time] = risk.ratio    
    }
    E1 = as.vector(rr)*init$E0
    Period <- init$Year
    #Simulate observed as NB(Eit, theta)
    theta = 1000
    YSIM <- lapply(1:nsim, function(i) rnegbin(E1, theta = theta))
    #Scale YSIM[i] to init$E0
    Ex <- scale.sim(YSIM, init, nsim, Time)
    vectors.sim <- list(Period = Period, Ex = Ex , E0_0 = init$E0, Y.vec=init$Y.vec)
    
    #SPACE-ONLY MODELS
    if(space==TRUE){
        #set up and run simulation models
        lassoresult.qp.st <- spacetime.lasso.sim(clusters, vectors.sim, Time, spacetime=FALSE, pois=FALSE, nsim, YSIM)
        lassoresult.p.st <- spacetime.lasso.sim(clusters, vectors.sim, Time, spacetime=FALSE, pois=TRUE, nsim, YSIM)
        
    }
    #SPACE-TIME MODELS = DEFAULT
    else{
        #set up and run simulation models
        lassoresult.qp.st <- spacetime.lasso.sim(clusters, vectors.sim, Time, spacetime=TRUE, pois=FALSE, nsim, YSIM)
        lassoresult.p.st <- spacetime.lasso.sim(clusters, vectors.sim, Time, spacetime=TRUE, pois=TRUE, nsim, YSIM)
        
    }
    
    #RR and Colors for Plotting
    riskratios.qp.st <- get.rr2(lassoresult.qp.st, vectors.sim,init, E1,Time, sim=TRUE)
    rrcolors.qp.st <- colormapping(riskratios.qp.st,Time)
    
    riskratios.p.st <- get.rr2(lassoresult.qp.st, vectors.sim,init, E1,Time, sim=TRUE)
    rrcolors.p.st <- colormapping(riskratios.p.st,Time)
    
    riskratios <- list(riskratios.qp.st = riskratios.qp.st, riskratios.p.st = riskratios.p.st)
    rrcolors <- list(rrcolors.qp.st = rrcolors.qp.st, rrcolors.p.st = rrcolors.p.st)
    
    #Detection
    ##QP - ST
    set <- detect.set(lassoresult.qp.st, vectors.sim, rr, Time, x, y, rMax, center, radius)
    incluster.qp.st <- detect.incluster(lassoresult.qp.st, vectors.sim, rr, set, timeperiod, Time, nsim, x, y, rMax, center, 
                                  radius, IC = "ic")
    detect.qp.st <- list(clust.diagnostics(incluster.qp.st , threshold[1]), clust.diagnostics(incluster.qp.st , threshold[2]))
    detect.out.qp.st <- (matrix(unlist(detect.qp.st),ncol=3, byrow=TRUE, 
            dimnames = list(c(paste0("incluster.any.", threshold[1]),
                              paste0("alldetect.",threshold[1]), 
                              paste0("potentialclusterdetect.",threshold[1]), 
                              paste0("trueclusterdetect.",threshold[1]),
                              paste0("incluster.any.",threshold[2]), paste0("alldetect.",threshold[2]),
                              paste0("potentialclusterdetect.",threshold[2]), 
                              paste0("trueclusterdetect.",threshold[2])),c("aic","aicc","bic"))))
    
    ##P - ST
    set <- detect.set(lassoresult.p.st, vectors.sim, rr, Time, x, y, rMax, center, radius)
    incluster.p.st <- detect.incluster(lassoresult.p.st, vectors.sim, rr, set, timeperiod, Time, nsim, x, y, rMax, center, 
                                  radius, IC = "ic")
    detect.p.st <- list(clust.diagnostics(incluster.p.st, threshold[1]), clust.diagnostics(incluster.p.st, threshold[2]))
    detect.out.p.st <- (matrix(unlist(detect.p.st),ncol=3, byrow=TRUE, 
                                dimnames = list(c(paste0("incluster.any.", threshold[1]),
                                                  paste0("alldetect.",threshold[1]), 
                                                  paste0("potentialclusterdetect.",threshold[1]), 
                                                  paste0("trueclusterdetect.",threshold[1]),
                                                  paste0("incluster.any.",threshold[2]), paste0("alldetect.",threshold[2]),
                                                  paste0("potentialclusterdetect.",threshold[2]), 
                                                  paste0("trueclusterdetect.",threshold[2])),c("aic","aicc","bic"))))
   
    return(list(lassoresult.qp.st = lassoresult.qp.st,
                lassoresult.p.st = lassoresult.p.st,
                riskratios = riskratios,
                rrcolors = rrcolors,
                rr.mat = rr,
                init.vec = vectors.sim,
                incluster.qp.st = incluster.qp.st,
                incluster.p.st = incluster.p.st,
                detect.qp.st = detect.qp.st,
                detect.p.st = detect.p.st,
                detect.out.p.st = detect.out.p.st,
                detect.out.qp.st = detect.out.qp.st))
}






#'findneighbors
#'uses global parameters. An object of class 'nb' must be passed here. Consider creating a neighborhood class via poly2nb.
#'
findneighbors <- function(nb, x, y, rMax,center, radius){
    clusters <- clusters.df(x,y,rMax, utm=TRUE, length(x))
    tmp <- clusters[clusters$center==center,]
    cluster <- tmp[(tmp$r <= radius),]
    
    last <- sort(cluster$last)
    myneigh <- NULL
    for(i in last){
        myneigh <- c(myneigh, unlist(nb[[i]]))
    }
    nbs <- setdiff(unique(myneigh),last)
    return(list(cluster = last, nbs = nbs))
}

scale <- function(init,Time,scaler=NULL,...){
    if(!is.null(scaler)){
        print(scaler)
        #scalenum == as.integer(scaler)
        #print(scalenum)
        std <- sapply(1:Time, function(i) 
            (((matrix(init$E0, ncol=Time)[,i])/scaler)*(sum(matrix(init$Y.vec, ncol=Time)[,i])/scaler))/(sum(matrix(init$E0, ncol=Time)[,i])/scaler))
        E0 <- as.vector(std)
    }
    else{
        std <- sapply(1:Time, function(i) 
            ((((matrix(init$E0, ncol=Time)[,i]))*(sum(matrix(init$Y.vec, ncol=Time)[,i])))/(sum(matrix(init$E0, ncol=Time)[,i]))))
        E0 <- as.vector(std)
    }
    return(E0)
}

scale.sim <- function(YSIM, init, nsim,Time,...){
    std <- lapply(1:nsim, function(i) sapply(1:Time, function(j) 
        (matrix(init$E0,ncol=Time)[,j])*(sum(matrix(YSIM[[i]],ncol=Time)[,j])/sum(matrix(init$E0,ncol=Time)[,j]))))
    E0 <- lapply(1:nsim, function(i) as.vector(std[[i]])) 
    return(E0)
}

get.rr2 <- function(lassoresult,vectors.sim,init,E1, Time, sim=TRUE,...){
    if(sim==TRUE){
        RRobs <- matrix(as.vector(E1)/as.vector(init$E0),ncol=Time)
    }
    if(sim==FALSE){
        RRobs <- matrix(as.vector(vectors$Y.vec)/as.vector(vectors$E0),ncol=Time)
    }
    return(list(RRbic=matrix(lassoresult$E.qbic,ncol=Time),
                RRaic=matrix(lassoresult$E.qaic,ncol=Time),
                RRaicc=matrix(lassoresult$E.qaicc,ncol=Time),
                RRobs= RRobs))
}

#PLOTTING
easyplot.st <- function(pdfname, res, mods, space=FALSE){
    if(space==TRUE){
        #QPois
        pdf_qp <- paste0(gsub(".pdf","", pdfname),mods[1], ".pdf")
        plotmap.s(pdf_qp, res, sub = res$rrcolors$rrcolors.qp.st)
        
        #Pois
        pdf_p <- paste0(gsub(".pdf","", pdfname),mods[2], ".pdf")
        plotmap.s(pdf_p, res, sub = res$rrcolors$rrcolors.p.st)
    }
    else{
        #QPois
        pdf_qp <- paste0(gsub(".pdf","", pdfname),mods[1], ".pdf")
        plotmap.st(pdf_qp, res, sub = res$rrcolors$rrcolors.qp.st)
        
        #Pois
        pdf_p <- paste0(gsub(".pdf","", pdfname),mods[2], ".pdf")
        plotmap.st(pdf_p, res, sub = res$rrcolors$rrcolors.p.st)

    }
    
}



plotmap.st <- function(pdfname,res, obs = NULL, sub=NULL){
    if(!is.null(obs)){
        firstrow = "Obs"
    }
    else{
        firstrow="Oracle"
    }
    if(!is.null(sub)){
        rrcolors <- sub
        print("ok")
    }
    else {
        rrcolors <- res$rrcolors
    }
    pdf(pdfname, height=11, width=10)
    #Maps of Observed Counts
    par(fig=c(0,.2,.6,1), mar=c(.5,0.5,0.5,0))
    plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
    #polygon(japan.poly2,col=res$rrcolors$colors.obs[,1],border=F)
    polygon(japan.poly2,col=rrcolors$colors.obs[,1],border=F)
    segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
    text(355,4120,paste0("Period 1 - ", firstrow),cex=1.00)
    
    par(fig=c(0.2,.4,.6,1), mar=c(.5,0.5,0.5,0), new=T)
    plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(japan.poly2,col=rrcolors$colors.obs[,2],border=F)
    segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
    text(355,4120,paste0("Period 2 - ", firstrow),cex=1.00)
    
    par(fig=c(0.4,.6,.6,1), mar=c(.5,0.5,0.5,0), new=T)
    plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(japan.poly2,col=rrcolors$colors.obs[,3],border=F)
    segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
    text(355,4120,paste0("Period 3 - ", firstrow),cex=1.00)
    
    par(fig=c(0.6,.8,.6,1), mar=c(.5,0.5,0.5,0), new=T)
    plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(japan.poly2,col=rrcolors$colors.obs[,4],border=F)
    segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
    text(355,4120,paste0("Period 4 - ", firstrow),cex=1.00)
    
    par(fig=c(0.8,1,.6,1), mar=c(.5,0.5,0.5,0), new=T)
    plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(japan.poly2,col=rrcolors$colors.obs[,5],border=F)
    segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
    text(355,4120,paste0("Period 5 - ", firstrow),cex=1.00)
    
    
    #Maps of AIC Path
    
    par(fig=c(0,.2,.4,.8), mar=c(.5,0.5,0.5,0), new=T)
    plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(japan.poly2,col=rrcolors$color.qaic[,1],border=F)
    segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
    text(355,4120,'Period 1 - QAIC',cex=1.00)
    
    par(fig=c(0.2,.4,.4,.8), mar=c(.5,0.5,0.5,0), new=T)   
    plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(japan.poly2,col=rrcolors$color.qaic[,2],border=F)
    segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
    text(355,4120,'Period 2 - QAIC',cex=1.00)
    
    par(fig=c(0.4,.6,.4,.8), mar=c(.5,0.5,0.5,0), new=T) 
    plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(japan.poly2,col=rrcolors$color.qaic[,3],border=F)
    segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
    text(355,4120,'Period 3 - QAIC',cex=1.00)
    
    par(fig=c(0.6,.8,.4,.8), mar=c(.5,0.5,0.5,0), new=T)
    plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(japan.poly2,col=rrcolors$color.qaic[,4],border=F)
    segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
    text(355,4120,'Period 4 - QAIC',cex=1.00)
    
    par(fig=c(0.8,1,.4,.8), mar=c(.5,0.5,0.5,0), new=T)
    plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(japan.poly2,col=rrcolors$color.qaic[,5],border=F)
    segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
    text(355,4120,'Period 5 - QAIC',cex=1.00)
    
    
    #Maps of AICc Path
    
    par(fig=c(0,.2,.2,.6), mar=c(.5,0.5,0.5,0), new=T)
    plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(japan.poly2,col=rrcolors$color.qaicc[,1],border=F)
    segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
    text(355,4120,'Period 1 - QAICc',cex=1.00)
    
    par(fig=c(0.2,.4,.2,.6), mar=c(.5,0.5,0.5,0), new=T) 
    plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(japan.poly2,col=rrcolors$color.qaicc[,2],border=F)
    segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
    text(355,4120,'Period 2 - QAICc',cex=1.00)
    
    par(fig=c(0.4,.6,.2,.6), mar=c(.5,0.5,0.5,0), new=T) 
    plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(japan.poly2,col=rrcolors$color.qaicc[,3],border=F)
    segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
    text(355,4120,'Period 3 - QAICc',cex=1.00)
    
    par(fig=c(0.6,.8,.2,.6), mar=c(.5,0.5,0.5,0), new=T)
    plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(japan.poly2,col=rrcolors$color.qaicc[,4],border=F)
    segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
    text(355,4120,'Period 4 - QAICc',cex=1.00)
    
    par(fig=c(0.8,1,.2,.6), mar=c(.5,0.5,0.5,0), new=T)
    plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(japan.poly2,col=rrcolors$color.qaicc[,5],border=F)
    segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
    text(355,4120,'Period 5 - QAICc',cex=1.00)
    
    
    #Maps of BIC Path
    
    par(fig=c(0,.2,0,.4), mar=c(.5,0.5,0.5,0), new=T)
    plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(japan.poly2,col=rrcolors$color.qbic[,1],border=F)
    segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
    text(355,4120,'Period 1 - QBIC',cex=1.00)
    
    par(fig=c(0.2,.4,0,.4), mar=c(.5,0.5,0.5,0), new=T)
    plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(japan.poly2,col=rrcolors$color.qbic[,2],border=F)
    segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
    text(355,4120,'Period 2 - QBIC',cex=1.00)
    
    par(fig=c(0.4,.6,0,.4), mar=c(.5,0.5,0.5,0), new=T)
    plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(japan.poly2,col=rrcolors$color.qbic[,3],border=F)
    segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
    text(355,4120,'Period 3 - QBIC',cex=1.00)
    
    
    par(fig=c(0.6,.8,0,.4), mar=c(.5,0.5,0.5,0), new=T)
    plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(japan.poly2,col=rrcolors$color.qbic[,4],border=F)
    segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
    text(355,4120,'Period 4 - QBIC',cex=1.00)
    
    
    par(fig=c(0.8,1,0,.4), mar=c(.5,0.5,0.5,0), new=T)
    plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(japan.poly2,col=rrcolors$color.qbic[,5],border=F)
    segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
    text(355,4120,'Period 5 - QBIC',cex=1.00)
    
    #Turn off pdf development
    dev.off()
}

plotmap.s <- function(pdfname,res, obs = NULL, sub=NULL){
    if(!is.null(obs)){
        firstrow = "Obs"
    }
    else{
        firstrow="Oracle"
    }
    if(!is.null(sub)){
        rrcolors <-  sub
        print("ok -spaceonly")
    }
    else{
        rrcolors <- res$rrcolors
    }
    pdf(pdfname, height=11, width=10)
    #Maps of Observed Counts
    par(fig=c(0,.2,.6,1), mar=c(.5,0.5,0.5,0))
    plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(japan.poly2,col=rrcolors$colors.obs[,1],border=F)
    segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
    text(355,4120,paste0(firstrow),cex=1.00)
    
    #Maps of AIC Path
    
    par(fig=c(0,.2,.4,.8), mar=c(.5,0.5,0.5,0), new=T)
    plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(japan.poly2,col=rrcolors$color.qaic[,1],border=F)
    segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
    text(355,4120,'QAIC',cex=1.00)
    
    #Maps of AICc Path
    
    par(fig=c(0,.2,.2,.6), mar=c(.5,0.5,0.5,0), new=T)
    plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(japan.poly2,col=rrcolors$color.qaicc[,1],border=F)
    segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
    text(355,4120,'QAICc',cex=1.00)
    
    
    #Maps of BIC Path
    
    par(fig=c(0,.2,0,.4), mar=c(.5,0.5,0.5,0), new=T)
    plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(japan.poly2,col=rrcolors$color.qbic[,1],border=F)
    segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
    text(355,4120,'QBIC',cex=1.00)
    
    
    #Turn off pdf development
    dev.off()
}


#Wrapper for plot.st





clust.diagnostics <- function(incluster, threshold){
    #thresholding of prop.alldetect
    alldetect.aic <- paste0((length(which(unlist(incluster$prop.alldetect.aic)> threshold))/nsim)*100, "%")
    alldetect.aicc <- paste0((length(which(unlist(incluster$prop.alldetect.aicc)> threshold))/nsim)*100, "%")
    alldetect.bic <- paste0((length(which(unlist(incluster$prop.alldetect.bic)> threshold))/nsim)*100, "%")
    
    #thresholding of prop.wasdetect
    potentialclusterdetect.aic <- paste0((length(which(unlist(incluster$prop.wasdetect.aic)> threshold))/nsim)*100, "%")
    potentialclusterdetect.aicc <- paste0((length(which(unlist(incluster$prop.wasdetect.aicc)> threshold))/nsim)*100, "%")
    potentialclusterdetect.bic <- paste0((length(which(unlist(incluster$prop.wasdetect.bic)> threshold))/nsim)*100, "%")
    
    #thresholding of prop.shoulddetect
    trueclusterdetect.aic <- paste0((length(which(unlist(incluster$prop.shoulddetect.aic)> threshold))/nsim)*100, "%")
    trueclusterdetect.aicc <- paste0((length(which(unlist(incluster$prop.shoulddetect.aicc)> threshold))/nsim)*100, "%")
    trueclusterdetect.bic <- paste0((length(which(unlist(incluster$prop.shoulddetect.bic)> threshold))/nsim)*100, "%")
    
    return(list(incluster.any.aic = incluster$incluster.any.aic, incluster.any.aicc = incluster$incluster.any.aicc,
           incluster.any.bic = incluster$incluster.any.bic,
           alldetect.aic = alldetect.aic, alldetect.aicc = alldetect.aicc, alldetect.bic = alldetect.bic,
           potentialclusterdetect.aic = potentialclusterdetect.aic, potentialclusterdetect.aicc = potentialclusterdetect.aicc,
           potentialclusterdetect.bic = potentialclusterdetect.bic,
           trueclusterdetect.aic = trueclusterdetect.aic, trueclusterdetect.aicc = trueclusterdetect.aicc,
           trueclusterdetect.bic = trueclusterdetect.bic
    ))
}

