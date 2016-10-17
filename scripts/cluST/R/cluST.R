#' cluST.R
#' @name cluST
#' @title Spatial and Spatial-Temporal Clustering using Lasso
#' @description This package will host the functions and code for a spatial and spatial-temporal clustering model using Lasso regression. We assume a Poisson distribution for counts. The model includes a spatially-structured random effect to account for heterogeneity in a Poisson model. Model selection is based on QAIC/QAICc/QBIC.
#' @details This package contains C++ code which is compiled at run-time. These functions are located in the "src" file.
#' Dependencies and inputs are documented here.
#' Data used for this analysis comes from MidWest Poverty data.
#' @references Xu, Jiale, Gangnon, Ronald: "Stagewise and Stepwise Methods for Space and Space-Time Cluster Detection"
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
#' clustersDF(x1,y1,rMax, utm=TRUE, length(x1))
#' clustersDF(lat, long, utm=FALSE, length(lat))

clustersDF <- function(xP,yP, r.max, utm=FALSE,n){
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

dpoisson <- function(y, lambda, log = FALSE) {
    if(log == FALSE) 
        return(lambda^y * exp(-lambda)/factorial(y))
    else
        return(y*ifelse(lambda==0,1,log(lambda))-lambda)
}

###DO I need to document cpp functions?

#' Creates a List Arranged by Time Period with Expected and Observed Counts and Time Period
#' 
#' @param year vector of periods or years in dataset. Should be imported as a factor.
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
setVectors <- function(year, expect, observed,Time, byrow=TRUE) {
    if (byrow==TRUE){
        E0=as.vector(matrix(expect, byrow=T, ncol=Time))
        Y.vec <- as.vector(matrix(observed,byrow=T, ncol=Time))
        Year <- as.vector(matrix(year, byrow=T, ncol=Time)) 
    }
    else {
        E0=as.vector(matrix(expect, ncol=Time))
        Y.vec <- as.vector(matrix(observed, ncol=Time))
        Year <- as.vector(matrix(year, ncol=Time))
    }
    return(list(
        E0 = E0,
        Y.vec = Y.vec,
        Year = Year))
}


#' spaceMat
#' 
#' This function creates a sparse matrix of 1's of all potential clusters for the Lasso algorithm to cycle over; this incorporates space
#' @param clusters clusters dataframe from (clustersDF function) that includes the center, x,y, r (radius), n (counter), and last (last observation in potential cluster)
#' @param numCenters the number of centers
#' @return returns sparse matrix of 1's
spaceMat <- function(clusters, numCenters){
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

#' timeMat
#' 
#' This function creates a sparse matrix of 1's of all of the potential time periods for the cluster to be in. The number of 
#' potential time periods is determined by [(Time*(Time-1)]/2.
#' @param Time Number of time periods in the data
#' @return returns sparse matrix of 1's as indicators of membership in the time cluster
timeMat <-function(Time){
    block <- Matrix(diag(1,Time),sparse=TRUE)
    master <- block
    for(i in 1:(Time-2)){
        diag(block[(i+1):Time,])<-1
        master <- cBind(master, block[,1:(Time-i)])        
    }
    master <- cBind(master, Matrix(rep(1,Time)))
    return(master)
}
    
    
#' spaceTimeMat
#' 
#' This function takes the Kronecker product of the space and time matrices to create the space-time matrix
#' @param timeMat Time matrix
#' @param spaceMat Space matrix
#' @return Returns sparse space time matrix          
spaceTimeMat <- function(clusters, numCenters, Time){
    space <- spaceMat(clusters, numCenters)
    time <- timeMat(Time)
    spaceTimeMat <- kronecker(time, space)
    return(spaceTimeMat)
}

#' myoverdisp
#' 
#' This function calculates the overdispersion parameter for the QIC 'c' overdispersion parameter.
#' @param 
#' @return returns sparse matrix of 1's
myoverdisp <- function(object) {
    with(object,sum((weights * residuals^2)[weights > 0])/df.residual)
}

#' spacetimeLasso
#' 
#' This function runs the Lasso regularization technique on our large sparse matric of potential space-time clusters.
#' @param potClus number of potential clusters. This will usually be the same as 'numCenters'
#' @param clusters clusters dataframe from (clustersDF function) that includes the center, x,y, r (radius), n (counter), and last (last observation in potential cluster)
#' @param numCenters the number of centers
#' @param vectors takes in the list of expected and observed counts from setVectors function
#' @param Time number of time periods in the dataset
#' @param spacetime indicator of whether the cluster detection method should be run on all space-time clusters(default) or on only the potential space clusters.
#' @return This function will return a list with the expected counts as selected by QBIC, QAIC, QAICc, a list of original expected counts (Ex),
#' a list of observed counts (Yx), the lasso object, a list of K values (number of unique values in each decision path), and n (length of unique centers in the clusters dataframe)
#' @export
spacetimeLasso <- function(potClus, clusters, numCenters, vectors, Time, spacetime=TRUE){
    n <- length(unique(clusters$center))
    potClus <- n
    numCenters <- n
    if(spacetime==TRUE){
        sparseMAT <- spaceTimeMat(clusters, numCenters, Time)
    }
    else{
        sparseMAT <- spaceMat(clusters, numCenters)
    }
    Ex <- vectors$E0
    Yx <- vectors$Y.vec
    lasso <- glmnet(sparseMAT, Yx, family=("poisson"), alpha=1, offset=log(Ex))
    coefs.lasso.all <- coef(lasso)
    intercept <- rep(1, dim(sparseMAT)[1])
    sparseMAT <- cBind(intercept, sparseMAT)
    xbetaPath<- sparseMAT%*%coefs.lasso.all
    mu <- sapply(1:length(lasso$lambda), function(i) Ex * exp(xbetaPath[,i]))    
    loglike <- sapply(1:length(lasso$lambda), function(i) sum(dpoisson(Yx, mu[,i],log=TRUE)))
    K <- sapply(1:length(lasso$lambda), function(i) length(unique(xbetaPath[,i])))
    
    offset_reg <- glm(Yx ~ offset(log(Ex)),family=poisson)
    overdisp <- myoverdisp(offset_reg)
    
    if(spacetime==TRUE){
        #QBIC
        PLL.qbic  <- (loglike/overdisp)-log(n*Time)/2*K
        qbicMax <- which.max(PLL.qbic)
        E.qbic <- mu[,qbicMax]
        
        #QAIC
        PLL.qaic = (loglike/overdisp) - K
        qaicMax <- which.max(PLL.qaic)
        E.qaic <- mu[,qaicMax]
        
        #QAICc
        PLL.qaicc=(loglike/overdisp)- ((K*n*Time)/(n*Time-K-1))
        qaiccMax <- which.max(PLL.qaicc)
        E.qaicc <- mu[,qaiccMax]
    }
    else{
        #BIC
        PLL.qbic  <- (loglike)-log(n*Time)/2*K
        qbicMax <- which.max(PLL.qbic)
        E.qbic <- mu[,qbicMax]
        
        #AIC
        PLL.qaic = (loglike) - K
        qaicMax <- which.max(PLL.qaic)
        E.qaic <- mu[,qaicMax]
        
        #AICc
        PLL.qaicc=(loglike)- ((K*n*Time)/(n*Time-K-1))
        qaiccMax <- which.max(PLL.qaicc)
        E.qaicc <- mu[,qaiccMax]
    }
    return(list(E.qbic, E.qaic, E.qaicc,Ex, Yx, lasso, K,n))    
}


#' spacetimeLasso.sim
#' 
#' This function runs the Lasso regularization technique on our large sparse matric of potential space-time clusters.It is specifically created to use with simulations.
#' @param potClus number of potential clusters. This will usually be the same as 'numCenters'
#' @param clusters clusters dataframe from (clustersDF function) that includes the center, x,y, r (radius), n (counter), and last (last observation in potential cluster)
#' @param numCenters the number of centers
#' @param vectors takes in the list of expected and observed counts from setVectors function
#' @param Time number of time periods in the dataset
#' @param spacetime indicator of whether the cluster detection method should be run on all space-time clusters(default) or on only the potential space clusters.
#' @return This function will return a list with the expected counts as selected by QBIC, QAIC, QAICc, a list of original expected counts (Ex),
#' a list of observed counts (Yx), the lasso object, a list of K values (number of unique values in each decision path), and n (length of unique centers in the clusters dataframe)
#' @export
spacetimeLasso.sim <- function(potClus, clusters, numCenters, vectors, Time, spacetime=TRUE,nsim,YSIM){
    n <- length(unique(clusters$center))
    potClus <- n
    numCenters <- n
    if(spacetime==TRUE){
        sparseMAT <- spaceTimeMat(clusters, numCenters, Time)
    }
    else{
        sparseMAT <- spaceMat(clusters, numCenters)
    }
    Ex <- vectors$E0
    Yx <- YSIM
    lasso <- lapply(1:nsim, function(i) glmnet(sparseMAT, Yx[,i], family=("poisson"), alpha=1, offset=log(Ex)))
    coefs.lasso.all <- lapply(1:nsim, function(i) coef(lasso[[i]]))
    intercept <- rep(1, dim(sparseMAT)[1])
    sparseMAT <- cBind(intercept, sparseMAT)
    xbetaPath<- lapply(1:nsim, function(i) sparseMAT%*%coefs.lasso.all[[i]])
    mu <- lapply(1:nsim, function(j) sapply(1:length(lasso[[j]]$lambda), function(i) Ex * exp(xbetaPath[[j]][,i])))    
    loglike <- lapply(1:nsim, function(k) sapply(1:length(lasso[[k]]$lambda), 
                                                 function(i) sum(dpoisson(Yx[,k], mu[[k]][,i],log=TRUE))))
    K <- lapply(1:nsim, function(j) sapply(1:length(lasso[[j]]$lambda), function(i) length(unique(xbetaPath[[j]][,i]))))
    offset_reg <- lapply(1:nsim, function(i) glm(Yx[,i] ~ offset(log(Ex)),family=poisson))
    overdisp <- lapply(1:nsim, function(i) myoverdisp(offset_reg[[i]]))
    mean_mu <- sapply(1:nsim, function(j) sapply(1:nrow(YSIM), function(i) mean(mu[[1]][i,j]:mu[[nsim]][i,j])))
    
    if(spacetime==TRUE){
        #QBIC
        PLL.qbic  <- lapply(1:nsim, function(i) (loglike[[i]]/overdisp[[i]])-log(n*Time)/2*K[[i]])
        mean_qbic <- which.max(unlist(lapply(1:nsim, function(i) mean(PLL.qbic[[1]][i]:PLL.qbic[[nsim]][i]))))
        E.qbic <- mean_mu[,mean_qbic]
        
        #QAIC
        PLL.qaic = lapply(1:nsim, function(i) (loglike[[i]]/overdisp[[i]]) - K[[i]])
        mean_qaic <- which.max(unlist(lapply(1:nsim, function(i) mean(PLL.qaic[[1]][i]:PLL.qaic[[nsim]][i]))))
        E.qaic <- mean_mu[,mean_qaic]
        
        #QAICc
        PLL.qaicc=lapply(1:nsim, function(i) (loglike[[i]]/overdisp[[i]])- ((K[[i]]*n*Time)/(n*Time-K[[i]]-1)))
        mean_qaicc <- which.max(unlist(lapply(1:nsim, function(i) mean(PLL.qaicc[[1]][i]:PLL.qaicc[[nsim]][i]))))
        E.qaicc <- mean_mu[,mean_qaicc]
    }
    else{
        #BIC
        PLL.qbic  <- lapply(1:nsim, function(i) (loglike[[i]])-log(n*Time)/2*K[[i]])
        mean_qbic <- which.max(unlist(lapply(1:nsim, function(i) mean(PLL.qbic[[1]][i]:PLL.qbic[[nsim]][i]))))
        E.qbic <- mean_mu[,mean_qbic]
        
        #AIC
        PLL.qaic = lapply(1:nsim, function(i) (loglike[[i]]) - K[[i]])
        mean_qaic <- which.max(unlist(lapply(1:nsim, function(i) mean(PLL.qaic[[1]][i]:PLL.qaic[[nsim]][i]))))
        E.qaic <- mean_mu[,mean_qaic]
        
        #AICc
        PLL.qaicc=lapply(1:nsim, function(i) (loglike[[i]])- ((K[[i]]*n*Time)/(n*Time-K[[i]]-1)))
        mean_qaicc <- which.max(unlist(lapply(1:nsim, function(i) mean(PLL.qaicc[[1]][i]:PLL.qaicc[[nsim]][i]))))
        E.qaicc <- mean_mu[,mean_qaicc]
    }
    mean_Yx <- round(rowSums(Yx)/nsim)
    return(list(E.qbic, E.qaic, E.qaicc,Ex, mean_Yx, Yx, lasso, K,n))    
}



#' setRR
#' 
#' This function will create vectors of the risk ratios as determined by observed counts, QBIC, QAIC, and QAICc, respectively.
#' @param lassoresult List of QBIC, QAIC, QAICc estimates from the mylasso function
#' @param vectors dataframe of initial vectors of the observed and expected counts
#' @param Time number of time period
#' @return This returns a list of the risk ratios (observed over expected) as determined by 1) pure observed/expected counts,
#' 2) observed based on QBIC path/expected; 3) observed based on QAIC path/expected; 4) observed based on QAICc path/expected.
#' @export
setRR <- function(lassoresult, vectors, Time, YSIM){
    if(missing(YSIM)){
        RRobs <- matrix(as.vector(vectors$Y.vec)/as.vector(vectors$E0),ncol=Time)
        RRbic <- matrix(lassoresult[[1]]/as.vector(vectors$E0),ncol=Time)
        RRaic <- matrix(lassoresult[[2]]/as.vector(vectors$E0),ncol=Time)
        RRaicc <- matrix(lassoresult[[3]]/as.vector(vectors$E0),ncol=Time)
    }
    else{
        RRobs <- matrix(as.vector(lassoresult[[5]])/as.vector(vectors$E0),ncol=Time)
        RRbic <- matrix(lassoresult[[1]]/as.vector(vectors$E0),ncol=Time)
        RRaic <- matrix(lassoresult[[2]]/as.vector(vectors$E0),ncol=Time)
        RRaicc <- matrix(lassoresult[[3]]/as.vector(vectors$E0),ncol=Time)    
    }
    return(list(RRobs=RRobs, RRbic=RRbic, RRaic=RRaic, RRaicc=RRaicc))  
}

#' redblue
#' 
#' This function establishes the spread of reds and blues for the risk ratios to be mapped to. Higher risk ratios will be deeper red colors and lower risk ratios will be deeper blue colors.
#' @param x this will be the risk ratios shrunk to be on the scale of half risk to twice the risk as end points.
#' @return colors
redblue=function(x) { 
    y=colorRamp(brewer.pal(11,"RdBu")[11:1])(x); rgb(y[,1],y[,2],y[,3],max=255) 
}


#' redblue
#' 
#' This function establishes the spread of reds and blues for the risk ratios to be mapped to. Higher risk ratios will be deeper red colors and lower risk ratios will be deeper blue colors.
#' @param x this will be the risk ratios shrunk to be on the scale of half risk to twice the risk as end points.
#' @return returns vectors ofcolors for each time period, where risk ratios have been constrained to be between half risk and twice the risk
#' @export
colormapping <- function(riskratios,Time) {
    color.obs <- sapply(1:Time, function(i) redblue(log(2*pmax(1/2,pmin(riskratios$RRobs[,i],2)))/log(4)))
    color.qbic <- sapply(1:Time, function(i) redblue(log(2*pmax(1/2,pmin(riskratios$RRbic[,i],2)))/log(4))) 
    color.qaic <- sapply(1:Time, function(i) redblue(log(2*pmax(1/2,pmin(riskratios$RRaic[,i],2)))/log(4)))
    color.qaicc <- sapply(1:Time, function(i) redblue(log(2*pmax(1/2,pmin(riskratios$RRaicc[,i],2)))/log(4)))
    return(list(colors.obs = color.obs, color.qbic = color.qbic, color.qaic = color.qaic, color.qaicc = color.qaicc)) 
}













