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
    print("Creating space-time matrix")
    if(spacetime==TRUE){
        sparseMAT <- spaceTimeMat(clusters, numCenters, Time)
    }
    else{
        sparseMAT <- spaceMat(clusters, numCenters)
    }
    print("Space-time matrix created")
    Ex <- vectors$E0
    Yx <- vectors$Y.vec
    print("Running Lasso - stay tuned")
    lasso <- glmnet(sparseMAT, Yx, family=("poisson"), alpha=1, offset=log(Ex))
    print("Lasso complete - extracting estimates and paths")
    coefs.lasso.all <- coef(lasso)
    intercept <- rep(1, dim(sparseMAT)[1])
    sparseMAT <- cBind(intercept, sparseMAT)
    xbetaPath<- sparseMAT%*%coefs.lasso.all
    mu <- sapply(1:length(lasso$lambda), function(i) Ex * exp(xbetaPath[,i]))    
    loglike <- sapply(1:length(lasso$lambda), function(i) sum(dpoisson(Yx, mu[,i],log=TRUE)))
    K <- sapply(1:length(lasso$lambda), function(i) length(unique(xbetaPath[,i])))
    
    offset_reg <- glm(Yx ~ offset(log(Ex)),family=poisson)
    overdisp <- myoverdisp(offset_reg)
    print("Selecting best paths")
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
    print("Returning results")
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
#' @param nsim number of simulations
#' @param YSIM vector of simulated observed counts
#' @param getRR default is TRUE. Determines how relative risk is calculated. If getRR=TRUE, then relative risk is first calculated in each simulation and averaged over simulations.
#' If getRR is FALSE, then relative risk will come from averaging the expected counts from each of the Lasso results by the number of simulations. This should then
#' by used in conjunction with the setRR function to get the relative risk. TODO integrate this into one flow.
#' @export
spacetimeLasso.sim <- function(potClus, clusters, numCenters, vectors, Time, spacetime=TRUE,nsim,YSIM, getRR=TRUE){
    n <- length(unique(clusters$center))
    potClus <- n
    numCenters <- n
    print("Creating space-time matrix")
    if(spacetime==TRUE){
        sparseMAT <- spaceTimeMat(clusters, numCenters, Time)
    }
    else{
        sparseMAT <- spaceMat(clusters, numCenters)
    }
    print("Space-time matrix created")
    Ex <- vectors[[2]]
    Yx <- YSIM
    Period <- vectors[[1]]
    print("Running Lasso - stay tuned")
    lasso <- lapply(1:nsim, function(i) glmnet(sparseMAT, Yx[,i], family=("poisson"), alpha=1, offset=log(Ex[[i]])))
    print("Lasso complete - extracting estimates and paths")
    coefs.lasso.all <- lapply(1:nsim, function(i) coef(lasso[[i]]))
    intercept <- rep(1, dim(sparseMAT)[1])
    sparseMAT <- cBind(intercept, sparseMAT)
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
    overdisp <- lapply(1:nsim, function(i) myoverdisp(offset_reg[[i]]))
    print("Selecting best paths")
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
            select_muRR.qbic <- lapply(1:nsim, function(i) select_mu.qbic[[i]]/vectors[[3]])
            E.qbic <- Reduce("+", select_muRR.qbic)/nsim   
        }
               
        #QAIC
        PLL.qaic = lapply(1:nsim, function(i) (loglike[[i]]/overdisp[[i]]) - K[[i]])
        select.qaic <- lapply(1:nsim, function(i) which.max(unlist(PLL.qaic[[i]])))
        if(getRR==FALSE){
            select_mu.qaic <- lapply(1:nsim, function(i) sapply(select.qaic[[i]], function(j) mu[[i]][,j]))
            E.qaic <- Reduce("+", select_mu.qaic)/nsim    
        }
        else{
            select_mu.qaic <- lapply(1:nsim, function(i) sapply(select.qaic[[i]], function(j) mu[[i]][,j]*Ex[[i]]))
            select_muRR.qaic <- lapply(1:nsim, function(i) select_mu.qaic[[i]]/vectors[[3]])
            E.qaic <- Reduce("+", select_muRR.qaic)/nsim  
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
            E.qaicc <- Reduce("+", select_muRR.qaicc)/nsim  
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
    print("Returning results")
    return(list(nsim = nsim, E.qbic = E.qbic, E.qaic = E.qaic, E.qaicc = E.qaicc,Ex = Ex,mu = mu, Yx = Yx, PLL.qbic = PLL.qbic, 
                PLL.qaic = PLL.qaic, PLL.qaicc = PLL.qaicc, select.qbic = select.qbic, select.qaic = select.qaic, 
                select.qaicc = select.qaicc, select_mu.qbic = select_mu.qbic, select_mu.qaic = select_mu.qaic, 
                select_mu.qaicc = select_mu.qaicc, xbetaPath = xbetaPath, coefs.lasso.all = coefs.lasso.all))    
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
setRR <- function(lassoresult, vectors, Time, sim=FALSE){
    if(sim==FALSE){
        RRobs <- matrix(as.vector(vectors$Y.vec)/as.vector(vectors$E0),ncol=Time)
        RRbic <- matrix(lassoresult$E.qbic/as.vector(vectors$E0),ncol=Time)
        RRaic <- matrix(lassoresult$E.qaic/as.vector(vectors$E0),ncol=Time)
        RRaicc <- matrix(lassoresult$E.qaicc/as.vector(vectors$E0),ncol=Time)
        print("Relative risks from observed data")
    }
    else{
        E0_avg <- Reduce("+", vectors$E0)/length(vectors$E0)
        RRobs <- matrix(as.vector(E0_avg)/as.vector(vectors$E0_fit),ncol=Time)
        RRbic <- matrix(lassoresult$E.qbic/as.vector(vectors$E0_fit),ncol=Time)
        RRaic <- matrix(lassoresult$E.qaic/as.vector(vectors$E0_fit),ncol=Time)
        RRaicc <- matrix(lassoresult$E.qaicc/as.vector(vectors$E0_fit),ncol=Time) 
        print("Relative risks from simulated data")
    }
    return(list(RRobs=RRobs, RRbic=RRbic, RRaic=RRaic, RRaicc=RRaicc))  
}


#' getRR
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
getRR <- function(lassoresult,vectors, Time, sim=TRUE){
    E0_avg <- Reduce("+", vectors$E0)/length(vectors$E0)
    RRobs <- matrix(as.vector(E0_avg)/as.vector(vectors$E0_fit),ncol=Time)
    print("Relative risk ratios from simulated data - average RR over nsim")
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
redblue=function(x) { 
    y=colorRamp(brewer.pal(11,"RdBu")[11:1])(x); rgb(y[,1],y[,2],y[,3],max=255) 
}


#' colormapping
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

#'scale
#'
#'This function standardizes the fitted values by the sum of the observed/the sum of the expected within each time period. Inputs should all be in vector form.
#'The function handles standardizing within time period as long.
#'@param Y.vec vector of observed (in vector format)
#'@param out.sim simulated data based on the simulate(out) step above #TODO integrate this into this step above and streamline into a single function
#'@param nsim number of simulations performed in out.sim
#'@param Time number of time periods
#'@export
scale <- function(Y.vec, out.sim, nsim,Time){
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
probmap <- function(lassoresult, vectors.sim, rr, nsim, Time, colormap=FALSE){
    prob.sim <- lapply(1:nsim, function(x) matrix(0, nrow(rr)*Time))
    indx <- which(rr !=1)
    rr.sim <- lapply(1:nsim, function(i) lassoresult$select_mu.qbic[[i]]/vectors.sim$E0_fit)
    alpha <- lapply(1:nsim, function(i) lapply(1:Time, function(k) 
        sort(table(matrix(rr.sim[[i]], ncol=Time)[,k]),decreasing=TRUE)[1]))
    for(j in 1:length(prob.sim)){
        for(i in 1:length(indx)){
            if (rr.sim[[j]][indx[i]] >= as.numeric(attributes(alpha[[j]][[1]])[[1]])) {
                prob.sim[[j]][indx[i]] <- 1
            }
            else {
                prob.sim[[j]][indx[i]] <- 0
            }
        }
    }
    myprobs <- Reduce("+", prob.sim)/nsim
    if (colormap==TRUE){
        color.probmap <- sapply(1:Time, function(i) redblue(log(2*pmax(1/2,pmin(matrix(myprobs,ncol=Time)[,i],2)))/log(4)))    
        return(list(probabilities = myprobs, colors = color.probmap))
    }
    else{
        return(probabilities = myprobs)   
    }
}    
