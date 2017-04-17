#'
#'Set initial conditions to enable detection diagnostic functions
#'
#'This function will extract and set the necessary vectors in order to detect what percent of clusters
#'in the simulation were detected (true positives), what percent of clusters were incorrectly in the background (false negatives),
#'and what percent were incorrectly identified as clusters (false positives). 
#'@param lassoresult List of QBIC, QAIC, QAICc estimates from the mylasso.sim function
#'@param vectors.sim  dataframe of initial vectors of the observed and expected counts that went into simulation function
#'@param rr risk ratio matrix that was used in the simulation
#'@param x x coordinates
#'@param y y coordinates
#'@param rMax set maximum radius for a cluster
#'@param center center of the simulation cluster
#'@param radius radius of the simulated cluster
#'@param nullmod indicator for whether or not the simulation being run is the null model or not. Default is not null.
#'@return returns a list of 1) indx_truth, which is a vector of indices where our cluster
#'should be, 2) indx, which contains the index of cluster as determined where the risk ratio in the rr matrix is not 1,
#'3) rr.simBIC, rr.simAIC, rr.simAICc - the risk ratios as determined by BIC, AIC, and AICc (respectively), and 4) alphaBIC, alphaAIC, alphaAICc - the calculated
#'background risk rates as determined by BIC, AIC, and AICc or QBIC, QAIC, QAICc. 
#'@export
#'@example
#'x=dframe2$utmx/1000
#'y=dframe2$utmy/1000
#'rMax=30 
#'Time=5  
#'nsim=100
#'center=50
#'radius= 18
#'timeperiod=5
#'risk.ratio=2
#'threshold = 0.5
#'res <- clust.sim.all(x,y,rMax,dframe$period, dframe$expdeath, dframe$death, Time,nsim,center, radius, risk.ratio, timeperiod, colors=TRUE, utm=TRUE, byrow=TRUE, threshold)
#'set <- detect_set(res$lassoresult, res$init.vec, res$rr.mat, Time=5, x, y, rMax, center, radius)
#'
#'#detect_set for null model
#'set <- detect_set(res$lassoresult, res$init.vec, res$rr.mat, Time=5,x, y, rMax, 0, radius, nullmod=TRUE)

detect_set <- function(lassoresult, vectors.sim, rr, Time, x, y, rMax, center, radius,nullmod,...){
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
        rr.simBIC = rr.simBIC,
        rr.simAIC = rr.simAIC,
        rr.simAICc = rr.simAICc,
        alphaBIC = alphaBIC,
        alphaAIC = alphaAIC,
        alphaAICc = alphaAICc))
}


