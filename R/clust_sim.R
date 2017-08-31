#'Simulation functions for space and space-time cluster detection with lasso. 
#'

#' vectors_space_sim
#' 
#' This function will collapse a space-time vector onto space only
#' @param x vector coordinates (unique regardless of time period)
#' @param Ex list of simulated and standardized expected counts
#' @param YSIM simulated observed
#' @param Time number of time periods
#' @param timeperiod explicit vector of timeperiods
#' @param init initial list of vectors, inherited from function setVectors.
#' @return returns space-time 
#' 
vectors_space_sim <- function(x,Ex, YSIM,Time, init){
    id <- rep(1:length(x), times = Time)
    if(length(id)!=length(as.vector(Ex[[1]]))) stop("Length of ID var not equal to number of observations")
    vectors.sim.s <- list(Period = rep("1", length(x)),
                          Ex = lapply(1:length(YSIM), function(i) tapply(as.vector(matrix(Ex[[i]], ncol=Time)), id, function(x) sum(x))),
                          E0_0 = tapply(as.vector(matrix(init$E0, ncol=Time)), id, function(x) sum(x)),
                          Y.vec = tapply(as.vector(matrix(init$Y.vec, ncol=Time)), id, function(x) round(sum(x))))
    YSIM.s <- lapply(1:length(YSIM), function(i) tapply(as.vector(matrix(YSIM[[i]], ncol=Time)), id, function(x) round(sum(x))))
    return(list(vectors.sim.s = vectors.sim.s, YSIM.s = YSIM.s))
}
    
    
#'
#'clust_sim
#'
#'This helper function runs both the space and space-time Lasso model simulations for all 4 models simulataneously: Quasi-Poisson vs. Poisson in both space and space-time.
#' This function is to be run on simulated data and all four models are run on the same simulated set. 
#'A separate function (clust.sim) can be used for running simulations on individual models and (clust) can be used for observed data.
#'@param clst list; output from toclust function. Must be of class clst.
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
#'@param utm default is true
#'@param byrow TRUE
#'@param threshold vector or value as threshold for cluster detection
#'@param space space and space-time. Default is to run all four models: quasi-poisson and poisson for both space and space-time. User can specify, space = space,
#'space = spacetime, or space = both.
#'@param theta default is 1000. Can add in overdispersion to simulated model by changing this value.
#'@param nullmod default is NULL; otherwise will run null model
#'@param floor floor default is TRUE. When TRUE, it limits phi (overdispersion parameter) to be greater or equal to 1. If FALSE, will allow for under dispersion.
#'@return returns list of lists
#'@export

clust_sim <- function(clst,covars, x,y, rMax, Time, nsim, center, radius, risk.ratio, 
                          timeperiod, utm=TRUE, byrow=TRUE, threshold, space = c("space", "spacetime", "both"), 
                      theta = NULL,nullmod=NULL, floor){
    expected <- clst$required_df$expected
    observed <- clst$required_df$observed
    period <- clst$required_df$timeperiod
    #initial user setting
    if(!is.null(clst$othercovariates_df)){
        covars <- clst$othercovariates_df
    }
    else{
        covars <- NULL
    }
    if(utm==FALSE){
        message("Coordinates are assumed to be in lat/long coordinates. For utm coordinates, please specify 'utm=TRUE'")
        utm=FALSE
    }
    else{
        utm=TRUE
    }
    if(byrow==FALSE){
        byrow=FALSE
    }
    else{
        byrow=TRUE
        message("Data assumed to be in panel data. To use vector data instead, please specify 'byrow=FALSE'")
    }
    if(floor==FALSE){
        floor <- FALSE
    }
    else{
        floor <- TRUE
    }
    if(!is.null(nullmod)){
        if(isTRUE(nullmod)) warning("Null mod has been set to true")
        nullmod=TRUE
        message("Running null models. For simulation models, leave nullmod NULL")
    }
    else{
        nullmod=NULL
    }
    if(!is.null(theta)){
        theta = theta
        message("Running model with user-specified theta value")
    }
    else{
        theta = 1000
        message("Running model with default theta value of 1000")
    }
    if(length(space) > 1) stop("You must select either `space`, `spacetime`, or `both`")
    space <- match.arg(space, several.ok = FALSE)
    switch(space, 
           # space = clustAll_sim.space(x, y, rMax,period, expected, observed, Time, nsim, center, radius, risk.ratio,
           #                             timeperiod,colors=NULL,utm, byrow, threshold, space=TRUE),
           # spacetime = clustAll_sim.spacetime(x, y, rMax,period, expected, observed, Time, nsim, center, radius, risk.ratio,
           #                                     timeperiod,colors=NULL,utm, byrow, threshold, space=FALSE),
           both = clustAll_sim(x, y, rMax,period, expected, observed, covars, Time, nsim, center, radius, risk.ratio,
                                     timeperiod,utm, byrow, threshold, theta, nullmod, floor))
}


#'clustAll_sim
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
#'@param utm utm TRUE/FALSE as to whether or not the x and y coordinates are in UTM (TRUE) or LAT/LONG(FALSE)
#'@param byrow  byrow default is True. If data should be imported by column then set to FALSE
#'@param threshold  vector or value as threshold for cluster detection
#'@param theta default is 1000. Can add in overdispersion to simulated model by changing this value.
#'@param nullmod if TRUE, then null models will be run. Otherwise, default is null.
#'@param floor floor default is TRUE. When TRUE, it limits phi (overdispersion parameter) to be greater or equal to 1. If FALSE, will allow for under dispersion.
#'@return returns list of lists
#'@export

clustAll_sim <- function(x, y, rMax, period, expected, observed, covars,Time, nsim, center, radius, risk.ratio, 
                               timeperiod,utm, byrow, threshold, theta = theta, nullmod=nullmod,floor=floor,...){
    message("Running both Space and Space-Time Models")
    
    #set up clusters and fitted values
    clusters <- clusters2df(x,y,rMax, utm=utm, length(x))
    n <- length(x)
    init <- setVectors(period, expected, observed, Time=Time, byrow=byrow)
    
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
        message(paste("Running model for period",timeperiod[1]))
    }
    E1 = as.vector(rr)*init$E0
    Period <- init$Year
    #Simulate observed as NB(Eit, theta)
    #theta = 1000
    YSIM <- lapply(1:nsim, function(i) MASS::rnegbin(E1, theta = theta))
    Ex <- scale_sim(YSIM, init, nsim, Time)
    #create vectors.sim for spacetime
    vectors.sim <- list(Period = Period, Ex = Ex , E0_0 = init$E0, Y.vec=init$Y.vec)
    
    #create vectors.sim for space-only
    spacevecs <- vectors_space_sim(x, Ex, YSIM, Time, init)
    vectors.sim.s <- spacevecs$vectors.sim.s 
    YSIM.s <- spacevecs$YSIM.s
    
    # #SPACE-ONLY MODELS
    #set up and run simulation models
    message("RUNNING: SPACE-ONLY QUASI-POISSON")
    lassoresult.qp.s <- spacetimeLasso_sim(clusters, vectors.sim.s, covars, 1, spacetime=FALSE, pois=FALSE, nsim, YSIM.s, floor)
    message("RUNNING: SPACE-ONLY POISSON")
    lassoresult.p.s <- spacetimeLasso_sim(clusters, vectors.sim.s, covars, 1, spacetime=FALSE, pois=TRUE, nsim, YSIM.s, floor)
    
    #SPACE-TIME MODELS 
    #set up and run simulation models
    message("RUNNING: SPACE-TIME QUASI-POISSON")
    lassoresult.qp.st <- spacetimeLasso_sim(clusters, vectors.sim, covars, Time, spacetime=TRUE, pois=FALSE, nsim, YSIM, floor)
    message("RUNNING: SPACE-TIME POISSON")
    lassoresult.p.st <- spacetimeLasso_sim(clusters, vectors.sim, covars, Time, spacetime=TRUE, pois=TRUE, nsim, YSIM, floor)
    
    message("All models ran successfully")
    
    #RR and Colors for Plotting
    ##SPACE-ONLY
    initial.s <- list(E0 = unlist(vectors.sim.s$E0_0))
    #id <- rep(1:length(x), times = length(timeperiod))
    id <- rep(1:length(x), times=Time)
    riskratios.qp.s <- get_rr(lassoresult.qp.s, vectors.sim.s,initial.s, 
                               tapply(as.vector(matrix(E1, ncol=Time)), id, function(x) mean(x)),
                               1, sim=TRUE)
    rrcolors.qp.s <- colormapping(riskratios.qp.s,1)
    
    riskratios.p.s <- get_rr(lassoresult.p.s, vectors.sim.s, initial.s, 
                              tapply(as.vector(matrix(E1, ncol=Time)), id, function(x) mean(x)),
                              1, sim=TRUE)
    rrcolors.p.s <- colormapping(riskratios.p.s,1)
    
    ##SPACE-TIME
    riskratios.qp.st <- get_rr(lassoresult.qp.st, vectors.sim,init, E1,Time, sim=TRUE)
    rrcolors.qp.st <- colormapping(riskratios.qp.st,Time)
    
    riskratios.p.st <- get_rr(lassoresult.p.st, vectors.sim,init, E1,Time, sim=TRUE)
    rrcolors.p.st <- colormapping(riskratios.p.st,Time)
    
    #COMBINE RISKRATIOS INTO LISTS
    riskratios <- list(riskratios.qp.s = riskratios.qp.s, riskratios.p.s = riskratios.p.s, 
                       riskratios.qp.st = riskratios.qp.st, riskratios.p.st = riskratios.p.st)
    rrcolors <- list(rrcolors.qp.s = rrcolors.qp.s, rrcolors.p.s = rrcolors.p.s,
                     rrcolors.qp.st = rrcolors.qp.st, rrcolors.p.st = rrcolors.p.st)
    
    #DETECTION
    ##QP - Space
    set <- detect_set(lassoresult.qp.s, vectors.sim.s, as.matrix(rr[,timeperiod[1]]), 1, x, y, rMax, center, radius, nullmod,nsim)
    incluster.qp.s <- detect_incluster(lassoresult.qp.s, vectors.sim.s, as.matrix(rr[,timeperiod[1]]), set, 1, 1, nsim, x, y, rMax, center, 
                                       radius, IC = "ic", under=FALSE, nullmod)
    detect.qp.s <- list(clustDiagnostics(incluster.qp.s, threshold[1], nullmod, nsim), clustDiagnostics(incluster.qp.s , threshold[2], nullmod,nsim))
    if(!is.null(nullmod)){
        detect.out.qp.s <- (matrix(unlist(detect.qp.s), ncol=3, byrow=TRUE,
                                   dimnames = list(c(
                                      # paste0("null.any.", threshold[1]),
                                       paste0("null.summary.mean.", threshold[1]),
                                       paste0("null.summary.median.", threshold[1]),
                                       paste0("null.summary.sd.", threshold[1]),
                                       paste0("prop.null.", threshold[1]),
                                              
                                     #  paste0("null.any.", threshold[2]),
                                       paste0("null.summary.mean.", threshold[2]),
                                       paste0("null.summary.median.", threshold[2]),
                                       paste0("null.summary.sd.", threshold[2]),
                                       paste0("prop.null.", threshold[2])),
                                       c("aic", "aicc", "bic")
                                   )))
    }
    else {
        detect.out.qp.s <- (matrix(unlist(detect.qp.s),ncol=3, byrow=TRUE,
                                                       dimnames = list(c(
                                                           paste0("incluster.any.", threshold[1]),
                                                           paste0("outcluster.", threshold[1]),
                                                           paste0("alldetect.",threshold[1]),
                                                           paste0("potentialclusterdetect.",threshold[1]),
                                                           paste0("trueclusterdetect.",threshold[1]),
                                                           paste0("alldetect.summary.mean.", threshold[1]),
                                                           paste0("alldetect.summary.median.", threshold[1]),
                                                           paste0("alldetect.summary.sd", threshold[1]),
                                                           paste0("potentialdetect.summary.mean.", threshold[1]),
                                                           paste0("potentialdetect.summary.median.", threshold[1]),
                                                           paste0("potentialdetect.summary.sd.", threshold[1]),
                                                           paste0("truedetect.summary.mean.", threshold[1]),
                                                           paste0("truedetect.summary.median.", threshold[1]),
                                                           paste0("truedetect.summary.sd.", threshold[1]),
                                                           
                                                           paste0("incluster.any.", threshold[2]),
                                                           paste0("outcluster.", threshold[2]),
                                                           paste0("alldetect.",threshold[2]),
                                                           paste0("potentialclusterdetect.",threshold[2]),
                                                           paste0("trueclusterdetect.",threshold[2]),
                                                           paste0("alldetect.summary.mean.", threshold[2]),
                                                           paste0("alldetect.summary.median.", threshold[2]),
                                                           paste0("alldetect.summary.sd", threshold[2]),
                                                           paste0("potentialdetect.summary.mean.", threshold[2]),
                                                           paste0("potentialdetect.summary.median.", threshold[2]),
                                                           paste0("potentialdetect.summary.sd.", threshold[2]),
                                                           paste0("truedetect.summary.mean.", threshold[2]),
                                                           paste0("truedetect.summary.median.", threshold[2]),
                                                           paste0("truedetect.summary.sd.", threshold[2])),
                                                           c("aic","aicc","bic"))))
    }
    
    ##P - Space
    set <- detect_set(lassoresult.p.s, vectors.sim.s, as.matrix(rr[,timeperiod[1]]), 1, x, y, rMax, center, radius, nullmod,nsim)
    incluster.p.s <- detect_incluster(lassoresult.p.s, vectors.sim.s, as.matrix(rr[,timeperiod[1]]), set, 1, 1, nsim, x, y, rMax, center, 
                                      radius, IC = "ic", under=FALSE,nullmod)
    detect.p.s <- list(clustDiagnostics(incluster.p.s, threshold[1], nullmod,nsim), clustDiagnostics(incluster.p.s, threshold[2], nullmod,nsim))
    if(!is.null(nullmod)){
        detect.out.p.s <- (matrix(unlist(detect.p.s), ncol=3, byrow=TRUE,
                                  dimnames = list(c(
                                #      paste0("null.any.", threshold[1]),
                                      paste0("null.summary.mean.", threshold[1]),
                                      paste0("null.summary.median.", threshold[1]),
                                      paste0("null.summary.sd.", threshold[1]),
                                      paste0("prop.null.", threshold[1]),
                                      
                                   #   paste0("null.any.", threshold[2]),
                                      paste0("null.summary.mean.", threshold[2]),
                                      paste0("null.summary.median.", threshold[2]),
                                      paste0("null.summary.sd.", threshold[2]),
                                      paste0("prop.null.", threshold[2])),
                                      c("aic", "aicc", "bic")
                                  )))
    }
    else {
        detect.out.p.s <- (matrix(unlist(detect.p.s),ncol=3, byrow=TRUE,
                                  dimnames = list(c(
                                      paste0("incluster.any.", threshold[1]),
                                      paste0("outcluster.", threshold[1]),
                                      paste0("alldetect.",threshold[1]),
                                      paste0("potentialclusterdetect.",threshold[1]),
                                      paste0("trueclusterdetect.",threshold[1]),
                                      paste0("alldetect.summary.mean.", threshold[1]),
                                      paste0("alldetect.summary.median.", threshold[1]),
                                      paste0("alldetect.summary.sd", threshold[1]),
                                      paste0("potentialdetect.summary.mean.", threshold[1]),
                                      paste0("potentialdetect.summary.median.", threshold[1]),
                                      paste0("potentialdetect.summary.sd.", threshold[1]),
                                      paste0("truedetect.summary.mean.", threshold[1]),
                                      paste0("truedetect.summary.median.", threshold[1]),
                                      paste0("truedetect.summary.sd.", threshold[1]),
                                      
                                      paste0("incluster.any.", threshold[2]),
                                      paste0("outcluster.", threshold[2]),
                                      paste0("alldetect.",threshold[2]),
                                      paste0("potentialclusterdetect.",threshold[2]),
                                      paste0("trueclusterdetect.",threshold[2]),
                                      paste0("alldetect.summary.mean.", threshold[2]),
                                      paste0("alldetect.summary.median.", threshold[2]),
                                      paste0("alldetect.summary.sd", threshold[2]),
                                      paste0("potentialdetect.summary.mean.", threshold[2]),
                                      paste0("potentialdetect.summary.median.", threshold[2]),
                                      paste0("potentialdetect.summary.sd.", threshold[2]),
                                      paste0("truedetect.summary.mean.", threshold[2]),
                                      paste0("truedetect.summary.median.", threshold[2]),
                                      paste0("truedetect.summary.sd.", threshold[2])),
                                      c("aic","aicc","bic"))))
    }

    ##QP - SPACETIME
    set <- detect_set(lassoresult.qp.st, vectors.sim, rr, Time, x, y, rMax, center, radius, nullmod,nsim)
    incluster.qp.st <- detect_incluster(lassoresult.qp.st, vectors.sim, rr, set, timeperiod, Time, nsim, x, y, rMax, center, 
                                        radius, IC = "ic", under=FALSE, nullmod)
    detect.qp.st <- list(clustDiagnostics(incluster.qp.st , threshold[1], nullmod,nsim), clustDiagnostics(incluster.qp.st , threshold[2], nullmod,nsim))
    
    if(!is.null(nullmod)){
        detect.out.qp.st <- (matrix(unlist(detect.qp.st), ncol=3, byrow=TRUE,
                                   dimnames = list(c(
                                   #    paste0("null.any.", threshold[1]),
                                       paste0("null.summary.mean.", threshold[1]),
                                       paste0("null.summary.median.", threshold[1]),
                                       paste0("null.summary.sd.", threshold[1]),
                                       paste0("prop.null.", threshold[1]),
                                       
                                   #    paste0("null.any.", threshold[2]),
                                       paste0("null.summary.mean.", threshold[2]),
                                       paste0("null.summary.median.", threshold[2]),
                                       paste0("null.summary.sd.", threshold[2]),
                                       paste0("prop.null.", threshold[2])),
                                       c("aic", "aicc", "bic")
                                   )))
    }
    else{
        detect.out.qp.st <- (matrix(unlist(detect.qp.st),ncol=3, byrow=TRUE,
                                   dimnames = list(c(
                                       paste0("incluster.any.", threshold[1]),
                                       paste0("outcluster.", threshold[1]),
                                       paste0("alldetect.",threshold[1]),
                                       paste0("potentialclusterdetect.",threshold[1]),
                                       paste0("trueclusterdetect.",threshold[1]),
                                       paste0("alldetect.summary.mean.", threshold[1]),
                                       paste0("alldetect.summary.median.", threshold[1]),
                                       paste0("alldetect.summary.sd", threshold[1]),
                                       paste0("potentialdetect.summary.mean.", threshold[1]),
                                       paste0("potentialdetect.summary.median.", threshold[1]),
                                       paste0("potentialdetect.summary.sd.", threshold[1]),
                                       paste0("truedetect.summary.mean.", threshold[1]),
                                       paste0("truedetect.summary.median.", threshold[1]),
                                       paste0("truedetect.summary.sd.", threshold[1]),
                                       
                                       paste0("incluster.any.", threshold[2]),
                                       paste0("outcluster.", threshold[2]),
                                       paste0("alldetect.",threshold[2]),
                                       paste0("potentialclusterdetect.",threshold[2]),
                                       paste0("trueclusterdetect.",threshold[2]),
                                       paste0("alldetect.summary.mean.", threshold[2]),
                                       paste0("alldetect.summary.median.", threshold[2]),
                                       paste0("alldetect.summary.sd", threshold[2]),
                                       paste0("potentialdetect.summary.mean.", threshold[2]),
                                       paste0("potentialdetect.summary.median.", threshold[2]),
                                       paste0("potentialdetect.summary.sd.", threshold[2]),
                                       paste0("truedetect.summary.mean.", threshold[2]),
                                       paste0("truedetect.summary.median.", threshold[2]),
                                       paste0("truedetect.summary.sd.", threshold[2])),
                                       c("aic","aicc","bic"))))
    }
    
    ##P - SPACETIME
    set <- detect_set(lassoresult.p.st, vectors.sim, rr, Time, x, y, rMax, center, radius, nullmod,nsim)
    incluster.p.st <- detect_incluster(lassoresult.p.st, vectors.sim, rr, set, timeperiod, Time, nsim, x, y, rMax, center, 
                                       radius, IC = "ic", under=FALSE, nullmod)
    detect.p.st <- list(clustDiagnostics(incluster.p.st, threshold[1], nullmod,nsim), clustDiagnostics(incluster.p.st, threshold[2], nullmod,nsim))
    
    if(!is.null(nullmod)){
        detect.out.p.st <- (matrix(unlist(detect.p.st), ncol=3, byrow=TRUE,
                                   dimnames = list(c(
                                     #  paste0("null.any.", threshold[1]),
                                       paste0("null.summary.mean.", threshold[1]),
                                       paste0("null.summary.median.", threshold[1]),
                                       paste0("null.summary.sd.", threshold[1]),
                                       paste0("prop.null.", threshold[1]),
                                       
                                      # paste0("null.any.", threshold[2]),
                                       paste0("null.summary.mean.", threshold[2]),
                                       paste0("null.summary.median.", threshold[2]),
                                       paste0("null.summary.sd.", threshold[2]),
                                       paste0("prop.null.", threshold[2])),
                                       c("aic", "aicc", "bic")
                                   )))
    }
    else {
        detect.out.p.st <- (matrix(unlist(detect.p.st),ncol=3, byrow=TRUE,
                                   dimnames = list(c(
                                       paste0("incluster.any.", threshold[1]),
                                       paste0("outcluster.", threshold[1]),
                                       paste0("alldetect.",threshold[1]),
                                       paste0("potentialclusterdetect.",threshold[1]),
                                       paste0("trueclusterdetect.",threshold[1]),
                                       paste0("alldetect.summary.mean.", threshold[1]),
                                       paste0("alldetect.summary.median.", threshold[1]),
                                       paste0("alldetect.summary.sd", threshold[1]),
                                       paste0("potentialdetect.summary.mean.", threshold[1]),
                                       paste0("potentialdetect.summary.median.", threshold[1]),
                                       paste0("potentialdetect.summary.sd.", threshold[1]),
                                       paste0("truedetect.summary.mean.", threshold[1]),
                                       paste0("truedetect.summary.median.", threshold[1]),
                                       paste0("truedetect.summary.sd.", threshold[1]),
                                       
                                       paste0("incluster.any.", threshold[2]),
                                       paste0("outcluster.", threshold[2]),
                                       paste0("alldetect.",threshold[2]),
                                       paste0("potentialclusterdetect.",threshold[2]),
                                       paste0("trueclusterdetect.",threshold[2]),
                                       paste0("alldetect.summary.mean.", threshold[2]),
                                       paste0("alldetect.summary.median.", threshold[2]),
                                       paste0("alldetect.summary.sd", threshold[2]),
                                       paste0("potentialdetect.summary.mean.", threshold[2]),
                                       paste0("potentialdetect.summary.median.", threshold[2]),
                                       paste0("potentialdetect.summary.sd.", threshold[2]),
                                       paste0("truedetect.summary.mean.", threshold[2]),
                                       paste0("truedetect.summary.median.", threshold[2]),
                                       paste0("truedetect.summary.sd.", threshold[2])),
                                       c("aic","aicc","bic"))))
    }
    
    
    #RETURN
    return(list(lassoresult.qp.st = lassoresult.qp.st,
                lassoresult.p.st = lassoresult.p.st,
                lassoresult.qp.s = lassoresult.qp.s,
                lassoresult.p.s = lassoresult.p.s,
                riskratios = riskratios,
                rrcolors = rrcolors,
                rr.mat = rr,
                init.vec = vectors.sim,
                init.vec.s = vectors.sim.s ,
                #return spacetime first
                incluster.qp.st = incluster.qp.st,
                incluster.p.st = incluster.p.st,
                detect.qp.st = detect.qp.st,
                detect.p.st = detect.p.st,
                detect.out.p.st = detect.out.p.st,
                detect.out.qp.st = detect.out.qp.st,
                #return space-only second
                incluster.qp.s = incluster.qp.s,
                incluster.p.s = incluster.p.s,
                detect.qp.s = detect.qp.s,
                detect.p.s = detect.p.s,
                detect.out.p.s = detect.out.p.s,
                detect.out.qp.s = detect.out.qp.s))
}

