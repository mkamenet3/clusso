#' @title
#' vectors_space_sim
#' @description 
#' Simulation functions for space and space-time cluster detection with lasso. This function will collapse a space-time vector onto space only
#' @param x vector coordinates (unique regardless of time period)
#' @param Ex list of simulated and standardized expected counts
#' @param YSIM simulated observed
#' @param Time number of time periods
#' @param init initial list of vectors, inherited from function setVectors.
#' @return returns space-time 
#' 
vectors_space_sim <- function(x,Ex, YSIM,Time, init){
    id <- rep(1:length(x), times = Time)
    if(length(id)!=length(as.vector(Ex[[1]]))){
        stop("Length of ID var not equal to number of observations")  
    } 
    if(!is.null(init$covars)){
        covars.sim.s <- sapply(1:ncol(init$covars), 
                           function(i) tapply(as.vector(matrix(init$covars[,i], ncol=Time)),id, 
                                              function(x) sum(x)))
    }
    else{
        covars.sim.s <- NULL
    }
    vectors.sim.s <- list(Period = rep("1", length(x)),
                          Ex = lapply(1:length(YSIM), function(i) tapply(as.vector(matrix(Ex[[i]], ncol=Time)), id, function(x) sum(x))),
                          E0_0 = tapply(as.vector(matrix(init$E0, ncol=Time)), id, function(x) sum(x)),
                          Y.vec = tapply(as.vector(matrix(init$Y.vec, ncol=Time)), id, function(x) round(sum(x))),
                          covars.s = as.data.frame(covars.sim.s))
    YSIM.s <- lapply(1:length(YSIM), function(i) tapply(as.vector(matrix(YSIM[[i]], ncol=Time)), id, function(x) round(sum(x))))
    return(list(vectors.sim.s = vectors.sim.s, YSIM.s = YSIM.s))
}
    

#'@title
#'clust_sim
#'@description 
#'This helper function runs both the space and space-time Lasso model simulations for all 4 models simulataneously: Quasi-Poisson vs. Poisson in both space and space-time.
#' This function is to be run on simulated data and all four models are run on the same simulated set. 
#'A separate function (clust.sim) can be used for running simulations on individual models and (clust) can be used for observed data.
#'@param clst list; output from toclust function. Must be of class clst.
#'@param x x coordinates (easting/latitude); if utm coordinates, scale to km.
#'@param y y coordinates (northing/longitude); if utm coordinates, scale to km.
#'@param rMax set max radius (in km)
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
#'@param overdispfloor overdispfloor default is TRUE. When TRUE, it limits phi (overdispersion parameter) to be greater or equal to 1. If FALSE, will allow for under dispersion.
#'@param collapsetime alternative definition for space-only model to instead collapse expected and observed counts across time. TODO
#'@return returns list of lists
#'@export
#'@examples
#'\donttest{
#'data(japanbreastcancer)
#'#Set Some Initial Conditions
#'x1=utmJapan$utmx/1000
#'y1=utmJapan$utmy/1000
#'#Set initial parameters
#'rMax <- 20 
#'Time=5
#'nsim=2 #number of simulations
#'center=1
#'radius=18
#'timeperiods = c(1:3)
#'risk.ratio=2
#'theta = 10000
#'overdispfloor = TRUE
#'nullmod <- TRUE
#'dframe <- dframe[,-1]
#'clst <- toclust(japanbreastcancer, expected = japanbreastcancer$expdeath, 
#'  observed = japanbreastcancer$death,timeperiod = japanbreastcancer$period, covars = FALSE)
#'res <- clust_sim(clst, x1,y1,rMax, Time, nsim,center, radius, risk.ratio, 
#'  timeperiod, utm=TRUE, byrow=TRUE, 
#'threshold, space= "both",theta = theta, nullmod = TRUE, overdispfloor)
#'}
clust_sim <- function(clst, x,y, rMax, Time, nsim, center, radius, risk.ratio, 
                          timeperiod, utm=TRUE, byrow=TRUE, threshold, space = c("space", "spacetime", "both"), 
                      theta = NULL,nullmod=NULL, overdispfloor,collapsetime=FALSE){
    if(is(clst, "clst")!=TRUE) stop("clst element not of class `clst`. This is required for the clust_sim function.")
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
    if(overdispfloor==FALSE){
        overdispfloor <- FALSE
    }
    else{
        overdispfloor <- TRUE
    }
    if(!is.null(nullmod)){
        if(isTRUE(nullmod)) warning("Null mod has been set to true")
        nullmod<-TRUE
        message("Running null models. For simulation models, leave nullmod NULL")
    }
    else{
        nullmod<-NULL
    }
    if(is.infinite(theta)){
        theta <- theta
        message("Running model with no overdispersion (pure Poisson model)")
    }
    else if(!is.null(theta) & !is.infinite(theta)){
        theta <- theta
        message("Running model with user-specified theta value")
    }
    else{
        theta <- 1000
        message("Running model with default theta value of 1000")
    }
    if(!is.null(threshold)){
        thresh <- threshold
    }
    else{
        thresh <- NULL
        message("No threshold diagnostics specified.")
    }
    if(length(space) > 1) stop("You must select either `space`, `spacetime`, or `both`")
    space <- match.arg(space, several.ok = FALSE)
    switch(space, 
           #TODO
           # space = 
           # spacetime = 
           both = clustAll_sim(x, y, rMax,period, expected, observed, covars, Time, nsim, center, radius, risk.ratio,
                                     timeperiod,utm, byrow, thresh, theta, nullmod, overdispfloor, collapsetime))
}


#'
#'@title
#'clustAll_sim
#' @description 
#'This function runs both the space and space-time Lasso model simulations for all 4 models simulataneously: Quasi-Poisson vs. Poisson in both space and space-time.
#' This function is to be run on simulated data and all four models are run on the same simulated set. 
#'A separate function (clust.sim) can be used for running simulations on individual models and (clust) can be used for observed data.
#'@param x x coordinates (easting/latitude); if utm coordinates, scale to km.
#'@param y y coordinates (northing/longitude); if utm coordinates, scale to km.
#'@param rMax set max radius (in km)
#'@param period vector of periods or years in dataset. Should be imported as a factor.
#'@param expected vector of expected counts. Expected counts must match up with the year and observed vectors.
#'@param observed vector of observed counts. Observed counts must match up with the year and expected vectors.
#'@param covars dataframe of covariates, if supplied to `toclust` function. 
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
#'@param thresh  vector or value as threshold for cluster detection
#'@param theta default is 1000. Can add in overdispersion to simulated model by changing this value.
#'@param nullmod if TRUE, then null models will be run. Otherwise, default is null.
#'@param overdispfloor overdispfloor default is TRUE. When TRUE, it limits phi (overdispersion parameter) to be greater or equal to 1. If FALSE, will allow for under dispersion.
#'@param collapsetime alternative definition for space-only model to instead collapse expected and observed counts across time. TODO
#'@inheritParams clust_sim
#'@return returns list of lists


clustAll_sim <- function(x, y, rMax, period, expected, observed, covars,Time, nsim, center, radius, risk.ratio, 
                               timeperiod,utm, byrow, thresh, theta = theta, nullmod=nullmod,
                         overdispfloor=overdispfloor, collapsetime=FALSE){
    message("Running both Space and Space-Time Models")
    
    #set up clusters and fitted values
    clusters <- clusters2df(x,y,rMax, utm=utm, length(x))
    n <- length(x)
    init <- setVectors(period, expected, observed, covars, Time, byrow)
    
    ##Multiple centers/centroids for clusters?
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
    
    #create rr (space-time) and rr.s (space-only) matrices of 1's
    rr = matrix(1, nrow=n, ncol=Time)
    rr.s = matrix(1, nrow=n, ncol=Time)
    
    ##Space-time
    #Create cluster across the time periods
    timelength <- length(timeperiod)
    if(timelength > 1){
        rr[cluster$last, timeperiod[1]:tail(timeperiod, n=1)] <- risk.ratio
        message(paste("Running model for periods",timeperiod[1],"through", tail(timeperiod, n=1)))
    }
    else if(timelength==1){
        rr[cluster$last, timeperiod:Time] = risk.ratio
        message(paste("Running model for period",timeperiod[1]))
    }
    #check for errors
    #if(isTRUE(all.equal(timeperiod,which(unique(rr, fromLast=TRUE)[1,]!=1)))==FALSE) stop("Timeperiods not equal to time elements in space-time rr matrix.")
    
    ##Space-only
    allTime <- 1:Time
    rr.s[cluster$last, allTime[1]:tail(allTime, n=1)] <- risk.ratio
    #Check for errors
    #if(isTRUE(all.equal(allTime,which(unique(rr.s, fromLast=TRUE)[1,]!=1)))==FALSE) stop("Timeperiods not equal to time elements in space-only rr.s matrix.")
    
    ##Expected Counts and simulations
    #Expected matrices
    E1 <- as.vector(rr)*init$E0
    E1.s <- as.vector(rr.s)*init$E0
    Period <- init$Year
    
    #Simulate observed as NB(Eit, theta)
    if(is.infinite(theta)){
        YSIM <- lapply(1:nsim, function(i) rpois(length(E1), lambda = E1))
        YSIM.s <- lapply(1:nsim, function(i) rpois(length(E1.s), lambda = E1.s))
    }
    else{
        YSIM <- lapply(1:nsim, function(i) MASS::rnegbin(E1, theta = theta))
        YSIM.s <- lapply(1:nsim, function(i) MASS::rnegbin(E1.s, theta = theta))
    }
    
    #Scale
    Ex <- scale_sim(YSIM, init, nsim, Time)
    Ex.s <- scale_sim(YSIM.s, init, nsim, Time)
    
    #create vectors.sim for spacetime
    vectors.sim <- list(Period = Period, Ex = Ex , E0_0 = init$E0, 
                        Y.vec=init$Y.vec, covars = covars)
    vectors.sim.s <- list(Period = Period, Ex = Ex.s , E0_0 = init$E0, 
                        Y.vec=init$Y.vec, covars = covars)
    
    #create vectors.sim for space-only ##MAYBE WILL RE_INTRODUCE THIS CODE AT SOME POINT
    # spacevecs <- vectors_space_sim(x, Ex, YSIM, Time, init)
    # vectors.sim.s <- spacevecs$vectors.sim.s 
    # YSIM.s <- spacevecs$YSIM.s
    
####################################################################################
####################################################################################
####################################################################################    
#TEST THIS CHUNK    
    n_uniq <- length(unique(clusters$center))
    potClus <- n_uniq
    numCenters <- n_uniq
    #CREATE sparseMAT and cache it for use throughout this function
    if(collapsetime==FALSE){
        sparseMAT <- spacetimeMat(clusters, numCenters, Time)
        SOAR::Store(sparseMAT)
        message("Space-time matrix created")
        
    }
    else{
        sparseMAT <- spaceMat(clusters, numCenters)
        SOAR::Store(sparseMAT)
        message("Spatial matrix created")
        if(nrow(covars)==0){
            covars<- NULL
        }
    }
#TEST THIS CHUNK        
####################################################################################
####################################################################################
####################################################################################    
    

    
    
    # #SPACE-ONLY MODELS
    #set up and run simulation models
    message("RUNNING: SPACE-ONLY QUASI-POISSON")
    lassoresult.qp.s <- spacetimeLasso_sim(sparseMAT, n_uniq, vectors.sim.s, Time, spacetime=TRUE, pois=FALSE, nsim, YSIM.s, overdispfloor)

    message("RUNNING: SPACE-ONLY POISSON")
    lassoresult.p.s <- spacetimeLasso_sim(sparseMAT, n_uniq, vectors.sim.s, Time, spacetime=TRUE, pois=TRUE, nsim, YSIM.s, overdispfloor)
    
    #SPACE-TIME MODELS 
    #set up and run simulation models
    message("RUNNING: SPACE-TIME QUASI-POISSON")
    lassoresult.qp.st <- spacetimeLasso_sim(sparseMAT,n_uniq,  vectors.sim, Time, spacetime=TRUE, pois=FALSE, nsim, YSIM, overdispfloor)
    message("RUNNING: SPACE-TIME POISSON")
    lassoresult.p.st <- spacetimeLasso_sim(sparseMAT, n_uniq, vectors.sim, Time, spacetime=TRUE, pois=TRUE, nsim, YSIM, overdispfloor)
    
    message("All models ran successfully")
    
    #RR and Colors for Plotting
    ##SPACE-ONLY
    initial.s <- list(E0 = unlist(vectors.sim.s$E0_0))
    #id <- rep(1:length(x), times=Time)
    ################################################################
    ###Space - QP
    ####RR
    riskratios.qp.s <- get_rr(lassoresult.qp.s, vectors.sim.s,initial.s, 
                              E1.s,
                              Time, sim=TRUE, cv= NULL)
    rrcolors.qp.s <- colormapping(riskratios.qp.s, Time = 5, cv=NULL, prob=FALSE)
    
    ####Probs and threshold probs
    pb.qp.s <- get_prob(lassoresult.qp.s, initial.s,
                        E1.s ,n, Time, nsim, thresh)
    probcolors.qp.s <- colormapping(pb.qp.s$probs, Time, cv = NULL, prob=TRUE)
    probcolors.qp.s.thresh1 <- colormapping(pb.qp.s$probs.thresh1, Time, cv = NULL, prob=TRUE)
    probcolors.qp.s.thresh2 <- colormapping(pb.qp.s$probs.thresh2, Time, cv = NULL, prob=TRUE)
    
    ###Space - P
    ####RR
    riskratios.p.s <- get_rr(lassoresult.p.s, vectors.sim.s, initial.s, 
                             E1.s,
                             Time, sim=TRUE, cv=NULL)
    rrcolors.p.s <- colormapping(riskratios.p.s,Time, cv=NULL, prob=FALSE)
    
    ####Probs and thresh probs
    pb.p.s <- get_prob(lassoresult.p.s, initial.s,
                       E1.s ,n, Time, nsim,thresh)
    probcolors.p.s <- colormapping(pb.p.s$probs, Time, cv = NULL, prob = TRUE)
    probcolors.p.s.thresh1 <- colormapping(pb.p.s$probs.thresh1, Time, cv = NULL, prob=TRUE)
    probcolors.p.s.thresh2 <- colormapping(pb.p.s$probs.thresh2, Time, cv = NULL, prob=TRUE)
     
    ################################
    
    ###Space-TIME - QP
    ####RR
    riskratios.qp.st <- get_rr(lassoresult.qp.st, vectors.sim,init, E1,Time, sim=TRUE, cv=NULL)
    rrcolors.qp.st <- colormapping(riskratios.qp.st,Time, cv=NULL, prob = FALSE)
    
    ####Probs and thresh probs
    pb.qp.st <- get_prob(lassoresult.qp.st,init, E1, n, Time, nsim,thresh)
    probcolors.qp.st <- colormapping(pb.qp.st$probs, Time, cv = NULL, prob = TRUE)
    probcolors.qp.st.thresh1 <- colormapping(pb.qp.st$probs.thresh1, Time, cv = NULL, prob = TRUE)
    probcolors.qp.st.thresh2 <- colormapping(pb.qp.st$probs.thresh2, Time, cv = NULL, prob = TRUE)
    
    
    ###Space-TIME - P
    ####RR
    riskratios.p.st <- get_rr(lassoresult.p.st, vectors.sim,init, E1,Time, sim=TRUE, cv=NULL)
    rrcolors.p.st <- colormapping(riskratios.p.st,Time, cv=NULL, prob = FALSE)
    
    ####Probs and thresh probs
    pb.p.st <- get_prob(lassoresult.p.st, init, E1, n, Time, nsim,thresh)
    probcolors.p.st <- colormapping(pb.p.st$probs, Time, cv = NULL, prob = TRUE)
    probcolors.p.st.thresh1 <- colormapping(pb.p.st$probs.thresh1, Time, cv = NULL, prob = TRUE)
    probcolors.p.st.thresh2 <- colormapping(pb.p.st$probs.thresh2, Time, cv = NULL, prob = TRUE)
    
    ################################################################
    
    #COMBINE RISKRATIOS INTO LISTS
    riskratios <- list(riskratios.qp.s = riskratios.qp.s, riskratios.p.s = riskratios.p.s, 
                       riskratios.qp.st = riskratios.qp.st, riskratios.p.st = riskratios.p.st)
    rrcolors <- list(rrcolors.qp.s = rrcolors.qp.s, rrcolors.p.s = rrcolors.p.s,
                     rrcolors.qp.st = rrcolors.qp.st, rrcolors.p.st = rrcolors.p.st)
    #COMBINE Probabilities INTO LISTS
    probrates <- list(pb.qp.s = pb.qp.s, pb.p.s = pb.p.s,
                      pb.qp.st = pb.qp.st, pb.p.st = pb.p.st)
    probcolors <- list(probcolors.qp.s = probcolors.qp.s, probcolors.p.s = probcolors.p.s,
                       probcolors.qp.st = probcolors.qp.st, probcolors.p.st = probcolors.p.st)
    probcolors.thresh1 <- list(probcolors.qp.s = probcolors.qp.s.thresh1, probcolors.p.s = probcolors.p.s.thresh1,
                               probcolors.qp.st = probcolors.qp.st.thresh1, probcolors.p.st = probcolors.p.st.thresh1)
    probcolors.thresh2 <- list(probcolors.qp.s = probcolors.qp.s.thresh2, probcolors.p.s = probcolors.p.s.thresh2,
                               probcolors.qp.st = probcolors.qp.st.thresh2, probcolors.p.st = probcolors.p.st.thresh2)
    
    
    #DETECTION
    ################################################################
    ##QP - Space
    set <- detect_set(lassoresult.qp.s, vectors.sim.s, rr.s, Time, x, y, rMax, center, radius, nullmod,nsim)
    incluster.qp.s <- detect_incluster(sparseMAT,lassoresult.qp.s, vectors.sim.s, rr.s, set, 1:Time, Time, nsim, under=FALSE, nullmod,risk.ratio,
                                       center, radius,x,y,rMax,thresh)
    if(!is.null(nullmod)){
        detect.out.qp.s <- (matrix(unlist(incluster.qp.s), ncol=3, byrow=TRUE,
                                   dimnames = list(c(
                                       paste0("prop.null.")),
                                       c("aic", "aicc", "bic")
                                   )))
        detect.out.thresh.qp.s <- NULL
    }
    else {
        detect.out.qp.s <- (matrix(unlist(incluster.qp.s[1:12]),ncol=3, byrow=TRUE,
                                                       dimnames = list(c(
                                                           paste0("incluster.centroid.", "nothresh"),
                                                           paste0("outcluster.centroid.", "nothresh"),
                                                           paste0("notinpotclus.","nothresh"),
                                                           paste0("inpotclus.","nothresh")),
                                                           c("aic","aicc","bic"))))
        detect.out.thresh.qp.s <- incluster.qp.s[13]
    }
    ################################################################
    ##P - Space
    set <- detect_set(lassoresult.p.s, vectors.sim.s, rr.s, Time, x, y, rMax, center, radius, nullmod,nsim)
    incluster.p.s <- detect_incluster(sparseMAT,lassoresult.p.s, vectors.sim.s, rr.s, set, 1:Time, Time, nsim, under=FALSE, nullmod,risk.ratio,
                                      center, radius,x,y,rMax,thresh)
    
    if(!is.null(nullmod)){
        detect.out.p.s <-  (matrix(unlist(incluster.p.s), ncol=3, byrow=TRUE,
                                   dimnames = list(c(
                                       paste0("prop.null.")),
                                       c("aic", "aicc", "bic")
                                   )))
        detect.out.thresh.p.s <- NULL
    }
    else {
        detect.out.p.s <- (matrix(unlist(incluster.p.s[1:12]),ncol=3, byrow=TRUE,
                                   dimnames = list(c(
                                       paste0("incluster.centroid.", "nothresh"),
                                       paste0("outcluster.centroid.", "nothresh"),
                                       paste0("notinpotclus.","nothresh"),
                                       paste0("inpotclus.","nothresh")),
                                       c("aic","aicc","bic"))))
        detect.out.thresh.p.s <- incluster.p.s[13]
                                      
    }
    ################################################################
    ##QP - SPACETIME
    set <- detect_set(lassoresult.qp.st, vectors.sim, rr, Time, x, y, rMax, center, radius, nullmod,nsim)
    incluster.qp.st <- detect_incluster(sparseMAT,lassoresult.qp.st, vectors.sim, rr, set, timeperiod, Time, nsim, under=FALSE, nullmod,risk.ratio,
                                        center, radius,x,y,rMax,thresh)
 
    
    if(!is.null(nullmod)){
        detect.out.qp.st <- (matrix(unlist(incluster.qp.st), ncol=3, byrow=TRUE,
                                    dimnames = list(c(
                                        paste0("prop.null.")),
                                        c("aic", "aicc", "bic")
                                    )))
        detect.out.thresh.qp.st <- NULL
    }
    else{
        detect.out.qp.st <- (matrix(unlist(incluster.qp.st[1:12]),ncol=3, byrow=TRUE,
                                   dimnames = list(c(
                                       paste0("incluster.centroid.", "nothresh"),
                                       paste0("outcluster.centroid.", "nothresh"),
                                       paste0("notinpotclus.","nothresh"),
                                       paste0("inpotclus.","nothresh")),
                                       c("aic","aicc","bic"))))
        detect.out.thresh.qp.st <- incluster.qp.st[13]
    }
    ################################################################
    ##P - SPACETIME
    set <- detect_set(lassoresult.p.st, vectors.sim, rr, Time, x, y, rMax, center, radius, nullmod,nsim)
    incluster.p.st <- detect_incluster(sparseMAT,lassoresult.p.st, vectors.sim, rr, set, timeperiod, Time, nsim, under=FALSE, nullmod,risk.ratio,
                                       center, radius,x,y,rMax,thresh)
    
    
    if(!is.null(nullmod)){
        detect.out.p.st <- (matrix(unlist(incluster.p.st), ncol=3, byrow=TRUE,
                                   dimnames = list(c(
                                       paste0("prop.null.")),
                                       c("aic", "aicc", "bic")
                                   )))
        detect.out.thresh.p.st <- NULL
    }
    else {
        detect.out.p.st <- (matrix(unlist(incluster.p.st[1:12]),ncol=3, byrow=TRUE,
                                   dimnames = list(c(
                                       paste0("incluster.centroid.", "nothresh"),
                                       paste0("outcluster.centroid.", "nothresh"),
                                       paste0("notinpotclus.","nothresh"),
                                       paste0("inpotclus.","nothresh")),
                                       c("aic","aicc","bic"))))
        detect.out.thresh.p.st <- incluster.p.st[13]
    }
    ################################################################
    SOAR::Remove(sparseMAT)
    
    #RETURN
    return(list(lassoresult.qp.st = lassoresult.qp.st,
                lassoresult.p.st = lassoresult.p.st,
                lassoresult.qp.s = lassoresult.qp.s,
                lassoresult.p.s = lassoresult.p.s,
                riskratios = riskratios,
                rrcolors = rrcolors,
                probrates = probrates,
                probcolors = probcolors,
                probcolors.thresh1 = probcolors.thresh1,
                probcolors.thresh2 = probcolors.thresh2,
                rr.mat = rr,
                rr.mat.s = rr.s,
                init.vec = vectors.sim,
                init.vec.s = vectors.sim.s ,
                #return spacetime first
                incluster.qp.st = incluster.qp.st,
                incluster.p.st = incluster.p.st,
                detect.out.thresh.qp.st = detect.out.thresh.qp.st,
                detect.out.thresh.p.st = detect.out.thresh.p.st,
                detect.out.p.st = detect.out.p.st,
                detect.out.qp.st = detect.out.qp.st,
                #return space-only second
                incluster.qp.s = incluster.qp.s,
                incluster.p.s = incluster.p.s,
                detect.out.thresh.qp.s = detect.out.thresh.qp.s,
                detect.out.thresh.p.s = detect.out.thresh.p.s,
                detect.out.p.s = detect.out.p.s,
                detect.out.qp.s = detect.out.qp.s))
}

