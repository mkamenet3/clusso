#'Simulation functions for space and space-time cluster detection with lasso. 
#'


#'clust.sim
#'
#'This function runs both the space and space-time Lasso model simulations. This function is to be run on simulated data for one model at a time (space or 
#'spacetime, Poisson or Quasi-Poisson). A separate function (clust) can be used for observed data. For both quasi-Poisson and Poisson use *clust.sim.all* and 
#'and specify space="space" for the two spatial models or "spacetime" for the two spatio-temporal models. If you specify space="both" it will run all four models. 
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
    clusters <- clusters2df(x,y,rMax, utm=TRUE, length(x))
    n <- length(x)
    init <- setVectors(period, expected, observed, Time=Time, byrow=row)
    
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
    lassoresult <- spacetimeLasso_sim(clusters, vectors.sim, Time, spacetime=spacetime, pois=pois, nsim, YSIM)
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




#'
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
#'@param space space and space-time. Default is to run all four models: quasi-poisson and poisson for both space and space-time. User can specify, space = space,
#'space = spacetime, or space = both.
#'@return
#'@details Optional functions include:
#'- 1) utm - default is FALSE. If you have utm coordinates, you want to change this to TRUE.
#'@export
#'TODO allow user to change theta parameter in simulation
clust.sim.all <- function(x, y, rMax, period, expected, observed, Time, nsim, center, radius, risk.ratio, 
                          timeperiod,colors=NULL,utm=TRUE, byrow=TRUE, threshold, space = c("space", "spacetime", "both"),...){
    #initial user setting
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
    space <- match.arg(space, several.ok = FALSE)
    if(length(space) > 1) stop("You must select either `space`, `spacetime`, or `both`")
    print(byrow, utm)
    switch(space, 
           space = clust.sim.all.space(x, y, rMax,period, expected, observed, Time, nsim, center, radius, risk.ratio, 
                                       timeperiod,colors=NULL,utm, byrow, threshold, space=TRUE),
           spacetime = clust.sim.all.spacetime(x, y, rMax,period, expected, observed, Time, nsim, center, radius, risk.ratio, 
                                               timeperiod,colors=NULL,utm, byrow, threshold, space=FALSE),
           both = clust.sim.all.both(x, y, rMax,period, expected, observed, Time, nsim, center, radius, risk.ratio, 
                                     timeperiod,colors=NULL,utm, byrow, threshold,...))
}

#'clust.sim.all.space
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
#'@param space space and space-time. Default is to run all four models: quasi-poisson and poisson for both space and space-time. User can specify, space = space,
#'space = spacetime, or space = both.
#'@return
#'@details Optional functions include:
#'- 1) utm - default is FALSE. If you have utm coordinates, you want to change this to TRUE.
#'@export
#'TODO allow user to change theta parameter in simulation
clust.sim.all.space <- function(x, y, rMax, period, expected, observed, Time, nsim, center, radius, risk.ratio, 
                                timeperiod,colors=NULL,utm, byrow, threshold, space = TRUE,...){
    #set up for space only
    timeperiod = 1
    Time = 1
    #Averaging over time for observed and expected and creating universal period "1"
    id <- rep(1:length(x), each = Time)
    observed  <- tapply(observed, id, function(x) round(mean(x)))
    expected <- tapply(expected, id, function(x) mean(x))
    period <- rep("1", length(x))
    message("Running Space-Only Model")
    
    #set up clusters and fitted values
    clusters <- clusters2df(x,y,rMax, utm, length(x))
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
    rr[cluster$last, timeperiod:Time] = risk.ratio    
    E1 = as.vector(rr)*init$E0
    Period <- init$Year
    #Simulate observed as NB(Eit, theta)
    theta = 1000
    YSIM <- lapply(1:nsim, function(i) rnegbin(E1, theta = theta))
    Ex <- scale.sim(YSIM, init, nsim, Time)
    vectors.sim <- list(Period = Period, Ex = Ex , E0_0 = init$E0, Y.vec=init$Y.vec)
    
    # #SPACE-ONLY MODELS
    #set up and run simulation models
    lassoresult.qp.s <- spacetimeLasso_sim(clusters, vectors.sim, Time, spacetime=FALSE, pois=FALSE, nsim, YSIM)
    lassoresult.p.s <- spacetimeLasso_sim(clusters, vectors.sim, Time, spacetime=FALSE, pois=TRUE, nsim, YSIM)
    
    #RR and Colors for Plotting
    riskratios.qp.s <- get.rr2(lassoresult.qp.s, vectors.sim,init, E1,Time, sim=TRUE)
    rrcolors.qp.s <- colormapping(riskratios.qp.s,Time)
    riskratios.p.s <- get.rr2(lassoresult.qp.s, vectors.sim,init, E1,Time, sim=TRUE)
    rrcolors.p.s <- colormapping(riskratios.p.s,Time)
    ##Combine these into a single list
    riskratios <- list(riskratios.qp.s = riskratios.qp.s, riskratios.p.s = riskratios.p.s)
    rrcolors <- list(rrcolors.qp.s = rrcolors.qp.s, rrcolors.p.s = rrcolors.p.s)
    
    #Detection
    ##QP - Space
    set <- detect_set(lassoresult.qp.s, vectors.sim, rr, Time, x, y, rMax, center, radius)
    incluster.qp.s <- detect.incluster(lassoresult.qp.s, vectors.sim, rr, set, timeperiod, Time, nsim, x, y, rMax, center, 
                                       radius, IC = "ic")
    detect.qp.s <- list(clust.diagnostics(incluster.qp.s, threshold[1]), clust.diagnostics(incluster.qp.s , threshold[2]))
    detect.out.qp.s <- (matrix(unlist(detect.qp.s),ncol=3, byrow=TRUE, 
                               dimnames = list(c(paste0("incluster.any.", threshold[1]),
                                                 paste0("alldetect.",threshold[1]), 
                                                 paste0("potentialclusterdetect.",threshold[1]), 
                                                 paste0("trueclusterdetect.",threshold[1]),
                                                 paste0("incluster.any.",threshold[2]), paste0("alldetect.",threshold[2]),
                                                 paste0("potentialclusterdetect.",threshold[2]), 
                                                 paste0("trueclusterdetect.",threshold[2])),c("aic","aicc","bic"))))
    
    ##P - Space
    set <- detect_set(lassoresult.p.s, vectors.sim, rr, Time, x, y, rMax, center, radius)
    incluster.p.s <- detect.incluster(lassoresult.p.s, vectors.sim, rr, set, timeperiod, Time, nsim, x, y, rMax, center, 
                                      radius, IC = "ic")
    detect.p.s <- list(clust.diagnostics(incluster.p.s, threshold[1]), clust.diagnostics(incluster.p.s, threshold[2]))
    detect.out.p.s <- (matrix(unlist(detect.p.s),ncol=3, byrow=TRUE, 
                              dimnames = list(c(paste0("incluster.any.", threshold[1]),
                                                paste0("alldetect.",threshold[1]), 
                                                paste0("potentialclusterdetect.",threshold[1]), 
                                                paste0("trueclusterdetect.",threshold[1]),
                                                paste0("incluster.any.",threshold[2]), paste0("alldetect.",threshold[2]),
                                                paste0("potentialclusterdetect.",threshold[2]), 
                                                paste0("trueclusterdetect.",threshold[2])),c("aic","aicc","bic"))))
    
    return(list(lassoresult.qp.s = lassoresult.qp.s,
                lassoresult.p.s = lassoresult.p.s,
                riskratios = riskratios,
                rrcolors = rrcolors,
                rr.mat = rr,
                init.vec = vectors.sim,
                incluster.qp.s = incluster.qp.s,
                incluster.p.s = incluster.p.s,
                detect.qp.s = detect.qp.s,
                detect.p.s = detect.p.s,
                detect.out.p.s = detect.out.p.s,
                detect.out.qp.s = detect.out.qp.s))
}

#'clust.sim.all.spacetime
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
#'@param space space and space-time. Default is to run all four models: quasi-poisson and poisson for both space and space-time. User can specify, space = space,
#'space = spacetime, or space = both.
#'@return
#'@details Optional functions include:
#'- 1) utm - default is FALSE. If you have utm coordinates, you want to change this to TRUE.
#'@export
#'TODO allow user to change theta parameter in simulation
clust.sim.all.spacetime <- function(x, y, rMax, period, expected, observed, Time, nsim, center, radius, risk.ratio, 
                                    timeperiod,colors=NULL,utm, byrow, threshold, space = FALSE,...){
    
    message("Running Space-Time Model")
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
    }
    E1 = as.vector(rr)*init$E0
    Period <- init$Year
    #Simulate observed as NB(Eit, theta)
    theta = 1000
    YSIM <- lapply(1:nsim, function(i) rnegbin(E1, theta = theta))
    Ex <- scale.sim(YSIM, init, nsim, Time)
    vectors.sim <- list(Period = Period, Ex = Ex , E0_0 = init$E0, Y.vec=init$Y.vec)
    print(str(vectors.sim))
    #SPACE-TIME MODELS 
    #set up and run simulation models
    lassoresult.qp.st <- spacetimeLasso_sim(clusters, vectors.sim, Time, spacetime=TRUE, pois=FALSE, nsim, YSIM)
    lassoresult.p.st <- spacetimeLasso_sim(clusters, vectors.sim, Time, spacetime=TRUE, pois=TRUE, nsim, YSIM)
    
    #RR and Colors for Plotting
    riskratios.qp.st <- get.rr2(lassoresult.qp.st, vectors.sim,init, E1,Time, sim=TRUE)
    rrcolors.qp.st <- colormapping(riskratios.qp.st,Time)
    
    riskratios.p.st <- get.rr2(lassoresult.qp.st, vectors.sim,init, E1,Time, sim=TRUE)
    rrcolors.p.st <- colormapping(riskratios.p.st,Time)
    
    riskratios <- list(riskratios.qp.st = riskratios.qp.st, riskratios.p.st = riskratios.p.st)
    rrcolors <- list(rrcolors.qp.st = rrcolors.qp.st, rrcolors.p.st = rrcolors.p.st)
    
    #Detection
    ##QP - ST
    set <- detect_set(lassoresult.qp.st, vectors.sim, rr, Time, x, y, rMax, center, radius)
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
    set <- detect_set(lassoresult.p.st, vectors.sim, rr, Time, x, y, rMax, center, radius)
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

#'clust.sim.all.both
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
#'@param space space and space-time. Default is to run all four models: quasi-poisson and poisson for both space and space-time. User can specify, space = space,
#'space = spacetime, or space = both.
#'@return
#'@details Optional functions include:
#'- 1) utm - default is FALSE. If you have utm coordinates, you want to change this to TRUE.
#'@export
#'TODO allow user to change theta parameter in simulation
clust.sim.all.both <- function(x, y, rMax, period, expected, observed, Time, nsim, center, radius, risk.ratio, 
                               timeperiod,colors=NULL,utm, byrow, threshold, ...){
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
    theta = 1000
    YSIM <- lapply(1:nsim, function(i) rnegbin(E1, theta = theta))
    Ex <- scale.sim(YSIM, init, nsim, Time)
    #create vectors.sim for spacetime
    vectors.sim <- list(Period = Period, Ex = Ex , E0_0 = init$E0, Y.vec=init$Y.vec)
    
    #create vectors.sim for space-only
    id <- rep(1:length(x), times = length(timeperiod))
    vectors.sim.s <- list(Period = rep("1", length(x)),
                          Ex = lapply(1:nsim, function(i) tapply(as.vector(matrix(Ex[[i]], ncol=Time)[, timeperiod]), id, function(x) mean(x))),
                          E0_0 = tapply(as.vector(matrix(init$E0, ncol=Time)[, timeperiod]), id, function(x) mean(x)),
                          Y.vec = tapply(as.vector(matrix(init$Y.vec, ncol=Time)[, timeperiod]), id, function(x) round(mean(x))))
    YSIM.s <- lapply(1:nsim, function(i) tapply(as.vector(matrix(YSIM[[i]], ncol=Time)[, timeperiod]), id, function(x) round(mean(x))))
    
    # #SPACE-ONLY MODELS
    #set up and run simulation models
    message("RUNNING: SPACE-ONLY QUASI-POISSON")
    lassoresult.qp.s <- spacetimeLasso_sim(clusters, vectors.sim.s, 1, spacetime=FALSE, pois=FALSE, nsim, YSIM.s)
    message("RUNNING: SPACE-ONLY POISSON")
    lassoresult.p.s <- spacetimeLasso_sim(clusters, vectors.sim.s, 1, spacetime=FALSE, pois=TRUE, nsim, YSIM.s)
    
    #SPACE-TIME MODELS 
    #set up and run simulation models
    message("RUNNING: SPACE-TIME QUASI-POISSON")
    lassoresult.qp.st <- spacetimeLasso_sim(clusters, vectors.sim, Time, spacetime=TRUE, pois=FALSE, nsim, YSIM)
    message("RUNNING: SPACE-TIME POISSON")
    lassoresult.p.st <- spacetimeLasso_sim(clusters, vectors.sim, Time, spacetime=TRUE, pois=TRUE, nsim, YSIM)
    
    #RR and Colors for Plotting
    ##SPACE-ONLY
    initial.s <- list(E0 = unlist(vectors.sim.s$E0_0))
    riskratios.qp.s <- get.rr2(lassoresult.qp.s, vectors.sim.s,initial.s, 
                               tapply(as.vector(matrix(E1, ncol=Time)[,timeperiod]), id, function(x) mean(x)),
                               1, sim=TRUE)
    rrcolors.qp.s <- colormapping(riskratios.qp.s,1)
    
    riskratios.p.s <- get.rr2(lassoresult.qp.s, vectors.sim.s, initial.s, 
                              tapply(as.vector(matrix(E1, ncol=Time)[,timeperiod]), id, function(x) mean(x)),
                              1, sim=TRUE)
    rrcolors.p.s <- colormapping(riskratios.p.s,1)
    
    ##SPACE-TIME
    riskratios.qp.st <- get.rr2(lassoresult.qp.st, vectors.sim,init, E1,Time, sim=TRUE)
    rrcolors.qp.st <- colormapping(riskratios.qp.st,Time)
    
    riskratios.p.st <- get.rr2(lassoresult.qp.st, vectors.sim,init, E1,Time, sim=TRUE)
    rrcolors.p.st <- colormapping(riskratios.p.st,Time)
    
    #COMBINE RISKRATIOS INTO LISTS
    riskratios <- list(riskratios.qp.s = riskratios.qp.s, riskratios.p.s = riskratios.p.s, 
                       riskratios.qp.st = riskratios.qp.st, riskratios.p.st = riskratios.p.st)
    rrcolors <- list(rrcolors.qp.s = rrcolors.qp.s, rrcolors.p.s = rrcolors.p.s,
                     rrcolors.qp.st = rrcolors.qp.st, rrcolors.p.st = rrcolors.p.st)
    
    #DETECTION
    ##QP - Space
    set <- detect_set(lassoresult.qp.s, vectors.sim.s, as.matrix(rr[,timeperiod[1]]), 1, x, y, rMax, center, radius)
    incluster.qp.s <- detect.incluster(lassoresult.qp.s, vectors.sim.s, as.matrix(rr[,timeperiod[1]]), set, 1, 1, nsim, x, y, rMax, center, 
                                       radius, IC = "ic")
    detect.qp.s <- list(clust.diagnostics(incluster.qp.s, threshold[1]), clust.diagnostics(incluster.qp.s , threshold[2]))
    detect.out.qp.s <- (matrix(unlist(detect.qp.s),ncol=3, byrow=TRUE, 
                               dimnames = list(c(paste0("incluster.any.", threshold[1]),
                                                 paste0("alldetect.",threshold[1]), 
                                                 paste0("potentialclusterdetect.",threshold[1]), 
                                                 paste0("trueclusterdetect.",threshold[1]),
                                                 paste0("incluster.any.",threshold[2]), paste0("alldetect.",threshold[2]),
                                                 paste0("potentialclusterdetect.",threshold[2]), 
                                                 paste0("trueclusterdetect.",threshold[2])),c("aic","aicc","bic"))))
    ##P - Space
    set <- detect_set(lassoresult.p.s, vectors.sim.s, as.matrix(rr[,timeperiod[1]]), 1, x, y, rMax, center, radius)
    incluster.p.s <- detect.incluster(lassoresult.p.s, vectors.sim.s, as.matrix(rr[,timeperiod[1]]), set, 1, 1, nsim, x, y, rMax, center, 
                                      radius, IC = "ic")
    detect.p.s <- list(clust.diagnostics(incluster.p.s, threshold[1]), clust.diagnostics(incluster.p.s, threshold[2]))
    detect.out.p.s <- (matrix(unlist(detect.p.s),ncol=3, byrow=TRUE, 
                              dimnames = list(c(paste0("incluster.any.", threshold[1]),
                                                paste0("alldetect.",threshold[1]), 
                                                paste0("potentialclusterdetect.",threshold[1]), 
                                                paste0("trueclusterdetect.",threshold[1]),
                                                paste0("incluster.any.",threshold[2]), paste0("alldetect.",threshold[2]),
                                                paste0("potentialclusterdetect.",threshold[2]), 
                                                paste0("trueclusterdetect.",threshold[2])),c("aic","aicc","bic"))))
    
    
    ##QP - SPACETIME
    set <- detect_set(lassoresult.qp.st, vectors.sim, rr, Time, x, y, rMax, center, radius)
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
    
    ##P - SPACETIME
    set <- detect_set(lassoresult.p.st, vectors.sim, rr, Time, x, y, rMax, center, radius)
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

