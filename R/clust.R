#'Detect a cluster in space or spacetime using Lasso 
#'
#'
#'
#'@title
#'vectors_space
#' 
#' @description 
#' This function will collapse a space-time vector onto space only
#' @param x vector coordinates (unique regardless of time period)
#' @param Ex list of simulated and standardized expected counts
#' @param Yx observed
#' @param Time number of time periods
#' @param init initial list of vectors, inherited from function setVectors.
#' @return returns space-time 
#' 
vectors_space <- function(x,Ex, Yx,Time, init){
    id <- rep(1:length(x), times = Time)
    if(length(id)!=length(as.vector(Ex))){
        stop("Length of ID var not equal to number of observations")
    }
    if(!is.null(init$covars)){
        covars.s <- sapply(1:ncol(init$covars), 
                       function(i) tapply(as.vector(matrix(init$covars[,i], ncol=Time)),id, 
                                          function(x) sum(x)))
    }
    else{
        covars.s <- NULL
    }
    vectors.s <- list(Period = rep("1", length(x)),
                      Ex = tapply(as.vector(matrix(Ex, ncol=Time)), id, function(x) sum(x)),
                      E0_0 = tapply(as.vector(matrix(init$E0, ncol=Time)), id, function(x) sum(x)),
                      Y.vec = tapply(as.vector(matrix(init$Y.vec, ncol=Time)), id, function(x) round(sum(x))),
                      covars.s = as.data.frame(covars.s))
    Yx.s <- tapply(as.vector(matrix(Yx, ncol=Time)), id, function(x) round(sum(x))) 
    res <- list(vectors.s = vectors.s, Yx.s = Yx.s)
    return(res)
}

#' Set up for clust
#'
#'@title
#'clust
#' @description
#' Runs helper function for both the space and space-time Lasso model on observed data. A separate function (clust_sim) can be used for simulating data and running diagnostics on simulations.
#'@param clst list; output from toclust function. Must be of class clst.
#'@param x x coordinates (easting/latitude); if utm coordinates, scale to km.
#'@param y y coordinates (northing/longitude); if utm coordinates, scale to km.
#'@param rMax set max radius (in km)
#'@param Time Number of time periods or years in your dataset. Must be declared as numeric.
#'@param utm default is TRUE. If FALSE, then will run long/lat data
#'@param byrow default is True. If data should be imported by column then set to FALSE
#'@param space space and space-time. Default is to run all four models: quasi-poisson and poisson for both space and space-time. User can specify, space = space,
#'space = spacetime, or space = both.
#'@param overdispfloor overdispfloor default is TRUE. When TRUE, it limits phi (overdispersion parameter) to be greater or equal to 1. If FALSE, will allow for under dispersion.
#'@param cv option for cross-validation
#'@export
#'@return returns list of cluster detection results ready to plot
#'@examples
#'\donttest{
#'data(japanbreastcancer)
#'#Set Some Initial Conditions
#'x1=utmJapan$utmx/1000
#'y1=utmJapan$utmy/1000
#'rMax <- 20 
#'Time=5
#'japanbreastcancer <- japanbreastcancer[,-1] #get rid of indicator column
#'clst <- toclust(japanbreastcancer, expected = japanbreastcancer$expdeath, 
#'  observed = japanbreastcancer$death,timeperiod = japanbreastcancer$period, covars = FALSE)
#'system.time(res <- clust(clst, x1,y1, rMax, Time, utm=TRUE, byrow=TRUE, 
#'  space="both", overdispfloor=TRUE, cv=NULL))
#'  }


clust <- function(clst, x,y,rMax, Time, utm=TRUE, byrow=TRUE,space = c("space","spacetime", "both"),overdispfloor=NULL, cv = NULL){
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
    if(!is.null(cv)){
        message("Running cross-validation method for path selection. For AIC, AICc, and BIC, set `cv = NULL`")
        cv = cv
    }
    else{
        cv = NULL
    }
    if(length(space) > 1) stop("You must select either `space`, `spacetime`, or `both`")
    space <- match.arg(space, several.ok = FALSE)
    switch(space, 
           #TODO
           # space = 
           # spacetime = 
           both = clustAll(x, y, rMax,period, expected, observed, covars, Time, utm, byrow, overdispfloor, cv))
}
    
#' Detect a cluster in space or spacetime using Lasso on observed data    
#' @title
#'clustAll
#' @description 
#'This function runs both the space and space-time Lasso model. This function is to be run on observed data. A separate function (clust) is the helper function which will have 
#'flexibility to specify the space or spacetime or both models to be run (TODO).
#'@param x x coordinates (easting/latitude); if utm coordinates, scale to km.
#'@param y y coordinates (northing/longitude); if utm coordinates, scale to km.
#'@param rMax set max radius (in km)
#'@param period vector of periods or years in dataset. Should be imported as a factor.
#'@param expected vector of expected counts. Expected counts must match up with the year and observed vectors.
#'@param observed vector of observed counts. Observed counts must match up with the year and expected vectors.
#'@param covars matrix of covariates.
#'@param Time Number of time periods or years in your dataset. Must be declared as numeric.
#'@param utm default is TRUE. If FALSE, then will run long/lat data
#'@param byrow default is True. If data should be imported by column then set to FALSE
#'@param overdispfloor overdispfloor default is TRUE. When TRUE, it limits phi (overdispersion parameter) to be greater or equal to 1. If FALSE, will allow for under dispersion.
#'@param cv option for cross-validation - numeric input specifies how many folds; default is 10
#'@inheritParams clust
#'@return list of output from detection
    
clustAll <- function(x,y,rMax, period, expected, observed, covars,Time, utm, byrow, overdispfloor, cv){    
    message("Running both Space and Space-Time Models")
    
    #set up clusters and fitted values
    clusters <- clusters2df(x,y,rMax, utm=utm, length(x))
    n <- length(x)
    init <- setVectors(period, expected, observed, covars, Time, byrow)
    E1 <- init$E0
    Ex <- scale(init, Time)
    Yx <- init$Y.vec
    #timeperiod <- 1:Time
    #set vectors
    vectors <- list(Period = init$Year, Ex=Ex, E0_0=init$E0, Y.vec=init$Y.vec, covars = covars)
    vectors.s <- list(Period = init$Year, Ex=Ex, E0_0=init$E0, Y.vec=init$Y.vec, covars = covars)
    #spacevecs <- vectors_space(x, Ex, Yx, Time,init)
    #vectors.s <- spacevecs$vectors.s
    
    #run lasso
    lassoresult.p.st <- spacetimeLasso(clusters, vectors,Time, spacetime=TRUE,pois=TRUE, overdispfloor, cv)
    lassoresult.qp.st <- spacetimeLasso(clusters, vectors, Time, spacetime=TRUE,pois=FALSE, overdispfloor, cv)
    lassoresult.p.s <- spacetimeLasso(clusters, vectors.s, Time, spacetime=TRUE,pois=TRUE, overdispfloor,cv)
    lassoresult.qp.s <- spacetimeLasso(clusters, vectors.s, Time, spacetime=TRUE,pois=FALSE, overdispfloor,cv)
    
    message("All models ran successfully")
    
    
    #space time
    ##risk ratios
    riskratios.p.st <- get_rr(lassoresult.p.st, vectors,init,E1,Time, sim=FALSE, cv)
    riskratios.qp.st <- get_rr(lassoresult.qp.st, vectors,init,E1,Time, sim=FALSE, cv)
    ##color mapping
    rrcolors.p.st <- colormapping(riskratios.p.st,Time, cv)
    rrcolors.qp.st <- colormapping(riskratios.qp.st,Time, cv)
    
    #space only
    ##risk ratios
    initial.s <- list(E0 = unlist(vectors.s$E0_0))
    id <- rep(1:length(x), times=Time)
    riskratios.p.s <- get_rr(lassoresult.p.s, vectors.s,initial.s,
                              tapply(as.vector(matrix(E1, ncol=Time)), id, function(x) mean(x)),
                              1,sim=FALSE, cv)
    riskratios.qp.s <- get_rr(lassoresult.qp.s,vectors.s,initial.s,
                               tapply(as.vector(matrix(E1, ncol=Time)), id, function(x) mean(x)),
                               1,sim=FALSE, cv)
    ##color mapping
    rrcolors.p.s <- colormapping(riskratios.p.s,1, cv)
    rrcolors.qp.s <- colormapping(riskratios.qp.s,1, cv)
    
    #COMBINE RISKRATIOS INTO LISTS
    riskratios <- list(riskratios.qp.s = riskratios.qp.s, riskratios.p.s = riskratios.p.s, 
                       riskratios.qp.st = riskratios.qp.st, riskratios.p.st = riskratios.p.st)
    rrcolors <- list(rrcolors.qp.s = rrcolors.qp.s, rrcolors.p.s = rrcolors.p.s,
                     rrcolors.qp.st = rrcolors.qp.st, rrcolors.p.st = rrcolors.p.st)
    
    return(list(lassoresult.p.st = lassoresult.p.st,
                lassoresult.qp.st = lassoresult.qp.st,
                lassoresult.p.s = lassoresult.p.s,
                lassoresult.qp.s = lassoresult.qp.s,
                riskratios = riskratios,
                rrcolors = rrcolors,
                init.vec = vectors,
                init.vec.s = vectors.s))
}




