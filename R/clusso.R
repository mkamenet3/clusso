#'Detect a cluster in space or spacetime using Lasso 
#'
#' Set up for clusso
#'
#'@title
#'clusso
#' @description
#' Runs helper function for both the space and space-time Lasso model on observed data. 
#'@param df name of dataframe.
#'@param expected Name of variable that contains the expected counts.
#'@param observed Name of variable that contains the observed counts.
#'@param timeperiod Name of variable that contains the timeperiod in which counts were observed (as factor). 
#'If spatial-only analysis, create a column that has a single value (ex: "Time1") and convert this to a factor.
#'@param covars are there additional covariates in the dataframe beyond the three required? If so, set to \code{TRUE}. Default is \code{FALSE}.
#'@param id Optional. If your dataframe contains an ID variable that should not be a covariate, set the name here.
#'@param x x coordinates (easting/latitude); if utm coordinates, scale to km.
#'@param y y coordinates (northing/longitude); if utm coordinates, scale to km.
#'@param rMax set max radius (in km).
#'@param utm default is \code{TRUE}. If \code{FALSE}, then coordinates will be interpreted as Longitude/Latitude and the haversine formula will be used to determine the distance between points.
#'@param longdat Is the data in panel/long format? Default is \code{TRUE}. For wide format, specify \code{FALSE} (TODO).
#'@param analysis A string specifying if the spatial (\code{"space"}), spatio-temporal (\code{"spacetime"}), or both spatial and spatio-temporal (\code{"both"}) analysis should be executed. Default is \code{"both"}. 
#'@param model A string specifying which model to use, Poisson or binomial. For Poisson, specify \code{"poisson"} and both the Poisson and quasi-Poisson model results are returned. For binomial, specify \code{"binomial"}.
#'@param maxclust Upper limit on the maximum number of clusters you expect to find in the region. This equivalent to setting \code{dfmax} in the lasso. If none supplied, default is \code{11}.
#'@param overdispfloor overdispfloor default is \code{TRUE}. When TRUE, it limits \eqn{\phi1} (overdispersion parameter) to be greater or equal to 1. If FALSE, will allow for under-dispersion in the model.
#'@param cv Numeric argument for the number of folds to use if using k-fold cross-validation. Default is \code{NULL}, indicating that cross-validation should not be performed in favor of \code{clusso}.
#'@param collapsetime Default is \code{FALSE}. Alternative definition for space-only model to instead collapse expected and observed counts across time. 
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
#'clst <- toclusso(japanbreastcancer, expected = expdeath, observed=death,timeperiod = period)
#'system.time(resreal <- clusso(clst, x,y, rMax, Time, utm=TRUE, analysis="both", maxclust=10))
#'  }

# clusso <- function(clst, x,y,rMax, Time, utm=TRUE, longdat=TRUE, analysis = c("space","spacetime", "both"),model = c("poisson", "binomial"),maxclust = 11,overdispfloor=TRUE, cv = NULL, collapsetime=FALSE){
clusso2 <- function(df, expected, observed, timeperiod,covars,id= NULL,x,y,rMax, utm=TRUE, longdat=TRUE, analysis = c("space","spacetime", "both"),model = c("poisson", "binomial"),maxclust = 11,overdispfloor=TRUE, cv = NULL, collapsetime=FALSE){
    expected <- eval(substitute(expected),df)
    observed <- eval(substitute(observed),df)
    timeperiod <- eval(substitute(timeperiod),df)
    if(!is.null(id)){
        id<- eval(substitute(id),df)
    }
    else{
        id <- NULL
    }
    if((missing(covars) | covars==FALSE)){
        covars <- FALSE
    }
    else{
        covars <- TRUE
    }
    clst <- toclusso(df, expected, observed,timeperiod, covars, id)
    #if(is(clst, "clst")!=TRUE) stop("clst element not of class `clst`. This is required for the clusso function.")
    expected <- clst$required_df$expected
    observed <- clst$required_df$observed
    period <- clst$required_df$timeperiod
    Time <- length(unique(period))
    #initial user setting
    if(!is.null(clst$othercovariates_df)){
        covars <- clst$othercovariates_df
    }
    else{
        covars <- NULL
    }
    
    if((missing(utm) | utm==TRUE)){
        utm=TRUE
    }
    else{
        message("Coordinates are assumed to be in lat/long coordinates. For UTM coordinates, please specify 'utm=TRUE' or leave empty for default (TRUE).")
        utm=FALSE
    }
    if((missing(longdat) | longdat==TRUE)){
        longdat=TRUE
    }
    else{
        longdat=FALSE
        #message("Data assumed to be in panel data. To use vector data instead, please specify 'longdat=FALSE'")
    }
    if(missing(maxclust)){
        maxclust = 11 + Time
        #print(maxclust)
    }
    else{
        maxclust = maxclust + Time
        #print(maxclust)
    }
    if((missing(overdispfloor) | overdispfloor==TRUE)){
        overdispfloor <- TRUE
    }
    else{
        overdispfloor <- FALSE
    }
    if(!is.null(cv)){
        message("Running cross-validation method for path selection. For AIC, AICc, and BIC, set `cv = NULL`")
        cv = cv
    }
    else{
        cv = NULL
    }
    if(collapsetime==TRUE){
        warning("You've selected to collapse time periods on space. This will results in averaging of counts by geographic unit.
                To allow time to vary across geographic unit, select `collapsetime=FALSE` (DEFAULT).")
    }
    else{
        collapsetime=FALSE
    }
    if(length(model) > 1) stop("You must select either `poisson` or `binomial`")
    model <- match.arg(model)
    if(length(analysis) > 1) stop("You must select either `space`, `spacetime`, or `both`")
    analysis <- match.arg(analysis, several.ok = FALSE)
    if(model=="poisson"){
        #print("Poisson models initiated")
        switch(analysis, 
               space = clussoPois(analysis="space",x, y, rMax,period, expected, observed, covars, 
                                  Time, utm, longdat, maxclust,overdispfloor, cv, collapsetime),
               spacetime = clussoPois(analysis="spacetime",x, y, rMax,period, expected, observed, covars, 
                                      Time, utm, longdat, maxclust,overdispfloor, cv, collapsetime),
               both = clussoPois(analysis="both",x, y, rMax,period, expected, observed, covars, 
                                 Time, utm, longdat, maxclust,overdispfloor, cv, collapsetime))
    }
    else if(model=="binomial"){
        #print("Binomial models initiated")
        switch(analysis, 
               space = clussoBinom(analysis="space",x, y, rMax,period, expected, observed, covars,
                                   Time, utm, longdat, maxclust, cv, collapsetime),
               spacetime = clussoBinom(analysis="spacetime",x, y, rMax,period, expected, observed, covars,
                                       Time, utm, longdat, maxclust, cv, collapsetime),
               both = clussoBinom(analysis="both",x, y, rMax,period, expected, observed, covars, 
                                  Time, utm, longdat, maxclust, cv, collapsetime))
    }
    else{
        warning("You have not specified a model. Please set 'model' argument to 'poisson' or 'binomial'.")
    }
    
    }

clusso <- function(clst, x,y,rMax, Time, utm=TRUE, longdat=TRUE, analysis = c("space","spacetime", "both"),model = c("poisson", "binomial"),maxclust = 11,overdispfloor=TRUE, cv = NULL, collapsetime=FALSE){
    if(is(clst, "clst")!=TRUE) stop("clst element not of class `clst`. This is required for the clusso function.")
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

    if((missing(utm) | utm==TRUE)){
        utm=TRUE
    }
    else{
        message("Coordinates are assumed to be in lat/long coordinates. For UTM coordinates, please specify 'utm=TRUE' or leave empty for default (TRUE).")
        utm=FALSE
    }
    if((missing(longdat) | longdat==TRUE)){
        longdat=TRUE
    }
    else{
        longdat=FALSE
        #message("Data assumed to be in panel data. To use vector data instead, please specify 'longdat=FALSE'")
    }
    if(missing(maxclust)){
        maxclust = 11 + Time
        #print(maxclust)
    }
    else{
        maxclust = maxclust + Time
        #print(maxclust)
    }
    if((missing(overdispfloor) | overdispfloor==TRUE)){
        overdispfloor <- TRUE
    }
    else{
        overdispfloor <- FALSE
    }
    if(!is.null(cv)){
        message("Running cross-validation method for path selection. For AIC, AICc, and BIC, set `cv = NULL`")
        cv = cv
    }
    else{
        cv = NULL
    }
    if(collapsetime==TRUE){
        warning("You've selected to collapse time periods on space. This will results in averaging of counts by geographic unit.
                To allow time to vary across geographic unit, select `collapsetime=FALSE` (DEFAULT).")
    }
    else{
        collapsetime=FALSE
    }
    if(length(model) > 1) stop("You must select either `poisson` or `binomial`")
    model <- match.arg(model)
    if(length(analysis) > 1) stop("You must select either `space`, `spacetime`, or `both`")
    analysis <- match.arg(analysis, several.ok = FALSE)
    if(model=="poisson"){
        print("Poisson models initiated")
        switch(analysis,
               space = clussoPois(analysis="space",x, y, rMax,period, expected, observed, covars,
                                  Time, utm, longdat, maxclust,overdispfloor, cv, collapsetime),
               spacetime = clussoPois(analysis="spacetime",x, y, rMax,period, expected, observed, covars,
                                      Time, utm, longdat, maxclust,overdispfloor, cv, collapsetime),
               both = clussoPois(analysis="both",x, y, rMax,period, expected, observed, covars,
                                 Time, utm, longdat, maxclust,overdispfloor, cv, collapsetime))
    }
    else if(model=="binomial"){
        print("Binomial models initiated")
        switch(analysis,
               space = clussoBinom(analysis="space",x, y, rMax,period, expected, observed, covars,
                                   Time, utm, longdat, maxclust, cv, collapsetime),
               spacetime = clussoBinom(analysis="spacetime",x, y, rMax,period, expected, observed, covars,
                                       Time, utm, longdat, maxclust, cv, collapsetime),
               both = clussoBinom(analysis="both",x, y, rMax,period, expected, observed, covars,
                                  Time, utm, longdat, maxclust, cv, collapsetime))
    }
    else{
        warning("You have not specified a model. Please set 'model' argument to 'poisson' or 'binomial'.")
    }

}

    


#' Detect a cluster in space or spacetime using the LASSO: Poisson model
#' @title
#'clussoPois
#' @description 
#'This function runs both the space and space-time LASSO poisson model. This function is to be run on observed data. A separate function (clusso) is the helper function which will have 
#'flexibility to specify the space or spacetime or both models to be run.
#'@param analysis A string specifying if the spatial (\code{"space"}), spatio-temporal (\code{"spacetime"}), or both spatial and spatio-temporal (\code{"both"}) analysis should be executed. Default is \code{"both"}. 
#'@param x x coordinates (easting/latitude); if utm coordinates, scale to km.
#'@param y y coordinates (northing/longitude); if utm coordinates, scale to km.
#'@param rMax set max radius (in km)
#'@param period vector of periods or years in data set. Should be imported as a factor.
#'@param expected vector of expected counts. Expected counts must match up with the year and observed vectors.
#'@param observed vector of observed counts. Observed counts must match up with the year and expected vectors.
#'@param covars matrix of covariates.
#'@param Time Number of time periods or years in your dataset. Must be declared as numeric.
#'@param utm default is \code{TRUE}. If \code{FALSE}, then coordinates will be interpreted as Longitude/Latitude and the haversine formula will be used to determine the distance between points.
#'@param longdat Is the data in panel/long format? Default is \code{TRUE}. For wide format, specify \code{FALSE} (TODO).
#'@param overdispfloor overdispfloor default is \code{TRUE}. When TRUE, it limits \eqn{\phi1} (overdispersion parameter) to be greater or equal to 1. If FALSE, will allow for under-dispersion in the model.
#'@param maxclust Upper limit on the maximum number of clusters you expect to find in the region. This equivalent to setting \code{dfmax} in the lasso. If none supplied, default is \code{11}.
#'@param cv Numeric argument for the number of folds to use if using k-fold cross-validation. Default is \code{NULL}, indicating that cross-validation should not be performed in favor of \code{clusso}.
#'@param collapsetime Default is \code{FALSE}. Alternative definition for space-only model to instead collapse expected and observed counts across time. 
#'@inheritParams clusso
#'@return list of lists output from detection

clussoPois <- function(analysis,x,y,rMax, period, expected, observed, covars,Time, utm, longdat, maxclust, overdispfloor, cv, collapsetime){
    model <- "poisson"
    if(analysis=="space"){
        analysis_name<-"spatial"
    }
    else if(analysis=="spacetime"){
        analysis_name <- "spatio-temporal" 
    }
    else {
        analysis_name <- "both spatial and spatio-temporal"
    }
    message(paste0("Running Poisson ", analysis_name," model(s)."))
    #set up clusters and fitted values
    clusters <- clusters2df(x,y,rMax, utm=utm, length(x))
    n <- length(x)
    init <- setVectors(period, expected, observed, covars, Time, longdat)
    E1 <- init$E0
    Ex <- clusso::scale(init, Time)
    Yx <- init$Y.vec
    #set vectors
    if(analysis=="spacetime"){
        vectors <- list(Period = init$Year, Ex=Ex, E0_0=init$E0, Y.vec=init$Y.vec, covars = covars)    
    }
    else if (analysis=="space"){
        vectors.s <- list(Period = init$Year, Ex=Ex, E0_0=init$E0, Y.vec=init$Y.vec, covars = covars)    
    }
    else{
        #both
        vectors <- list(Period = init$Year, Ex=Ex, E0_0=init$E0, Y.vec=init$Y.vec, covars = covars)    
        vectors.s <- list(Period = init$Year, Ex=Ex, E0_0=init$E0, Y.vec=init$Y.vec, covars = covars)    
    }
    #create sparseMAT once and cache it   
    n_uniq <- length(unique(clusters$center))
    potClus <- n
    numCenters <- n
    #CREATE sparseMAT and cache it for use throughout this function
    if(collapsetime==FALSE){
        sparseMAT <- spacetimeMat(clusters, numCenters, Time)
        #Create time matrix - not lasso'd
        time_period <- factor(rep(1:Time, each=n_uniq))
        timeMat <- Matrix::Matrix(model.matrix(~ time_period - 1), sparse=TRUE)
        #add this to sparsemat
        sparseMAT <- cbind(sparseMAT, timeMat)
        SOAR::Store(sparseMAT)
        #message("Space-time matrix created")
    }
    else {
        sparseMAT <- spaceMat(clusters, numCenters)
        SOAR::Store(sparseMAT)
        #message("Creating space-only matrix")
        if(nrow(covars)==0){
            covars <- NULL
        }
    }
    #RUN LASSO
    if(analysis=="spacetime"){
        lassoresult.qp.st <- spacetimeLasso(model, sparseMAT, n_uniq, vectors,Time, 
                                            quasi=TRUE, maxclust, overdispfloor, cv)
        lassoresult.p.st <- spacetimeLasso(model, sparseMAT, n_uniq, vectors, Time, 
                                           quasi=FALSE, maxclust, overdispfloor, cv)

        message("All models ran successfully")
        #space time
        ##risk ratios
        riskratios.p.st <- get_rr(lassoresult.p.st, vectors,init,E1,Time, sim=FALSE, cv)
        riskratios.qp.st <- get_rr(lassoresult.qp.st, vectors,init,E1,Time, sim=FALSE, cv)
        ##color mapping
        rrcolors.p.st <- colormapping(riskratios.p.st,Time, cv, prob=FALSE)
        rrcolors.qp.st <- colormapping(riskratios.qp.st,Time, cv, prob = FALSE)
        #COMBINE RISKRATIOS INTO LISTS
        riskratios <- list(riskratios.qp.st = riskratios.qp.st, riskratios.p.st = riskratios.p.st)
        rrcolors <- list(rrcolors.qp.st = rrcolors.qp.st, rrcolors.p.st = rrcolors.p.st)
        return(list(lassoresult.p.st = lassoresult.p.st,
                    lassoresult.qp.st = lassoresult.qp.st,
                    riskratios = riskratios,
                    rrcolors = rrcolors,
                    init.vec = vectors))
    }
    else if (analysis=="space"){
        lassoresult.qp.s <- spacetimeLasso(model, sparseMAT, n_uniq, vectors.s, Time, 
                                           quasi=TRUE, maxclust, overdispfloor,cv)
        lassoresult.p.s <- spacetimeLasso(model, sparseMAT, n_uniq, vectors.s, Time, 
                                          quasi=FALSE, maxclust, overdispfloor,cv)
        message("All models ran successfully")
        #space only
        ##risk ratios
        initial.s <- list(E0 = unlist(vectors.s$E0_0))
        riskratios.p.s <- get_rr(lassoresult.p.s, vectors.s,initial.s,
                                 as.vector(matrix(E1, ncol=Time)),
                                 Time,sim=FALSE, cv)
        riskratios.qp.s <- get_rr(lassoresult.qp.s,vectors.s,initial.s,
                                  as.vector(matrix(E1, ncol=Time)),
                                  Time,sim=FALSE, cv)
        ##color mapping
        rrcolors.p.s <- colormapping(riskratios.p.s,Time, cv, prob=FALSE)
        rrcolors.qp.s <- colormapping(riskratios.qp.s,Time, cv, prob=FALSE)
        #COMBINE RISKRATIOS INTO LISTS
        riskratios <- list(riskratios.qp.s = riskratios.qp.s, riskratios.p.s = riskratios.p.s)
        rrcolors <- list(rrcolors.qp.s = rrcolors.qp.s, rrcolors.p.s = rrcolors.p.s)
        return(list(lassoresult.p.s = lassoresult.p.s,
                    lassoresult.qp.s = lassoresult.qp.s,
                    riskratios = riskratios,
                    rrcolors = rrcolors,
                    init.vec.s = vectors.s))
    }
    else{
        #both
        lassoresult.qp.st <- spacetimeLasso(model, sparseMAT, n_uniq, vectors,Time, 
                                            quasi=TRUE, maxclust, overdispfloor, cv)
        lassoresult.p.st <- spacetimeLasso(model, sparseMAT, n_uniq, vectors, Time, 
                                           quasi=FALSE, maxclust, overdispfloor, cv)
        lassoresult.qp.s <- spacetimeLasso(model, sparseMAT, n_uniq, vectors.s, Time, 
                                           quasi=TRUE, maxclust, overdispfloor,cv)
        lassoresult.p.s <- spacetimeLasso(model,sparseMAT, n_uniq, vectors.s, Time, 
                                          quasi=FALSE, maxclust, overdispfloor,cv)
        message("All models ran successfully")
        #space time
        ##risk ratios
        riskratios.p.st <- get_rr(lassoresult.p.st, vectors,init,E1,Time, sim=FALSE, cv)
        riskratios.qp.st <- get_rr(lassoresult.qp.st, vectors,init,E1,Time, sim=FALSE, cv)
        ##color mapping
        rrcolors.p.st <- colormapping(riskratios.p.st,Time, cv, prob=FALSE)
        rrcolors.qp.st <- colormapping(riskratios.qp.st,Time, cv, prob = FALSE)
        #space only
        ##risk ratios
        initial.s <- list(E0 = unlist(vectors.s$E0_0))
        riskratios.p.s <- get_rr(lassoresult.p.s, vectors.s,initial.s,
                                 as.vector(matrix(E1, ncol=Time)),
                                 Time,sim=FALSE, cv)
        riskratios.qp.s <- get_rr(lassoresult.qp.s,vectors.s,initial.s,
                                  as.vector(matrix(E1, ncol=Time)),
                                  Time,sim=FALSE, cv)
        ##color mapping
        rrcolors.p.s <- colormapping(riskratios.p.s,Time, cv, prob=FALSE)
        rrcolors.qp.s <- colormapping(riskratios.qp.s,Time, cv, prob=FALSE)
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
    SOAR::Remove(sparseMAT)
    #add warning about maxclust
    numclust.AIC <- c(lassoresult.p.s$numclust.qaic, lassoresult.p.st$numclust.qaic, 
                      lassoresult.qp.s$numclust.qaic, lassoresult.qp.st$numclust.qaic)
    numclust.AICc <- c(lassoresult.p.s$numclust.qaicc, lassoresult.p.st$numclust.qaicc, 
                       lassoresult.qp.s$numclust.qaicc, lassoresult.qp.st$numclust.qaicc)
    numclust.BIC <- c(lassoresult.p.s$numclust.qbic, lassoresult.p.st$numclust.qbic, 
                      lassoresult.qp.s$numclust.qbic, lassoresult.qp.st$numclust.qbic)
    numclust.cv <-  c(lassoresult.p.s$numclust.cv, lassoresult.p.st$numclust.cv, 
                      lassoresult.qp.s$numclust.cv, lassoresult.qp.st$numclust.cv)
    nclusters <- c(numclust.AIC, numclust.AICc, numclust.BIC,numclust.cv)
    if(any(nclusters==maxclust)){
        warning("Number of clusters detected by at least one criterion equals maxclust. Consider increasing maxclust.")
    }
}


#' Detect a cluster in space or spacetime using the LASSO: Binomial model
#' @title
#'clussoBinom
#' @description 
#'This function runs both the space and space-time LASSO binomial model. This function is to be run on observed data. A separate function (clusso) is the helper function which will have 
#'flexibility to specify the space or spacetime or both models to be run .
#'@param analysis A string specifying if the spatial (\code{"space"}), spatio-temporal (\code{"spacetime"}), or both spatial and spatio-temporal (\code{"both"}) analysis should be executed. Default is \code{"both"}. 
#'@param x x coordinates (easting/latitude); if utm coordinates, scale to km.
#'@param y y coordinates (northing/longitude); if utm coordinates, scale to km.
#'@param rMax set max radius (in km)
#'@param period vector of periods or years in data set. Should be imported as a factor.
#'@param numcases vector of number of cases (or successes).
#'@param n Total number of both cases and controls (or number of trials)
#'@param covars matrix of covariates.
#'@param Time Number of time periods or years in your dataset. Must be declared as numeric.
#'@param utm default is \code{TRUE}. If \code{FALSE}, then coordinates will be interpreted as Longitude/Latitude and the haversine formula will be used to determine the distance between points.
#'@param longdat Is the data in panel/long format? Default is \code{TRUE}. For wide format, specify \code{FALSE} (TODO).
#'@param maxclust Upper limit on the maximum number of clusters you expect to find in the region. This equivalent to setting \code{dfmax} in the lasso. If none supplied, default is \code{11}.
#'@param cv Numeric argument for the number of folds to use if using k-fold cross-validation. Default is \code{NULL}, indicating that cross-validation should not be performed in favor of \code{clusso}.
#'@param collapsetime Default is \code{FALSE}. Alternative definition for space-only model to instead collapse expected and observed counts across time. 
#'@inheritParams clusso
#'@return list of lists output from detection

clussoBinom <- function(analysis,x,y,rMax, period, expected, observed, covars,Time, utm, longdat, maxclust,cv, collapsetime){  
    model <- "binomial"
    if(analysis=="space"){
        analysis_name<-"spatial"
    }
    else if(analysis=="spacetime"){
        analysis_name <- "spatio-temporal" 
    }
    else {
        analysis_name <- "both spatial and spatio-temporal"
    }
    message(paste0("Running binomial ", analysis_name," model(s)."))
    #set up clusters and fitted values
    clusters <- clusters2df(x,y,rMax, utm=utm, length(x))
    n_uniq <- length(unique(clusters$center))
    init <- setVectors(period, expected, observed, covars, Time, longdat) #TODO but ntrials=E0 and Y.vec=cases
    Yx <- init$Y.vec #ncases
    Ex <- init$E0 #ntrials
    #set vectors
    if(analysis=="spacetime"){
        vectors <- list(Period = init$Year, Ex=Ex, E0_0=init$E0, Y.vec=init$Y.vec, covars = covars)    
    }
    else if (analysis=="space"){
        vectors.s <- list(Period = init$Year, Ex=Ex, E0_0=init$E0, Y.vec=init$Y.vec, covars = covars)    
    }
    else{
        #both
        vectors <- list(Period = init$Year, Ex=Ex, E0_0=init$E0, Y.vec=init$Y.vec, covars = covars)    
        vectors.s <- list(Period = init$Year, Ex=Ex, E0_0=init$E0, Y.vec=init$Y.vec, covars = covars)    
    }
    #create sparseMAT once and cache it   
    potClus <- n_uniq
    numCenters <- n_uniq
    #CREATE sparseMAT and cache it for use throughout this function
    if(collapsetime==FALSE){
        sparseMAT <- spacetimeMat(clusters, numCenters, Time)
        #Create time matrix - not lasso'd
        time_period <- factor(rep(1:Time, each=n_uniq))
        timeMat <- Matrix::Matrix(model.matrix(~ time_period - 1), sparse=TRUE)
        #add this to sparsemat
        sparseMAT <- cbind(sparseMAT, timeMat)
        SOAR::Store(sparseMAT)
        #message("Space-time matrix created")
    }
    else {
        sparseMAT <- spaceMat(clusters, numCenters)
        SOAR::Store(sparseMAT)
        #message("Creating space-only matrix")
        if(nrow(covars)==0){
            covars <- NULL
        }
    }
    #RUN LASSO
    if(analysis=="spacetime"){
        lassoresult.p.st <- spacetimeLasso(model, sparseMAT, n_uniq, vectors, Time, 
                                           quasi=FALSE, maxclust, overdispfloor, cv)
        lassoresult.qp.st <- spacetimeLasso(model, sparseMAT, n_uniq, vectors, Time,
                                            quasi=TRUE, maxclust, overdispfloor, cv)
        message("All models ran successfully")
        #space time
        ##risk ratios
        riskratios.p.st <- get_rr(lassoresult.p.st, vectors,init,E1,Time, sim=FALSE, cv)
        riskratios.qp.st <- get_rr(lassoresult.qp.st, vectors,init,E1,Time, sim=FALSE, cv)
        ##color mapping
        rrcolors.p.st <- colormapping(riskratios.p.st,Time, cv, prob=FALSE)
        rrcolors.qp.st <- colormapping(riskratios.qp.st,Time, cv, prob = FALSE)
        #COMBINE RISKRATIOS INTO LISTS
        riskratios <- list(riskratios.qp.st = riskratios.qp.st, riskratios.p.st = riskratios.p.st)
        rrcolors <- list(rrcolors.qp.st = rrcolors.qp.st, rrcolors.p.st = rrcolors.p.st)
        return(list(lassoresult.p.st = lassoresult.p.st,
                    lassoresult.qp.st = lassoresult.qp.st,
                    riskratios = riskratios,
                    rrcolors = rrcolors,
                    init.vec = vectors))
    }
    else if (analysis=="space"){
        lassoresult.p.s <- spacetimeLasso(model, sparseMAT, n_uniq, vectors.s, Time, 
                                          quasi=FALSE, maxclust, overdispfloor,cv)
        lassoresult.qp.s <- spacetimeLasso(model, sparseMAT, n_uniq, vectors.s, Time, 
                                           quasi=TRUE, maxclust, overdispfloor,cv)
        message("All models ran successfully")
        #space only
        ##risk ratios
        initial.s <- list(E0 = unlist(vectors.s$E0_0))
        riskratios.p.s <- get_rr(lassoresult.p.s, vectors.s,initial.s,
                                 as.vector(matrix(E1, ncol=Time)),
                                 Time,sim=FALSE, cv)
        riskratios.qp.s <- get_rr(lassoresult.qp.s,vectors.s,initial.s,
                                  as.vector(matrix(E1, ncol=Time)),
                                  Time,sim=FALSE, cv)
        ##color mapping
        rrcolors.p.s <- colormapping(riskratios.p.s,Time, cv, prob=FALSE)
        rrcolors.qp.s <- colormapping(riskratios.qp.s,Time, cv, prob=FALSE)
        #COMBINE RISKRATIOS INTO LISTS
        riskratios <- list(riskratios.qp.s = riskratios.qp.s, riskratios.p.s = riskratios.p.s)
        rrcolors <- list(rrcolors.qp.s = rrcolors.qp.s, rrcolors.p.s = rrcolors.p.s)
        return(list(lassoresult.p.s = lassoresult.p.s,
                    lassoresult.qp.s = lassoresult.qp.s,
                    riskratios = riskratios,
                    rrcolors = rrcolors,
                    init.vec.s = vectors.s))
    }
    else{
        #both

        lassoresult.p.st <- spacetimeLasso(model, sparseMAT, n_uniq, vectors,Time, 
                                           quasi=FALSE, maxclust, overdispfloor, cv)
        lassoresult.qp.st <- spacetimeLasso(model, sparseMAT, n_uniq, vectors, Time, 
                                            quasi=TRUE, maxclust, overdispfloor, cv)
        lassoresult.p.s <- spacetimeLasso(model, sparseMAT, n_uniq, vectors.s, Time, 
                                          quasi=FALSE, maxclust, overdispfloor,cv)
        lassoresult.qp.s <- spacetimeLasso(model, sparseMAT, n_uniq, vectors.s, Time, 
                                           quasi=TRUE, maxclust, overdispfloor,cv)
        message("All models ran successfully")
        #space time
        ##risk ratios
        riskratios.p.st <- get_rr(lassoresult.p.st, vectors,init,E1,Time, sim=FALSE, cv)
        riskratios.qp.st <- get_rr(lassoresult.qp.st, vectors,init,E1,Time, sim=FALSE, cv)
        ##color mapping
        rrcolors.p.st <- colormapping(riskratios.p.st,Time, cv, prob=FALSE)
        rrcolors.qp.st <- colormapping(riskratios.qp.st,Time, cv, prob = FALSE)
        #space only
        ##risk ratios
        initial.s <- list(E0 = unlist(vectors.s$E0_0))
        riskratios.p.s <- get_rr(lassoresult.p.s, vectors.s,initial.s,
                                 as.vector(matrix(E1, ncol=Time)),
                                 Time,sim=FALSE, cv)
        riskratios.qp.s <- get_rr(lassoresult.qp.s,vectors.s,initial.s,
                                  as.vector(matrix(E1, ncol=Time)),
                                  Time,sim=FALSE, cv)
        ##color mapping
        rrcolors.p.s <- colormapping(riskratios.p.s,Time, cv, prob=FALSE)
        rrcolors.qp.s <- colormapping(riskratios.qp.s,Time, cv, prob=FALSE)
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
    SOAR::Remove(sparseMAT)
    #add warning about maxclust
    numclust.AIC <- c(lassoresult.p.s$numclust.qaic, lassoresult.p.st$numclust.qaic, 
                      lassoresult.qp.s$numclust.qaic, lassoresult.qp.st$numclust.qaic)
    numclust.AICc <- c(lassoresult.p.s$numclust.qaicc, lassoresult.p.st$numclust.qaicc, 
                       lassoresult.qp.s$numclust.qaicc, lassoresult.qp.st$numclust.qaicc)
    numclust.BIC <- c(lassoresult.p.s$numclust.qbic, lassoresult.p.st$numclust.qbic, 
                      lassoresult.qp.s$numclust.qbic, lassoresult.qp.st$numclust.qbic)
    numclust.cv <-  c(lassoresult.p.s$numclust.cv, lassoresult.p.st$numclust.cv, 
                      lassoresult.qp.s$numclust.cv, lassoresult.qp.st$numclust.cv)
    nclusters <- c(numclust.AIC, numclust.AICc, numclust.BIC,numclust.cv)
    if(any(nclusters==maxclust)){
        warning("Number of clusters detected by at least one criterion equals maxclust. Consider increasing maxclust.")
    }
    
    
}



