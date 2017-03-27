#'Detect a cluster in space or spacetime using Lasso 
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
#'@param pois default is FALSE (default is to run Quasi-Poisson)
#'@param utm default is TRUE. If FALSE, then will run long/lat data
#'@param byrow default is True. If data should be imported by column then set to FALSE
#'@return
#'@details Optional functions include:
#'- 1) utm - default is FALSE. If you have utm coordinates, you want to change this to TRUE.
#'@examples 
#'res <- clust(x,y,rMax, dframe$period, dframe$expdeath, dframe$death, Time, spacetime=TRUE, pois=FALSE, utm=TRUE, byrow=TRUE)
#'@export

clust <- function(x, y, rMax, period, expected, observed, Time, spacetime=TRUE, pois = FALSE,utm=TRUE, byrow=TRUE,...){
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
    clusters <- clusters2df(x,y,rMax, utm=utm, length(x))
    n <- length(x)
    #init <- setVectors(dframe$period, dframe$expdeath, dframe$death, Time=Time, byrow=row)
    init <- setVectors(period, expected, observed, Time, row)
    E0 <- scale(init, Time, scaler = 1)
    vectors <- list(Period = init$Year, E0=E0, E0_0=init$E0, Y.vec=init$Y.vec)
    lassoresult <- spacetimeLasso(clusters, vectors, Time, spacetime=spacetime,pois=pois)
    riskratios <- get.rr(lassoresult, vectors, Time, sim=FALSE)
    rrcolors <- colormapping(riskratios,Time)
    return(list(lassoresult = lassoresult,
                riskratios = riskratios,
                rrcolors = rrcolors,
                init.vec = init))
}


