#' Create the clusters dataframe
#' 
#' @param xP vector of x coordinates (easting/latitude); if utm coordinates, scale to km.
#' @param yP vector of y coordinates (northing/longitude); if utm coordinates, scale to km.
#' @param r.max Maximum radius for potential clusters (in km).
#' @param utm TRUE/FALSE as to whether or not the x and y coordinates are in UTM (TRUE) or LAT/LONG(FALSE).
#' @param n Number of coordinate pairs/number of centers.
#' @return This function returns a dataframe that contains. 
#' @examples
#' rMax = 20
#' #utm example
#' x_utm <- c(399786.6, 360917.0, 385175.1, 371603.4, 388154.2, 375023.3)
#' y_utm <- c(4047756, 4023885, 4025749, 4018172, 4047900, 4068053)
#' #lat/long example
#' x_latlon <- c(36.569996, 36.350001, 36.370002, 36.299997, 36.570002, 36.749997)
#' y_latlon<- c(7.88, 7.45, 7.72, 7.57, 7.57, 7.60)
#' clusters2df(x_utm, y_utm, rMax, utm=TRUE, length(x_utm))
#' clusters2df(x_latlon, y_latlon, rMax, utm=FALSE, length(x_latlon))

clusters2df <- function(xP,yP, r.max, utm=FALSE,n){
    message("Creating radius-based potential clusters")
    indR = (1:n)[!duplicated(cbind(xP,yP))] 
    if(utm==FALSE){
        if(mean(nchar(vapply(strsplit(as.character(xP), "[.]"),"[", 1, FUN.VALUE=character(1))))>2 & 
           mean(nchar(vapply(strsplit(as.character(yP), "[.]"),"[", 1, FUN.VALUE=character(1))))>2)
            message("Your coordinates may be in UTM due to character length. Please double check.")
        tmpR <- (as.matrix(geosphere::distm(cbind(xP, yP), fun=geosphere::distHaversine))[indR,])/1000    
    } 
    else{
        if(mean(nchar(vapply(strsplit(as.character(xP), "[.]"),"[", 1, FUN.VALUE=character(1))))<=2 & 
           mean(nchar(vapply(strsplit(as.character(yP), "[.]"),"[", 1, FUN.VALUE=character(1))))<=2)
            message("Your coordinates may be Lat/Long due to character length. Please double check.")
        tmpR = as.matrix(dist(cbind(xP,yP)))[indR,]
    }
    lastR = apply(tmpR, 1, function(x,r) order(x)[1:sum(x<=r)],r=r.max)
    ncR = unlist(lapply(lastR, length))
    lastR = unlist(lastR)
    rR=unlist(apply(tmpR,1, function(x,r) { sort(x[x<=r]) },r=r.max))
    
    clustersR=data.frame(center=rep(indR,ncR),
                         x=xP[rep(indR,ncR)],y=yP[rep(indR,ncR)],
                         r=rR, 
                         n=unlist(lapply(ncR,seq)),
                         last=lastR)    
    return(clustersR)
}



