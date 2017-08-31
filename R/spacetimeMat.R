#' Spatio-temporal Sparse Matrix
#' 
#' This function takes the Kronecker product of the space and time matrices to create the space-time matrix
#' @param clusters the potential clusters dataframe created by *clusters2df*
#' @param numCenters number of centroids in space matrix
#' @param Time number of time periods
#' @return Returns sparse space time matrix     
#' @export     
#' @examples 
#' spacetimeMat(clusters, numCenters, Time)
spacetimeMat <- function(clusters, numCenters, Time){
    space <- spaceMat(clusters, numCenters)
    time <- timeMat(Time)
    spacetimeMatrix <- kronecker(time, space)
    return(spacetimeMatrix)
}



#' timeMat
#' 
#' This function creates a sparse matrix of 1's of all of the potential time periods for the cluster to be in. The number of 
#' potential time periods is determined by [(Time*(Time-1)]/2.
#' @param Time Number of time periods in the data
#' @return returns sparse matrix of 1's as indicators of membership in the time cluster
#' @export
timeMat <-function(Time){
    block <- Matrix::Matrix(diag(1,Time),sparse=TRUE)
    master <- block
    for(i in 1:(Time-2)){
        diag(block[(i+1):Time,])<-1
        master <- Matrix::cBind(master, block[,1:(Time-i)])        
    }
    master <- Matrix::cBind(master, Matrix::Matrix(rep(1,Time)))
    return(master)
}