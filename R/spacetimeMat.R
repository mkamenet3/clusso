#' Spatio-temporal Sparse Matrix
#' 
#' This function takes the Kronecker product of the space and time matrices to create the space-time matrix.
#' @param clusters Clusters dataframe that includes the center, x,y, r (radius), n (counter), and last (last observation in potential cluster).
#' @param numCenters Number of geographic regions/centroids.
#' @param Time Number of timeperiods in the dataset. This is calculated based on the number of unique timeperiods (factor levels) supplied to \code{clusso}.
#' @return Returns sparse space time matrix.     
#' @export     
spacetimeMat <- function(clusters, numCenters, Time){
    space <- spaceMat(clusters, numCenters)
    time <- timeMat(Time)
    spacetimeMatrix <- kronecker(time, space)
    return(spacetimeMatrix)
}

#' timeMat
#' 
#' This function creates a sparse matrix of 1's indicating  all of the time periods a potential cluster could be in. The number of 
#' potential time periods is determined by \code{[(Time*(Time-1))]/2}.
#' @param Time Number of timeperiods in the dataset. This is calculated based on the number of unique timeperiods (factor levels) supplied to \code{clusso}.
#' @return Returns sparse matrix of 1's as indicators of membership in the time dimension.
#' @export
timeMat <-function(Time){
    block <- Matrix::Matrix(diag(1,Time),sparse=TRUE)
    master <- block
    if(Time==2){
        master <- cbind(master, Matrix::Matrix(rep(1,Time), sparse=TRUE))
    }
    else if(Time==1){
       master <- Matrix::Matrix(diag(1,Time),sparse=TRUE)
    }
    else{
        for(i in 1:(Time-2)){
            diag(block[(i+1):Time,]) <-1
            master <- cbind(master, block[,1:(Time-i)])        
        }
        master <- cbind(master, Matrix::Matrix(rep(1,Time), sparse=TRUE))
    }
    return(master)
}
