#' Creating sparse spatial matrix based on clusters dataframe
#' 
#' This function creates a sparse matrix of 1's of all potential clusters for the Lasso algorithm to cycle over; this will take in the potential clusters dataframe 
#' created by ```clusters2df()``` function
#' @param clusters clusters dataframe from (cluster.df function) that includes the center, x,y, r (radius), n (counter), and last (last observation in potential cluster)
#' @param numCenters the number of centers/centroids
#' @return returns sparse matrix of 1's

spaceMat <- function(clusters, numCenters){
    potClus <- numCenters
    mymat <- NULL
    for(i in 1:nrow(clusters)){
        myvec <- list(as(max_colCpp(numCenters, i, clusters$n, clusters$last), "sparseVector")) 
       # myvec <- list(sparseVector(max_colCpp(numCenters, i, clusters$n, clusters$last))) 
        mymat <- c(mymat,myvec) 
    }
    xx <- NULL
    jj <- NULL
    ii <- NULL
    for(k in 1:length(mymat)){
        xx<- c(xx, mymat[[k]]@x)
        jj <- c(jj, mymat[[k]]@i)
        ii <- c(ii, rep(k,length(mymat[[k]]@x)))
    }
    return(Matrix::t(Matrix::sparseMatrix(i = ii, j = jj, x =xx, dims = c(length(mymat), numCenters))))
}
