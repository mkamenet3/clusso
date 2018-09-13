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
