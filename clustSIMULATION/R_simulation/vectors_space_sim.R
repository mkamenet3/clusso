#'Collapse space-tie onto space by summing counts over time by spatial unit - simulation.
#'
#' @title
#' vectors_space_sim
#' @description 
#' Simulation functions for space and space-time cluster detection with lasso. This function will collapse a space-time vector onto space only. Function used for simulation.
#' @param x vector coordinates (unique regardless of time period)
#' @param Ex list of simulated and standardized expected counts
#' @param YSIM simulated observed
#' @param Time number of time periods
#' @param init initial list of vectors, inherited from function setVectors.
#' @return returns returns spatial units with counts aggregated over time periods.
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
