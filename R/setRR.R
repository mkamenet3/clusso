#' These functions will extract the risk ratios from the output
#' 
#' set_rr
#' 
#' This function will create vectors of the risk ratios as determined by observed counts, QBIC, QAIC, and QAICc, respectively.
#' @param lassoresult List of QBIC, QAIC, QAICc estimates from the mylasso function
#' @param vectors dataframe of initial vectors of the observed and expected counts
#' @param Time number of time period
#' @param sim default is FALSE. If True, instead of returning observed over expected counts it will return the oracle.
#' @return This returns a list of the risk ratios (observed over expected) as determined by 1) pure observed/expected counts,
#' 2) observed based on QBIC path/expected; 3) observed based on QAIC path/expected; 4) observed based on QAICc path/expected.
#' @export
set_rr<- function(lassoresult, vectors, Time, sim=FALSE){
    if(sim==FALSE){
        RRobs <- matrix(as.vector(vectors$Y.vec)/as.vector(vectors$E0),ncol=Time)
        RRbic <- matrix(lassoresult$E.qbic/as.vector(vectors$E0),ncol=Time)
        RRaic <- matrix(lassoresult$E.qaic/as.vector(vectors$E0),ncol=Time)
        RRaicc <- matrix(lassoresult$E.qaicc/as.vector(vectors$E0),ncol=Time)
        message("Relative risks from observed data")
    }
    else{
        E0_avg <- Reduce("+", vectors$E0)/length(vectors$E0)
        RRobs <- matrix(as.vector(E0_avg)/as.vector(vectors$E0_fit),ncol=Time)
        RRbic <- matrix(lassoresult$E.qbic/as.vector(vectors$E0_fit),ncol=Time)
        RRaic <- matrix(lassoresult$E.qaic/as.vector(vectors$E0_fit),ncol=Time)
        RRaicc <- matrix(lassoresult$E.qaicc/as.vector(vectors$E0_fit),ncol=Time) 
        message("Relative risks from simulated data")
    }
    return(list(RRobs=RRobs, RRbic=RRbic, RRaic=RRaic, RRaicc=RRaicc))  
}


#' get_rr
#' 
#' This function will extract the relative risk ratios as determined by observed counts, QBIC, QAIC, and QAICc, respectively. This method determines the 
#' relative risk by averageing the relative risk in each simulation over the number of simulations. This is in contrast to $setRR$ where the expected counts
#' are first average over the number of simulations and then we determine the relative risk to the initial expected counts. Both methods return the same results,
#' however this method returns a smoother background ratio that does not vary across time periods as much
#' @param lassoresult List of QBIC, QAIC, QAICc estimates from the mylasso function
#' @param vectors.sim dataframe of initial vectors of the observed and expected counts
#' @param init initial settings
#' @param E1 initial expected counts
#' @param Time number of time period
#' @param sim default is TRUE. Signals whether these values are to be extracted from a simulated model or real data.
#' @return This returns a list of the risk ratios (observed over expected) as determined by 1) pure observed/expected counts,
#' 2) observed based on QBIC path/expected; 3) observed based on QAIC path/expected; 4) observed based on QAICc path/expected.
#' @export
#' 
get_rr <- function(lassoresult,vectors.sim,init,E1, Time, sim=TRUE){
    if(sim==TRUE){
        RRobs <- matrix(as.vector(E1)/as.vector(init$E0),ncol=Time)
    }
    if(sim==FALSE){
        RRobs <- matrix(as.vector(vectors.sim$Y.vec)/as.vector(vectors.sim$E0),ncol=Time)
        message("Returning results for real (non-sim) data")
    }
    return(list(RRbic=matrix(lassoresult$E.qbic,ncol=Time),
                RRaic=matrix(lassoresult$E.qaic,ncol=Time),
                RRaicc=matrix(lassoresult$E.qaicc,ncol=Time),
                RRobs= RRobs))
}



