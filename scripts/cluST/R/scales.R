#' Re-scaling data to get SMR (standardized mortality rate)
#'
#'This function standardizes the fitted values by the sum of the observed/the sum of the expected within each time period. Inputs should all be in vector form.
#'The function handles standardizing within time period as long.
#'@param YSIM list of simulated observed
#'@param init list of initial settings in spacetimeLasso 
#'@param Time number of time periods
#'@export
#'
#'@examples
#'
#'#Simulation
#'#set up
#'Time <- 2; theta = 1000; nsim = 2
#'period <- rep(seq(1,2),5)
#'expected <- rnegbin(n = 10,mu = 15,theta = 1000)
#'observed <- rnegbin(E, theta=1000)
#'#set init
#'init <- setVectors(period, expected, observed, Time, byrow=TRUE)
#'ysim <- lapply(1:nsim, function(i) rnegbin(expected, theta = theta))
#'scale.sim(ysim, init, nsim, Time))
#'
scale <- function(init,Time,scaler=NULL,...){
    # if(!is.null(scaler)){
    #     print(scaler)
    #     std <- sapply(1:Time, function(i) 
    #         (((matrix(init$E0, ncol=Time)[,i])/scaler)*(sum(matrix(init$Y.vec, ncol=Time)[,i])/scaler))/(sum(matrix(init$E0, ncol=Time)[,i])/scaler))
    #     E0 <- as.vector(std)
    # }
    # else{
    std <- sapply(1:Time, function(i) 
        ((((matrix(init$E0, ncol=Time)[,i]))*(sum(matrix(init$Y.vec, ncol=Time)[,i])))/(sum(matrix(init$E0, ncol=Time)[,i]))))
    E0 <- as.vector(std)
    # }
    return(E0)
}

#scaling for simulation
scale.sim <- function(YSIM, init, nsim,Time,...){
    std <- lapply(1:nsim, function(i) sapply(1:Time, function(j) 
        (matrix(init$E0,ncol=Time)[,j])*(sum(matrix(YSIM[[i]],ncol=Time)[,j])/sum(matrix(init$E0,ncol=Time)[,j]))))
    E0 <- lapply(1:nsim, function(i) as.vector(std[[i]])) 
    return(E0)
}

