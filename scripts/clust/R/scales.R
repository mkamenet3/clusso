#' Re-scaling data to get SMR (standardized mortality rate)
#'
#'This function standardizes the fitted values by the sum of the observed/the sum of the expected within each time period. Inputs should all be in vector form.
#'The function handles standardizing within time period as long.
#'@param init list of initial settings in spacetimeLasso 
#'@param Time number of time periods
#'@export
#'
#'@examples
#'
#'#set up
#'set.seed(2)
#'Time <- 2
#'theta = 1000
#'nsim = 2
#'period <- rep(seq(1,2),5)
#'expected <- rnegbin(n = 10,mu = 15,theta = 1000)
#'observed <- rnegbin(expected, theta=1000)
#'#set init
#'init <- setVectors(period, expected, observed, Time, byrow=TRUE)
#'scale(init, Time)
scale <- function(init,Time){
    if(isTRUE(length(init$E0)%%Time ==0) == FALSE) {stop("Length of Expected vector (E0) not divisible by number of time periods")}
    if(isTRUE(length(init$E0)!=length(init$Y.vec))){stop("Lengths of Expected (E0) and Observed (Y.vec) not equivalent")}
    std <- sapply(1:Time, function(i) 
        ((((matrix(init$E0, ncol=Time)[,i]))*(sum(matrix(init$Y.vec, ncol=Time)[,i])))/(sum(matrix(init$E0, ncol=Time)[,i]))))
    E0 <- as.vector(std)
    return(E0)
}

#'This function standardizes the fitted values by the sum of the *simulated* observed/the sum of the expected within each time period. Inputs should all be in vector form.
#'The function handles standardizing within time period as long. This function should be used with simulation study.
#'@param YSIM list of simulated observed
#'@param init list of initial settings in spacetimeLasso 
#'@param nsim number of simulations
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
#'scale_sim(ysim, init, nsim, Time))

scale_sim <- function(YSIM, init, nsim,Time){
    if(isTRUE(length(init$E0)!= mean(unlist(lapply(YSIM, function(x) length(x)))))) {
        stop("Length of Expected vector (E0) not divisible by number of observed YSIM[[i]] in at least one simulation")
        }
    std <- lapply(1:nsim, function(i) sapply(1:Time, function(j) 
        (matrix(init$E0,ncol=Time)[,j])*(sum(matrix(YSIM[[i]],ncol=Time)[,j])/sum(matrix(init$E0,ncol=Time)[,j]))))
    E0 <- lapply(1:nsim, function(i) as.vector(std[[i]])) 
    return(E0)
}

