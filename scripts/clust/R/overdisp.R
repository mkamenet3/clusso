#' Calculate overdispersion for quasi-likelihood models
#' 
#' This function calculates the overdispersion parameter for the QIC 'c' overdispersion parameter.
#' @param offset_reg object of class glm
#' @param sim for a simulation run this takes on TRUE, else it should be FALSE. Default is set to TRUE.
#' @param floor Default is set to TRUE so model will not estimate any under-dispersion where phi < 1. The floor is set for phi to be 1. In a case where this 
#' @return returns phi overdispersion parameter
#' @examples
#' offset <- glm(Yx ~ 1 + as.factor(vectors$Period) + offset(log(Ex)), family=poisson)
#' overdisp.est <- overdisp(offset_reg)
overdisp <- function(offset_reg, sim = TRUE, floor = TRUE) {
    if(sim==TRUE){
        stopifnot(inherits(offset_reg[[1]], c("glm", "lm")))
        phi <- max(unlist(lapply(1:nsim, function(i) deviance(offset_reg[[i]])/df.residual(offset_reg[[i]]))))
    }
    else{
        stopifnot(inherits(offset_reg, c("glm", "lm")))
        phi <- max(unlist(deviance(offset_reg)/df.residual(offset_reg)))
    }
    if(floor == TRUE & phi < 1){
        message(paste("Underdispersion detected (", phi,"). Setting phi to 1"))
        phi <- 1
    }
    if(floor == FALSE & phi < 1){
        message(paste("Underdispersion detect (", phi,").\n 'Floor' argument was set to FALSE to underdispersed model will be run"))
    }
    return(phi)
}
