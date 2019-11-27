#' Calculate overdispersion for quasi-likelihood models
#' 
#' This function calculates the overdispersion parameter for the QIC 'c' overdispersion parameter.
#' @param offset_reg Object of class glm.
#' @param sim For simulation run this takes on \code{TRUE}, else it should be \code{FALSE}. Default is set to \code{TRUE}.
#' @param overdispfloor  Default is \code{TRUE}. When \code{TRUE}, it limits \eqn{\phi} (overdispersion parameter) to be greater or equal to 1. If \code{FALSE}, it will allow for under-dispersion in the model.
#' @return returns \eqn{\phi} overdispersion parameter.

overdisp <- function(offset_reg, sim = TRUE, overdispfloor = TRUE) {
    if(sim==TRUE){
        stopifnot(inherits(offset_reg[[1]], c("glm", "lm")))
        phi <- max(unlist(lapply(1:length(offset_reg), function(i) deviance(offset_reg[[i]])/df.residual(offset_reg[[i]]))))
    }
    else{
        stopifnot(inherits(offset_reg, c("glm", "lm")))
        phi <- max(unlist(deviance(offset_reg)/df.residual(offset_reg)))
    }
    if(overdispfloor == TRUE & phi < 1){
        message(paste("Underdispersion detected (", phi,"). Setting phi to 1"))
        phi <- 1
    }
    if(overdispfloor == FALSE & phi < 1){
        message(paste("Underdispersion detect (", phi,").\n 'Floor' argument was set to FALSE to underdispersed model will be run"))
    }
    return(phi)
}
