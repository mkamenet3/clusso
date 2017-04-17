#' Calculate overdispersion for quasi-likelihood models
#' 
#' This function calculates the overdispersion parameter for the QIC 'c' overdispersion parameter.
#' @param 
#' @return returns sparse matrix of 1's
#' @example 
#' offset_reg <- glm(Yx ~ 1 + as.factor(vectors$Period) + offset(log(Ex)),family=poisson)
#' overdisp.est <- overdisp(offset_reg)
overdisp <- function(object) {
    with(object,sum((weights * residuals^2)[weights > 0])/df.residual)
}
