#' Binomial Function
#' 
#' @title dbinom
#' @param linpred
#' @return Returns vector of Binomial log-likelihood for each proposed path of the LASSO tuning parameter values.
#' 
#' @example 
#' set.seed(1)
#' 
dbinom <- function(x, n, phat){
    #x <- observed$observedsum#number of successes
    #n <- 5#number of trials
    #phat<- exp(linpred)/(1+exp(linpred))
    loglik_i <- x*log(phat) + (n-x)*log(1-phat)
    return(sum(loglik_i))
}

invlogit <- function(linpred){
    phat<- exp(linpred)/(1+exp(linpred))
    return(phat)
}
# ccs <- c(0,0,1,0,1,0,1,1,1,1,0,1,0,1,1,0,1,0,1,1,
#          0,1,1,1,1,0,1,0,0,0,1,0,1,0,1,1,0,0,1,1)
# region <- rep(c("R1","R2", "R3","R4"), each=10)
# tp <- rep(c("1","2"), times=20)
# df <- cbind.data.frame(ccs=ccs, region=region, time=tp)
# tapply(ccs, region, table)
# 
# 
# expected <- calcexpected(df$ccs, df$time, df$region, family="binomial")
# observed <- df %>%
#     group_by(region,time) %>%
#     summarise(observedsum = sum(ccs))
# 
# mats <- matrix(rep(0,80), ncol=10)
# mats[,1] <- rep(1,8)
# mats[,2] <- c(rep(1,7),0)
# mats[,3] <- c(rep(1,6),0,0)
# mats[,4] <- c(rep(1,5),0,0,0)
# mats[,5] <- c(rep(1,4),0,0,0,0)
# mats[,6] <- c(rep(0,4),rep(1,4))
# mats[,7] <- c(rep(0,3),rep(1,5))
# mats[,8] <- c(rep(0,2),rep(1,6))
# mats[,9] <- c(rep(0,1),rep(1,7))
# mats[,10] <- c(rep(0,2), rep(1,4), rep(0,2))
# 
# 
# lasso<- glmnet::glmnet(mats, cbind(observed$observedsum,5-observed$observedsum), family = "binomial", alpha=1, 
#                        standardize = FALSE, intercept = FALSE)
# 
# coefs.lasso.all <- coef(lasso)[-1,]
# xbetaPath<- mats%*%coefs.lasso.all
# loglikes <- sapply(1:ncol(xbetaPath), function(i) dbinom(xbetaPath[,i]))
