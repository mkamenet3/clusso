
#'@title
#'setVectors
#'@description
#'  Creates a List Arranged by Time Period with Expected and Observed Counts and Time Period
#' 
#' @param period Vector of timeperiods in the data set. If this is variable is not a factor in the dataframe, then it will be automatically converted to one by \code{clusso()} with a warning message. Periods must match up with observed and expected vectors.
#' @param expect Vector of expected counts. Expected counts must match up with period and observed vectors.
#' @param observed Vector of observed counts. Observed counts must match up with period and expected vectors.
#' @param covars Dataframe of covariates. \code{NULL} if no covariates supplied. 
#' @param Time Number of timeperiods in the dataset. This is calculated based on the number of unique timeperiods (factor levels) supplied to \code{clusso}.
#' @param byrow Set to \code{TRUE}. Data from the dataset should be imported by row, meaning that data must be sorted by geographic identifier, with repeated measurements for each geographic unit repeated. 
#' @return This function returns a list of expected and observed counts along with the period. 
#' @export
#' @examples
#' period1 <- rep(c("1","2"),times=5)
#' expected <- MASS::rnegbin(n = 10,mu = 15,theta = 1000)
#' observed <- MASS::rnegbin(expected, theta=1000)
#' Time = 2
#' covars <- NULL
#' setVectors(period1, expected, observed, covars, Time)

setVectors <- function(period, expect, observed, covars, Time, byrow=TRUE) {
    if (byrow==TRUE){
        if(period[1] == period[2]) warning("Please check the format of the data, you may must sort data by geographic identifier and then time period. It appears that the data is sorted by time period only.")
        E0=as.vector(matrix(expect, byrow=TRUE, ncol=Time))
        Y.vec <- as.vector(matrix(observed,byrow=TRUE, ncol=Time))
        Year <- as.vector(matrix(period, byrow=TRUE, ncol=Time)) 
        if(!is.null(covars)){
            covars_df <- sapply(covars, function(x) matrix(x, byrow=TRUE, ncol=Time))    
        }
        else{
            covars_df <- NULL
        }
    }
    else {
        E0=as.vector(matrix(expect, ncol=Time))
        Y.vec <- as.vector(matrix(observed, ncol=Time))
        Year <- as.vector(matrix(period, ncol=Time))
        if(!is.null(covars)){
            covars_df <- sapply(covars, function(x) matrix(x, ncol=Time))
        }
        else{
            covars_df <-NULL
        }
    }
    return(list(
        E0 = E0,
        Y.vec = Y.vec,
        Year = Year,
        covars = covars_df))
}
