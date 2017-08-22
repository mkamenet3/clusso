#' Creates a List Arranged by Time Period with Expected and Observed Counts and Time Period
#' 
#' @param period vector of periods or years in dataset. Should be imported as a factor.
#' @param expect vector of expected counts. Expected counts must match up with the year and observed vectors.
#' @param observed vector of observed counts. Observed counts must match up with the year and expected vectors.
#' @param Time Number of time periods or years in your dataset. Must be declared as numeric.
#' @param byrow default is set to TRUE. Data from the dataset should be imported by row. This is most often the case
#' when you have a dataframe ordered by an identifier and then the period/time frame within that id listed chronologically (in panel format by identifier).
#' If you are simulating data and have each observed/expected vector separate and create the period vector with repetitions of each time
#' period by group, this should be set to false.
#' @return This function returns a list of expected and observed counts along with the period. 
#' @export
#' @examples
#' period1 <- c(rep("1",5),rep("2",5))
#' period2 <- rep(seq(1,2),5)
#' expected <- rnegbin(n = 10,mu = 15,theta = 1000)
#' observed <- rnegbin(expected, theta=1000)
#' Time = 2
#' setVectors(period1, expected, observed, Time, byrow=TRUE))
#' setVectors(period2, expected, observed, Time, byrow=FALSE)



setVectors <- function(period, expect, observed,Time, byrow=TRUE) {
    if (byrow==TRUE){
        if(period[1] == period[2]) warning("Please check the format of the data, you may want byrow=FALSE. It appears that the time periods appear sequentially")
        E0=as.vector(matrix(expect, byrow=T, ncol=Time))
        Y.vec <- as.vector(matrix(observed,byrow=T, ncol=Time))
        Year <- as.vector(matrix(period, byrow=T, ncol=Time)) 
    }
    else {
        if(period[1] != period[2]) warning("Please check the format of the data, you may want byrow=TRUE. It appears that the time periods do not appear sequentially")
        E0=as.vector(matrix(expect, ncol=Time))
        Y.vec <- as.vector(matrix(observed, ncol=Time))
        Year <- as.vector(matrix(period, ncol=Time))
    }
    return(list(
        E0 = E0,
        Y.vec = Y.vec,
        Year = Year))
}
