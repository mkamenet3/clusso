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
#' setVectors(period, expected, observed, Time, byrow=TRUE)
#' 
setVectors <- function(period, expect, observed,Time, byrow=TRUE,...) {
    if (byrow==TRUE){
        E0=as.vector(matrix(expect, byrow=T, ncol=Time))
        Y.vec <- as.vector(matrix(observed,byrow=T, ncol=Time))
        Year <- as.vector(matrix(period, byrow=T, ncol=Time)) 
    }
    else {
        E0=as.vector(matrix(expect, ncol=Time))
        Y.vec <- as.vector(matrix(observed, ncol=Time))
        Year <- as.vector(matrix(period, ncol=Time))
    }
    return(list(
        E0 = E0,
        Y.vec = Y.vec,
        Year = Year))
}
