#'Create clst object for creating potential spatial and spatio-temporal clusters. 
#'
#'
#'
#'@title
#'toclust
#' 
#' @description 
#' This function will create a clst object which will contain the expected, observed, and time period information necessary to run further cluster functions.
#' If dataframe is fed in, assumes panel format. Also used for converting covariates to proper object.
#' @param df dataframe of variables
#' @return clst object
#' 
toclust<- function(df){
    if(inherits(df,"data.frame") == FALSE){
        stop(message("Input must be a dataframe with clearly labeled covariates"))
    }
    if(all(unique(is.na(geoglist))!=FALSE)){
        stop("Length of expected, observed, and timeperiods vectors does not match")
    }
    clustobj <- df
    #class(clustobj) <- "clst"
    return(clustobj)
}