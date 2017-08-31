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
#' @param expected vector or column in a dataframe of expected counts. Format must be supplied as df$var
#' @param observed vector or column in a dataframe of observed counts. Format must be supplied as df$var
#' @param timeperiod vector or column in a dataframe of timeperiod (should be converted to a factor beforehand). Format must be supplied as df$var
#' @param covars are there additional covariates in the dataframe beyond the three required? Default is FALSE
#' @return clst object
#' 
toclust<- function(df, expected, observed, timeperiod, covars=FALSE){
    cl <- match.call()
    if(inherits(df,"data.frame") == FALSE){
        stop("Input must be a dataframe with clearly labeled covariates")
    }
    if(is.null(expected) | is.null(observed) | is.null(timeperiod)){
        stop("Must supply expected, observed, and timeperiod data for clust to run.")
    }
    if(inherits(timeperiod, "factor") == FALSE){
        timeperiod <- as.factor(timeperiod)
        warning(paste("timeperiod argument was not supplied as a factor. I've converted it to have these levels:", 
                      paste(as.character(unique(levels(timeperiod))), collapse = ","),
                      ". Please check that this is correct before proceeding."))
    }
    if(length(expected) != length(observed) | length(expected) != length(timeperiod) | length(observed)!=length(timeperiod)){
        stop("Lengths of at least one of the three required parameters (expected, observed, timeperiod) are not equal")
    }
    #print(expected)
    
    #print(unlist(strsplit(as.character(cl[[3]]),"[$]"))[3])
    requiredcolNames <- c(unlist(strsplit(as.character(cl[[3]]),"[$]"))[3],
                          unlist(strsplit(as.character(cl[[4]]),"[$]"))[3],
                          unlist(strsplit(as.character(cl[[5]]),"[$]"))[3])
    reqix <- which(names(df) %in% requiredcolNames)
    required_df <- cbind.data.frame(expected, observed,timeperiod)
    if(covars==TRUE){
        othercovariates_df <- df[,-reqix]
    }
    else{
        othercovariates_df <- NULL
    }
    res <- list(required_df = required_df,
                othercovariates_df = othercovariates_df)
    return(res)
}