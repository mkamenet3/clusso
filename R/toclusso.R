#' @title
#'toclusso
#' 
#' Creates \code{clst} object.
#' @description 
#' Creates \code{clst} object for creating potential spatial and spatio-temporal clusters. This function creates a \code{clst} 
#' object which will contain the expected, observed, and timeperiod information necessary to run \pkg{clusso}.
#' If dataframe is fed in, assumes panel format - see \code{vignette} for details. 
#' @param df name of dataframe.
#' @param expected Name of variable that contains the expected counts.
#' @param observed Name of variable that contains the observed counts.
#' @param timeperiod Name of variable that contains the timeperiod in which counts were observed (as factor). 
#' If spatial-only analysis, create a column that has a single value (ex: "Time1") and convert this to a factor.
#' @param covars are there additional covariates in the dataframe beyond the three required? If so, set to \code{TRUE}. Default is \code{FALSE}.
#' @return \code{clst} object
#'@examples
#' @export
#'\donttest{
#'data(japanbreastcancer)
#'clst <- toclusso(japanbreastcancer, expected = expdeath, observed=death,timeperiod = period, covars = FALSE)  
#'}

toclusso <- function(df, expected, observed, timeperiod, covars=FALSE){
    if((missing(covars) | covars==FALSE)){
        covars <- FALSE
    }
    else{
        covars <- TRUE
    }
    #print(as.character(expected))
    #a <- as.character(expected)
    #print(a)
    #print(str(expected))
    # cl <- match.call()
    # print(cl)
    # print(str(cl))
    expected <- eval(substitute(expected),df)
    observed <- eval(substitute(observed),df)
    timeperiod <- eval(substitute(timeperiod),df)
    
    
    #print(expected)
    if(inherits(df,"data.frame") == FALSE){
        stop("Input must be a dataframe with clearly labeled covariates")
    }
    if(is.null(expected) | is.null(observed) | is.null(timeperiod)){
        stop("Must supply expected, observed, and timeperiod data for clusso() to run.")
    }
    if(inherits(timeperiod, "factor") == FALSE){
        timeperiod <- as.factor(timeperiod)
        warning(paste("timeperiod argument was not supplied as a factor. I've converted it to have these levels:", 
                      paste(as.character(unique(levels(timeperiod))), collapse = ","),
                      ". Please check that this is correct before proceeding."))
    }
    if(length(expected) != length(observed) | length(expected) != length(timeperiod) | length(observed)!=length(timeperiod)){
        stop("Lengths of at least one of the three required parameters (expected, observed, timeperiod) are not equal. Please check your data.")
    }
    
    requiredcolNames <- names(df)[1:3]
    
   # requiredcolNames <- c("expected", "observed", "timeperiod")
        # requiredcolNames <- c(unlist(strsplit(as.character(cl[[3]]),"[$]"))[3],
    #                       unlist(strsplit(as.character(cl[[4]]),"[$]"))[3],
    #                       unlist(strsplit(as.character(cl[[5]]),"[$]"))[3])
    #message(requiredcolNames)
    
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
    class(res) <- "clst"
    return(res)
}