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
#' @param id Optional. If your dataframe contains an ID variable that should not be a covariate, set the name here.
#' @return \code{clst} object
#'@examples
#' @export
#'\donttest{
#'data(japanbreastcancer)
#'clst <- toclusso(japanbreastcancer, expected = expdeath, observed=death,timeperiod = period, covars = FALSE)  
#'}

toclusso <- function(df, expected, observed, timeperiod, covars=FALSE, id=NULL){
    if((missing(covars) | covars==FALSE)){
        covars <- FALSE
    }
    else{
        covars <- TRUE
    }
    expect <- eval(substitute(expected),df)
    observe <- eval(substitute(observed),df)
    period <- eval(substitute(timeperiod),df)
    
    if(inherits(df,"data.frame") == FALSE){
        stop("Input must be a dataframe with clearly labeled covariates")
    }
    if(is.null(expect) | is.null(observe) | is.null(period)){
        stop("Must supply expected, observed, and timeperiod data for clusso() to run.")
    }
    if(inherits(period, "factor") == FALSE){
        period <- as.factor(period)
        warning(paste("timeperiod argument was not supplied as a factor. I've converted it to have these levels:", 
                      paste(as.character(unique(levels(period))), collapse = ","),
                      ". Please check that this is correct before proceeding."))
    }
    if(length(expect) != length(observe) | length(expect) != length(period) | length(observe)!=length(period)){
        stop("Lengths of at least one of the three required parameters (expected, observed, timeperiod) are not equal. Please check your data.")
    }
    
    # if(!is.null(id)){
    #     ids <- substitute(id)
    #     print('a')
    # }
    # else{
    #     ids <- NULL
    #     print('b')
    # }
    # print("id:", id)
    requiredcolNames <- c(substitute(expected),
      substitute(observed),
      substitute(timeperiod),
      substitute(id))

    reqix <- which(names(df) %in% requiredcolNames)
    #print(requiredcolNames)
    #print(names(df))
    #print(reqix)
    required_df <- cbind.data.frame(expect, observe,period)
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