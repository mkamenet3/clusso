#' @title
#'toclusso
#' 
#' Creates \code{clst} object.
#' @description 
#' Creates \code{clst} object for creating potential spatial and spatio-temporal clusters. This function creates a \code{clst} 
#' object which will contain the expected, observed, and timeperiod information necessary to run \pkg{clusso}.
#' If dataframe is fed in, assumes panel format - see \code{vignette} for details. 
#' @param df Name of dataframe.
#' @param expected Name of variable that contains the expected counts.
#' @param observed Name of variable that contains the observed counts.
#' @param timeperiod Name of variable that contains the timeperiod in which counts were observed (as factor). 
#' If spatial-only analysis, create a column that has a single value (ex: "Time1") and convert this to a factor.
#' @param covars Are there additional covariates in the dataframe beyond the three required? If so, set to \code{TRUE}. Default is \code{FALSE}.
#' @param id Optional. If your dataframe contains an ID variable that should not be a covariate, set the name here.
#' @param requiredcolNames Vector of required variable names (not covariates) as extracted from \code{clusso()} arguments: \code{expected}, \code{observed}, and \code{timeperiod}.
#' @return \code{clst} object; list of lists. First element of list (called required_df) contains the expected, observed, and time period vectors in a dataframe. The second element (called othercovariates_df) is a dataframe containing covariates that should not be penalized in the model. The second element is \code{NULL} if no other covariates are adjusted for in the model.
#' @export
toclusso <- function(df, expected, observed, timeperiod, covars, id=NULL, requiredcolNames){
    if((missing(covars) | covars==FALSE)){
        covars <- FALSE
    }
    else{
        covars <- TRUE
    }

    if(is.null(expected) | is.null(observed) | is.null(timeperiod)){
        stop("Must supply expected, observed, and timeperiod data for clusso() to run.")
    }
    if(inherits(timeperiod, "factor") == FALSE){
        timeperiod <- as.factor(timeperiod)
        warning(paste("Argument timeperiod was not supplied as a factor. I've converted it to have these levels:", 
                      paste(as.character(unique(levels(timeperiod))), collapse = ","),
                      ". Please check that this is correct before proceeding."))
    }
    if(length(expected) != length(observed) | length(expected) != length(timeperiod) | length(observed)!=length(timeperiod)){
        stop("Lengths of at least one of the three required parameters (expected, observed, timeperiod) are not equal. Please check your data.")
    }
    reqix <- which(names(df) %in% requiredcolNames)
    required_df <- cbind.data.frame(expected = expected, 
                                    observed = observed,
                                    timeperiod = timeperiod)
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