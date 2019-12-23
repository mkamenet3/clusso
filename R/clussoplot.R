#'Plot coefficient paths from clusso
#'
#'@title
#'clussoplot
#'@description Plot the coefficient paths from the LASSO-identified potential clusters.
#'@param outclusso outclusso Object with output from \code{clusso}.
#'@param analysis A string specifying if the spatial (\code{"space"}), spatio-temporal (\code{"spacetime"}), or both spatial and spatio-temporal (\code{"both"}) analysis should be executed.  
#'@param model A string specifying which model to use, Poisson or binomial. For Poisson, specify \code{"poisson"} and both the Poisson and quasi-Poisson model results are returned. For binomial, specify \code{"binomial"} and both the binomial and quasi-binomial model results are returned.
#'@param Time Number of time periods in the analysis.
#'@param collapsetime Default is \code{FALSE}. Alternative definition for space-only model to instead collapse expected and observed counts across time by summing the counts.
#'@param cv Boolean. Takes \code{TRUE} if cross-validation was used. Default is \code{FALSE}.
#'@export
#'@import data.table
#'@return Returns a ggplot.
#'@examples
#'\donttest{
#'clussoplot(resreal, analysis="both", model="poisson",Time=5)
#'}
clussoplot <- function(outclusso, analysis=c("space","spacetime","both"), model = c("poisson", "binomial"), Time, collapsetime=FALSE,cv=FALSE){
    analysis <- match.arg(analysis, several.ok = FALSE)
    if(is.null(collapsetime)| missing(collapsetime)){
        collapsetime <- FALSE
    }
    else if(collapsetime==TRUE){
        collapsetime <- TRUE
    }
    if (length(model)>1){
        stop("You must select either `poisson` or `binomial`")
    }
    else if(model=="poisson"){
        model <- c("Poisson", "Quasi-Poisson")
    }
    else if(model=="binomial"){
        model <- c("Binomial", "Quasi-Binomial")
    }
    else {stop("Unknown model type. If you think this was by error, please submit an issue.")}
    if(cv==TRUE){
        maxdim <- dim(outclusso$lassoresult.p.st$lasso$glmnet.fit$beta)[1]
        switch(analysis,
               space = clussoplotCV(outclusso, analysistype=c("p.s", "qp.s"), model, Time, maxdim, collapsetime),
               spacetime = clussoplotCV(outclusso, analysistype = c("p.st", "qp.st"), model,Time, maxdim, collapsetime),
               both = clussoplotCV(outclusso, analysistype = c("p.s", "qp.s","p.st", "qp.st"), model,Time, maxdim, collapsetime))    
    }
    else{
        #dims
        maxdim <-dim(outclusso$lassoresult.qp.st$lasso$beta)[1]
        switch(analysis,
               space = clussoplotIC(outclusso, analysistype=c("p.s", "qp.s"), model, Time, maxdim, collapsetime),
               spacetime = clussoplotIC(outclusso, analysistype = c("p.st", "qp.st"), model, Time, maxdim, collapsetime),
               both = clussoplotIC(outclusso, analysistype = c("p.s", "qp.s","p.st", "qp.st"), model, Time, maxdim, collapsetime))    
    }
    
}

#'Helper: Plot coefficient paths (based on information criteria) from clusso
#'
#'@title
#'clussoplotIC 
#'@description Plot the coefficient paths from the LASSO-identified potential clusters (based on information criteria).
#'@param outclusso outclusso Object with output from \code{clusso}.
#'@param analysistype A string specifying if the spatial (\code{"space"}), spatio-temporal (\code{"spacetime"}), or both spatial and spatio-temporal (\code{"both"}) analysis should be executed.  
#'@param model A string specifying which model to use, Poisson or binomial. For Poisson, specify \code{"poisson"} and both the Poisson and quasi-Poisson model results are returned. For binomial, specify \code{"binomial"} and both the binomial and quasi-binomial model results are returned.
#'@param Time Number of time periods in the analysis.
#'@param maxdim Maximum number of potential clusters.
#'@param collapsetime Default is \code{FALSE}. Alternative definition for space-only model to instead collapse expected and observed counts across time by summing the counts.
#'@import data.table
#'@return Returns plots based on information criteria.
clussoplotIC <- function(outclusso, analysistype, model,Time, maxdim, collapsetime){

    for (i in 1:length(analysistype)){
        #Create labels for plots
        #labtype <- ifelse(substr(analysistype[i],1,1)=="p","Poisson", "Quasi-Poisson")
        labtype <- ifelse(substr(analysistype[i],1,1)=="p",model[1], model[2])
        dimtype <- ifelse(substr(analysistype[i],(nchar(analysistype[i])+1)-1,
                                 nchar(analysistype[i]))=="s","Space", "Space-Time")
        prefix <- paste0("outclusso$lassoresult.",analysistype[i])
        lams <- log(eval(parse(text=paste0(prefix,"$lasso$lambda"))))
        changepoints_ix <- which(diff(eval(parse(text=paste0(prefix,"$lasso$df"))))!=0)
        coefsix <- lapply(1:length(changepoints_ix), 
                          function(i) which(Matrix::t(eval(parse(text=paste0(prefix,"$lasso$beta"))))[changepoints_ix[i],]!=0))
        
        for (i in 1:length(changepoints_ix)){
            names(coefsix[[i]])<- unlist(coefsix[[i]]) 
        }
        
        
        
        coefs <- lapply(1:length(changepoints_ix),
                        function(i) Matrix::t(eval(parse(text=paste0(prefix,"$lasso$beta"))))[changepoints_ix[[i]],coefsix[[i]]])
        for (i in 1:length(changepoints_ix)){
            names(coefs[[i]])<- names(coefsix[[i]])
        }
        outclussodframe <- data.table::rbindlist(lapply(coefs, function(x) data.table::data.table(t(x))),
                                                 fill = TRUE)
        outclussodframe[is.na(outclussodframe)] <- 0
        outclussodframe$lams <- lams[changepoints_ix]
        outclussodframe$k <- eval(parse(text=paste0(prefix, "$lasso$df")))[changepoints_ix]-Time
        #convert to long and exclude unpenalized time
        s <- var <- NULL
        
        if(collapsetime==TRUE){
            outclusso_long <- tidyr::gather(outclussodframe, s, var , -c("lams", "k"), factor_key = TRUE) 
        }
        else{
            outclusso_long <- tidyr::gather(outclussodframe, s, var , -c("lams", "k"), factor_key = TRUE) %>%
                dplyr::filter(!(s %in% (maxdim-Time):maxdim))
        }
      
        #extract nclusters identified by AIC, AICc, and BIC
        numclust.qaic <- eval(parse(text=paste0(prefix,"$numclust.qaic")))
        numclust.qaicc <- eval(parse(text=paste0(prefix,"$numclust.qaicc")))
        numclust.qbic <- eval(parse(text=paste0(prefix,"$numclust.qbic")))
        kaic <- outclusso_long$lams[which(outclusso_long$k==numclust.qaic)][1]
        kaicc <- outclusso_long$lams[which(outclusso_long$k==numclust.qaicc)][1]
        kbic <- outclusso_long$lams[which(outclusso_long$k==numclust.qbic)][1]

        #PLOT!
        p <- ggplot2::ggplot(outclusso_long,ggplot2::aes(x=lams, y=var, color=s)) +
            ggplot2::geom_line(size=1.5) +
            ggplot2::theme_bw() +
            ggplot2::ylab("Coefficients") +
            ggplot2::xlab(parse(text=paste0('"(log)"', ' ~ lambda '))) +
            ggplot2::ggtitle(paste0("Solution Paths - Potential Clusters: ",labtype,", ",  dimtype)) +
            ggplot2::geom_hline(yintercept=0, lwd=1.5) +
            ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size=14),
                  legend.title = ggplot2::element_blank(),
                  legend.position="bottom",
                  legend.text=ggplot2::element_text(size=14),
                  text = ggplot2::element_text(size=14)) +
            ggplot2::scale_x_reverse() +
            ggplot2::geom_vline(xintercept=kaic, linetype=2, lwd=1.5) +
            ggplot2::geom_vline(xintercept=kaicc, linetype=2, lwd=1.5) +
            ggplot2::geom_vline(xintercept=kbic, linetype=2, lwd=1.5) +
            ggplot2::annotate("text", x=kaic, y= 0.11,hjust=-0.25,vjust=0.8,label = paste0("(Q)AIC\n  k=", numclust.qaic)) +
            ggplot2::annotate("text", x=kaicc, y= -0.03,hjust=-0.25,vjust=0.8,label = paste0("(Q)AICc\n  k=", numclust.qaicc)) +
            ggplot2::annotate("text", x=kbic, y= 0.11,hjust=-0.25,vjust=0.8,label = paste0("(Q)BIC\n  k=", numclust.qbic)) 
        print(p)
    }
}


#'Helper: Plot coefficient paths (based on cross-validation) from clusso
#'
#'@title
#'clussoplotCV
#'@description Plot the coefficient paths from the LASSO-identified potential clusters (based on cross-validation).
#'@param outclusso outclusso Object with output from \code{clusso}.
#'@param analysistype A string specifying if the spatial (\code{"space"}), spatio-temporal (\code{"spacetime"}), or both spatial and spatio-temporal (\code{"both"}) analysis should be executed.   
#'@param model A string specifying which model to use, Poisson or binomial. For Poisson, specify \code{"poisson"} and both the Poisson and quasi-Poisson model results are returned. For binomial, specify \code{"binomial"} and both the binomial and quasi-binomial model results are returned.
#'@param Time Number of time periods in the analysis.
#'@param maxdim maximum number of potential clusters.
#'@import data.table
#'@return Returns plots based on cross-validation.
clussoplotCV <- function(outclusso, analysistype,model, Time, maxdim){
    for (i in 1:length(analysistype)){
        #Create labels for plots
        #labtype <- ifelse(substr(analysistype[i],1,1)=="p","Poisson", "Quasi-Poisson")
        labtype <- ifelse(substr(analysistype[i],1,1)=="p",model[1], model[2])
        dimtype <- ifelse(substr(analysistype[i],(nchar(analysistype[i])+1)-1,
                                 nchar(analysistype[i]))=="s","Space", "Space-Time")
        prefix <- paste0("outclusso$lassoresult.",analysistype[i])
        lams <- log(eval(parse(text=paste0(prefix,"$lasso$glmnet.fit$lambda"))))
        changepoints_ix <- which(diff(eval(parse(text=paste0(prefix,"$lasso$glmnet.fit$df"))))!=0)
        coefsix <- lapply(1:length(changepoints_ix), 
                          function(i) which(Matrix::t(eval(parse(text=paste0(prefix,"$lasso$glmnet.fit$beta"))))[changepoints_ix[i],]!=0))
        
        for (i in 1:length(changepoints_ix)){
            names(coefsix[[i]])<- unlist(coefsix[[i]]) 
        }
        
        
        
        coefs <- lapply(1:length(changepoints_ix),
                        function(i) Matrix::t(eval(parse(text=paste0(prefix,"$lasso$glmnet.fit$beta"))))[changepoints_ix[[i]],coefsix[[i]]])
        for (i in 1:length(changepoints_ix)){
            names(coefs[[i]])<- names(coefsix[[i]])
        }
        outclussodframe <- data.table::rbindlist(lapply(coefs, function(x) data.table::data.table(t(x))),
                                                 fill = TRUE)
        outclussodframe[is.na(outclussodframe)] <- 0
        outclussodframe$lams <- lams[changepoints_ix]
        outclussodframe$k <- eval(parse(text=paste0(prefix, "$lasso$glmnet.fit$df")))[changepoints_ix]-Time
        #convert to long and exclude unpenalized time
        s <- var <- NULL
        outclusso_long <- tidyr::gather(outclussodframe, s, var, -c("lams", "k"), factor_key = TRUE) %>%
            dplyr::filter(!(s %in% (maxdim-Time):maxdim))
        numclust.cv <- eval(parse(text=paste0(prefix,"$numclust.cv")))
        kcv <- outclusso_long$lams[which(outclusso_long$k==numclust.cv)][1]
        if(is.na(kcv)){
            message(paste0("In analysis ",labtype,", ",
                           dimtype,
                           " , number of detected clusters exceeds maxclust. Please increase maxclust and re-run clusso." ))
        }
        
        #PLOT!
        p <- ggplot2::ggplot(outclusso_long,ggplot2::aes(x=lams, y=var, color=s)) +
            ggplot2::geom_line(size=1.5) +
            ggplot2::theme_bw() +
            ggplot2::ylab("Coefficients") +
            ggplot2::xlab(parse(text=paste0('"(log)"', ' ~ lambda '))) +
            ggplot2::ggtitle(paste0("Solution Paths - Potential Clusters: ",labtype,", ",  dimtype)) +
            ggplot2::geom_hline(yintercept=0, lwd=1.5) +
            ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size=14),
                  legend.title = ggplot2::element_blank(),
                  legend.position="bottom",
                  legend.text=ggplot2::element_text(size=14),
                  text = ggplot2::element_text(size=14)) +
            ggplot2::scale_x_reverse() +
            ggplot2::geom_vline(xintercept=kcv, linetype=2, lwd=1.5) +
            #ggplot2::geom_vline(xintercept=kaicc, linetype=2, lwd=1.5) +
            #ggplot2::geom_vline(xintercept=kbic, linetype=2, lwd=1.5) +
            ggplot2::annotate("text", x=kcv, y= 0.11,hjust=-0.25,vjust=0.8,label = paste0("CV\n  k=", numclust.cv)) 
            #ggplot2::annotate("text", x=kaicc, y= -0.03,hjust=-0.25,vjust=0.8,label = paste0("(Q)AICc\n  k=", numclust.qaicc)) +
            #ggplot2::annotate("text", x=kbic, y= 0.11,hjust=-0.25,vjust=0.8,label = paste0("(Q)BIC\n  k=", numclust.qbic)) 
        print(p)
    }
}
