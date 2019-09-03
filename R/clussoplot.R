#'Plot coefficient paths from clusso
#'
#'@title
#'clussoplot
#'@description Plot the coefficient paths from the LASSO-identified potential clusters.
#'@param outclusso outclusso Object with output from \code{clusso}.
#'@param analysis A string specifying if the spatial (\code{"space"}), spatio-temporal (\code{"spacetime"}), or both spatial and spatio-temporal (\code{"both"}) analysis  should be executed. 
#'@param Time Number of time periods in the analysis
#'@param cv Boolean. Takes TRUE if cross-validation was used. Default is FALSE.
#'@export
#'@return Returns a plot
#'@examples
#'\donttest{
#'clussoplot(resreal, analysis="both", Time=5)
#'}
clussoplot <- function(outclusso, analysis=c("space","spacetime","both"), Time, cv=FALSE){
    analysis <- match.arg(analysis, several.ok = FALSE)
    if(cv==TRUE){
        maxdim <- dim(rescv$lassoresult.p.st$lasso$glmnet.fit$beta)[1]
        switch(analysis,
               space = clussoplotCV(outclusso, analysistype=c("p.s", "qp.s"), Time, maxdim),
               spacetime = clussoplotCV(outclusso, analysistype = c("p.st", "qp.st"), Time, maxdim),
               both = clussoplotCV(outclusso, analysistype = c("p.s", "qp.s","p.st", "qp.st"), Time, maxdim))    
    }
    else{
        #dims
        maxdim <-dim(outclusso$lassoresult.qp.st$lasso$beta)[1]
        switch(analysis,
               space = clussoplotIC(outclusso, analysistype=c("p.s", "qp.s"), Time, maxdim),
               spacetime = clussoplotIC(outclusso, analysistype = c("p.st", "qp.st"), Time, maxdim),
               both = clussoplotIC(outclusso, analysistype = c("p.s", "qp.s","p.st", "qp.st"), Time, maxdim))    
    }
    
}

#'Helper: Plot coefficient paths (based on information criteria) from clusso
#'
#'@title
#'clussoplotIC 
#'@description Plot the coefficient paths from the LASSO-identified potential clusters (based on information criteria).
#'@param outclusso outclusso Object with output from \code{clusso}.
#'@param analysis A string specifying if the spatial (\code{"space"}), spatio-temporal (\code{"spacetime"}), or both spatial and spatio-temporal (\code{"both"}) analysis should be executed. 
#'@param Time Number of time periods in the analysis
#'@param maxdim maximum number of potential clusters
#'@return Returns plots
clussoplotIC <- function(outclusso, analysistype, Time, maxdim){
    for (i in 1:length(analysistype)){
        #Create labels for plots
        labtype <- ifelse(substr(analysistype[i],1,1)=="p","Poisson", "Quasi-Poisson")
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
        outclusso_long <- tidyr::gather(outclussodframe, s, var, -c("lams", "k"), factor_key = TRUE) %>%
            dplyr::filter(!(s %in% (maxdim-Time):maxdim))
        #extract nclusters identified by AIC, AICc, and BIC
        numclust.qaic <- eval(parse(text=paste0(prefix,"$numclust.qaic")))
        numclust.qaicc <- eval(parse(text=paste0(prefix,"$numclust.qaicc")))
        numclust.qbic <- eval(parse(text=paste0(prefix,"$numclust.qbic")))
        kaic <- outclusso_long$lams[which(outclusso_long$k==numclust.qaic)][1]
        kaicc <- outclusso_long$lams[which(outclusso_long$k==numclust.qaicc)][1]
        kbic <- outclusso_long$lams[which(outclusso_long$k==numclust.qbic)][1]

        #PLOT!
        p <- ggplot2::ggplot(outclusso_long,aes(x=lams, y=var, color=s)) +
            geom_line(size=1.5) +
            theme_bw() +
            ylab("Coefficients") +
            xlab(parse(text=paste0('"(log)"', ' ~ lambda '))) +
            ggtitle(paste0("Solution Paths - Potential Clusters: ",labtype,", ",  dimtype)) +
            geom_hline(yintercept=0, lwd=1.5) +
            theme(plot.title = element_text(hjust = 0.5, size=14),
                  legend.title = element_blank(),
                  legend.position="bottom",
                  legend.text=element_text(size=14),
                  text = element_text(size=14)) +
            scale_x_reverse() +
            geom_vline(xintercept=kaic, linetype=2, lwd=1.5) +
            geom_vline(xintercept=kaicc, linetype=2, lwd=1.5) +
            geom_vline(xintercept=kbic, linetype=2, lwd=1.5) +
            annotate("text", x=kaic, y= 0.11,hjust=-0.25,vjust=0.8,label = paste0("(Q)AIC\n  k=", numclust.qaic)) +
            annotate("text", x=kaicc, y= -0.03,hjust=-0.25,vjust=0.8,label = paste0("(Q)AICc\n  k=", numclust.qaicc)) +
            annotate("text", x=kbic, y= 0.11,hjust=-0.25,vjust=0.8,label = paste0("(Q)BIC\n  k=", numclust.qbic)) 
        print(p)
    }
}


#'Helper: Plot coefficient paths (based on cross-validation) from clusso
#'
#'@title
#'clussoplotCV
#'@description Plot the coefficient paths from the LASSO-identified potential clusters (based on cross-validation).
#'@param outclusso outclusso Object with output from \code{clusso}.
#'@param analysis A string specifying if the spatial (\code{"space"}), spatio-temporal (\code{"spacetime"}), or both spatial and spatio-temporal (\code{"both"}) analysis should be executed. 
#'@param Time Number of time periods in the analysis
#'@param maxdim maximum number of potential clusters
#'@return Returns plots
clussoplotCV <- function(outclusso, analysistype, Time, maxdim){
    for (i in 1:length(analysistype)){
        #Create labels for plots
        labtype <- ifelse(substr(analysistype[i],1,1)=="p","Poisson", "Quasi-Poisson")
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
        p <- ggplot2::ggplot(outclusso_long,aes(x=lams, y=var, color=s)) +
            geom_line(size=1.5) +
            theme_bw() +
            ylab("Coefficients") +
            xlab(parse(text=paste0('"(log)"', ' ~ lambda '))) +
            ggtitle(paste0("Solution Paths - Potential Clusters: ",labtype,", ",  dimtype)) +
            geom_hline(yintercept=0, lwd=1.5) +
            theme(plot.title = element_text(hjust = 0.5, size=14),
                  legend.title = element_blank(),
                  legend.position="bottom",
                  legend.text=element_text(size=14),
                  text = element_text(size=14)) +
            scale_x_reverse() +
            geom_vline(xintercept=kcv, linetype=2, lwd=1.5) +
            #geom_vline(xintercept=kaicc, linetype=2, lwd=1.5) +
            #geom_vline(xintercept=kbic, linetype=2, lwd=1.5) +
            annotate("text", x=kcv, y= 0.11,hjust=-0.25,vjust=0.8,label = paste0("CV\n  k=", numclust.cv)) 
            #annotate("text", x=kaicc, y= -0.03,hjust=-0.25,vjust=0.8,label = paste0("(Q)AICc\n  k=", numclust.qaicc)) +
            #annotate("text", x=kbic, y= 0.11,hjust=-0.25,vjust=0.8,label = paste0("(Q)BIC\n  k=", numclust.qbic)) 
        print(p)
    }
}
