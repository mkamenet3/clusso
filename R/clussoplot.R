#'Plot coefficient paths from clusso
#'
#'@title
#'clussoplot
#'@description Plot the coefficient paths from the LASSO-identified potential clusters.
#'@param
#'
#'
#'
#'@export
#'@return Returns a plot
#'@examples





# 
# library(data.table)
# library(Matrix)
# library(dplyr)
# library(tidyr)
# library(ggplot2)

clussoplot <- function(outclusso, analysis=c("space","spacetime","both"), Time){
    #dims
    maxdim <-dim(outclusso$lassoresult.qp.st$lasso$beta)[1]
    analysis <- match.arg(analysis, several.ok = FALSE)
    switch(analysis,
           space = clussoplotMaster(outclusso, analysistype=c("p.s", "qp.s"), Time, maxdim),
           spacetime = clussoplotMaster(),
           both = clussoplotMaster())
}

clussoplotMaster <- function(outclusso, analysistype, Time, maxdim){
    for (i in 1:length(analysistype)){
        print(analysistype[i])
        labtype <- ifelse(substr(analysistype[i],1,1)=="p","Poisson", "Quasi-Poisson")
        print(labtype)
        prefix <- paste0("outclusso$lassoresult.",analysistype[i])
        lams <- log(eval(parse(text=paste0(prefix,"$lasso$lambda"))))
        
        
        #lams <- log(outclusso$lassoresult.qp.st$lasso$lambda)
        changepoints_ix <- which(diff(eval(parse(text=paste0(prefix,"$lasso$df"))))!=0)
        coefsix <- lapply(1:length(changepoints_ix), 
                          function(i) which( t(eval(parse(text=paste0(prefix,"$lasso$beta"))))[changepoints_ix[i],]!=0))
        
        for (i in 1:length(changepoints_ix)){
            names(coefsix[[i]])<- unlist(coefsix[[i]]) 
        }
        
        
        
        coefs <- lapply(1:length(changepoints_ix),
                        function(i) t(eval(parse(text=paste0(prefix,"$lasso$beta"))))[changepoints_ix[[i]],coefsix[[i]]])
        for (i in 1:length(changepoints_ix)){
            names(coefs[[i]])<- names(coefsix[[i]])
        }
        outclussodframe <- data.table::rbindlist(lapply(coefs, function(x) data.table(t(x))),
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
        #Create labels for plots
        #if(substr(analysistype[i],1,1)=="p"){labtype="Poisson"}
        # labtype <- ifelse(substr(eval(parse(text=analysistype[i])),1,1)=="p","Poisson", "Quasi-Poisson")
        # print(labtype)
        #if(substr(analysistype[i],1,1)=="q"){labtype="Quasi-Poisson"}
        
        #PLOT!
        p <- ggplot2::ggplot(outclusso_long,aes(x=lams, y=var, color=s)) +
            geom_line(size=1.5) +
            theme_bw() +
            ylab("Coefficients") +
            xlab(parse(text=paste0('"(log)"', ' ~ lambda '))) +
            ggtitle(paste0("Solution Paths - Potential Clusters: ",labtype)) +
            geom_hline(yintercept=0, lwd=1.5) +
            theme(plot.title = element_text(hjust = 0.5, size=18),
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
