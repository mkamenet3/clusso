#' Auxiliary Plotting Functions for Japanese Breast Cancer Data
#' 
#' Wrapper function easyplot
#' @param prefect prefects dataframe
#' @param polygons polygons dataframe
#' @param pdfname pdfname of what the output pdf will be called
#' @param res resultant list from clust_ function
#' @param mods string vector of which models you ran
#' @param space If space="space", then the Quasi-Poisson and Poisson spatial only models will be run; if space="spacetime" then the Quasi-Poisson and Poisson
#' spatio-temporal models will be run; if space="both" then all four models will be run
#' @param obs default is NULL. If not null, then will add "observed" instead of "oracle" label to plot for comparison map.
#' @export
#' @examples
#' pdfname <- paste0("figures/simulations/sim","_","center","_",center,"radius",radius,"_","start","_",as.numeric(paste(timeperiod, collapse = "")),"_","rr","_",gsub("[.]","",risk.ratio),".pdf")
#' easyplot(pdfname, res, mods, space="both")


easyplot <- function(prefect, polygons, pdfname, res, mods, space=c("space", "spacetime", "both"), obs=NULL){
    if(is.null(space)){ stop("You must specify `space`, `spacetime` or `both`")}
    space <- match.arg(space, several.ok = FALSE)
    pdf_qp.s <- paste0(gsub(".pdf","", pdfname),mods[1],"space" ,".pdf")
    pdf_p.s <- paste0(gsub(".pdf","", pdfname),mods[2], "space",".pdf")
    pdf_qp.st <- paste0(gsub(".pdf","", pdfname),mods[1], "ST",".pdf")
    pdf_p.st <- paste0(gsub(".pdf","", pdfname),mods[2],"ST", ".pdf")
    switch(space, 
           space = {
               plotmap_S(prefect, polygons, pdf_qp.s, res, obs, sub = res$rrcolors$rrcolors.qp.s) 
               plotmap_S(prefect, polygons, pdf_p.s, res, obs, sub = res$rrcolors$rrcolors.p.s)},
           spacetime = {
               plotmap_ST(prefect, polygons, pdf_qp.st, res, obs, sub = res$rrcolors$rrcolors.qp.st)
               plotmap_ST(prefect, polygons, pdf_p.st, res, obs, sub = res$rrcolors$rrcolors.p.st)},
           both = {
               plotmap_S(prefect, polygons, pdf_qp.s, res, obs, sub = res$rrcolors$rrcolors.qp.s)
               plotmap_S(prefect, polygons, pdf_p.s, res, obs, sub = res$rrcolors$rrcolors.p.s)
               plotmap_ST(prefect, polygons, pdf_qp.st, res, obs, sub = res$rrcolors$rrcolors.qp.st)
               plotmap_ST(prefect, polygons, pdf_p.st, res, obs, sub = res$rrcolors$rrcolors.p.st)})
}


#' #Space-time plotting
#' @param prefect prefects dataframe
#' @param polygons polygons dataframe
#' @param pdfname pdfname of what the output pdf will be called
#' @param res resultant list from clust_ function
#' @param obs if observed is to be plotted or oracle from simulation
#' @param sub optional parameter if you want to just use the function for plotting different vectors
#' 
plotmap_ST <- function(prefect, polygons, pdfname,res, obs, sub){
    if(!is.null(obs)){
        firstrow = "Obs"
    }
    else{
        firstrow="Oracle"
    }
    if(!is.null(sub)){
        rrcolors <- sub
    }
    else {
        rrcolors <- res$rrcolors
    }
    pdf(pdfname, height=11, width=10)
    #Maps of Observed Counts
    par(fig=c(0,.2,.6,1), mar=c(.5,0.5,0.5,0))
    plot(polygons,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(polygons,col=rrcolors$colors.obs[,1],border=F)
    segments(prefect$x1,prefect$y1,prefect$x2,prefect$y2)
    text(355,4120,paste0("Period 1 - ", firstrow),cex=1.00)
    
    par(fig=c(0.2,.4,.6,1), mar=c(.5,0.5,0.5,0), new=T)
    plot(polygons,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(polygons,col=rrcolors$colors.obs[,2],border=F)
    segments(prefect$x1,prefect$y1,prefect$x2,prefect$y2)
    text(355,4120,paste0("Period 2 - ", firstrow),cex=1.00)
    
    par(fig=c(0.4,.6,.6,1), mar=c(.5,0.5,0.5,0), new=T)
    plot(polygons,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(polygons,col=rrcolors$colors.obs[,3],border=F)
    segments(prefect$x1,prefect$y1,prefect$x2,prefect$y2)
    text(355,4120,paste0("Period 3 - ", firstrow),cex=1.00)
    
    par(fig=c(0.6,.8,.6,1), mar=c(.5,0.5,0.5,0), new=T)
    plot(polygons,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(polygons,col=rrcolors$colors.obs[,4],border=F)
    segments(prefect$x1,prefect$y1,prefect$x2,prefect$y2)
    text(355,4120,paste0("Period 4 - ", firstrow),cex=1.00)
    
    par(fig=c(0.8,1,.6,1), mar=c(.5,0.5,0.5,0), new=T)
    plot(polygons,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(polygons,col=rrcolors$colors.obs[,5],border=F)
    segments(prefect$x1,prefect$y1,prefect$x2,prefect$y2)
    text(355,4120,paste0("Period 5 - ", firstrow),cex=1.00)
    
    
    #Maps of AIC Path
    
    par(fig=c(0,.2,.4,.8), mar=c(.5,0.5,0.5,0), new=T)
    plot(polygons,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(polygons,col=rrcolors$color.qaic[,1],border=F)
    segments(prefect$x1,prefect$y1,prefect$x2,prefect$y2)
    text(355,4120,'Period 1 - QAIC',cex=1.00)
    
    par(fig=c(0.2,.4,.4,.8), mar=c(.5,0.5,0.5,0), new=T)   
    plot(polygons,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(polygons,col=rrcolors$color.qaic[,2],border=F)
    segments(prefect$x1,prefect$y1,prefect$x2,prefect$y2)
    text(355,4120,'Period 2 - QAIC',cex=1.00)
    
    par(fig=c(0.4,.6,.4,.8), mar=c(.5,0.5,0.5,0), new=T) 
    plot(polygons,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(polygons,col=rrcolors$color.qaic[,3],border=F)
    segments(prefect$x1,prefect$y1,prefect$x2,prefect$y2)
    text(355,4120,'Period 3 - QAIC',cex=1.00)
    
    par(fig=c(0.6,.8,.4,.8), mar=c(.5,0.5,0.5,0), new=T)
    plot(polygons,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(polygons,col=rrcolors$color.qaic[,4],border=F)
    segments(prefect$x1,prefect$y1,prefect$x2,prefect$y2)
    text(355,4120,'Period 4 - QAIC',cex=1.00)
    
    par(fig=c(0.8,1,.4,.8), mar=c(.5,0.5,0.5,0), new=T)
    plot(polygons,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(polygons,col=rrcolors$color.qaic[,5],border=F)
    segments(prefect$x1,prefect$y1,prefect$x2,prefect$y2)
    text(355,4120,'Period 5 - QAIC',cex=1.00)
    
    
    #Maps of AICc Path
    
    par(fig=c(0,.2,.2,.6), mar=c(.5,0.5,0.5,0), new=T)
    plot(polygons,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(polygons,col=rrcolors$color.qaicc[,1],border=F)
    segments(prefect$x1,prefect$y1,prefect$x2,prefect$y2)
    text(355,4120,'Period 1 - QAICc',cex=1.00)
    
    par(fig=c(0.2,.4,.2,.6), mar=c(.5,0.5,0.5,0), new=T) 
    plot(polygons,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(polygons,col=rrcolors$color.qaicc[,2],border=F)
    segments(prefect$x1,prefect$y1,prefect$x2,prefect$y2)
    text(355,4120,'Period 2 - QAICc',cex=1.00)
    
    par(fig=c(0.4,.6,.2,.6), mar=c(.5,0.5,0.5,0), new=T) 
    plot(polygons,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(polygons,col=rrcolors$color.qaicc[,3],border=F)
    segments(prefect$x1,prefect$y1,prefect$x2,prefect$y2)
    text(355,4120,'Period 3 - QAICc',cex=1.00)
    
    par(fig=c(0.6,.8,.2,.6), mar=c(.5,0.5,0.5,0), new=T)
    plot(polygons,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(polygons,col=rrcolors$color.qaicc[,4],border=F)
    segments(prefect$x1,prefect$y1,prefect$x2,prefect$y2)
    text(355,4120,'Period 4 - QAICc',cex=1.00)
    
    par(fig=c(0.8,1,.2,.6), mar=c(.5,0.5,0.5,0), new=T)
    plot(polygons,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(polygons,col=rrcolors$color.qaicc[,5],border=F)
    segments(prefect$x1,prefect$y1,prefect$x2,prefect$y2)
    text(355,4120,'Period 5 - QAICc',cex=1.00)
    
    
    #Maps of BIC Path
    
    par(fig=c(0,.2,0,.4), mar=c(.5,0.5,0.5,0), new=T)
    plot(polygons,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(polygons,col=rrcolors$color.qbic[,1],border=F)
    segments(prefect$x1,prefect$y1,prefect$x2,prefect$y2)
    text(355,4120,'Period 1 - QBIC',cex=1.00)
    
    par(fig=c(0.2,.4,0,.4), mar=c(.5,0.5,0.5,0), new=T)
    plot(polygons,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(polygons,col=rrcolors$color.qbic[,2],border=F)
    segments(prefect$x1,prefect$y1,prefect$x2,prefect$y2)
    text(355,4120,'Period 2 - QBIC',cex=1.00)
    
    par(fig=c(0.4,.6,0,.4), mar=c(.5,0.5,0.5,0), new=T)
    plot(polygons,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(polygons,col=rrcolors$color.qbic[,3],border=F)
    segments(prefect$x1,prefect$y1,prefect$x2,prefect$y2)
    text(355,4120,'Period 3 - QBIC',cex=1.00)
    
    
    par(fig=c(0.6,.8,0,.4), mar=c(.5,0.5,0.5,0), new=T)
    plot(polygons,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(polygons,col=rrcolors$color.qbic[,4],border=F)
    segments(prefect$x1,prefect$y1,prefect$x2,prefect$y2)
    text(355,4120,'Period 4 - QBIC',cex=1.00)
    
    
    par(fig=c(0.8,1,0,.4), mar=c(.5,0.5,0.5,0), new=T)
    plot(polygons,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(polygons,col=rrcolors$color.qbic[,5],border=F)
    segments(prefect$x1,prefect$y1,prefect$x2,prefect$y2)
    text(355,4120,'Period 5 - QBIC',cex=1.00)
    
    #Turn off pdf development
    dev.off()
}

#' #Space-only plotting
#' @param prefect prefects dataframe
#' @param polygons polygons dataframe
#' @param pdfname pdfname of what the output pdf will be called
#' @param res resultant list from clust_ function
#' @param obs if observed is to be plotted or oracle from simulation
#' @param sub optional parameter if you want to just use the function for plotting different vectors
#' 
plotmap_S <- function(prefect, polygons, pdfname,res, obs, sub){
    if(!is.null(obs)){
        firstrow = "Obs"
    }
    else{
        firstrow="Oracle"
    }
    if(!is.null(sub)){
        rrcolors <-  sub
        }
    else{
        rrcolors <- res$rrcolors
    }
    pdf(pdfname, height=11, width=10)
    #Maps of Observed Counts
    par(fig=c(0,.2,.6,1), mar=c(.5,0.5,0.5,0))
    plot(polygons,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(polygons,col=rrcolors$colors.obs[,1],border=F)
    segments(prefect$x1,prefect$y1,prefect$x2,prefect$y2)
    text(355,4120,paste0(firstrow),cex=1.00)
    
    #Maps of AIC Path
    
    par(fig=c(0,.2,.4,.8), mar=c(.5,0.5,0.5,0), new=T)
    plot(polygons,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(polygons,col=rrcolors$color.qaic[,1],border=F)
    segments(prefect$x1,prefect$y1,prefect$x2,prefect$y2)
    text(355,4120,'QAIC',cex=1.00)
    
    #Maps of AICc Path
    
    par(fig=c(0,.2,.2,.6), mar=c(.5,0.5,0.5,0), new=T)
    plot(polygons,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(polygons,col=rrcolors$color.qaicc[,1],border=F)
    segments(prefect$x1,prefect$y1,prefect$x2,prefect$y2)
    text(355,4120,'QAICc',cex=1.00)
    
    
    #Maps of BIC Path
    
    par(fig=c(0,.2,0,.4), mar=c(.5,0.5,0.5,0), new=T)
    plot(polygons,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(polygons,col=rrcolors$color.qbic[,1],border=F)
    segments(prefect$x1,prefect$y1,prefect$x2,prefect$y2)
    text(355,4120,'QBIC',cex=1.00)
    
    
    #Turn off pdf development
    dev.off()
}
