#' Auxiliary Plotting Functions for Japanese Breast Cancer Data
#' 
#' Wrapper function easyplot
#' @param prefect prefects dataframe
#' @param polygons polygons dataframe
#' @param pdfname pdfname of what the output pdf will be called
#' @param rescols result colors to plot
#' @param mods string vector of which models you ran
#' @param space If space="space", then the Quasi-Poisson and Poisson spatial only models will be run; if space="spacetime" then the Quasi-Poisson and Poisson
#' spatio-temporal models will be run; if space="both" then all four models will be run
#' @param probmap TRUE or FALSE if probabilities should be mapped or relative risks
#' @param cv default is NULL. If not null, then results from cross-validation will be plotted (for real data example)
#' @param obs default is NULL. If not null, then will add "observed" instead of "oracle" label to plot for comparison map.
#' @param rr TRUE or FALSE. TRUE indicates risk ratios should be plotted, FALSE indicates probabilities should be plotted.
#' @export
#' @examples 
#' \donttest{
#' mods <- c("QuasiPoisson", "Poisson")
#' pdfname1 <- paste0("rrmaps.pdf")
#' pdfname2 <- paste0("probabilitymaps.pdf")
#' easyplot(prefect = japan.prefect2 , polygons = japan.poly2 ,pdfname1,
#'  res$rrcolors, mods, space="both", probmap=FALSE, obs = NULL, rr=TRUE)
#' easyplot(prefect = japan.prefect2 , polygons = japan.poly2 ,pdfname2, 
#'  res$probcolors, mods, space="both", probmap=TRUE, obs = NULL,rr=FALSE)
#' }

easyplot <- function(prefect, polygons, pdfname, rescols, mods, space=c("space", "spacetime", "both", "BIC"), probmap, cv=NULL,obs=NULL,rr){
    if(is.null(space)){ stop("You must specify `space`, `spacetime` or `both`")}
    space <- match.arg(space, several.ok = FALSE)
    pdf_qp.s <- paste0(gsub(".pdf","", pdfname),mods[1],"space" ,".pdf")
    pdf_p.s <- paste0(gsub(".pdf","", pdfname),mods[2], "space",".pdf")
    pdf_qp.st <- paste0(gsub(".pdf","", pdfname),mods[1], "ST",".pdf")
    pdf_p.st <- paste0(gsub(".pdf","", pdfname),mods[2],"ST", ".pdf")
    
    if(probmap==FALSE & is.null(cv)){
        message("RR Maps")
        switch(space, 
               space = {
                   plotmap_S(prefect, polygons, pdf_qp.s, rescols$rrcolors.qp.s, obs)
                   plotmap_S(prefect, polygons, pdf_p.s, rescols$rrcolors.p.s, obs)},
               spacetime = {
                   plotmap_ST(prefect, polygons, pdf_qp.st, rescols$rrcolors.qp.st, obs)
                   plotmap_ST(prefect, polygons, pdf_p.st, rescols$rrcolors.p.st, obs)},
               both = {
                   plotmap_ST(prefect, polygons, pdf_qp.s, rescols$rrcolors.qp.s, obs, rr)
                   plotmap_ST(prefect, polygons, pdf_p.s, rescols$rrcolors.p.s, obs, rr)
                   plotmap_ST(prefect, polygons, pdf_qp.st, rescols$rrcolors.qp.st, obs, rr)
                   plotmap_ST(prefect, polygons, pdf_p.st, rescols$rrcolors.p.st, obs, rr)},
               BIC = {
                   plotmap_bic(prefect, polygons, pdfname, rescols, rr)
               })
    }
    if(probmap==TRUE & is.null(cv)){
        message("Probability Maps")
        switch(space, 
               space = {
                   plotmap_S(prefect, polygons, pdf_qp.s, rescols$probcolors.qp.s, obs) 
                   plotmap_S(prefect, polygons, pdf_p.s, rescols$probcolors.p.s, obs)},
               spacetime = {
                   plotmap_ST(prefect, polygons, pdf_qp.st, rescols$probcolors.qp.st, obs)
                   plotmap_ST(prefect, polygons, pdf_p.st, rescols$probcolors.p.s, obs)},
               both = {
                   plotmap_ST(prefect, polygons, pdf_qp.s, rescols$probcolors.qp.s, obs, rr)
                   plotmap_ST(prefect, polygons, pdf_p.s, rescols$probcolors.p.s, obs, rr)
                   plotmap_ST(prefect, polygons, pdf_qp.st, rescols$probcolors.qp.st, obs, rr)
                   plotmap_ST(prefect, polygons, pdf_p.st, rescols$probcolors.qp.st, obs, rr)})
    }
    if(!is.null(cv)){
        message("CV Maps")
        switch(space, 
               # TODO
               # space = {
               #    # plotmap_S(prefect, polygons, pdf_qp.s, res, obs, sub = res$probcolors$probcolors.qp.s) 
               #     plotmap_S(prefect, polygons, pdf_p.s, res, obs, sub = res$probcolors$probcolors.p.s)},
               # spacetime = {
               #     plotmap_ST(prefect, polygons, pdf_qp.st, res, obs, sub = res$probcolors$probcolors.qp.st)
               #     plotmap_ST(prefect, polygons, pdf_p.st, res, obs, sub = res$probcolors$probcolors.p.s)},
               both = {
                   plotmap_ST_cv(prefect, polygons, pdf_qp.s,  rescols$rrcolors.qp.s, obs)
                   plotmap_ST_cv(prefect, polygons, pdf_p.s, rescols$rrcolors.p.s, obs)
                   plotmap_ST_cv(prefect, polygons, pdf_qp.st, rescols$rrcolors.qp.st, obs)
                   plotmap_ST_cv(prefect, polygons, pdf_p.st, rescols$rrcolors.p.st, obs)})
    }
    
}


#' Space-time plotting
#' @param prefect prefects dataframe
#' @param polygons polygons dataframe
#' @param pdfname pdfname of what the output pdf will be called
#' @param res resultant list from clust_ function
#' @param obs if observed is to be plotted or oracle from simulation
#' 
plotmap_ST_cv <- function(prefect, polygons, pdfname,res, obs){
    if(!is.null(obs)){
        firstrow = "Observed"
    }
    else{
        firstrow="Oracle"
    }
    pdf(pdfname, height=11, width=10)
    #Maps of Observed Counts
    par(fig=c(0,.2,.6,1), mar=c(.5,0.5,0.5,0))
    plot(polygons,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(polygons,col=res$colors.obs[,1],border=F)
    segments(prefect$x1,prefect$y1,prefect$x2,prefect$y2)
    text(355,4120,paste0("Period 1 - ", firstrow),cex=1.00)
    
    par(fig=c(0.2,.4,.6,1), mar=c(.5,0.5,0.5,0), new=T)
    plot(polygons,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(polygons,col=res$colors.obs[,2],border=F)
    segments(prefect$x1,prefect$y1,prefect$x2,prefect$y2)
    text(355,4120,paste0("Period 2 - ", firstrow),cex=1.00)
    
    par(fig=c(0.4,.6,.6,1), mar=c(.5,0.5,0.5,0), new=T)
    plot(polygons,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(polygons,col=res$colors.obs[,3],border=F)
    segments(prefect$x1,prefect$y1,prefect$x2,prefect$y2)
    text(355,4120,paste0("Period 3 - ", firstrow),cex=1.00)
    
    par(fig=c(0.6,.8,.6,1), mar=c(.5,0.5,0.5,0), new=T)
    plot(polygons,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(polygons,col=res$colors.obs[,4],border=F)
    segments(prefect$x1,prefect$y1,prefect$x2,prefect$y2)
    text(355,4120,paste0("Period 4 - ", firstrow),cex=1.00)
    
    par(fig=c(0.8,1,.6,1), mar=c(.5,0.5,0.5,0), new=T)
    plot(polygons,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(polygons,col=res$colors.obs[,5],border=F)
    segments(prefect$x1,prefect$y1,prefect$x2,prefect$y2)
    text(355,4120,paste0("Period 5 - ", firstrow),cex=1.00)
    
    
    #Maps of CV Path
    
    par(fig=c(0,.2,.4,.8), mar=c(.5,0.5,0.5,0), new=T)
    plot(polygons,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(polygons,col=res$color.cv[,1],border=F)
    segments(prefect$x1,prefect$y1,prefect$x2,prefect$y2)
    text(355,4120,'Period 1 - CV',cex=1.00)
    
    par(fig=c(0.2,.4,.4,.8), mar=c(.5,0.5,0.5,0), new=T)   
    plot(polygons,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(polygons,col=res$color.cv[,2],border=F)
    segments(prefect$x1,prefect$y1,prefect$x2,prefect$y2)
    text(355,4120,'Period 2 - CV',cex=1.00)
    
    par(fig=c(0.4,.6,.4,.8), mar=c(.5,0.5,0.5,0), new=T) 
    plot(polygons,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(polygons,col=res$color.cv[,3],border=F)
    segments(prefect$x1,prefect$y1,prefect$x2,prefect$y2)
    text(355,4120,'Period 3 - CV',cex=1.00)
    
    par(fig=c(0.6,.8,.4,.8), mar=c(.5,0.5,0.5,0), new=T)
    plot(polygons,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(polygons,col=res$color.cv[,4],border=F)
    segments(prefect$x1,prefect$y1,prefect$x2,prefect$y2)
    text(355,4120,'Period 4 - CV',cex=1.00)
    
    par(fig=c(0.8,1,.4,.8), mar=c(.5,0.5,0.5,0), new=T)
    plot(polygons,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(polygons,col=res$color.cv[,5],border=F)
    segments(prefect$x1,prefect$y1,prefect$x2,prefect$y2)
    text(355,4120,'Period 5 - CV',cex=1.00)
    
    #legend
    par(fig=c(.35,.75,0,.1), new=T)
    plot(1, xlim=c(0.6,1.5), ylim=c(0.2,1), axes=F, type='n',  xlab="", ylab="")
    rect(seq(.6,1.4,length=50)[-50],.5,seq(.65,1.4,length=50)[-1],.62,col=redblue(0:50/50),border=F)
    text(seq(.6,1.4,length=5),rep(.45,5),seq(0,2,length.out=5),srt=330,adj=0)
    
    
    #Turn off pdf development
    dev.off()
}

#' Space-only plotting
#' @param prefect prefects dataframe
#' @param polygons polygons dataframe
#' @param pdfname pdfname of what the output pdf will be called
#' @param res resultant list from clust_ function
#' @param obs if observed is to be plotted or oracle from simulation
plotmap_S_cv <- function(prefect, polygons, pdfname,res, obs){
    if(!is.null(obs)){
        firstrow = "Observed"
    }
    else{
        firstrow="Oracle"
    }
    pdf(pdfname, height=11, width=10)
    #Maps of Observed Counts
    par(fig=c(0,.2,.6,1), mar=c(.5,0.5,0.5,0))
    plot(polygons,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(polygons,col=res$colors.obs[,1],border=F)
    segments(prefect$x1,prefect$y1,prefect$x2,prefect$y2)
    text(355,4120,paste0(firstrow),cex=1.00)
    
    #Maps of AIC Path
    
    par(fig=c(0,.2,.4,.8), mar=c(.5,0.5,0.5,0), new=T)
    plot(polygons,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(polygons,col=res$color.cv[,1],border=F)
    segments(prefect$x1,prefect$y1,prefect$x2,prefect$y2)
    text(355,4120,'CV',cex=1.00)
    
    #legend
    par(fig=c(.35,.75,0,.1), new=T)
    plot(1, xlim=c(0.6,1.5), ylim=c(0.2,1), axes=F, type='n',  xlab="", ylab="")
    rect(seq(.6,1.4,length=50)[-50],.5,seq(.65,1.4,length=50)[-1],.62,col=redblue(0:50/50),border=F)
    text(seq(.6,1.4,length=5),rep(.45,5),seq(0,2,length.out=5),srt=330,adj=0)
   
    #Turn off pdf development
    dev.off()
}


#' Space-time plotting
#' @param prefect prefects dataframe
#' @param polygons polygons dataframe
#' @param pdfname pdfname of what the output pdf will be called
#' @param res resultant list from clust_ function
#' @param obs if observed is to be plotted or oracle from simulation
#' @param rr if FALSE, will print probability map legend, if TRUE will print legend for risk ratios (redblue scheme)
plotmap_ST <- function(prefect, polygons, pdfname,res, obs,rr){
    if(!is.null(obs)){
        firstrow = "Observed"
    }
    else{
        firstrow="Oracle"
    }
    pdf(pdfname, height=11, width=10)
    #Maps of Observed Counts
    par(fig=c(0,.2,.6,1), mar=c(.5,0.5,0.5,0))
    plot(polygons,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(polygons,col=res$colors.obs[,1],border=F)
    segments(prefect$x1,prefect$y1,prefect$x2,prefect$y2)
    text(355,4120,paste0("Period 1 - ", firstrow),cex=1.00)
    
    par(fig=c(0.2,.4,.6,1), mar=c(.5,0.5,0.5,0), new=T)
    plot(polygons,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(polygons,col=res$colors.obs[,2],border=F)
    segments(prefect$x1,prefect$y1,prefect$x2,prefect$y2)
    text(355,4120,paste0("Period 2 - ", firstrow),cex=1.00)
    
    par(fig=c(0.4,.6,.6,1), mar=c(.5,0.5,0.5,0), new=T)
    plot(polygons,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(polygons,col=res$colors.obs[,3],border=F)
    segments(prefect$x1,prefect$y1,prefect$x2,prefect$y2)
    text(355,4120,paste0("Period 3 - ", firstrow),cex=1.00)
    
    par(fig=c(0.6,.8,.6,1), mar=c(.5,0.5,0.5,0), new=T)
    plot(polygons,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(polygons,col=res$colors.obs[,4],border=F)
    segments(prefect$x1,prefect$y1,prefect$x2,prefect$y2)
    text(355,4120,paste0("Period 4 - ", firstrow),cex=1.00)
    
    par(fig=c(0.8,1,.6,1), mar=c(.5,0.5,0.5,0), new=T)
    plot(polygons,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(polygons,col=res$colors.obs[,5],border=F)
    segments(prefect$x1,prefect$y1,prefect$x2,prefect$y2)
    text(355,4120,paste0("Period 5 - ", firstrow),cex=1.00)
    
    
    #Maps of AIC Path
    
    par(fig=c(0,.2,.4,.8), mar=c(.5,0.5,0.5,0), new=T)
    plot(polygons,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(polygons,col=res$color.qaic[,1],border=F)
    segments(prefect$x1,prefect$y1,prefect$x2,prefect$y2)
    text(355,4120,'Period 1 - QAIC',cex=1.00)
    
    par(fig=c(0.2,.4,.4,.8), mar=c(.5,0.5,0.5,0), new=T)   
    plot(polygons,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(polygons,col=res$color.qaic[,2],border=F)
    segments(prefect$x1,prefect$y1,prefect$x2,prefect$y2)
    text(355,4120,'Period 2 - QAIC',cex=1.00)
    
    par(fig=c(0.4,.6,.4,.8), mar=c(.5,0.5,0.5,0), new=T) 
    plot(polygons,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(polygons,col=res$color.qaic[,3],border=F)
    segments(prefect$x1,prefect$y1,prefect$x2,prefect$y2)
    text(355,4120,'Period 3 - QAIC',cex=1.00)
    
    par(fig=c(0.6,.8,.4,.8), mar=c(.5,0.5,0.5,0), new=T)
    plot(polygons,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(polygons,col=res$color.qaic[,4],border=F)
    segments(prefect$x1,prefect$y1,prefect$x2,prefect$y2)
    text(355,4120,'Period 4 - QAIC',cex=1.00)
    
    par(fig=c(0.8,1,.4,.8), mar=c(.5,0.5,0.5,0), new=T)
    plot(polygons,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(polygons,col=res$color.qaic[,5],border=F)
    segments(prefect$x1,prefect$y1,prefect$x2,prefect$y2)
    text(355,4120,'Period 5 - QAIC',cex=1.00)
    
    
    #Maps of AICc Path
    
    par(fig=c(0,.2,.2,.6), mar=c(.5,0.5,0.5,0), new=T)
    plot(polygons,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(polygons,col=res$color.qaicc[,1],border=F)
    segments(prefect$x1,prefect$y1,prefect$x2,prefect$y2)
    text(355,4120,'Period 1 - QAICc',cex=1.00)
    
    par(fig=c(0.2,.4,.2,.6), mar=c(.5,0.5,0.5,0), new=T) 
    plot(polygons,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(polygons,col=res$color.qaicc[,2],border=F)
    segments(prefect$x1,prefect$y1,prefect$x2,prefect$y2)
    text(355,4120,'Period 2 - QAICc',cex=1.00)
    
    par(fig=c(0.4,.6,.2,.6), mar=c(.5,0.5,0.5,0), new=T) 
    plot(polygons,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(polygons,col=res$color.qaicc[,3],border=F)
    segments(prefect$x1,prefect$y1,prefect$x2,prefect$y2)
    text(355,4120,'Period 3 - QAICc',cex=1.00)
    
    par(fig=c(0.6,.8,.2,.6), mar=c(.5,0.5,0.5,0), new=T)
    plot(polygons,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(polygons,col=res$color.qaicc[,4],border=F)
    segments(prefect$x1,prefect$y1,prefect$x2,prefect$y2)
    text(355,4120,'Period 4 - QAICc',cex=1.00)
    
    par(fig=c(0.8,1,.2,.6), mar=c(.5,0.5,0.5,0), new=T)
    plot(polygons,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(polygons,col=res$color.qaicc[,5],border=F)
    segments(prefect$x1,prefect$y1,prefect$x2,prefect$y2)
    text(355,4120,'Period 5 - QAICc',cex=1.00)
    
    
    #Maps of BIC Path
    
    par(fig=c(0,.2,0,.4), mar=c(.5,0.5,0.5,0), new=T)
    plot(polygons,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(polygons,col=res$color.qbic[,1],border=F)
    segments(prefect$x1,prefect$y1,prefect$x2,prefect$y2)
    text(355,4120,'Period 1 - QBIC',cex=1.00)
    
    par(fig=c(0.2,.4,0,.4), mar=c(.5,0.5,0.5,0), new=T)
    plot(polygons,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(polygons,col=res$color.qbic[,2],border=F)
    segments(prefect$x1,prefect$y1,prefect$x2,prefect$y2)
    text(355,4120,'Period 2 - QBIC',cex=1.00)
    
    par(fig=c(0.4,.6,0,.4), mar=c(.5,0.5,0.5,0), new=T)
    plot(polygons,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(polygons,col=res$color.qbic[,3],border=F)
    segments(prefect$x1,prefect$y1,prefect$x2,prefect$y2)
    text(355,4120,'Period 3 - QBIC',cex=1.00)
    
    
    par(fig=c(0.6,.8,0,.4), mar=c(.5,0.5,0.5,0), new=T)
    plot(polygons,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(polygons,col=res$color.qbic[,4],border=F)
    segments(prefect$x1,prefect$y1,prefect$x2,prefect$y2)
    text(355,4120,'Period 4 - QBIC',cex=1.00)
    
    
    par(fig=c(0.8,1,0,.4), mar=c(.5,0.5,0.5,0), new=T)
    plot(polygons,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(polygons,col=res$color.qbic[,5],border=F)
    segments(prefect$x1,prefect$y1,prefect$x2,prefect$y2)
    text(355,4120,'Period 5 - QBIC',cex=1.00)
    
    
    #legend
    if(rr==TRUE) {
        par(fig=c(.35,.75,0,.1), new=T)
        plot(1, xlim=c(0.6,1.5), ylim=c(0.2,1), axes=F, type='n',  xlab="", ylab="")
        rect(seq(.6,1.4,length=50)[-50],.5,seq(.65,1.4,length=50)[-1],.62,col=redblue(0:50/50),border=F)
        #text(seq(.6,1.4,length=5),rep(.45,5),seq(0,2,length.out=5),srt=330,adj=0)
        text(seq(0.6, 1.4, length = 5), rep(0.45, 5), seq(0.5, 1.5, length.out = 5), srt = 330, adj = 0)
    }
    else{
        par(fig=c(.35,.75,0,.1), new=T)
        plot(1, xlim=c(0.6,1.5), ylim=c(0.2,1), axes=F, type='n',  xlab="", ylab="")
        rect(seq(.6,1.4,length=50)[-50],.5,seq(.65,1.4,length=50)[-1],.62,col=greys(0:50/50),border=F)
        #text(seq(.6,1.4,length=6),rep(.45,6),c('0.0','0.2','0.4','0.6','0.8','1.0'),srt=330,adj=0)
        text(seq(0.6, 1.4, length = 5), rep(0.45, 5), seq(0.5, 1.5, length.out = 5), srt = 330, adj = 0)
    }
    
    #Turn off pdf development
    dev.off()
}

#' Space-only plotting
#' @param prefect prefects dataframe
#' @param polygons polygons dataframe
#' @param pdfname pdfname of what the output pdf will be called
#' @param res resultant list from clust_ function
#' @param obs if observed is to be plotted or oracle from simulation
#' @param rr is TRUE, result in legend for redblue (risk ratios). Default is grey scale legend
plotmap_S <- function(prefect, polygons, pdfname,res, obs, rr){
    if(!is.null(obs)){
        firstrow = "Observed"
    }
    else{
        firstrow="Oracle"
    }
    pdf(pdfname, height=11, width=10)
    #Maps of Observed Counts
    par(fig=c(0,.2,.6,1), mar=c(.5,0.5,0.5,0))
    plot(polygons,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(polygons,col=res$colors.obs[,1],border=F)
    segments(prefect$x1,prefect$y1,prefect$x2,prefect$y2)
    text(355,4120,paste0(firstrow),cex=1.00)
    
    #Maps of AIC Path
    
    par(fig=c(0,.2,.4,.8), mar=c(.5,0.5,0.5,0), new=T)
    plot(polygons,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(polygons,col=res$color.qaic[,1],border=F)
    segments(prefect$x1,prefect$y1,prefect$x2,prefect$y2)
    text(355,4120,'QAIC',cex=1.00)
    
    #Maps of AICc Path
    
    par(fig=c(0,.2,.2,.6), mar=c(.5,0.5,0.5,0), new=T)
    plot(polygons,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(polygons,col=res$color.qaicc[,1],border=F)
    segments(prefect$x1,prefect$y1,prefect$x2,prefect$y2)
    text(355,4120,'QAICc',cex=1.00)
    
    
    #Maps of BIC Path
    
    par(fig=c(0,.2,0,.4), mar=c(.5,0.5,0.5,0), new=T)
    plot(polygons,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(polygons,col=res$color.qbic[,1],border=F)
    segments(prefect$x1,prefect$y1,prefect$x2,prefect$y2)
    text(355,4120,'QBIC',cex=1.00)
    
    #legend
    if(rr==TRUE) {
        par(fig=c(.35,.75,0,.1), new=T)
        plot(1, xlim=c(0.6,1.5), ylim=c(0.2,1), axes=F, type='n',  xlab="", ylab="")
        rect(seq(.6,1.4,length=50)[-50],.5,seq(.65,1.4,length=50)[-1],.62,col=redblue(0:50/50),border=F)
        #text(seq(.6,1.4,length=5),rep(.45,5),seq(0,2,length.out=5),srt=330,adj=0)
        text(seq(0.6, 1.4, length = 5), rep(0.45, 5), seq(0.5, 1.5, length.out = 5), srt = 330, adj = 0)
    }
    else{
        par(fig=c(.35,.75,0,.1), new=T)
        plot(1, xlim=c(0.6,1.5), ylim=c(0.2,1), axes=F, type='n',  xlab="", ylab="")
        rect(seq(.6,1.4,length=50)[-50],.5,seq(.65,1.4,length=50)[-1],.62,col=greys(0:50/50),border=F)
        #text(seq(.6,1.4,length=6),rep(.45,6),c('0.0','0.2','0.4','0.6','0.8','1.0'),srt=330,adj=0)
        text(seq(0.6, 1.4, length = 5), rep(0.45, 5), seq(0.5, 1.5, length.out = 5), srt = 330, adj = 0)
    }
    
    #Turn off pdf development
    dev.off()
}


######################################################################
#' Space-time plotting
#' @param prefect prefects dataframe
#' @param polygons polygons dataframe
#' @param pdfname pdfname of what the output pdf will be called
#' @param res resultant list from clust_ function
#' @param rr if FALSE, will print probability map legend, if TRUE will print legend for risk ratios (redblue scheme)
plotmap_bic <- function(prefect, polygons, pdfname,res,rr){
    #Observed plot
    firstrow = "Observed"
    #Observed
    pdf(paste0("observed",pdfname), height=11, width=10)
    #Maps of Observed Counts
    par(fig=c(0,.2,0.668,1), mar=c(.5,0.5,0.5,0))
    plot(polygons,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(polygons,col=res$rrcolors$rrcolors.qp.st$colors.obs[,1],border=F)
    segments(prefect$x1,prefect$y1,prefect$x2,prefect$y2)
    text(355,4120,paste0("Period 1 - ", firstrow),cex=0.75)
    
    par(fig=c(0.2,.4,0.668,1), mar=c(.5,0.5,0.5,0), new=T)
    plot(polygons,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(polygons,col=res$rrcolors$rrcolors.qp.st$colors.obs[,2],border=F)
    segments(prefect$x1,prefect$y1,prefect$x2,prefect$y2)
    text(355,4120,paste0("Period 2 - ", firstrow),cex=0.75)
    
    par(fig=c(0.4,.6,0.668,1), mar=c(.5,0.5,0.5,0), new=T)
    plot(polygons,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(polygons,col=res$rrcolors$rrcolors.qp.st$colors.obs[,3],border=F)
    segments(prefect$x1,prefect$y1,prefect$x2,prefect$y2)
    text(355,4120,paste0("Period 3 - ", firstrow),cex=0.75)
    
    par(fig=c(0.6,.8,0.668,1), mar=c(.5,0.5,0.5,0), new=T)
    plot(polygons,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(polygons,col=res$rrcolors$rrcolors.qp.st$colors.obs[,4],border=F)
    segments(prefect$x1,prefect$y1,prefect$x2,prefect$y2)
    text(355,4120,paste0("Period 4 - ", firstrow),cex=0.75)
    
    par(fig=c(0.8,1,0.668,1), mar=c(.5,0.5,0.5,0), new=T)
    plot(polygons,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(polygons,col=res$rrcolors$rrcolors.qp.st$colors.obs[,5],border=F)
    segments(prefect$x1,prefect$y1,prefect$x2,prefect$y2)
    text(355,4120,paste0("Period 5 - ", firstrow),cex=0.75)
    #observed legend
    par(fig=c(0.4,0.6,0,0.1), new=T)
    plot(1, xlim=c(0.6,1.5), ylim=c(0.1,1), axes=F, type='n',  xlab="", ylab="")
    rect(seq(.6,1.4,length=50)[-50],.5,seq(.65,1.4,length=50)[-1],.62,col=redblue(0:50/50),border=F)
    text(seq(.6,1.4,length=5),rep(.45,5),seq(0,2,length.out=5),srt=330,adj=0)
    dev.off()
    
    
    #BIC AIC SELECTIONS MAPPINGS
    pdf(pdfname, height=11, width=10)
    #Map QP-ST
    print(pdfname)
    par(fig=c(0,.2,0.501,0.835), mar=c(.5,0.5,0.5,0), new=T)
    plot(polygons,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(polygons,col=res$rrcolors$rrcolors.qp.st$color.qbic[,1],border=F)
    segments(prefect$x1,prefect$y1,prefect$x2,prefect$y2)
    text(355,4120,'Period 1: BIC, quasi-Poisson',cex=0.75)
    
    par(fig=c(0.2,.4,0.501,0.835), mar=c(.5,0.5,0.5,0), new=T)   
    plot(polygons,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(polygons,col=res$rrcolors$rrcolors.qp.st$color.qbic[,2],border=F)
    segments(prefect$x1,prefect$y1,prefect$x2,prefect$y2)
    text(355,4120,'Period 2: BIC, quasi-Poisson',cex=0.75)
    
    par(fig=c(0.4,.6,0.501,0.835), mar=c(.5,0.5,0.5,0), new=T) 
    plot(polygons,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(polygons,col=res$rrcolors$rrcolors.qp.st$color.qbic[,3],border=F)
    segments(prefect$x1,prefect$y1,prefect$x2,prefect$y2)
    text(355,4120,'Period 3: BIC, quasi-Poisson',cex=0.75)
    
    par(fig=c(0.6,.8,0.501,0.835), mar=c(.5,0.5,0.5,0), new=T)
    plot(polygons,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(polygons,col=res$rrcolors$rrcolors.qp.st$color.qbic[,4],border=F)
    segments(prefect$x1,prefect$y1,prefect$x2,prefect$y2)
    text(355,4120,'Period 4: BIC, quasi-Poisson',cex=0.75)
    
    par(fig=c(0.8,1,0.501,0.835), mar=c(.5,0.5,0.5,0), new=T)
    plot(polygons,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(polygons,col=res$rrcolors$rrcolors.qp.st$color.qbic[,5],border=F)
    segments(prefect$x1,prefect$y1,prefect$x2,prefect$y2)
    text(355,4120,'Period 5: BIC, quasi-Poisson',cex=0.75)
    
    
    #Maps of P,ST
    
    par(fig=c(0,.2,0.334,0.668), mar=c(.5,0.5,0.5,0), new=T)
    plot(polygons,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(polygons,col=res$rrcolors$rrcolors.p.st$color.qaic[,1],border=F)
    segments(prefect$x1,prefect$y1,prefect$x2,prefect$y2)
    text(355,4120,'Period 1: AIC, quasi-Poisson',cex=0.75)
    
    par(fig=c(0.2,.4,0.334,0.668), mar=c(.5,0.5,0.5,0), new=T) 
    plot(polygons,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(polygons,col=res$rrcolors$rrcolors.p.st$color.qaic[,2],border=F)
    segments(prefect$x1,prefect$y1,prefect$x2,prefect$y2)
    text(355,4120,'Period 2: AIC, quasi-Poisson',cex=0.75)
    
    par(fig=c(0.4,.6,0.334,0.668), mar=c(.5,0.5,0.5,0), new=T) 
    plot(polygons,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(polygons,col=res$rrcolors$rrcolors.p.st$color.qbic[,3],border=F)
    segments(prefect$x1,prefect$y1,prefect$x2,prefect$y2)
    text(355,4120,'Period 3: AIC, quasi-Poisson',cex=0.75)
    
    par(fig=c(0.6,.8,0.334,0.668), mar=c(.5,0.5,0.5,0), new=T)
    plot(polygons,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(polygons,col=res$rrcolors$rrcolors.p.st$color.qaic[,4],border=F)
    segments(prefect$x1,prefect$y1,prefect$x2,prefect$y2)
    text(355,4120,'Period 4: AIC, quasi-Poisson',cex=0.75)
    
    par(fig=c(0.8,1,0.334,0.668), mar=c(.5,0.5,0.5,0), new=T)
    plot(polygons,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(polygons,col=res$rrcolors$rrcolors.p.st$color.qaic[,5],border=F)
    segments(prefect$x1,prefect$y1,prefect$x2,prefect$y2)
    text(355,4120,'Period 5: AIC, quasi-Poisson',cex=0.75)
    
    
    #Maps of QP,S
    
    par(fig=c(0,.2,0.167,0.501), mar=c(.5,0.5,0.5,0), new=T)
    plot(polygons,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(polygons,col=res$rrcolors$rrcolors.p.st$color.qbic[,1],border=F)
    segments(prefect$x1,prefect$y1,prefect$x2,prefect$y2)
    text(355,4120,'Period 1: BIC, Poisson',cex=0.75)
    
    par(fig=c(0.2,.4,0.167,0.501), mar=c(.5,0.5,0.5,0), new=T)
    plot(polygons,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(polygons,col=res$rrcolors$rrcolors.p.st$color.qbic[,2],border=F)
    segments(prefect$x1,prefect$y1,prefect$x2,prefect$y2)
    text(355,4120,'Period 2: BIC, Poisson',cex=0.75)
    
    par(fig=c(0.4,.6,0.167,0.501), mar=c(.5,0.5,0.5,0), new=T)
    plot(polygons,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(polygons,col=res$rrcolors$rrcolors.p.st$color.qbic[,3],border=F)
    segments(prefect$x1,prefect$y1,prefect$x2,prefect$y2)
    text(355,4120,'Period 3: BIC, Poisson',cex=0.75)
    
    
    par(fig=c(0.6,.8,0.167,0.501), mar=c(.5,0.5,0.5,0), new=T)
    plot(polygons,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(polygons,col=res$rrcolors$rrcolors.p.st$color.qbic[,4],border=F)
    segments(prefect$x1,prefect$y1,prefect$x2,prefect$y2)
    text(355,4120,'Period 4: BIC, Poisson',cex=0.75)
    
    
    par(fig=c(0.8,1,0.167,0.501), mar=c(.5,0.5,0.5,0), new=T)
    plot(polygons,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(polygons,col=res$rrcolors$rrcolors.p.st$color.qbic[,5],border=F)
    segments(prefect$x1,prefect$y1,prefect$x2,prefect$y2)
    text(355,4120,'Period 5: BIC, Poisson',cex=0.75)
    
    
    #Maps of P,S
    
    par(fig=c(0,.2,0,0.334), mar=c(.5,0.5,0.5,0), new=T)
    plot(polygons,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(polygons,col=res$rrcolors$rrcolors.p.st$color.qaic[,1],border=F)
    segments(prefect$x1,prefect$y1,prefect$x2,prefect$y2)
    text(355,4120,'Period 1: AIC, Poisson',cex=0.75)
    
    par(fig=c(0.2,.4,0,0.334), mar=c(.5,0.5,0.5,0), new=T)
    plot(polygons,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(polygons,col=res$rrcolors$rrcolors.p.st$color.qaic[,2],border=F)
    segments(prefect$x1,prefect$y1,prefect$x2,prefect$y2)
    text(355,4120,'Period 2: AIC, Poisson',cex=0.75)
    
    par(fig=c(0.4,.6,0,0.334), mar=c(.5,0.5,0.5,0), new=T)
    plot(polygons,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(polygons,col=res$rrcolors$rrcolors.p.st$color.qaic[,3],border=F)
    segments(prefect$x1,prefect$y1,prefect$x2,prefect$y2)
    text(355,4120,'Period 3: AIC, Poisson',cex=0.75)
    
    
    par(fig=c(0.6,.8,0,0.334), mar=c(.5,0.5,0.5,0), new=T)
    plot(polygons,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(polygons,col=res$rrcolors$rrcolors.p.st$color.qaic[,4],border=F)
    segments(prefect$x1,prefect$y1,prefect$x2,prefect$y2)
    text(355,4120,'Period 4: AIC, Poisson',cex=0.75)
    
    
    par(fig=c(0.8,1,0,0.334), mar=c(.5,0.5,0.5,0), new=T)
    plot(polygons,type='n',asp=1,axes=F,xlab='',ylab='')
    polygon(polygons,col=res$rrcolors$rrcolors.p.st$color.qaic[,5],border=F)
    segments(prefect$x1,prefect$y1,prefect$x2,prefect$y2)
    text(355,4120,'Period 5: AIC, Poisson',cex=0.75)
    #legend
    par(fig=c(.35,.75,0,.1), new=T)
    plot(1, xlim=c(0.6,1.5), ylim=c(0.2,1), axes=F, type='n',  xlab="", ylab="")
    rect(seq(.6,1.4,length=50)[-50],.5,seq(.65,1.4,length=50)[-1],.62,col=redblue(0:50/50),border=F)
    #text(seq(.6,1.4,length=5),rep(.45,5),seq(0,2,length.out=5),srt=330,adj=0)
    text(seq(0.6, 1.4, length = 5), rep(0.45, 5), seq(0.5, 1.5, length.out = 5), srt = 330, adj = 0)
    
    #Turn off pdf development
    dev.off()
}
