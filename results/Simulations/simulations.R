########################################################################################################
########################################################################################################
########################################################################################################
#Space and Space-Time Cluster Detection Using the Lasso 
##Running simulations and testing functions
##This should later be integrated into vignettes/examples

#Maria Kamenetsky
#1-15-2017
########################################################################################################
########################################################################################################
########################################################################################################

##############################################
##############################################
#Set-Up
##############################################
##############################################

####################################################
#Source Files, Scripts, and Packages
####################################################
#import libraries
library("Rcpp")
library("MASS")
library("RColorBrewer")
library("geosphere")
library(Matrix)
library(glmnet)
library(maps)
library(truncnorm)

#Source .cpp files
sourceCpp("scripts/cluST/src/maxcol.cpp")
sourceCpp("scripts/cluST/src/st_matCpp.cpp")
sourceCpp("scripts/cluST/src/prod_yx.cpp")


#temporarily source my cluST.R file
source("scripts/cluST//R//cluST.R")

########################################################################################################

####################################################
#Import Data
####################################################
set.seed(1152017)

dframe1 <- read.csv("data/JBC/jap.breast.F.9.10.11.csv")
dframe2 <- read.csv("data//JBC//utmJapan.csv")
dframe3 <- aggregate(dframe1, by=list(as.factor(rep(1:(nrow(dframe1)/4),each=4))), FUN="sum")
dframe=data.frame(id=as.factor(dframe3$id/4),period=as.factor(dframe3$year),death=dframe3$death,expdeath=dframe3$expdeath)
levels(dframe$period) <- c("1","2","3","4","5")

dframe.poly2 <- read.csv("data/JBC/japan_poly2.csv")
japan.poly2 <- dframe.poly2[,2:3]
dframe.prefect2 <- read.csv("data/JBC/japan_prefect2.csv")
japan.prefect2 <- dframe.prefect2[,2:5]


########################################################################################################
########################################################################################################
########################################################################################################

####################################################
#rmax 30, time = 3,4,5, rr.ratio = 1.5, center = 200, r = 18km
####################################################
#Set Some Initial Conditions
x=dframe2$utmx/1000
y=dframe2$utmy/1000
rMax=30 
Time=5  
nsim=100
center=200
radius=18
period_start=3
risk.ratio=1.5




res <- clust.sim(x,y,rMax,dframe$period, dframe$expdeath, dframe$death, Time, spacetime=TRUE, 
                  nsim,center, radius, risk.ratio, period_start, colors=TRUE, utm=TRUE, byrow=TRUE)

#make maps
pdfname <- paste0("figures/simulations/japan","_","rMax","_",rMax,"nsim","_",nsim, "center","_",center,"radius","_",radius, "start","_",period_start, "rr","_",risk.ratio,".pdf")
pdf(pdfname, height=11, width=10)

#Maps of Observed Counts
par(fig=c(0,.2,.6,1), mar=c(.5,0.5,0.5,0))
plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
polygon(japan.poly2,col=res$rrcolors$colors.obs[,1],border=F)
segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
text(355,4120,'Period 1 - Obs',cex=1.00)

par(fig=c(0.2,.4,.6,1), mar=c(.5,0.5,0.5,0), new=T)
plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
polygon(japan.poly2,col=res$rrcolors$colors.obs[,2],border=F)
segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
text(355,4120,'Period 2 - Obs',cex=1.00)

par(fig=c(0.4,.6,.6,1), mar=c(.5,0.5,0.5,0), new=T)
plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
polygon(japan.poly2,col=res$rrcolors$colors.obs[,3],border=F)
segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
text(355,4120,'Period 3 - Obs',cex=1.00)

par(fig=c(0.6,.8,.6,1), mar=c(.5,0.5,0.5,0), new=T)
plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
polygon(japan.poly2,col=res$rrcolors$colors.obs[,4],border=F)
segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
text(355,4120,'Period 4 - Obs',cex=1.00)

par(fig=c(0.8,1,.6,1), mar=c(.5,0.5,0.5,0), new=T)
plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
polygon(japan.poly2,col=res$rrcolors$colors.obs[,5],border=F)
segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
text(355,4120,'Period 5 - Obs',cex=1.00)


#Maps of AIC Path

par(fig=c(0,.2,.4,.8), mar=c(.5,0.5,0.5,0), new=T)
plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
polygon(japan.poly2,col=res$rrcolors$color.qaic[,1],border=F)
segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
text(355,4120,'Period 1 - QAIC',cex=1.00)

par(fig=c(0.2,.4,.4,.8), mar=c(.5,0.5,0.5,0), new=T)   
plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
polygon(japan.poly2,col=res$rrcolors$color.qaic[,2],border=F)
segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
text(355,4120,'Period 2 - QAIC',cex=1.00)

par(fig=c(0.4,.6,.4,.8), mar=c(.5,0.5,0.5,0), new=T) 
plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
polygon(japan.poly2,col=res$rrcolors$color.qaic[,3],border=F)
segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
text(355,4120,'Period 3 - QAIC',cex=1.00)

par(fig=c(0.6,.8,.4,.8), mar=c(.5,0.5,0.5,0), new=T)
plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
polygon(japan.poly2,col=res$rrcolors$color.qaic[,4],border=F)
segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
text(355,4120,'Period 4 - QAIC',cex=1.00)

par(fig=c(0.8,1,.4,.8), mar=c(.5,0.5,0.5,0), new=T)
plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
polygon(japan.poly2,col=res$rrcolors$color.qaic[,5],border=F)
segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
text(355,4120,'Period 5 - QAIC',cex=1.00)


#Maps of AICc Path

par(fig=c(0,.2,.2,.6), mar=c(.5,0.5,0.5,0), new=T)
plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
polygon(japan.poly2,col=res$rrcolors$color.qaicc[,1],border=F)
segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
text(355,4120,'Period 1 - QAICc',cex=1.00)

par(fig=c(0.2,.4,.2,.6), mar=c(.5,0.5,0.5,0), new=T) 
plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
polygon(japan.poly2,col=res$rrcolors$color.qaicc[,2],border=F)
segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
text(355,4120,'Period 2 - QAICc',cex=1.00)

par(fig=c(0.4,.6,.2,.6), mar=c(.5,0.5,0.5,0), new=T) 
plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
polygon(japan.poly2,col=res$rrcolors$color.qaicc[,3],border=F)
segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
text(355,4120,'Period 3 - QAICc',cex=1.00)

par(fig=c(0.6,.8,.2,.6), mar=c(.5,0.5,0.5,0), new=T)
plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
polygon(japan.poly2,col=res$rrcolors$color.qaicc[,4],border=F)
segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
text(355,4120,'Period 4 - QAICc',cex=1.00)

par(fig=c(0.8,1,.2,.6), mar=c(.5,0.5,0.5,0), new=T)
plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
polygon(japan.poly2,col=res$rrcolors$color.qaicc[,5],border=F)
segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
text(355,4120,'Period 5 - QAICc',cex=1.00)


#Maps of BIC Path

par(fig=c(0,.2,0,.4), mar=c(.5,0.5,0.5,0), new=T)
plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
polygon(japan.poly2,col=res$rrcolors$color.qbic[,1],border=F)
segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
text(355,4120,'Period 1 - QBI',cex=1.00)

par(fig=c(0.2,.4,0,.4), mar=c(.5,0.5,0.5,0), new=T)
plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
polygon(japan.poly2,col=res$rrcolors$color.qbic[,2],border=F)
segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
text(355,4120,'Period 2 - QBIC',cex=1.00)

par(fig=c(0.4,.6,0,.4), mar=c(.5,0.5,0.5,0), new=T)
plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
polygon(japan.poly2,col=res$rrcolors$color.qbic[,3],border=F)
segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
text(355,4120,'Period 3 - QBIC',cex=1.00)


par(fig=c(0.6,.8,0,.4), mar=c(.5,0.5,0.5,0), new=T)
plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
polygon(japan.poly2,col=res$rrcolors$color.qbic[,4],border=F)
segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
text(355,4120,'Period 4 - QBIC',cex=1.00)


par(fig=c(0.8,1,0,.4), mar=c(.5,0.5,0.5,0), new=T)
plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
polygon(japan.poly2,col=res$rrcolors$color.qbic[,5],border=F)
segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
text(355,4120,'Period 5 - QBIC',cex=1.00)

#Turn off pdf development
dev.off()

########################################################################################################
########################################################################################################
########################################################################################################

####################################################
#rmax 30, time = 2,3, rr.ratio = 1.25, center = 200, r = 11km
####################################################
#Set Some Initial Conditions
x=dframe2$utmx/1000
y=dframe2$utmy/1000
rMax=30 
Time= 5
nsim=100
center=200
radius=11
period_start= c(2,3)
risk.ratio=1.25




res <- clust.sim(x,y,rMax,dframe$period, dframe$expdeath, dframe$death, Time, spacetime=TRUE, 
                 nsim,center, radius, risk.ratio, period_start, colors=TRUE, utm=TRUE, byrow=TRUE)

#make maps
pdfname <- paste0("figures/simulations/japan","_","rMax","_",rMax,"nsim","_",nsim, "center","_",center,"radius","_",radius, "start","_",period_start, "rr","_",risk.ratio,".pdf")
pdf(pdfname, height=11, width=10)

#Maps of Observed Counts
par(fig=c(0,.2,.6,1), mar=c(.5,0.5,0.5,0))
plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
polygon(japan.poly2,col=res$rrcolors$colors.obs[,1],border=F)
segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
text(355,4120,'Period 1 - Obs',cex=1.00)

par(fig=c(0.2,.4,.6,1), mar=c(.5,0.5,0.5,0), new=T)
plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
polygon(japan.poly2,col=res$rrcolors$colors.obs[,2],border=F)
segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
text(355,4120,'Period 2 - Obs',cex=1.00)

par(fig=c(0.4,.6,.6,1), mar=c(.5,0.5,0.5,0), new=T)
plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
polygon(japan.poly2,col=res$rrcolors$colors.obs[,3],border=F)
segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
text(355,4120,'Period 3 - Obs',cex=1.00)

par(fig=c(0.6,.8,.6,1), mar=c(.5,0.5,0.5,0), new=T)
plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
polygon(japan.poly2,col=res$rrcolors$colors.obs[,4],border=F)
segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
text(355,4120,'Period 4 - Obs',cex=1.00)

par(fig=c(0.8,1,.6,1), mar=c(.5,0.5,0.5,0), new=T)
plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
polygon(japan.poly2,col=res$rrcolors$colors.obs[,5],border=F)
segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
text(355,4120,'Period 5 - Obs',cex=1.00)


#Maps of AIC Path

par(fig=c(0,.2,.4,.8), mar=c(.5,0.5,0.5,0), new=T)
plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
polygon(japan.poly2,col=res$rrcolors$color.qaic[,1],border=F)
segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
text(355,4120,'Period 1 - QAIC',cex=1.00)

par(fig=c(0.2,.4,.4,.8), mar=c(.5,0.5,0.5,0), new=T)   
plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
polygon(japan.poly2,col=res$rrcolors$color.qaic[,2],border=F)
segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
text(355,4120,'Period 2 - QAIC',cex=1.00)

par(fig=c(0.4,.6,.4,.8), mar=c(.5,0.5,0.5,0), new=T) 
plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
polygon(japan.poly2,col=res$rrcolors$color.qaic[,3],border=F)
segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
text(355,4120,'Period 3 - QAIC',cex=1.00)

par(fig=c(0.6,.8,.4,.8), mar=c(.5,0.5,0.5,0), new=T)
plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
polygon(japan.poly2,col=res$rrcolors$color.qaic[,4],border=F)
segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
text(355,4120,'Period 4 - QAIC',cex=1.00)

par(fig=c(0.8,1,.4,.8), mar=c(.5,0.5,0.5,0), new=T)
plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
polygon(japan.poly2,col=res$rrcolors$color.qaic[,5],border=F)
segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
text(355,4120,'Period 5 - QAIC',cex=1.00)


#Maps of AICc Path

par(fig=c(0,.2,.2,.6), mar=c(.5,0.5,0.5,0), new=T)
plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
polygon(japan.poly2,col=res$rrcolors$color.qaicc[,1],border=F)
segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
text(355,4120,'Period 1 - QAICc',cex=1.00)

par(fig=c(0.2,.4,.2,.6), mar=c(.5,0.5,0.5,0), new=T) 
plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
polygon(japan.poly2,col=res$rrcolors$color.qaicc[,2],border=F)
segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
text(355,4120,'Period 2 - QAICc',cex=1.00)

par(fig=c(0.4,.6,.2,.6), mar=c(.5,0.5,0.5,0), new=T) 
plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
polygon(japan.poly2,col=res$rrcolors$color.qaicc[,3],border=F)
segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
text(355,4120,'Period 3 - QAICc',cex=1.00)

par(fig=c(0.6,.8,.2,.6), mar=c(.5,0.5,0.5,0), new=T)
plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
polygon(japan.poly2,col=res$rrcolors$color.qaicc[,4],border=F)
segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
text(355,4120,'Period 4 - QAICc',cex=1.00)

par(fig=c(0.8,1,.2,.6), mar=c(.5,0.5,0.5,0), new=T)
plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
polygon(japan.poly2,col=res$rrcolors$color.qaicc[,5],border=F)
segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
text(355,4120,'Period 5 - QAICc',cex=1.00)


#Maps of BIC Path

par(fig=c(0,.2,0,.4), mar=c(.5,0.5,0.5,0), new=T)
plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
polygon(japan.poly2,col=res$rrcolors$color.qbic[,1],border=F)
segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
text(355,4120,'Period 1 - QBI',cex=1.00)

par(fig=c(0.2,.4,0,.4), mar=c(.5,0.5,0.5,0), new=T)
plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
polygon(japan.poly2,col=res$rrcolors$color.qbic[,2],border=F)
segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
text(355,4120,'Period 2 - QBIC',cex=1.00)

par(fig=c(0.4,.6,0,.4), mar=c(.5,0.5,0.5,0), new=T)
plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
polygon(japan.poly2,col=res$rrcolors$color.qbic[,3],border=F)
segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
text(355,4120,'Period 3 - QBIC',cex=1.00)


par(fig=c(0.6,.8,0,.4), mar=c(.5,0.5,0.5,0), new=T)
plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
polygon(japan.poly2,col=res$rrcolors$color.qbic[,4],border=F)
segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
text(355,4120,'Period 4 - QBIC',cex=1.00)


par(fig=c(0.8,1,0,.4), mar=c(.5,0.5,0.5,0), new=T)
plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
polygon(japan.poly2,col=res$rrcolors$color.qbic[,5],border=F)
segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
text(355,4120,'Period 5 - QBIC',cex=1.00)

#Turn off pdf development
dev.off()



########################################################################################################
########################################################################################################
########################################################################################################

####################################################
#rmax 30, time = 4, rr.ratio = 1.5, center = 100, r = 18km, spaceonly
####################################################
#Set Some Initial Conditions
x=dframe2$utmx/1000
y=dframe2$utmy/1000
rMax=30 
Time= 1
nsim=100
center=100
radius=18
period_start= 1
risk.ratio=1.5

#create average dataframe for spaceonly

death <- with(dframe, tapply(death, id, function(x) round(mean(x))))
expdeath <- with(dframe, tapply(expdeath, id, function(x) mean(x)))
df <- cbind.data.frame(id = unique(dframe$id), period = rep("1", length(unique(dframe$id))), death = death, expdeath=expdeath)

res <- clust.sim(x,y,rMax,df$period, df$expdeath, df$death, Time, spacetime=FALSE, 
                 nsim,center, radius, risk.ratio, period_start, colors=TRUE, utm=TRUE, byrow=TRUE)

#make maps
pdfname <- paste0("figures/simulations/japan","_","rMax","_",rMax,"nsim","_",nsim, "center","_",center,"radius","_",radius, "start","_",period_start, "rr","_",risk.ratio,"_spaceonly",".pdf")
pdf(pdfname, height=11, width=10)

#Maps of Observed Counts
par(fig=c(0,.2,.6,1), mar=c(.5,0.5,0.5,0))
plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
polygon(japan.poly2,col=res$rrcolors$colors.obs[,1],border=F)
segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
text(355,4120,'Period 1 - Obs',cex=1.00)

#Maps of AIC Path

par(fig=c(0,.2,.4,.8), mar=c(.5,0.5,0.5,0), new=T)
plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
polygon(japan.poly2,col=res$rrcolors$color.qaic[,1],border=F)
segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
text(355,4120,'Period 1 - QAIC',cex=1.00)

#Maps of AICc Path

par(fig=c(0,.2,.2,.6), mar=c(.5,0.5,0.5,0), new=T)
plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
polygon(japan.poly2,col=res$rrcolors$color.qaicc[,1],border=F)
segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
text(355,4120,'Period 1 - QAICc',cex=1.00)


#Maps of BIC Path

par(fig=c(0,.2,0,.4), mar=c(.5,0.5,0.5,0), new=T)
plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
polygon(japan.poly2,col=res$rrcolors$color.qbic[,1],border=F)
segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
text(355,4120,'Period 1 - QBIC',cex=1.00)


#Turn off pdf development
dev.off()


########################################################################################################
########################################################################################################
########################################################################################################

####################################################
#rmax 30, time = 1,2, rr.ratio = 1.1, center = 75, r = 18km
####################################################
#Set Some Initial Conditions
x=dframe2$utmx/1000
y=dframe2$utmy/1000
rMax=30 
Time=5
nsim=100
center=75
radius=18
period_start=c(1,2)
risk.ratio=1.1




result <- clust.sim(x,y,rMax,dframe$period, dframe$expdeath, dframe$death, Time, spacetime=TRUE, 
                 nsim,center, radius, risk.ratio, period_start, colors=TRUE, utm=TRUE, byrow=TRUE)



###did it find at least 1 cell in the cluster?

set <- detect.set(result$lassoresult, result$init.vec, result$rr.mat, Time=5)
index_truth <- sapply(1:Time, function(i) which( result$rr.mat[,i] !=1))
index <- which( result$rr.mat !=1)

incluster <-detect.incluster.aic(result$lassoresult, result$init.vec, result$rr.mat, set, period=c(1,2),Time=5)





# 
# 
# #diagnostics
# set <- detect.set(res$lassoresult, res$init.vec, res$rr.mat, Time=5)
# 
# det <- detect(res$lassoresult, res$init.vec, res$rr.mat, period_start=1,period_end= 2, multi_period = FALSE, IC="aic", Time)
# 
# detin <- detect.incluster(res$lassoresult, res$init.vec, res$rr.mat, set,period_start=1, period_end=2, multi_period=FALSE, IC = "ic", Time=5) #this works
# detfalse <- detect.falsecluster(result$lassoresult, result$init.vec, result$rr.mat, set,period_start=1, period_end=2, multi_period=FALSE, IC = "aic", Time=5)
# detfalse <- detect.falsecluster(result$lassoresult, result$init.vec, result$rr.mat, set,period_start=1, period_end=2, multi_period=FALSE, IC = "ic", Time=5)
# 
# 
# detbk <- detect.inbackground(result$lassoresult, result$init.vec, result$rr.mat, set,period_start=1, period_end=2, multi_period=FALSE, IC = "bic", Time=5)


#TODO in background - clarify
#borders - best by plot?
#AIC, BIC should work
#Pois - space, quasi-pois - ST













