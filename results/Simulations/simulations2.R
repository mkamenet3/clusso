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
library(spdep)
library(sp)

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

####################################################
#Create Spatial Polygons and Create Matrix of Neighborhoods
####################################################
ind <- which(is.na(dframe.poly2[,2]))
names <- as.factor(seq(1,n,1))
coord.system <- '+proj=utm'
nrepeats <- rep(1, sum(is.na(dframe.poly2[, 2])))
poly <- as.matrix(dframe.poly2[, 2:3])
area.names <- as.character(1:length(nrepeats))
na.index <- which(is.na(poly[, 1]))
n <- length(nrepeats)
list.polygon <- NULL
list.polygon <- list(Polygon(poly[1:(na.index[1] - 1), ], 
                             hole = FALSE))
for (i in 1:(length(na.index) - 1)) {
    list.polygon <- c(list.polygon, list(Polygon(poly[(na.index[i] + 
                                                           1):(na.index[i + 1] - 1), ], hole = FALSE)))
}
list.polygons <- NULL
start <- 1
for (i in 1:length(nrepeats)) {
    end <- start + nrepeats[i] - 1
    temp.polygon <- NULL
    for (j in start:end) {
        print(c(i,j))
        temp.polygon <- c(temp.polygon, list(list.polygon[[j]]))
    }
    list.polygons <- c(list.polygons, list(Polygons(temp.polygon, 
                                                    ID = area.names[i])))
    start <- end + 1
}
mypoly <- SpatialPolygons(list.polygons)
plot(mypoly)
nb <- poly2nb(mypoly)


########################################################################################################
########################################################################################################
########################################################################################################

####################################################
#NULL MODEL - Quasi-Poisson
####################################################
#Set Some Initial Conditions
x=dframe2$utmx/1000
y=dframe2$utmy/1000
rMax=30 
Time=5  
nsim=100
center=200
radius=18
period_start=1
period_end = 5
risk.ratio=1




res <- clust.sim(x,y,rMax,dframe$period, dframe$expdeath, dframe$death, Time, spacetime=TRUE, pois=FALSE,
                 nsim,center, radius, risk.ratio, period_start, colors=TRUE, utm=TRUE, byrow=TRUE)

#save results
filename <- paste0("SimulationOutput/sim","_","center","_",center,"radius","_",radius, "start","_",period_start, "rr","_",risk.ratio,"NULL",".RData")
save(res, file = filename)



###did it find at least 1 cell in the cluster?
set <- detect.set(res$lassoresult, res$init.vec, res$rr.mat, Time=5, nb, x, y, rMax, 0, radius)
(incluster <- detect.incluster(res$lassoresult, res$init.vec, res$rr.mat, set, period_start = 1, 
                               period_end = 5, multi_period = TRUE,Time=5, nsim,nb, x, y, rMax, 0, radius, IC = "ic"))
#make maps
pdfname <- paste0("figures/simulations/japan","_","rMax","_",rMax,"nsim","_",nsim, "center","_",center,"radius","_",radius, "start","_", "rr","_",risk.ratio,"NULL",".pdf")
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
text(355,4120,'Period 1 - QBIC',cex=1.00)

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


####################################################
#NULL MODEL - POISSON
####################################################
#Set Some Initial Conditions
x=dframe2$utmx/1000
y=dframe2$utmy/1000
rMax=30 
Time=5  
nsim=100
center=200
radius=18
period_start=1
period_end = 5
risk.ratio=1




res <- clust.sim(x,y,rMax,dframe$period, dframe$expdeath, dframe$death, Time, spacetime=TRUE, pois=TRUE,
                 nsim,center, radius, risk.ratio, period_start, colors=TRUE, utm=TRUE, byrow=TRUE)

#save results
filename <- paste0("SimulationOutput/sim","_","center","_",center,"radius","_",radius, "start","_",period_start, "rr","_",risk.ratio,"NULLPOISSON",".RData")
save(res, file = filename)



###did it find at least 1 cell in the cluster?
set <- detect.set(res$lassoresult, res$init.vec, res$rr.mat, Time=5, nb, x, y, rMax, 0, radius)
(incluster <- detect.incluster(res$lassoresult, res$init.vec, res$rr.mat, set, period_start = 1, 
                               period_end = 5, multi_period = TRUE,Time=5, nsim,nb, x, y, rMax, 0, radius, IC = "ic"))
#make maps
pdfname <- paste0("figures/simulations/japan","_","rMax","_",rMax,"nsim","_",nsim, "center","_",center,"radius","_",radius, "start","_", "rr","_",risk.ratio,"NULLPOISSON",".pdf")
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
text(355,4120,'Period 1 - QBIC',cex=1.00)

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
#NULL spaceonly - Quasi-Poisson
####################################################
#Set Some Initial Conditions
x=dframe2$utmx/1000
y=dframe2$utmy/1000
rMax=30 
Time=1
nsim=100
center=200
radius=18
period_start=1
risk.ratio=1


#create average dataframe for spaceonly

death <- with(dframe, tapply(death, id, function(x) round(mean(x))))
expdeath <- with(dframe, tapply(expdeath, id, function(x) mean(x)))
df <- cbind.data.frame(id = unique(dframe$id), period = rep("1", length(unique(dframe$id))), death = death, expdeath=expdeath)

res <- clust.sim(x,y,rMax,df$period, df$expdeath, df$death, Time, spacetime=FALSE, pois=FALSE, 
                 nsim,center, radius, risk.ratio, period_start, colors=TRUE, utm=TRUE, byrow=TRUE)
#save results
filename <- paste0("SimulationOutput/sim","_","center","_",center,"radius","_",radius, "start","_","rr","_",risk.ratio,"spaceonlyQPOIS",".RData")
save(res, file = filename)

###did it find at least 1 cell in the cluster?
set <- detect.set(res$lassoresult, res$init.vec, res$rr.mat, Time=1, nb, x, y, rMax, 0, radius)
(incluster <- detect.incluster(res$lassoresult, res$init.vec, res$rr.mat, set, period_start = 1, 
                               period_end = 1, multi_period = FALSE,Time=1, nsim,nb, x, y, rMax, 0, radius, IC = "ic", space=FALSE))



#make maps
pdfname <- paste0("figures/simulations/sim","_","center","_",center,"radius","_",radius, "start","_","rr","_",risk.ratio,"spaceonlyQPOIS",".pdf")
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

####################################################
#NULL spaceonly - Poisson
####################################################
#Set Some Initial Conditions
x=dframe2$utmx/1000
y=dframe2$utmy/1000
rMax=30 
Time=1
nsim=100
center=200
radius=18
period_start=1
risk.ratio=1


#create average dataframe for spaceonly

death <- with(dframe, tapply(death, id, function(x) round(mean(x))))
expdeath <- with(dframe, tapply(expdeath, id, function(x) mean(x)))
df <- cbind.data.frame(id = unique(dframe$id), period = rep("1", length(unique(dframe$id))), death = death, expdeath=expdeath)

res <- clust.sim(x,y,rMax,df$period, df$expdeath, df$death, Time, spacetime=FALSE, pois=TRUE, 
                 nsim,center, radius, risk.ratio, period_start, colors=TRUE, utm=TRUE, byrow=TRUE)
#save results
filename <- paste0("SimulationOutput/sim","_","center","_",center,"radius","_",radius, "start","_","rr","_",risk.ratio,"spaceonlyPOISSONONLY",".RData")
save(res, file = filename)

###did it find at least 1 cell in the cluster?
set <- detect.set(res$lassoresult, res$init.vec, res$rr.mat, Time=1, nb, x, y, rMax, 0, radius)
(incluster <- detect.incluster(res$lassoresult, res$init.vec, res$rr.mat, set, period_start = 1, 
                               period_end = 1, multi_period = FALSE,Time=1, nsim,nb, x, y, rMax, 0, radius, IC = "ic", space=FALSE, nullmod=TRUE))



#make maps
pdfname <- paste0("figures/simulations/sim","_","center","_",center,"radius","_",radius, "start","_","rr","_",risk.ratio,"spaceonlyPOISSONONLY",".pdf")
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
#END NULL MODELS, BEGIN CLUSTER SIMULATIONS
########################################################################################################
########################################################################################################
########################################################################################################

####################################################
#Under SIM
# sum(as.vector(matrix(res$init.vec$E0[[1]],ncol=5)[c(21,9,23),c(4:5)]))
# [1] 29.21828
# sd(as.vector(matrix(res$init.vec$E0[[1]],ncol=5)[c(21,9,23),c(4:5)]))
# [1] 3.246741

#Under Null
# sum(as.vector(matrix(res$init.vec$E0[[1]],ncol=5)[c(21,9,23),c(4:5)]))
# [1] 26.58643
# sd(as.vector(matrix(res$init.vec$E0[[1]],ncol=5)[c(21,9,23),c(4:5)]))
# [1] 2.954276

##CALC SAYS SAMPLE SIZE OF 4 WILL GIVE US POWER OF 50%
####################################################
#Set Some Initial Conditions
x=dframe2$utmx/1000
y=dframe2$utmy/1000
rMax=30 
Time=5  
nsim=100
center=21
radius=18
period_start=4
risk.ratio=1.2




res <- clust.sim(x,y,rMax,dframe$period, dframe$expdeath, dframe$death, Time, spacetime=TRUE, pois=FALSE,
                  nsim,center, radius, risk.ratio, period_start, colors=TRUE, utm=TRUE, byrow=TRUE)

#save results
filename <- paste0("SimulationOutput/sim","_","center","_",center,"radius","_",radius, "start","_","rr","_","2",".RData")
save(res, file = filename)

###did it find at least 1 cell in the cluster?
set <- detect.set(res$lassoresult, res$init.vec, res$rr.mat, Time=5, nb, x, y, rMax, center, radius)
(incluster <- detect.incluster(res$lassoresult, res$init.vec, res$rr.mat, set, period_start = 4, 
                               period_end = 5, multi_period = TRUE,Time=5, nsim,nb, x, y, rMax, center, radius, IC = "ic"))
#make maps
pdfname <- paste0("figures/simulations/sim","_","center","_",center,"radius","_",radius, "start","_","rr","_","2_differensum_times",".pdf")
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
text(355,4120,'Period 1 - QBIC',cex=1.00)

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




res <- clust.sim(x,y,rMax,dframe$period, dframe$expdeath, dframe$death, Time, spacetime=TRUE, pois=FALSE
                 nsim,center, radius, risk.ratio, period_start, colors=TRUE, utm=TRUE, byrow=TRUE)
#save results
filename <- paste0("SimulationOutput/sim","_","center","_",center,"radius","_",radius, "start","_","rr","_","1_25",".RData")
save(res, file = filename)

###did it find at least 1 cell in the cluster?
set <- detect.set(res$lassoresult, res$init.vec, res$rr.mat, Time=5, nb, x, y, rMax, center, radius)
(incluster <- detect.incluster(res$lassoresult, res$init.vec, res$rr.mat, set, period_start = 2, 
                              period_end = 3, multi_period = TRUE,Time=5, nsim,nb, x, y, rMax, center, radius, IC = "ic"))


#make maps
pdfname <- paste0("figures/simulations/japan","_","rMax","_",rMax,"nsim","_",nsim, "center","_",center,"radius","_",radius, "start","_",period_start, "rr","_","1_25",".pdf")
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
text(355,4120,'Period 1 - QBIC',cex=1.00)

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
#rmax 30, time = spaceonly, rr.ratio = 1.5, center = 100, r = 18km, spaceonly
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

res <- clust.sim(x,y,rMax,df$period, df$expdeath, df$death, Time, spacetime=FALSE, pois=FALSE,
                 nsim,center, radius, risk.ratio, period_start, colors=TRUE, utm=TRUE, byrow=TRUE)
#save results
filename <- paste0("SimulationOutput/sim","_","center","_",center,"radius","_",radius, "start","_","rr","_",risk.ratio,"spaceonly",".RData")
save(res, file = filename)

###did it find at least 1 cell in the cluster?
set <- detect.set(res$lassoresult, res$init.vec, res$rr.mat, Time=1, nb, x, y, rMax, center, radius)
(incluster <- detect.incluster(res$lassoresult, res$init.vec, res$rr.mat, set, period_start = 1, 
                               period_end = 1, multi_period = FALSE,Time=1, nsim,nb, x, y, rMax, center, radius, IC = "ic", space=TRUE))



#make maps
pdfname <- paste0("figures/simulations/sim","_","center","_",center,"radius","_",radius, "start","_","rr","_",risk.ratio,"spaceonly",".pdf")
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




res <- clust.sim(x,y,rMax,dframe$period, dframe$expdeath, dframe$death, Time, spacetime=TRUE, pois=FALSE,
                 nsim,center, radius, risk.ratio, period_start, colors=TRUE, utm=TRUE, byrow=TRUE)

#save results
filename <- paste0("SimulationOutput/sim","_","center","_",center,"radius","_",radius, "start","_","rr","_","1_1",".RData")
save(res, file = filename)

###did it find at least 1 cell in the cluster?
set <- detect.set(res$lassoresult, res$init.vec, res$rr.mat, Time=5, nb, x, y, rMax, center, radius)
(incluster <- detect.incluster(res$lassoresult, res$init.vec, res$rr.mat, set, period_start = 1, 
                               period_end = 2, multi_period = TRUE,Time=5, nsim,nb, x, y, rMax, center, radius, IC = "ic"))


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
text(355,4120,'Period 1 - QBIC',cex=1.00)

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
#rmax 30, time = 3,4,5, rr.ratio = 2, center = 50, r = 18km
####################################################
#Set Some Initial Conditions
x=dframe2$utmx/1000
y=dframe2$utmy/1000
rMax=30 
Time=5
nsim=100
center=50
radius= 9
period_start=c(1:4)
risk.ratio=2




res <- clust.sim(x,y,rMax,dframe$period, dframe$expdeath, dframe$death, Time, spacetime=TRUE, pois=FALSE,
                 nsim,center, radius, risk.ratio, period_start, colors=TRUE, utm=TRUE, byrow=TRUE)

#save results
filename <- paste0("SimulationOutput/sim","_","center","_",center,"radius","_",radius, "start","_","rr","_",risk.ratio,".RData")
save(res, file = filename)

###did it find at least 1 cell in the cluster?
set <- detect.set(res$lassoresult, res$init.vec, res$rr.mat, Time=5, nb, x, y, rMax, center, radius)
(incluster <- detect.incluster(res$lassoresult, res$init.vec, res$rr.mat, set, period_start = 1, 
                               period_end = 4, multi_period = TRUE,Time=5, nsim,nb, x, y, rMax, center, radius, IC = "ic"))

#make maps
pdfname <- paste0("figures/Simulations/sim","_","center","_",center,"radius","_",radius, "start","_","rr","_",risk.ratio,".pdf")
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
text(355,4120,'Period 1 - QBIC',cex=1.00)

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
#rmax 30, time = 3,4,5, rr.ratio = 1.5, center = 50, r = 18km
####################################################
#Set Some Initial Conditions
x=dframe2$utmx/1000
y=dframe2$utmy/1000
rMax=30 
Time=5
nsim=100
center=25
radius=18
period_start=c(1:5)
risk.ratio=1.1


res <- clust.sim(x,y,rMax,dframe$period, dframe$expdeath, dframe$death, Time, spacetime=TRUE, pois=FALSE,
                 nsim,center, radius, risk.ratio, period_start, colors=TRUE, utm=TRUE, byrow=TRUE)

#save results
filename <- paste0("SimulationOutput/sim","_","center","_",center,"radius","_",radius, "start","_","rr","_","1_1",".RData")
save(res, file = filename)

###did it find at least 1 cell in the cluster?
set <- detect.set(res$lassoresult, res$init.vec, res$rr.mat, Time=5, nb, x, y, rMax, center, radius)
(incluster <- detect.incluster(res$lassoresult, res$init.vec, res$rr.mat, set, period_start = 1, 
                               period_end = 5, multi_period = TRUE,Time=5, nsim,nb, x, y, rMax, center, radius, IC = "ic"))

#make maps
pdfname <- paste0("figures/simulations/sim","_","center","_",center,"radius","_",radius, "start","_","rr","_","1_1",".pdf")
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
text(355,4120,'Period 1 - QBIC',cex=1.00)

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
# 
# ####################################################
# #rmax 30, time = 5, rr.ratio = 1.25, center = 50, r = 18km
# ####################################################
# #Set Some Initial Conditions
# x=dframe2$utmx/1000
# y=dframe2$utmy/1000
# rMax=30 
# Time=5
# nsim=100
# center=50
# radius=18
# period_start=c(5:5)
# risk.ratio=1.25
# 
# 
# 
# 
# res <- clust.sim(x,y,rMax,dframe$period, dframe$expdeath, dframe$death, Time, spacetime=TRUE, 
#                  nsim,center, radius, risk.ratio, period_start, colors=TRUE, utm=TRUE, byrow=TRUE)
# 
# #save results
# filename <- paste0("SimulationOutput/japan","_","center","_",center,"radius","_",radius, "start","_",period_start, "rr","_","1_25",".RData")
# save(res, file = filename)
# 
# ###did it find at least 1 cell in the cluster?
# set <- detect.set(res$lassoresult, res$init.vec, res$rr.mat, Time=5)
# incluster7 <- detect.incluster.ic(res$lassoresult, res$init.vec, res$rr.mat, set, period=c(5:5),Time=5)
# 
# #make maps
# pdfname <- paste0("figures/simulations/japan","_","rMax","_",rMax,"nsim","_",nsim, "center","_",center,"radius","_",radius, "start","_",period_start, "rr","_","1_25",".pdf")
# pdf(pdfname, height=11, width=10)
# 
# #Maps of Observed Counts
# par(fig=c(0,.2,.6,1), mar=c(.5,0.5,0.5,0))
# plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
# polygon(japan.poly2,col=res$rrcolors$colors.obs[,1],border=F)
# segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
# text(355,4120,'Period 1 - Obs',cex=1.00)
# 
# par(fig=c(0.2,.4,.6,1), mar=c(.5,0.5,0.5,0), new=T)
# plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
# polygon(japan.poly2,col=res$rrcolors$colors.obs[,2],border=F)
# segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
# text(355,4120,'Period 2 - Obs',cex=1.00)
# 
# par(fig=c(0.4,.6,.6,1), mar=c(.5,0.5,0.5,0), new=T)
# plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
# polygon(japan.poly2,col=res$rrcolors$colors.obs[,3],border=F)
# segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
# text(355,4120,'Period 3 - Obs',cex=1.00)
# 
# par(fig=c(0.6,.8,.6,1), mar=c(.5,0.5,0.5,0), new=T)
# plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
# polygon(japan.poly2,col=res$rrcolors$colors.obs[,4],border=F)
# segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
# text(355,4120,'Period 4 - Obs',cex=1.00)
# 
# par(fig=c(0.8,1,.6,1), mar=c(.5,0.5,0.5,0), new=T)
# plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
# polygon(japan.poly2,col=res$rrcolors$colors.obs[,5],border=F)
# segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
# text(355,4120,'Period 5 - Obs',cex=1.00)
# 
# 
# #Maps of AIC Path
# 
# par(fig=c(0,.2,.4,.8), mar=c(.5,0.5,0.5,0), new=T)
# plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
# polygon(japan.poly2,col=res$rrcolors$color.qaic[,1],border=F)
# segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
# text(355,4120,'Period 1 - QAIC',cex=1.00)
# 
# par(fig=c(0.2,.4,.4,.8), mar=c(.5,0.5,0.5,0), new=T)   
# plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
# polygon(japan.poly2,col=res$rrcolors$color.qaic[,2],border=F)
# segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
# text(355,4120,'Period 2 - QAIC',cex=1.00)
# 
# par(fig=c(0.4,.6,.4,.8), mar=c(.5,0.5,0.5,0), new=T) 
# plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
# polygon(japan.poly2,col=res$rrcolors$color.qaic[,3],border=F)
# segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
# text(355,4120,'Period 3 - QAIC',cex=1.00)
# 
# par(fig=c(0.6,.8,.4,.8), mar=c(.5,0.5,0.5,0), new=T)
# plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
# polygon(japan.poly2,col=res$rrcolors$color.qaic[,4],border=F)
# segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
# text(355,4120,'Period 4 - QAIC',cex=1.00)
# 
# par(fig=c(0.8,1,.4,.8), mar=c(.5,0.5,0.5,0), new=T)
# plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
# polygon(japan.poly2,col=res$rrcolors$color.qaic[,5],border=F)
# segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
# text(355,4120,'Period 5 - QAIC',cex=1.00)
# 
# 
# #Maps of AICc Path
# 
# par(fig=c(0,.2,.2,.6), mar=c(.5,0.5,0.5,0), new=T)
# plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
# polygon(japan.poly2,col=res$rrcolors$color.qaicc[,1],border=F)
# segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
# text(355,4120,'Period 1 - QAICc',cex=1.00)
# 
# par(fig=c(0.2,.4,.2,.6), mar=c(.5,0.5,0.5,0), new=T) 
# plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
# polygon(japan.poly2,col=res$rrcolors$color.qaicc[,2],border=F)
# segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
# text(355,4120,'Period 2 - QAICc',cex=1.00)
# 
# par(fig=c(0.4,.6,.2,.6), mar=c(.5,0.5,0.5,0), new=T) 
# plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
# polygon(japan.poly2,col=res$rrcolors$color.qaicc[,3],border=F)
# segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
# text(355,4120,'Period 3 - QAICc',cex=1.00)
# 
# par(fig=c(0.6,.8,.2,.6), mar=c(.5,0.5,0.5,0), new=T)
# plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
# polygon(japan.poly2,col=res$rrcolors$color.qaicc[,4],border=F)
# segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
# text(355,4120,'Period 4 - QAICc',cex=1.00)
# 
# par(fig=c(0.8,1,.2,.6), mar=c(.5,0.5,0.5,0), new=T)
# plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
# polygon(japan.poly2,col=res$rrcolors$color.qaicc[,5],border=F)
# segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
# text(355,4120,'Period 5 - QAICc',cex=1.00)
# 
# 
# #Maps of BIC Path
# 
# par(fig=c(0,.2,0,.4), mar=c(.5,0.5,0.5,0), new=T)
# plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
# polygon(japan.poly2,col=res$rrcolors$color.qbic[,1],border=F)
# segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
# text(355,4120,'Period 1 - QBIC',cex=1.00)
# 
# par(fig=c(0.2,.4,0,.4), mar=c(.5,0.5,0.5,0), new=T)
# plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
# polygon(japan.poly2,col=res$rrcolors$color.qbic[,2],border=F)
# segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
# text(355,4120,'Period 2 - QBIC',cex=1.00)
# 
# par(fig=c(0.4,.6,0,.4), mar=c(.5,0.5,0.5,0), new=T)
# plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
# polygon(japan.poly2,col=res$rrcolors$color.qbic[,3],border=F)
# segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
# text(355,4120,'Period 3 - QBIC',cex=1.00)
# 
# 
# par(fig=c(0.6,.8,0,.4), mar=c(.5,0.5,0.5,0), new=T)
# plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
# polygon(japan.poly2,col=res$rrcolors$color.qbic[,4],border=F)
# segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
# text(355,4120,'Period 4 - QBIC',cex=1.00)
# 
# 
# par(fig=c(0.8,1,0,.4), mar=c(.5,0.5,0.5,0), new=T)
# plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
# polygon(japan.poly2,col=res$rrcolors$color.qbic[,5],border=F)
# segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
# text(355,4120,'Period 5 - QBIC',cex=1.00)
# 
# #Turn off pdf development
# dev.off()

# 
# ########################################################################################################
# ########################################################################################################
# ########################################################################################################
# 
# ####################################################
# #rmax 30, time = 4,5, rr.ratio = 1.5, center = 50, r = 9.5km
# ####################################################
# #Set Some Initial Conditions
# x=dframe2$utmx/1000
# y=dframe2$utmy/1000
# rMax=30 
# Time=5
# nsim=100
# center=50
# radius=9.5
# period_start=c(4,5)
# risk.ratio=1.5
# 
# 
# 
# 
# res <- clust.sim(x,y,rMax,dframe$period, dframe$expdeath, dframe$death, Time, spacetime=TRUE, 
#                  nsim,center, radius, risk.ratio, period_start, colors=TRUE, utm=TRUE, byrow=TRUE)
# 
# #save results
# filename <- paste0("SimulationOutput/japan","_","center","_",center,"radius","_","9_5", "start","_",period_start, "rr","_","1_5",".RData")
# save(res, file = filename)
# 
# ###did it find at least 1 cell in the cluster?
# set <- detect.set(result$lassoresult, result$init.vec, result$rr.mat, Time=5)
# incluster8 <- detect.incluster.ic(result$lassoresult, result$init.vec, result$rr.mat, set, period=c(4:5),Time=5)
# 
# #make maps
# pdfname <- paste0("figures/simulations/japan","_","rMax","_",rMax,"nsim","_",nsim, "center","_",center,"radius","_","9_5", "start","_",period_start, "rr","_","1_5",".pdf")
# pdf(pdfname, height=11, width=10)
# 
# #Maps of Observed Counts
# par(fig=c(0,.2,.6,1), mar=c(.5,0.5,0.5,0))
# plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
# polygon(japan.poly2,col=res$rrcolors$colors.obs[,1],border=F)
# segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
# text(355,4120,'Period 1 - Obs',cex=1.00)
# 
# par(fig=c(0.2,.4,.6,1), mar=c(.5,0.5,0.5,0), new=T)
# plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
# polygon(japan.poly2,col=res$rrcolors$colors.obs[,2],border=F)
# segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
# text(355,4120,'Period 2 - Obs',cex=1.00)
# 
# par(fig=c(0.4,.6,.6,1), mar=c(.5,0.5,0.5,0), new=T)
# plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
# polygon(japan.poly2,col=res$rrcolors$colors.obs[,3],border=F)
# segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
# text(355,4120,'Period 3 - Obs',cex=1.00)
# 
# par(fig=c(0.6,.8,.6,1), mar=c(.5,0.5,0.5,0), new=T)
# plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
# polygon(japan.poly2,col=res$rrcolors$colors.obs[,4],border=F)
# segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
# text(355,4120,'Period 4 - Obs',cex=1.00)
# 
# par(fig=c(0.8,1,.6,1), mar=c(.5,0.5,0.5,0), new=T)
# plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
# polygon(japan.poly2,col=res$rrcolors$colors.obs[,5],border=F)
# segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
# text(355,4120,'Period 5 - Obs',cex=1.00)
# 
# 
# #Maps of AIC Path
# 
# par(fig=c(0,.2,.4,.8), mar=c(.5,0.5,0.5,0), new=T)
# plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
# polygon(japan.poly2,col=res$rrcolors$color.qaic[,1],border=F)
# segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
# text(355,4120,'Period 1 - QAIC',cex=1.00)
# 
# par(fig=c(0.2,.4,.4,.8), mar=c(.5,0.5,0.5,0), new=T)   
# plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
# polygon(japan.poly2,col=res$rrcolors$color.qaic[,2],border=F)
# segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
# text(355,4120,'Period 2 - QAIC',cex=1.00)
# 
# par(fig=c(0.4,.6,.4,.8), mar=c(.5,0.5,0.5,0), new=T) 
# plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
# polygon(japan.poly2,col=res$rrcolors$color.qaic[,3],border=F)
# segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
# text(355,4120,'Period 3 - QAIC',cex=1.00)
# 
# par(fig=c(0.6,.8,.4,.8), mar=c(.5,0.5,0.5,0), new=T)
# plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
# polygon(japan.poly2,col=res$rrcolors$color.qaic[,4],border=F)
# segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
# text(355,4120,'Period 4 - QAIC',cex=1.00)
# 
# par(fig=c(0.8,1,.4,.8), mar=c(.5,0.5,0.5,0), new=T)
# plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
# polygon(japan.poly2,col=res$rrcolors$color.qaic[,5],border=F)
# segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
# text(355,4120,'Period 5 - QAIC',cex=1.00)
# 
# 
# #Maps of AICc Path
# 
# par(fig=c(0,.2,.2,.6), mar=c(.5,0.5,0.5,0), new=T)
# plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
# polygon(japan.poly2,col=res$rrcolors$color.qaicc[,1],border=F)
# segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
# text(355,4120,'Period 1 - QAICc',cex=1.00)
# 
# par(fig=c(0.2,.4,.2,.6), mar=c(.5,0.5,0.5,0), new=T) 
# plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
# polygon(japan.poly2,col=res$rrcolors$color.qaicc[,2],border=F)
# segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
# text(355,4120,'Period 2 - QAICc',cex=1.00)
# 
# par(fig=c(0.4,.6,.2,.6), mar=c(.5,0.5,0.5,0), new=T) 
# plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
# polygon(japan.poly2,col=res$rrcolors$color.qaicc[,3],border=F)
# segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
# text(355,4120,'Period 3 - QAICc',cex=1.00)
# 
# par(fig=c(0.6,.8,.2,.6), mar=c(.5,0.5,0.5,0), new=T)
# plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
# polygon(japan.poly2,col=res$rrcolors$color.qaicc[,4],border=F)
# segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
# text(355,4120,'Period 4 - QAICc',cex=1.00)
# 
# par(fig=c(0.8,1,.2,.6), mar=c(.5,0.5,0.5,0), new=T)
# plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
# polygon(japan.poly2,col=res$rrcolors$color.qaicc[,5],border=F)
# segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
# text(355,4120,'Period 5 - QAICc',cex=1.00)
# 
# 
# #Maps of BIC Path
# 
# par(fig=c(0,.2,0,.4), mar=c(.5,0.5,0.5,0), new=T)
# plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
# polygon(japan.poly2,col=res$rrcolors$color.qbic[,1],border=F)
# segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
# text(355,4120,'Period 1 - QBIC',cex=1.00)
# 
# par(fig=c(0.2,.4,0,.4), mar=c(.5,0.5,0.5,0), new=T)
# plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
# polygon(japan.poly2,col=res$rrcolors$color.qbic[,2],border=F)
# segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
# text(355,4120,'Period 2 - QBIC',cex=1.00)
# 
# par(fig=c(0.4,.6,0,.4), mar=c(.5,0.5,0.5,0), new=T)
# plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
# polygon(japan.poly2,col=res$rrcolors$color.qbic[,3],border=F)
# segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
# text(355,4120,'Period 3 - QBIC',cex=1.00)
# 
# 
# par(fig=c(0.6,.8,0,.4), mar=c(.5,0.5,0.5,0), new=T)
# plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
# polygon(japan.poly2,col=res$rrcolors$color.qbic[,4],border=F)
# segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
# text(355,4120,'Period 4 - QBIC',cex=1.00)
# 
# 
# par(fig=c(0.8,1,0,.4), mar=c(.5,0.5,0.5,0), new=T)
# plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
# polygon(japan.poly2,col=res$rrcolors$color.qbic[,5],border=F)
# segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
# text(355,4120,'Period 5 - QBIC',cex=1.00)
# 
# #Turn off pdf development
# dev.off()



########################################################################################################
########################################################################################################
########################################################################################################

####################################################
#rmax 30, time = 2,5, rr.ratio = 1.25, center = 75, r = 11km
####################################################
#Set Some Initial Conditions
x=dframe2$utmx/1000
y=dframe2$utmy/1000
rMax=30 
Time=5
nsim=100
center=150
radius=11
period_start=c(2,5)
risk.ratio=1.25




res <- clust.sim(x,y,rMax,dframe$period, dframe$expdeath, dframe$death, Time, spacetime=TRUE, pois=FALSE,
                 nsim,center, radius, risk.ratio, period_start, colors=TRUE, utm=TRUE, byrow=TRUE)

#save results
filename <- paste0("SimulationOutput/sim","_","center","_",center,"radius","_",radius, "start","_","rr","_","1_25",".RData")
save(res, file = filename)

###did it find at least 1 cell in the cluster?
set <- detect.set(res$lassoresult, res$init.vec, res$rr.mat, Time=5, nb, x, y, rMax, center, radius)
(incluster <- detect.incluster(res$lassoresult, res$init.vec, res$rr.mat, set, period_start = 2, 
                               period_end = 5, multi_period = FALSE,Time=5, nsim,nb, x, y, rMax, center, radius, IC = "ic"))

#make maps
pdfname <- paste0("figures/simulations/sim","_","center","_",center,"radius","_",radius, "start","_","rr","_","1_25",".pdf")
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
text(355,4120,'Period 1 - QBIC',cex=1.00)

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
#rmax 30, time = 3,4,5, rr.ratio = 0.75 center = 30, r = 18km
####################################################
#Set Some Initial Conditions
x=dframe2$utmx/1000
y=dframe2$utmy/1000
rMax=30 
Time=5
nsim=100
center=30
radius=18
period_start=c(3:5)
risk.ratio=0.75




res <- clust.sim(x,y,rMax,dframe$period, dframe$expdeath, dframe$death, Time, spacetime=TRUE, pois=TRUE,
                 nsim,center, radius, risk.ratio, period_start, colors=TRUE, utm=TRUE, byrow=TRUE)

#save results
filename <- paste0("SimulationOutput/sim","_","center","_",center,"radius","_",radius, "start","_","rr","_","075",".RData")
save(res, file = filename)

###did it find at least 1 cell in the cluster?
set <- detect.set(res$lassoresult, res$init.vec, res$rr.mat, Time=5, nb, x, y, rMax, center, radius)
(incluster <- detect.incluster(res$lassoresult, res$init.vec, res$rr.mat, set, period_start = 3, 
                               period_end = 5, multi_period = TRUE,Time=5, nsim,nb, x, y, rMax, center, radius, IC = "ic", under=TRUE))

#make maps
pdfname <- paste0("figures/simulations/sim","_","center","_",center,"radius","_",radius, "start","_","rr","_","075",".pdf")
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
text(355,4120,'Period 1 - QBIC',cex=1.00)

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
#rmax 30, time = 5, rr.ratio = 0.75, center = 50, r = 18km
####################################################
#Set Some Initial Conditions
x=dframe2$utmx/1000
y=dframe2$utmy/1000
rMax=30 
Time=5
nsim=100
center=175
radius=18
period_start=c(4,5)
risk.ratio=2




res <- clust.sim(x,y,rMax,dframe$period, dframe$expdeath, dframe$death, Time, spacetime=TRUE,
                 nsim,center, radius, risk.ratio, period_start, colors=TRUE, utm=TRUE, byrow=TRUE)

#save results
filename <- paste0("SimulationOutput/sim","_","center","_",center,"radius","_",radius, "start","_","rr","_","2",".RData")
save(res, file = filename)

###did it find at least 1 cell in the cluster?
set <- detect.set(res$lassoresult, res$init.vec, res$rr.mat, Time=5, nb, x, y, rMax, center, radius)
(incluster <- detect.incluster(res$lassoresult, res$init.vec, res$rr.mat, set, period_start = 4, 
                               period_end = 5, multi_period = TRUE,Time=5, nsim,nb, x, y, rMax, center, radius, IC = "ic"))

#make maps
pdfname <- paste0("figures/simulations/sim","_","center","_",center,"radius","_",radius, "start","_","rr","_","2",".pdf")
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
text(355,4120,'Period 1 - QBIC',cex=1.00)

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
#rmax 30, time = 5, rr.ratio = 0.75, center = 50, r = 18km
####################################################
#Set Some Initial Conditions
x=dframe2$utmx/1000
y=dframe2$utmy/1000
rMax=30 
Time=5
nsim=100
center=160
radius= 5
period_start=c(1,4)
risk.ratio=2




res <- clust.sim(x,y,rMax,dframe$period, dframe$expdeath, dframe$death, Time, spacetime=TRUE, pois=TRUE,
                 nsim,center, radius, risk.ratio, period_start, colors=TRUE, utm=TRUE, byrow=TRUE)

#save results
filename <- paste0("SimulationOutput/sim","_","center","_",center,"radius","_",radius, "start","_","rr","_","2",".RData")
save(res, file = filename)

###did it find at least 1 cell in the cluster?
set <- detect.set(res$lassoresult, res$init.vec, res$rr.mat, Time=5, nb, x, y, rMax, center, radius)
(incluster <- detect.incluster(res$lassoresult, res$init.vec, res$rr.mat, set, period_start = 1, 
                               period_end = 4, multi_period = FALSE,Time=5, nsim,nb, x, y, rMax, center, radius, IC = "ic"))

#make maps
pdfname <- paste0("figures/simulations/sim","_","center","_",center,"radius","_",radius, "start","_","rr","_","2",".pdf")
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
text(355,4120,'Period 1 - QBIC',cex=1.00)

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
#rmax 30, time = 5, rr.ratio = 0.75, center = 50, r = 18km - POISSON ONLY
####################################################
#Set Some Initial Conditions
x=dframe2$utmx/1000
y=dframe2$utmy/1000
rMax=30 
Time=5
nsim=100
center=160
radius= 9
period_start=c(2:4)
risk.ratio=1.25




res <- clust.sim(x,y,rMax,dframe$period, dframe$expdeath, dframe$death, Time, spacetime=TRUE, pois=TRUE,
                 nsim,center, radius, risk.ratio, period_start, colors=TRUE, utm=TRUE, byrow=TRUE)

#save results
filename <- paste0("SimulationOutput/sim","_","center","_",center,"radius","_",radius, "start","_","rr","_","125_poisson",".RData")
save(res, file = filename)

###did it find at least 1 cell in the cluster?
set <- detect.set(res$lassoresult, res$init.vec, res$rr.mat, Time=5, nb, x, y, rMax, center, radius)
(incluster <- detect.incluster(res$lassoresult, res$init.vec, res$rr.mat, set, period_start = 2, 
                               period_end = 4, multi_period = TRUE,Time=5, nsim,nb, x, y, rMax, center, radius, IC = "ic"))

#make maps
pdfname <- paste0("figures/simulations/sim","_","center","_",center,"radius","_",radius, "start","_","rr","_","125_poisson",".pdf")
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
text(355,4120,'Period 1 - QBIC',cex=1.00)

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
#rmax 30, time = 5, rr.ratio = 2, center = 35, 160 , r = 18km - 2 clusters
####################################################
#Set Some Initial Conditions
x=dframe2$utmx/1000
y=dframe2$utmy/1000
rMax=30 
Time=5
nsim=100
center=c(35, 160)
radius= 9
period_start=c(2:4)
risk.ratio=2




res <- clust.sim(x,y,rMax,dframe$period, dframe$expdeath, dframe$death, Time, spacetime=TRUE, pois=FALSE,
                 nsim,center, radius, risk.ratio, period_start, colors=TRUE, utm=TRUE, byrow=TRUE)

#save results
filename <- paste0("SimulationOutput/sim","_","center","_",center,"radius","_",radius, "start","_","rr","_","15","twocenters",".RData")
save(res, file = filename)

###did it find at least 1 cell in the cluster?
set <- detect.set(res$lassoresult, res$init.vec, res$rr.mat, Time=5, nb, x, y, rMax, center, radius)
(incluster <- detect.incluster(res$lassoresult, res$init.vec, res$rr.mat, set, period_start = 2, 
                               period_end = 4, multi_period = TRUE,Time=5, nsim,nb, x, y, rMax, center, radius, IC = "ic"))

#make maps
pdfname <- paste0("figures/simulations/sim","_","center","_",center,"radius","_",radius, "start","_","rr","_","15","2centers",".pdf")
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
text(355,4120,'Period 1 - QBIC',cex=1.00)

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


