########################################################################################################
########################################################################################################
########################################################################################################
#Space-ONLY Analysis of Japanese Breast Cancer Data using cluST.R
#Maria Kamenetsky
#9-26-16
########################################################################################################
########################################################################################################
########################################################################################################

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
sourceCpp("../../scripts/cluST/src/maxcol.cpp")
sourceCpp("../../scripts/cluST/src/st_matCppTest.cpp")
sourceCpp("../../scripts/cluST/src/prod_yxTest.cpp")


#temporarily source my cluST.R file
source("../../scripts/cluST//R//cluST.R")



####################################################
#LOAD JBC Data and Set Up
####################################################

#Load Data
dframe1 <- read.csv("../../data/JBC/jap.breast.F.9.10.11.csv")
dframe2 <- read.csv("../../data//JBC//utmJapan.csv")
dframe3 <- aggregate(dframe1, by=list(as.factor(rep(1:(nrow(dframe1)/4),each=4))), FUN="sum")
dframe=data.frame(id=as.factor(dframe3$id/4),period=as.factor(dframe3$year),death=dframe3$death,expdeath=dframe3$expdeath)
levels(dframe$period) <- c("1","2","3","4","5")

#Limit to 1 Time Period for Space Only Model
dframe <- dframe[dframe$period=="1",]

#Set Some Initial Conditions
x1=dframe2$utmx/1000
y1=dframe2$utmy/1000
rMax <- 30 
Time=1

#Create Potential Clusters Dataframe
clusters <- clustersDF(x1,y1,rMax, utm=TRUE, length(x1))

#Set initial expected and observed
JBCinit <- setVectors(dframe$period, dframe$expdeath, dframe$death, Time=1, byrow=TRUE)

#Adjust for observed given expected counts as coming from negative binomial distribution
outinit <- glm.nb(JBCinit$Y.vec ~1)
out <- glm.nb(JBCinit$Y.vec ~ 1 + offset(log(JBCinit$E0)), init.theta = outinit$theta, 
              link=log,control=glm.control(maxit=10000))

#Set initial expected to the fitted values
E0 <- out$fitted


####################################################
#RUN Model
####################################################

#set initial conditions for function (all centers can be potential origin of the cluster)
potentialClus <- max(clusters$center)
numberCenters <- max(clusters$center)

JBCresults <- spacetimeLasso(potentialClus, clusters, numberCenters, JBCinit, Time, spacetime=FALSE)

####################################################
#Set Risk Ratio Vectors Based on QIC
####################################################
rr <- setRR(JBCresults, JBCinit, Time)

####################################################
#Map RR to Colors
####################################################
rrcolors <- colormapping(rr, Time)


####################################################
#Map Colors to Maps
####################################################

#Import Datasets with Map Polygons
dframe.poly2 <- read.csv("../../data/JBC/japan_poly2.csv")
japan.poly2 <- dframe.poly2[,2:3]
dframe.prefect2 <- read.csv("../../data/JBC/japan_prefect2.csv")
japan.prefect2 <- dframe.prefect2[,2:5]

#Create Empty PDF to Map Onto
pdf("../../figures/JBC/japan_map_spaceonly.pdf", height=11, width=10)

#Maps of Observed Counts
par(fig=c(0,.2,.6,1), mar=c(.5,0.5,0.5,0))
plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
polygon(japan.poly2,col=rrcolors$colors.obs[,1],border=F)
segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
text(355,4120,'Period 1 - Obs',cex=1.00)

#Maps of AIC Path

par(fig=c(0,.2,.4,.8), mar=c(.5,0.5,0.5,0), new=T)
plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
polygon(japan.poly2,col=rrcolors$color.qaic[,1],border=F)
segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
text(355,4120,'Period 1 - QAIC',cex=1.00)

#Maps of AICc Path

par(fig=c(0,.2,.2,.6), mar=c(.5,0.5,0.5,0), new=T)
plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
polygon(japan.poly2,col=rrcolors$color.qaicc[,1],border=F)
segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
text(355,4120,'Period 1 - QAICc',cex=1.00)


#Maps of BIC Path

par(fig=c(0,.2,0,.4), mar=c(.5,0.5,0.5,0), new=T)
plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
polygon(japan.poly2,col=rrcolors$color.qbic[,1],border=F)
segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
text(355,4120,'Period 1 - QBI',cex=1.00)


#Turn off pdf development
dev.off()






