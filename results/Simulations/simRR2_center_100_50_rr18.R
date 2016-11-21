########################################################################################################
########################################################################################################
########################################################################################################
#Simulations of Space-time Analysis Based on JBC Data using cluST.R
#Maria Kamenetsky
#10-3-16
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
sourceCpp("scripts/cluST/src/maxcol.cpp")
sourceCpp("scripts/cluST/src/st_matCpp.cpp")
sourceCpp("scripts/cluST/src/prod_yx.cpp")


#temporarily source my cluST.R file
source("scripts/cluST//R//cluST.R")



####################################################
#Create Some Fake Data
####################################################
set.seed(1032016)

#Load Data
dframe1 <- read.csv("data/JBC/jap.breast.F.9.10.11.csv")
dframe2 <- read.csv("data//JBC//utmJapan.csv")
dframe3 <- aggregate(dframe1, by=list(as.factor(rep(1:(nrow(dframe1)/4),each=4))), FUN="sum")
dframe=data.frame(id=as.factor(dframe3$id/4),period=as.factor(dframe3$year),death=dframe3$death,expdeath=dframe3$expdeath)
levels(dframe$period) <- c("1","2","3","4","5")

#Set Some Initial Conditions
x1=dframe2$utmx/1000
y1=dframe2$utmy/1000
rMax <- 30 
Time=5  

#Set number of simulations
nsim=100

#Number of time periods and centers
n <- 208
Time <- 5

#Set theta parameter
thetainit = 1000

#Set parameters of your fake cluster
center <- c(100,50)
r_list <- 18
cluster_end <- 3
rr.ratio<- 2

#Create Potential Clusters Dataframe
clusters <- clustersDF(x1,y1,rMax, utm=TRUE, length(x1))


#Set initial expected and observed
JBCinit <- setVectors(dframe$period, dframe$expdeath, dframe$death, Time=5, byrow=TRUE)

#Adjust for observed given expected counts as coming from negative binomial distribution
out <- glm.nb(JBCinit$Y.vec ~ 1 + as.factor(JBCinit$Year)  + offset(log(JBCinit$E0)), init.theta = thetainit, 
              link=log,control=glm.control(maxit=10000))
E0_fit <- out$fitted.values

#Set initial expected to the fitted values
E0_0 <- JBCinit$E0

####################################################
#Create Fake Clusters
####################################################
tmp <- clusters[clusters$center==center[1] | clusters$center==center[2],]
cluster <- tmp[(tmp$r <= r_list),]
rr = matrix(1, nrow=n, ncol=Time)
rr[cluster$last, cluster_end:Time] = rr.ratio
expect_fake <- as.vector(rr)*E0_0


#Set Vectors
Y.vec <- JBCinit$Y.vec
Period <- JBCinit$Year


####################################################
#Sim Response
####################################################
YSIM <- simulate(out, nsim=nsim)

####################################################
#Refit E0 for each simulation and scale
####################################################
#Adjust for observed given expected counts as coming from negative binomial distribution
out.sim <- lapply(1:nsim, function(i) glm.nb(YSIM[,i] ~ 1 + as.factor(JBCinit$Year)  + offset(log(expect_fake)), init.theta = thetainit, 
                                             link=log,control=glm.control(maxit=10000)))

#Set initial expected to the fitted values and standardize
E0 <- scale(Y.vec, out.sim, nsim, Time)

JBCinit.sim <- list(Period = Period, E0=E0, E0_fit=E0_fit, Y.vec=Y.vec)

####################################################
#Set up and Run Model
####################################################
potentialClusters <- max(clusters$center)
numCenters <- max(clusters$center)

JBCresults.sim <- spacetimeLasso.sim(potentialClusters, clusters, numCenters,
                           JBCinit.sim, Time, spacetime=TRUE, nsim, YSIM)

save(JBCresults.sim, file="SimulationOutput//simR2_center100_50_r18.RData")
####################################################
#Risk Ratios
####################################################
##Calculate average observed for simulated
##RR calculations
riskratios <- getRR(JBCresults.sim, JBCinit.sim, Time, sim=TRUE)
rrcolors <- colormapping(riskratios,Time)

##Probability map calculations
probcolors <- probmap(JBCresults.sim, JBCinit.sim, rr, nsim,Time, colormap=TRUE)



####################################################
#Make Maps
####################################################

#Import Datasets with Map Polygons
dframe.poly2 <- read.csv("data/JBC/japan_poly2.csv")
japan.poly2 <- dframe.poly2[,2:3]
dframe.prefect2 <- read.csv("data/JBC/japan_prefect2.csv")
japan.prefect2 <- dframe.prefect2[,2:5]

#Create Empty PDF to Map Onto
pdf("figures/simulations/japan_map_R2_2clus.pdf", height=11, width=10)

#Maps of Observed Counts
par(fig=c(0,.2,.6,1), mar=c(.5,0.5,0.5,0))
plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
polygon(japan.poly2,col=rrcolors$colors.obs[,1],border=F)
segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
text(355,4120,'Period 1 - Obs',cex=1.00)

par(fig=c(0.2,.4,.6,1), mar=c(.5,0.5,0.5,0), new=T)
plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
polygon(japan.poly2,col=rrcolors$colors.obs[,2],border=F)
segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
text(355,4120,'Period 2 - Obs',cex=1.00)

par(fig=c(0.4,.6,.6,1), mar=c(.5,0.5,0.5,0), new=T)
plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
polygon(japan.poly2,col=rrcolors$colors.obs[,3],border=F)
segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
text(355,4120,'Period 3 - Obs',cex=1.00)

par(fig=c(0.6,.8,.6,1), mar=c(.5,0.5,0.5,0), new=T)
plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
polygon(japan.poly2,col=rrcolors$colors.obs[,4],border=F)
segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
text(355,4120,'Period 4 - Obs',cex=1.00)

par(fig=c(0.8,1,.6,1), mar=c(.5,0.5,0.5,0), new=T)
plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
polygon(japan.poly2,col=rrcolors$colors.obs[,5],border=F)
segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
text(355,4120,'Period 5 - Obs',cex=1.00)


#Maps of AIC Path

par(fig=c(0,.2,.4,.8), mar=c(.5,0.5,0.5,0), new=T)
plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
polygon(japan.poly2,col=rrcolors$color.qaic[,1],border=F)
segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
text(355,4120,'Period 1 - QAIC',cex=1.00)

par(fig=c(0.2,.4,.4,.8), mar=c(.5,0.5,0.5,0), new=T)   
plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
polygon(japan.poly2,col=rrcolors$color.qaic[,2],border=F)
segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
text(355,4120,'Period 2 - QAIC',cex=1.00)

par(fig=c(0.4,.6,.4,.8), mar=c(.5,0.5,0.5,0), new=T) 
plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
polygon(japan.poly2,col=rrcolors$color.qaic[,3],border=F)
segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
text(355,4120,'Period 3 - QAIC',cex=1.00)

par(fig=c(0.6,.8,.4,.8), mar=c(.5,0.5,0.5,0), new=T)
plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
polygon(japan.poly2,col=rrcolors$color.qaic[,4],border=F)
segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
text(355,4120,'Period 4 - QAIC',cex=1.00)

par(fig=c(0.8,1,.4,.8), mar=c(.5,0.5,0.5,0), new=T)
plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
polygon(japan.poly2,col=rrcolors$color.qaic[,5],border=F)
segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
text(355,4120,'Period 5 - QAIC',cex=1.00)


#Maps of AICc Path

par(fig=c(0,.2,.2,.6), mar=c(.5,0.5,0.5,0), new=T)
plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
polygon(japan.poly2,col=rrcolors$color.qaicc[,1],border=F)
segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
text(355,4120,'Period 1 - QAICc',cex=1.00)

par(fig=c(0.2,.4,.2,.6), mar=c(.5,0.5,0.5,0), new=T) 
plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
polygon(japan.poly2,col=rrcolors$color.qaicc[,2],border=F)
segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
text(355,4120,'Period 2 - QAICc',cex=1.00)

par(fig=c(0.4,.6,.2,.6), mar=c(.5,0.5,0.5,0), new=T) 
plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
polygon(japan.poly2,col=rrcolors$color.qaicc[,3],border=F)
segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
text(355,4120,'Period 3 - QAICc',cex=1.00)

par(fig=c(0.6,.8,.2,.6), mar=c(.5,0.5,0.5,0), new=T)
plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
polygon(japan.poly2,col=rrcolors$color.qaicc[,4],border=F)
segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
text(355,4120,'Period 4 - QAICc',cex=1.00)

par(fig=c(0.8,1,.2,.6), mar=c(.5,0.5,0.5,0), new=T)
plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
polygon(japan.poly2,col=rrcolors$color.qaicc[,5],border=F)
segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
text(355,4120,'Period 5 - QAICc',cex=1.00)


#Maps of BIC Path

par(fig=c(0,.2,0,.4), mar=c(.5,0.5,0.5,0), new=T)
plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
polygon(japan.poly2,col=rrcolors$color.qbic[,1],border=F)
segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
text(355,4120,'Period 1 - QBI',cex=1.00)

par(fig=c(0.2,.4,0,.4), mar=c(.5,0.5,0.5,0), new=T)
plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
polygon(japan.poly2,col=rrcolors$color.qbic[,2],border=F)
segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
text(355,4120,'Period 2 - QBIC',cex=1.00)

par(fig=c(0.4,.6,0,.4), mar=c(.5,0.5,0.5,0), new=T)
plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
polygon(japan.poly2,col=rrcolors$color.qbic[,3],border=F)
segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
text(355,4120,'Period 3 - QBIC',cex=1.00)


par(fig=c(0.6,.8,0,.4), mar=c(.5,0.5,0.5,0), new=T)
plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
polygon(japan.poly2,col=rrcolors$color.qbic[,4],border=F)
segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
text(355,4120,'Period 4 - QBIC',cex=1.00)


par(fig=c(0.8,1,0,.4), mar=c(.5,0.5,0.5,0), new=T)
plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
polygon(japan.poly2,col=rrcolors$color.qbic[,5],border=F)
segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
text(355,4120,'Period 5 - QBIC',cex=1.00)

#Turn off pdf development
dev.off()



####################################################
#Make Maps - Probability Map
####################################################
#Create Empty PDF to Map Onto
pdf("figures/simulations/japan_probmap_R2_2clus100_50_rr18.pdf", height=11, width=10)
#Maps of Observed Counts
par(fig=c(0,.2,.6,1), mar=c(.5,0.5,0.5,0))
plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
polygon(japan.poly2,col=probcolors$color.probmap[,1],border=F)
segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
text(355,4120,'Period 1 - Obs',cex=1.00)

par(fig=c(0.2,.4,.6,1), mar=c(.5,0.5,0.5,0), new=T)
plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
polygon(japan.poly2,col=probcolors$color.probmap[,2],border=F)
segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
text(355,4120,'Period 2 - Obs',cex=1.00)

par(fig=c(0.4,.6,.6,1), mar=c(.5,0.5,0.5,0), new=T)
plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
polygon(japan.poly2,col=probcolors$color.probmap[,3],border=F)
segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
text(355,4120,'Period 3 - Obs',cex=1.00)

par(fig=c(0.6,.8,.6,1), mar=c(.5,0.5,0.5,0), new=T)
plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
polygon(japan.poly2,col=probcolors$color.probmap[,4],border=F)
segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
text(355,4120,'Period 4 - Obs',cex=1.00)

par(fig=c(0.8,1,.6,1), mar=c(.5,0.5,0.5,0), new=T)
plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
polygon(japan.poly2,col=probcolors$color.probmap[,5],border=F)
segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
text(355,4120,'Period 5 - Obs',cex=1.00)

#Probability Maps AIC

par(fig=c(0,.2,.4,.8), mar=c(.5,0.5,0.5,0), new=T)
plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
polygon(japan.poly2,col=probcolors$color.probmapAIC[,1],border=F)
segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
text(355,4120,'Period 1 - QAIC',cex=1.00)

par(fig=c(0.2,.4,.4,.8), mar=c(.5,0.5,0.5,0), new=T)   
plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
polygon(japan.poly2,col=probcolors$color.probmapAIC[,2],border=F)
segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
text(355,4120,'Period 2 - QAIC',cex=1.00)

par(fig=c(0.4,.6,.4,.8), mar=c(.5,0.5,0.5,0), new=T) 
plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
polygon(japan.poly2,col=probcolors$color.probmapAIC[,3],border=F)
segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
text(355,4120,'Period 3 - QAIC',cex=1.00)

par(fig=c(0.6,.8,.4,.8), mar=c(.5,0.5,0.5,0), new=T)
plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
polygon(japan.poly2,col=probcolors$color.probmapAIC[,4],border=F)
segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
text(355,4120,'Period 4 - QAIC',cex=1.00)

par(fig=c(0.8,1,.4,.8), mar=c(.5,0.5,0.5,0), new=T)
plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
polygon(japan.poly2,col=probcolors$color.probmapAIC[,5],border=F)
segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
text(355,4120,'Period 5 - QAIC',cex=1.00)

#Probability Map AICc


par(fig=c(0,.2,.2,.6), mar=c(.5,0.5,0.5,0), new=T)
plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
polygon(japan.poly2,col=probcolors$color.probmapAICc[,1],border=F)
segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
text(355,4120,'Period 1 - QAICc',cex=1.00)

par(fig=c(0.2,.4,.2,.6), mar=c(.5,0.5,0.5,0), new=T) 
plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
polygon(japan.poly2,col=probcolors$color.probmapAICc[,2],border=F)
segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
text(355,4120,'Period 2 - QAICc',cex=1.00)

par(fig=c(0.4,.6,.2,.6), mar=c(.5,0.5,0.5,0), new=T) 
plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
polygon(japan.poly2,col=probcolors$color.probmapAICc[,3],border=F)
segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
text(355,4120,'Period 3 - QAICc',cex=1.00)

par(fig=c(0.6,.8,.2,.6), mar=c(.5,0.5,0.5,0), new=T)
plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
polygon(japan.poly2,col=probcolors$color.probmapAICc[,4],border=F)
segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
text(355,4120,'Period 4 - QAICc',cex=1.00)

par(fig=c(0.8,1,.2,.6), mar=c(.5,0.5,0.5,0), new=T)
plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
polygon(japan.poly2,col=probcolors$color.probmapAICc[,5],border=F)
segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
text(355,4120,'Period 5 - QAICc',cex=1.00)

#Probability Map of BIC

par(fig=c(0,.2,0,.4), mar=c(.5,0.5,0.5,0), new=T)
plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
polygon(japan.poly2,col=probcolors$color.probmapBIC[,1],border=F)
segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
text(355,4120,'Period 1 - QBI',cex=1.00)

par(fig=c(0.2,.4,0,.4), mar=c(.5,0.5,0.5,0), new=T)
plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
polygon(japan.poly2,col=probcolors$color.probmapBIC[,2],border=F)
segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
text(355,4120,'Period 2 - QBIC',cex=1.00)

par(fig=c(0.4,.6,0,.4), mar=c(.5,0.5,0.5,0), new=T)
plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
polygon(japan.poly2,col=probcolors$color.probmapBIC[,3],border=F)
segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
text(355,4120,'Period 3 - QBIC',cex=1.00)


par(fig=c(0.6,.8,0,.4), mar=c(.5,0.5,0.5,0), new=T)
plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
polygon(japan.poly2,col=probcolors$color.probmapBIC[,4],border=F)
segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
text(355,4120,'Period 4 - QBIC',cex=1.00)


par(fig=c(0.8,1,0,.4), mar=c(.5,0.5,0.5,0), new=T)
plot(japan.poly2,type='n',asp=1,axes=F,xlab='',ylab='')
polygon(japan.poly2,col=probcolors$color.probmapBIC[,5],border=F)
segments(japan.prefect2$x1,japan.prefect2$y1,japan.prefect2$x2,japan.prefect2$y2)
text(355,4120,'Period 5 - QBIC',cex=1.00)


#Turn off pdf development
dev.off()