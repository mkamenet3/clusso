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
sink("simulations.txt")
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
n <- 208
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
#nsim=2
center=200
radius=18
timeperiod = c(1:5)
risk.ratio=1

res <- clust.sim(x,y,rMax,dframe$period, dframe$expdeath, dframe$death, Time, spacetime=TRUE, pois=FALSE,
                 nsim,center, radius, risk.ratio, timeperiod, colors=TRUE, utm=TRUE, byrow=TRUE)

#save results
filename <- paste0("SimulationOutput/sim","_","center","_",center,"radius",radius,"_", "start",
                   "_",as.numeric(paste(timeperiod, collapse = "")),"_","rr","_",gsub("[.]","",risk.ratio),"NULL",".RData")
save(res, file = filename)

#Detection
set <- detect.set(res$lassoresult, res$init.vec, res$rr.mat, Time=5, nb, x, y, rMax, 0, radius)
(incluster <- detect.incluster(res$lassoresult, res$init.vec, res$rr.mat, set, period_start = timeperiod[1], 
                               period_end = tail(timeperiod, n=1), multi_period = TRUE,Time=5, nsim,nb, x, y, rMax, 0, radius, IC = "ic"))


#make maps
pdfname <- paste0("figures/simulations/sim","_","center","_",center,"radius",radius,"_", "start", "_",as.numeric(paste(timeperiod, collapse = "")),
                  "_","rr","_",gsub("[.]","",risk.ratio),"NULL",".pdf")
plotmap.st(pdfname, res)


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
timeperiod = c(1:5)
risk.ratio=1

res <- clust.sim(x,y,rMax,dframe$period, dframe$expdeath, dframe$death, Time, spacetime=TRUE, pois=TRUE,
                 nsim,center, radius, risk.ratio, timeperiod, colors=TRUE, utm=TRUE, byrow=TRUE)

#save results
filename <- paste0("SimulationOutput/sim","_","center","_",center,"radius",radius,"_", "start",
                   "_",as.numeric(paste(timeperiod, collapse = "")),"_","rr","_",gsub("[.]","",risk.ratio),"NULLPOISSON",".RData")
save(res, file = filename)

#Detection
set <- detect.set(res$lassoresult, res$init.vec, res$rr.mat, Time=5, nb, x, y, rMax, 0, radius)
(incluster <- detect.incluster(res$lassoresult, res$init.vec, res$rr.mat, set, period_start = timeperiod[1], 
                               period_end = tail(timeperiod, n=1), multi_period = TRUE,Time=5, 
                               nsim,nb, x, y, rMax, 0, radius, IC = "ic", nullmod=TRUE))
#make maps
pdfname <- paste0("figures/simulations/sim","_","center","_",center,"radius",radius,"_", "start", "_",as.numeric(paste(timeperiod, collapse = "")),
                  "_","rr","_",gsub("[.]","",risk.ratio),"NULLPOISSON",".pdf")
plotmap.st(pdfname, res)


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
timeperiod=1
risk.ratio=1


#create average dataframe for spaceonly

death <- with(dframe, tapply(death, id, function(x) round(mean(x))))
expdeath <- with(dframe, tapply(expdeath, id, function(x) mean(x)))
df <- cbind.data.frame(id = unique(dframe$id), period = rep("1", length(unique(dframe$id))), death = death, expdeath=expdeath)

res <- clust.sim(x,y,rMax,df$period, df$expdeath, df$death, Time, spacetime=FALSE, pois=FALSE, 
                 nsim,center, radius, risk.ratio, period_start, colors=TRUE, utm=TRUE, byrow=TRUE)
#save results
filename <- paste0("SimulationOutput/sim","_","center","_",center,"radius",radius,"_", "start",
                   "_",as.numeric(paste(timeperiod, collapse = "")),"_","rr","_",gsub("[.]","",risk.ratio),"spaceonlyQPOIS",".RData")
save(res, file = filename)

#Detection
set <- detect.set(res$lassoresult, res$init.vec, res$rr.mat, Time=1, nb, x, y, rMax, 0, radius)
(incluster <- detect.incluster(res$lassoresult, res$init.vec, res$rr.mat, set, period_start = 1, 
                               period_end = 1, multi_period = FALSE,Time=1, nsim,nb, x, y, rMax, 0, radius, IC = "ic", space=FALSE, nullmod=TRUE))

#make maps
pdfname <- paste0("figures/simulations/sim","_","center","_",center,"radius",radius,"_", "start", "_",as.numeric(paste(timeperiod, collapse = "")),
                  "_","rr","_",gsub("[.]","",risk.ratio),"spaceonlyQPOIS",".pdf")
plotmap.s(pdfname, res)

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
timeperiod=1
risk.ratio=1


#create average dataframe for spaceonly

death <- with(dframe, tapply(death, id, function(x) round(mean(x))))
expdeath <- with(dframe, tapply(expdeath, id, function(x) mean(x)))
df <- cbind.data.frame(id = unique(dframe$id), period = rep("1", length(unique(dframe$id))), death = death, expdeath=expdeath)

res <- clust.sim(x,y,rMax,df$period, df$expdeath, df$death, Time, spacetime=FALSE, pois=TRUE, 
                 nsim,center, radius, risk.ratio, period_start, colors=TRUE, utm=TRUE, byrow=TRUE)
#save results
filename <- paste0("SimulationOutput/sim","_","center","_",center,"radius",radius,"_", "start",
                   "_",as.numeric(paste(timeperiod, collapse = "")),"_","rr","_",gsub("[.]","",risk.ratio),"spaceonlyPOISSONONLY",".RData")
save(res, file = filename)

#Detection
set <- detect.set(res$lassoresult, res$init.vec, res$rr.mat, Time=1, nb, x, y, rMax, 0, radius)
(incluster <- detect.incluster(res$lassoresult, res$init.vec, res$rr.mat, set, period_start = 1, 
                               period_end = 1, multi_period = FALSE,Time=1, nsim,nb, x, y, rMax, 
                               0, radius, IC = "ic", space=FALSE, nullmod=TRUE))

#make maps
pdfname <- paste0("figures/simulations/sim","_","center","_",center,"radius",radius,"_", "start", "_",as.numeric(paste(timeperiod, collapse = "")),
                   "_","rr","_",gsub("[.]","",risk.ratio),"spaceonlyPOISSONONLY",".pdf")
plotmap.s(pdfname, res)


########################################################################################################
########################################################################################################
########################################################################################################
#END NULL MODELS, BEGIN CLUSTER SIMULATIONS
########################################################################################################
########################################################################################################
########################################################################################################
#Center = 50

####################################################

####################################################
#Set Some Initial Conditions
x=dframe2$utmx/1000
y=dframe2$utmy/1000
rMax=30 
Time=5  
nsim=100
center=50
radius= 18
timeperiod=c(1:5)
risk.ratio=2

res <- clust.sim(x,y,rMax,dframe$period, dframe$expdeath, dframe$death, Time, spacetime=TRUE, pois=FALSE,
                  nsim,center, radius, risk.ratio, timeperiod, colors=TRUE, utm=TRUE, byrow=TRUE)

#save results
filename <- paste0("SimulationOutput/sim","_","center","_",center,"radius",radius,"_", "start",
                   "_",as.numeric(paste(timeperiod, collapse = "")),"_","rr","_",gsub("[.]","",risk.ratio),".RData")
save(res, file = filename)

#Detection
set <- detect.set(res$lassoresult, res$init.vec, res$rr.mat, Time=5, nb, x, y, rMax, center, radius)
(incluster <- detect.incluster(res$lassoresult, res$init.vec, res$rr.mat, set, period_start = timeperiod[1], 
                               period_end = timeperiod[2], multi_period = TRUE,Time=5, nsim,nb, x, y, rMax, center, radius, IC = "ic"))
#make maps
pdfname <- paste0("figures/simulations/sim","_","center","_",center,"radius",radius,"_", "start",
                  "_",as.numeric(paste(timeperiod, collapse = "")),"_","rr","_",gsub("[.]","",risk.ratio),".pdf")
plotmap.st(pdfname, res)


########################################################################################################
########################################################################################################
########################################################################################################


#Set Some Initial Conditions
x=dframe2$utmx/1000
y=dframe2$utmy/1000
rMax=30 
Time=5  
nsim=100
center=50
radius= 18
timeperiod=c(3:5)
risk.ratio=2

res <- clust.sim(x,y,rMax,dframe$period, dframe$expdeath, dframe$death, Time, spacetime=TRUE, pois=FALSE,
                 nsim,center, radius, risk.ratio, timeperiod, colors=TRUE, utm=TRUE, byrow=TRUE)

#save results
filename <- paste0("SimulationOutput/sim","_","center","_",center,"radius",radius,"_", "start",
                   "_",as.numeric(paste(timeperiod, collapse = "")),"_","rr","_",gsub("[.]","",risk.ratio),".RData")
save(res, file = filename)

#Detection
set <- detect.set(res$lassoresult, res$init.vec, res$rr.mat, Time=5, nb, x, y, rMax, center, radius)
(incluster <- detect.incluster(res$lassoresult, res$init.vec, res$rr.mat, set, period_start = timeperiod[1], 
                               period_end = timeperiod[2], multi_period = TRUE,Time=5, nsim,nb, x, y, rMax, center, radius, IC = "ic"))
#make maps
pdfname <- paste0("figures/simulations/sim","_","center","_",center,"radius",radius,"_", "start",
                  "_",as.numeric(paste(timeperiod, collapse = "")),"_","rr","_",gsub("[.]","",risk.ratio),".pdf")
plotmap.st(pdfname, res)

########################################################################################################
########################################################################################################
########################################################################################################

#Set Some Initial Conditions
x=dframe2$utmx/1000
y=dframe2$utmy/1000
rMax=30 
Time=5  
nsim=100
center=50
radius= 18
timeperiod=5
risk.ratio=2

res <- clust.sim(x,y,rMax,dframe$period, dframe$expdeath, dframe$death, Time, spacetime=TRUE, pois=FALSE,
                 nsim,center, radius, risk.ratio, timeperiod, colors=TRUE, utm=TRUE, byrow=TRUE)

#save results
filename <- paste0("SimulationOutput/sim","_","center","_",center,"radius",radius,"_", "start",
                   "_",as.numeric(paste(timeperiod, collapse = "")),"_","rr","_",gsub("[.]","",risk.ratio),".RData")
save(res, file = filename)

#Detection
set <- detect.set(res$lassoresult, res$init.vec, res$rr.mat, Time=5, nb, x, y, rMax, center, radius)
(incluster <- detect.incluster(res$lassoresult, res$init.vec, res$rr.mat, set, period_start = timeperiod[1], 
                               period_end = timeperiod[2], multi_period = TRUE,Time=5, nsim,nb, x, y, rMax, center, radius, IC = "ic"))
#make maps
pdfname <- paste0("figures/simulations/sim","_","center","_",center,"radius",radius,"_", "start",
                  "_",as.numeric(paste(timeperiod, collapse = "")),"_","rr","_",gsub("[.]","",risk.ratio),".pdf")
plotmap.st(pdfname, res)

########################################################################################################
########################################################################################################
########################################################################################################

#Set Some Initial Conditions
x=dframe2$utmx/1000
y=dframe2$utmy/1000
rMax=30 
Time=5  
nsim=100
center=50
radius= 18
timeperiod=c(2,4)
risk.ratio=2

res <- clust.sim(x,y,rMax,dframe$period, dframe$expdeath, dframe$death, Time, spacetime=TRUE, pois=FALSE,
                 nsim,center, radius, risk.ratio, timeperiod, colors=TRUE, utm=TRUE, byrow=TRUE)

#save results
filename <- paste0("SimulationOutput/sim","_","center","_",center,"radius",radius,"_", "start",
                   "_",as.numeric(paste(timeperiod, collapse = "")),"_","rr","_",gsub("[.]","",risk.ratio),".RData")
save(res, file = filename)

#Detection
set <- detect.set(res$lassoresult, res$init.vec, res$rr.mat, Time=5, nb, x, y, rMax, center, radius)
(incluster <- detect.incluster(res$lassoresult, res$init.vec, res$rr.mat, set, period_start = timeperiod[1], 
                               period_end = timeperiod[2], multi_period = FALSE,Time=5, nsim,nb, x, y, rMax, center, radius, IC = "ic"))
#make maps
pdfname <- paste0("figures/simulations/sim","_","center","_",center,"radius",radius,"_", "start",
                  "_",as.numeric(paste(timeperiod, collapse = "")),"_","rr","_",gsub("[.]","",risk.ratio),".pdf")
plotmap.st(pdfname, res)

########################################################################################################
########################################################################################################
########################################################################################################

#Set Some Initial Conditions
x=dframe2$utmx/1000
y=dframe2$utmy/1000
rMax=30 
Time=5  
nsim=100
center=50
radius= 11
timeperiod=c(3:5)
risk.ratio=2

res <- clust.sim(x,y,rMax,dframe$period, dframe$expdeath, dframe$death, Time, spacetime=TRUE, pois=FALSE,
                 nsim,center, radius, risk.ratio, timeperiod, colors=TRUE, utm=TRUE, byrow=TRUE)

#save results
filename <- paste0("SimulationOutput/sim","_","center","_",center,"radius",radius,"_", "start",
                   "_",as.numeric(paste(timeperiod, collapse = "")),"_","rr","_",gsub("[.]","",risk.ratio),".RData")
save(res, file = filename)

#Detection
set <- detect.set(res$lassoresult, res$init.vec, res$rr.mat, Time=5, nb, x, y, rMax, center, radius)
(incluster <- detect.incluster(res$lassoresult, res$init.vec, res$rr.mat, set, period_start = timeperiod[1], 
                               period_end = timeperiod[2], multi_period = TRUE,Time=5, nsim,nb, x, y, rMax, center, radius, IC = "ic"))
#make maps
pdfname <- paste0("figures/simulations/sim","_","center","_",center,"radius",radius,"_", "start",
                  "_",as.numeric(paste(timeperiod, collapse = "")),"_","rr","_",gsub("[.]","",risk.ratio),".pdf")
plotmap.st(pdfname, res)

########################################################################################################
########################################################################################################
########################################################################################################

#Set Some Initial Conditions
x=dframe2$utmx/1000
y=dframe2$utmy/1000
rMax=30 
Time=5  
nsim=100
center=50
radius= 9
timeperiod=c(3:5)
risk.ratio=2

res <- clust.sim(x,y,rMax,dframe$period, dframe$expdeath, dframe$death, Time, spacetime=TRUE, pois=FALSE,
                 nsim,center, radius, risk.ratio, timeperiod, colors=TRUE, utm=TRUE, byrow=TRUE)

#save results
filename <- paste0("SimulationOutput/sim","_","center","_",center,"radius",radius,"_", "start",
                   "_",as.numeric(paste(timeperiod, collapse = "")),"_","rr","_",gsub("[.]","",risk.ratio),".RData")
save(res, file = filename)

#Detection
set <- detect.set(res$lassoresult, res$init.vec, res$rr.mat, Time=5, nb, x, y, rMax, center, radius)
(incluster <- detect.incluster(res$lassoresult, res$init.vec, res$rr.mat, set, period_start = timeperiod[1], 
                               period_end = timeperiod[2], multi_period = TRUE,Time=5, nsim,nb, x, y, rMax, center, radius, IC = "ic"))
#make maps
pdfname <- paste0("figures/simulations/sim","_","center","_",center,"radius",radius,"_", "start",
                  "_",as.numeric(paste(timeperiod, collapse = "")),"_","rr","_",gsub("[.]","",risk.ratio),".pdf")
plotmap.st(pdfname, res)

########################################################################################################
########################################################################################################
########################################################################################################

#Set Some Initial Conditions
x=dframe2$utmx/1000
y=dframe2$utmy/1000
rMax=30 
Time=5  
nsim=100
center=50
radius= 11
timeperiod=c(3:5)
risk.ratio=1.5

res <- clust.sim(x,y,rMax,dframe$period, dframe$expdeath, dframe$death, Time, spacetime=TRUE, pois=FALSE,
                 nsim,center, radius, risk.ratio, timeperiod, colors=TRUE, utm=TRUE, byrow=TRUE)

#save results
filename <- paste0("SimulationOutput/sim","_","center","_",center,"radius",radius,"_", "start",
                   "_",as.numeric(paste(timeperiod, collapse = "")),"_","rr","_",gsub("[.]","",risk.ratio),".RData")
save(res, file = filename)

#Detection
set <- detect.set(res$lassoresult, res$init.vec, res$rr.mat, Time=5, nb, x, y, rMax, center, radius)
(incluster <- detect.incluster(res$lassoresult, res$init.vec, res$rr.mat, set, period_start = timeperiod[1], 
                               period_end = timeperiod[2], multi_period = TRUE,Time=5, nsim,nb, x, y, rMax, center, radius, IC = "ic"))
#make maps
pdfname <- paste0("figures/simulations/sim","_","center","_",center,"radius",radius,"_", "start",
                  "_",as.numeric(paste(timeperiod, collapse = "")),"_","rr","_",gsub("[.]","",risk.ratio),".pdf")
plotmap.st(pdfname, res)

########################################################################################################
########################################################################################################
########################################################################################################

#Set Some Initial Conditions
x=dframe2$utmx/1000
y=dframe2$utmy/1000
rMax=30 
Time=5  
nsim=100
center=50
radius= 11
timeperiod=c(3:5)
risk.ratio=1.25

res <- clust.sim(x,y,rMax,dframe$period, dframe$expdeath, dframe$death, Time, spacetime=TRUE, pois=FALSE,
                 nsim,center, radius, risk.ratio, timeperiod, colors=TRUE, utm=TRUE, byrow=TRUE)

#save results
filename <- paste0("SimulationOutput/sim","_","center","_",center,"radius",radius,"_", "start",
                   "_",as.numeric(paste(timeperiod, collapse = "")),"_","rr","_",gsub("[.]","",risk.ratio),".RData")
save(res, file = filename)

#Detection
set <- detect.set(res$lassoresult, res$init.vec, res$rr.mat, Time=5, nb, x, y, rMax, center, radius)
(incluster <- detect.incluster(res$lassoresult, res$init.vec, res$rr.mat, set, period_start = timeperiod[1], 
                               period_end = timeperiod[2], multi_period = TRUE,Time=5, nsim,nb, x, y, rMax, center, radius, IC = "ic"))
#make maps
pdfname <- paste0("figures/simulations/sim","_","center","_",center,"radius",radius,"_", "start",
                  "_",as.numeric(paste(timeperiod, collapse = "")),"_","rr","_",gsub("[.]","",risk.ratio),".pdf")
plotmap.st(pdfname, res)

########################################################################################################
########################################################################################################
########################################################################################################

#Set Some Initial Conditions
x=dframe2$utmx/1000
y=dframe2$utmy/1000
rMax=30 
Time=5  
nsim=100
center=50
radius= 11
timeperiod=c(3:5)
risk.ratio=1.1

res <- clust.sim(x,y,rMax,dframe$period, dframe$expdeath, dframe$death, Time, spacetime=TRUE, pois=FALSE,
                 nsim,center, radius, risk.ratio, timeperiod, colors=TRUE, utm=TRUE, byrow=TRUE)

#save results
filename <- paste0("SimulationOutput/sim","_","center","_",center,"radius",radius,"_", "start",
                   "_",as.numeric(paste(timeperiod, collapse = "")),"_","rr","_",gsub("[.]","",risk.ratio),".RData")
save(res, file = filename)

#Detection
set <- detect.set(res$lassoresult, res$init.vec, res$rr.mat, Time=5, nb, x, y, rMax, center, radius)
(incluster <- detect.incluster(res$lassoresult, res$init.vec, res$rr.mat, set, period_start = timeperiod[1], 
                               period_end = timeperiod[2], multi_period = TRUE,Time=5, nsim,nb, x, y, rMax, center, radius, IC = "ic"))
#make maps
pdfname <- paste0("figures/simulations/sim","_","center","_",center,"radius",radius,"_", "start",
                  "_",as.numeric(paste(timeperiod, collapse = "")),"_","rr","_",gsub("[.]","",risk.ratio),".pdf")
plotmap.st(pdfname, res)

########################################################################################################
########################################################################################################
########################################################################################################

#Set Some Initial Conditions
x=dframe2$utmx/1000
y=dframe2$utmy/1000
rMax=30 
Time=5  
nsim=100
center=50
radius= 9
timeperiod=c(1:5)
risk.ratio=1.1

res <- clust.sim(x,y,rMax,dframe$period, dframe$expdeath, dframe$death, Time, spacetime=TRUE, pois=FALSE,
                 nsim,center, radius, risk.ratio, timeperiod, colors=TRUE, utm=TRUE, byrow=TRUE)

#save results
filename <- paste0("SimulationOutput/sim","_","center","_",center,"radius",radius,"_", "start",
                   "_",as.numeric(paste(timeperiod, collapse = "")),"_","rr","_",gsub("[.]","",risk.ratio),".RData")
save(res, file = filename)

#Detection
set <- detect.set(res$lassoresult, res$init.vec, res$rr.mat, Time=5, nb, x, y, rMax, center, radius)
(incluster <- detect.incluster(res$lassoresult, res$init.vec, res$rr.mat, set, period_start = timeperiod[1], 
                               period_end = timeperiod[2], multi_period = TRUE,Time=5, nsim,nb, x, y, rMax, center, radius, IC = "ic"))
#make maps
pdfname <- paste0("figures/simulations/sim","_","center","_",center,"radius",radius,"_", "start",
                  "_",as.numeric(paste(timeperiod, collapse = "")),"_","rr","_",gsub("[.]","",risk.ratio),".pdf")
plotmap.st(pdfname, res)

########################################################################################################
########################################################################################################
########################################################################################################

#Set Some Initial Conditions
x=dframe2$utmx/1000
y=dframe2$utmy/1000
rMax=30 
Time=5  
nsim=100
center=50
radius= 9
timeperiod=c(1:5)
risk.ratio=1.1

res <- clust.sim(x,y,rMax,dframe$period, dframe$expdeath, dframe$death, Time, spacetime=TRUE, pois=FALSE,
                 nsim,center, radius, risk.ratio, timeperiod, colors=TRUE, utm=TRUE, byrow=TRUE)

#save results
filename <- paste0("SimulationOutput/sim","_","center","_",center,"radius",radius,"_", "start",
                   "_",as.numeric(paste(timeperiod, collapse = "")),"_","rr","_",gsub("[.]","",risk.ratio),".RData")
save(res, file = filename)

#Detection
set <- detect.set(res$lassoresult, res$init.vec, res$rr.mat, Time=5, nb, x, y, rMax, center, radius)
(incluster <- detect.incluster(res$lassoresult, res$init.vec, res$rr.mat, set, period_start = timeperiod[1], 
                               period_end = timeperiod[2], multi_period = TRUE,Time=5, nsim,nb, x, y, rMax, center, radius, IC = "ic"))
#make maps
pdfname <- paste0("figures/simulations/sim","_","center","_",center,"radius",radius,"_", "start",
                  "_",as.numeric(paste(timeperiod, collapse = "")),"_","rr","_",gsub("[.]","",risk.ratio),".pdf")
plotmap.st(pdfname, res)

sink()