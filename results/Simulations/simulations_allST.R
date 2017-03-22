########################################################################################################
########################################################################################################
########################################################################################################
#Space and Space-Time Cluster Detection Using the Lasso 
##Running simulations and testing functions - running all models on same simulated values for comparison
##This should later be integrated into vignettes/examples

#Maria Kamenetsky
#3-11-2017
########################################################################################################
########################################################################################################
########################################################################################################

##############################################
##############################################
#Set-Up
##############################################
##############################################
sink("simulations_compareST.txt")
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


########################################################################################################
########################################################################################################
########################################################################################################
#Global conditions

threshold <- c(0.9, 0.5)
mods <- c("QuasiPoisson", "Poisson")



#
#########################################################################################################
#########################################################################################################
#########################################################################################################
##END NULL MODELS, BEGIN CLUSTER SIMULATIONS
#########################################################################################################
#########################################################################################################
#########################################################################################################
##Center = 50
#
#####################################################
#TESTER - 1 SIM TO MAKE SURE IT WORKS
####################################################
#Set Some Initial Conditions
x=dframe2$utmx/1000
y=dframe2$utmy/1000
rMax=30
Time=5
#nsim=100
nsim=1
center=50
radius= 18
timeperiod=c(1:5)
risk.ratio=2


res <- clust.sim.all(x,y,rMax,dframe$period, dframe$expdeath, dframe$death, Time,
                 nsim,center, radius, risk.ratio, timeperiod, colors=TRUE, utm=TRUE, byrow=TRUE, threshold)

print(res$detect.out.qp.st)
print(res$detect.out.p.st)


#save results
filename <- paste0("SimulationOutput/sim","_","center","_",center,"radius",radius,"_", "start",
                   "_",as.numeric(paste(timeperiod, collapse = "")),"_","rr","_",gsub("[.]","",risk.ratio),".RData")
save(res, file = filename)


#make maps
pdfname <- paste0("figures/simulations/sim","_","center","_",center,"radius",radius,"_", "start",
                  "_",as.numeric(paste(timeperiod, collapse = "")),"_","rr","_",gsub("[.]","",risk.ratio),".pdf")
easyplot.st(pdfname, res, mods)


#
#########################################################################################################
#########################################################################################################
#########################################################################################################
#
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


res <- clust.sim.all(x,y,rMax,dframe$period, dframe$expdeath, dframe$death, Time,
                     nsim,center, radius, risk.ratio, timeperiod, colors=TRUE, utm=TRUE, byrow=TRUE, threshold)

print(res$detect.out.qp.st)
print(res$detect.out.p.st)


#save results
filename <- paste0("SimulationOutput/sim","_","center","_",center,"radius",radius,"_", "start",
                   "_",as.numeric(paste(timeperiod, collapse = "")),"_","rr","_",gsub("[.]","",risk.ratio),".RData")
save(res, file = filename)


#make maps
pdfname <- paste0("figures/simulations/sim","_","center","_",center,"radius",radius,"_", "start",
                  "_",as.numeric(paste(timeperiod, collapse = "")),"_","rr","_",gsub("[.]","",risk.ratio),".pdf")
easyplot.st(pdfname, res, mods)

#
#########################################################################################################
#########################################################################################################
#########################################################################################################
#
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


res <- clust.sim.all(x,y,rMax,dframe$period, dframe$expdeath, dframe$death, Time,
                     nsim,center, radius, risk.ratio, timeperiod, colors=TRUE, utm=TRUE, byrow=TRUE, threshold)

print(res$detect.out.qp.st)
print(res$detect.out.p.st)


#save results
filename <- paste0("SimulationOutput/sim","_","center","_",center,"radius",radius,"_", "start",
                   "_",as.numeric(paste(timeperiod, collapse = "")),"_","rr","_",gsub("[.]","",risk.ratio),".RData")
save(res, file = filename)


#make maps
pdfname <- paste0("figures/simulations/sim","_","center","_",center,"radius",radius,"_", "start",
                  "_",as.numeric(paste(timeperiod, collapse = "")),"_","rr","_",gsub("[.]","",risk.ratio),".pdf")
easyplot.st(pdfname, res, mods)


#########################################################################################################
#########################################################################################################
#########################################################################################################
#
#Set Some Initial Conditions
x=dframe2$utmx/1000
y=dframe2$utmy/1000
rMax=30
Time=5
nsim=100
center=50
radius= 18
timeperiod=c(3:5)
risk.ratio=1.5


res <- clust.sim.all(x,y,rMax,dframe$period, dframe$expdeath, dframe$death, Time,
                     nsim,center, radius, risk.ratio, timeperiod, colors=TRUE, utm=TRUE, byrow=TRUE, threshold)

print(res$detect.out.qp.st)
print(res$detect.out.p.st)


#save results
filename <- paste0("SimulationOutput/sim","_","center","_",center,"radius",radius,"_", "start",
                   "_",as.numeric(paste(timeperiod, collapse = "")),"_","rr","_",gsub("[.]","",risk.ratio),".RData")
save(res, file = filename)


#make maps
pdfname <- paste0("figures/simulations/sim","_","center","_",center,"radius",radius,"_", "start",
                  "_",as.numeric(paste(timeperiod, collapse = "")),"_","rr","_",gsub("[.]","",risk.ratio),".pdf")
easyplot.st(pdfname, res, mods)


#########################################################################################################
#########################################################################################################
#########################################################################################################
#
#Set Some Initial Conditions
x=dframe2$utmx/1000
y=dframe2$utmy/1000
rMax=30
Time=5
nsim=100
center=50
radius= 18
timeperiod=c(3:5)
risk.ratio=1.1


res <- clust.sim.all(x,y,rMax,dframe$period, dframe$expdeath, dframe$death, Time,
                     nsim,center, radius, risk.ratio, timeperiod, colors=TRUE, utm=TRUE, byrow=TRUE, threshold)

print(res$detect.out.qp.st)
print(res$detect.out.p.st)


#save results
filename <- paste0("SimulationOutput/sim","_","center","_",center,"radius",radius,"_", "start",
                   "_",as.numeric(paste(timeperiod, collapse = "")),"_","rr","_",gsub("[.]","",risk.ratio),".RData")
save(res, file = filename)


#make maps
pdfname <- paste0("figures/simulations/sim","_","center","_",center,"radius",radius,"_", "start",
                  "_",as.numeric(paste(timeperiod, collapse = "")),"_","rr","_",gsub("[.]","",risk.ratio),".pdf")
easyplot.st(pdfname, res, mods)


#########################################################################################################
#########################################################################################################
#########################################################################################################
#
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


res <- clust.sim.all(x,y,rMax,dframe$period, dframe$expdeath, dframe$death, Time,
                     nsim,center, radius, risk.ratio, timeperiod, colors=TRUE, utm=TRUE, byrow=TRUE, threshold)

print(res$detect.out.qp.st)
print(res$detect.out.p.st)


#save results
filename <- paste0("SimulationOutput/sim","_","center","_",center,"radius",radius,"_", "start",
                   "_",as.numeric(paste(timeperiod, collapse = "")),"_","rr","_",gsub("[.]","",risk.ratio),".RData")
save(res, file = filename)


#make maps
pdfname <- paste0("figures/simulations/sim","_","center","_",center,"radius",radius,"_", "start",
                  "_",as.numeric(paste(timeperiod, collapse = "")),"_","rr","_",gsub("[.]","",risk.ratio),".pdf")
easyplot.st(pdfname, res, mods)


#########################################################################################################
#########################################################################################################
#########################################################################################################
#
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


res <- clust.sim.all(x,y,rMax,dframe$period, dframe$expdeath, dframe$death, Time,
                     nsim,center, radius, risk.ratio, timeperiod, colors=TRUE, utm=TRUE, byrow=TRUE, threshold)

print(res$detect.out.qp.st)
print(res$detect.out.p.st)


#save results
filename <- paste0("SimulationOutput/sim","_","center","_",center,"radius",radius,"_", "start",
                   "_",as.numeric(paste(timeperiod, collapse = "")),"_","rr","_",gsub("[.]","",risk.ratio),".RData")
save(res, file = filename)


#make maps
pdfname <- paste0("figures/simulations/sim","_","center","_",center,"radius",radius,"_", "start",
                  "_",as.numeric(paste(timeperiod, collapse = "")),"_","rr","_",gsub("[.]","",risk.ratio),".pdf")
easyplot.st(pdfname, res, mods)

#########################################################################################################
#########################################################################################################
#########################################################################################################
#
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


res <- clust.sim.all(x,y,rMax,dframe$period, dframe$expdeath, dframe$death, Time,
                     nsim,center, radius, risk.ratio, timeperiod, colors=TRUE, utm=TRUE, byrow=TRUE, threshold)

print(res$detect.out.qp.st)
print(res$detect.out.p.st)


#save results
filename <- paste0("SimulationOutput/sim","_","center","_",center,"radius",radius,"_", "start",
                   "_",as.numeric(paste(timeperiod, collapse = "")),"_","rr","_",gsub("[.]","",risk.ratio),".RData")
save(res, file = filename)


#make maps
pdfname <- paste0("figures/simulations/sim","_","center","_",center,"radius",radius,"_", "start",
                  "_",as.numeric(paste(timeperiod, collapse = "")),"_","rr","_",gsub("[.]","",risk.ratio),".pdf")
easyplot.st(pdfname, res, mods)


#########################################################################################################
#########################################################################################################
#########################################################################################################
#
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


res <- clust.sim.all(x,y,rMax,dframe$period, dframe$expdeath, dframe$death, Time,
                     nsim,center, radius, risk.ratio, timeperiod, colors=TRUE, utm=TRUE, byrow=TRUE, threshold)

print(res$detect.out.qp.st)
print(res$detect.out.p.st)


#save results
filename <- paste0("SimulationOutput/sim","_","center","_",center,"radius",radius,"_", "start",
                   "_",as.numeric(paste(timeperiod, collapse = "")),"_","rr","_",gsub("[.]","",risk.ratio),".RData")
save(res, file = filename)


#make maps
pdfname <- paste0("figures/simulations/sim","_","center","_",center,"radius",radius,"_", "start",
                  "_",as.numeric(paste(timeperiod, collapse = "")),"_","rr","_",gsub("[.]","",risk.ratio),".pdf")
easyplot.st(pdfname, res, mods)

#########################################################################################################
#########################################################################################################
#########################################################################################################
#
#Set Some Initial Conditions
x=dframe2$utmx/1000
y=dframe2$utmy/1000
rMax=30
Time=5
nsim=100
center=50
radius= 9
timeperiod=c(3:5)
risk.ratio=1.5


res <- clust.sim.all(x,y,rMax,dframe$period, dframe$expdeath, dframe$death, Time,
                     nsim,center, radius, risk.ratio, timeperiod, colors=TRUE, utm=TRUE, byrow=TRUE, threshold)

print(res$detect.out.qp.st)
print(res$detect.out.p.st)


#save results
filename <- paste0("SimulationOutput/sim","_","center","_",center,"radius",radius,"_", "start",
                   "_",as.numeric(paste(timeperiod, collapse = "")),"_","rr","_",gsub("[.]","",risk.ratio),".RData")
save(res, file = filename)


#make maps
pdfname <- paste0("figures/simulations/sim","_","center","_",center,"radius",radius,"_", "start",
                  "_",as.numeric(paste(timeperiod, collapse = "")),"_","rr","_",gsub("[.]","",risk.ratio),".pdf")
easyplot.st(pdfname, res, mods)


#########################################################################################################
#########################################################################################################
#########################################################################################################
#
#Set Some Initial Conditions
x=dframe2$utmx/1000
y=dframe2$utmy/1000
rMax=30
Time=5
nsim=100
center=50
radius= 9
timeperiod=c(3:5)
risk.ratio=1.1


res <- clust.sim.all(x,y,rMax,dframe$period, dframe$expdeath, dframe$death, Time,
                     nsim,center, radius, risk.ratio, timeperiod, colors=TRUE, utm=TRUE, byrow=TRUE, threshold)

print(res$detect.out.qp.st)
print(res$detect.out.p.st)


#save results
filename <- paste0("SimulationOutput/sim","_","center","_",center,"radius",radius,"_", "start",
                   "_",as.numeric(paste(timeperiod, collapse = "")),"_","rr","_",gsub("[.]","",risk.ratio),".RData")
save(res, file = filename)


#make maps
pdfname <- paste0("figures/simulations/sim","_","center","_",center,"radius",radius,"_", "start",
                  "_",as.numeric(paste(timeperiod, collapse = "")),"_","rr","_",gsub("[.]","",risk.ratio),".pdf")
easyplot.st(pdfname, res, mods)

###########################################################################################
#VARY THE CENTER
###########################################################################################

#
#########################################################################################################
#########################################################################################################
#########################################################################################################
#
#Set Some Initial Conditions
x=dframe2$utmx/1000
y=dframe2$utmy/1000
rMax=30
Time=5
nsim=100
center=150
radius= 18
timeperiod=c(1:5)
risk.ratio=2


res <- clust.sim.all(x,y,rMax,dframe$period, dframe$expdeath, dframe$death, Time,
                     nsim,center, radius, risk.ratio, timeperiod, colors=TRUE, utm=TRUE, byrow=TRUE, threshold)

print(res$detect.out.qp.st)
print(res$detect.out.p.st)


#save results
filename <- paste0("SimulationOutput/sim","_","center","_",center,"radius",radius,"_", "start",
                   "_",as.numeric(paste(timeperiod, collapse = "")),"_","rr","_",gsub("[.]","",risk.ratio),".RData")
save(res, file = filename)


#make maps
pdfname <- paste0("figures/simulations/sim","_","center","_",center,"radius",radius,"_", "start",
                  "_",as.numeric(paste(timeperiod, collapse = "")),"_","rr","_",gsub("[.]","",risk.ratio),".pdf")
easyplot.st(pdfname, res, mods)

#
#########################################################################################################
#########################################################################################################
#########################################################################################################
#
#Set Some Initial Conditions
x=dframe2$utmx/1000
y=dframe2$utmy/1000
rMax=30
Time=5
nsim=100
center=150
radius= 18
timeperiod=c(3:5)
risk.ratio=2


res <- clust.sim.all(x,y,rMax,dframe$period, dframe$expdeath, dframe$death, Time,
                     nsim,center, radius, risk.ratio, timeperiod, colors=TRUE, utm=TRUE, byrow=TRUE, threshold)

print(res$detect.out.qp.st)
print(res$detect.out.p.st)


#save results
filename <- paste0("SimulationOutput/sim","_","center","_",center,"radius",radius,"_", "start",
                   "_",as.numeric(paste(timeperiod, collapse = "")),"_","rr","_",gsub("[.]","",risk.ratio),".RData")
save(res, file = filename)


#make maps
pdfname <- paste0("figures/simulations/sim","_","center","_",center,"radius",radius,"_", "start",
                  "_",as.numeric(paste(timeperiod, collapse = "")),"_","rr","_",gsub("[.]","",risk.ratio),".pdf")
easyplot.st(pdfname, res, mods)


#########################################################################################################
#########################################################################################################
#########################################################################################################
#
#Set Some Initial Conditions
x=dframe2$utmx/1000
y=dframe2$utmy/1000
rMax=30
Time=5
nsim=100
center=150
radius= 18
timeperiod=c(3:5)
risk.ratio=1.5


res <- clust.sim.all(x,y,rMax,dframe$period, dframe$expdeath, dframe$death, Time,
                     nsim,center, radius, risk.ratio, timeperiod, colors=TRUE, utm=TRUE, byrow=TRUE, threshold)

print(res$detect.out.qp.st)
print(res$detect.out.p.st)


#save results
filename <- paste0("SimulationOutput/sim","_","center","_",center,"radius",radius,"_", "start",
                   "_",as.numeric(paste(timeperiod, collapse = "")),"_","rr","_",gsub("[.]","",risk.ratio),".RData")
save(res, file = filename)


#make maps
pdfname <- paste0("figures/simulations/sim","_","center","_",center,"radius",radius,"_", "start",
                  "_",as.numeric(paste(timeperiod, collapse = "")),"_","rr","_",gsub("[.]","",risk.ratio),".pdf")
easyplot.st(pdfname, res, mods)


#########################################################################################################
#########################################################################################################
#########################################################################################################
#
#Set Some Initial Conditions
x=dframe2$utmx/1000
y=dframe2$utmy/1000
rMax=30
Time=5
nsim=100
center=150
radius= 18
timeperiod=c(3:5)
risk.ratio=1.1


res <- clust.sim.all(x,y,rMax,dframe$period, dframe$expdeath, dframe$death, Time,
                     nsim,center, radius, risk.ratio, timeperiod, colors=TRUE, utm=TRUE, byrow=TRUE, threshold)

print(res$detect.out.qp.st)
print(res$detect.out.p.st)


#save results
filename <- paste0("SimulationOutput/sim","_","center","_",center,"radius",radius,"_", "start",
                   "_",as.numeric(paste(timeperiod, collapse = "")),"_","rr","_",gsub("[.]","",risk.ratio),".RData")
save(res, file = filename)


#make maps
pdfname <- paste0("figures/simulations/sim","_","center","_",center,"radius",radius,"_", "start",
                  "_",as.numeric(paste(timeperiod, collapse = "")),"_","rr","_",gsub("[.]","",risk.ratio),".pdf")
easyplot.st(pdfname, res, mods)


#########################################################################################################
#########################################################################################################
#########################################################################################################
#
#Set Some Initial Conditions
x=dframe2$utmx/1000
y=dframe2$utmy/1000
rMax=30
Time=5
nsim=100
center=150
radius= 11
timeperiod=c(3:5)
risk.ratio=2


res <- clust.sim.all(x,y,rMax,dframe$period, dframe$expdeath, dframe$death, Time,
                     nsim,center, radius, risk.ratio, timeperiod, colors=TRUE, utm=TRUE, byrow=TRUE, threshold)

print(res$detect.out.qp.st)
print(res$detect.out.p.st)


#save results
filename <- paste0("SimulationOutput/sim","_","center","_",center,"radius",radius,"_", "start",
                   "_",as.numeric(paste(timeperiod, collapse = "")),"_","rr","_",gsub("[.]","",risk.ratio),".RData")
save(res, file = filename)


#make maps
pdfname <- paste0("figures/simulations/sim","_","center","_",center,"radius",radius,"_", "start",
                  "_",as.numeric(paste(timeperiod, collapse = "")),"_","rr","_",gsub("[.]","",risk.ratio),".pdf")
easyplot.st(pdfname, res, mods)


#########################################################################################################
#########################################################################################################
#########################################################################################################
#
#Set Some Initial Conditions
x=dframe2$utmx/1000
y=dframe2$utmy/1000
rMax=30
Time=5
nsim=100
center=150
radius= 11
timeperiod=c(3:5)
risk.ratio=1.5


res <- clust.sim.all(x,y,rMax,dframe$period, dframe$expdeath, dframe$death, Time,
                     nsim,center, radius, risk.ratio, timeperiod, colors=TRUE, utm=TRUE, byrow=TRUE, threshold)

print(res$detect.out.qp.st)
print(res$detect.out.p.st)


#save results
filename <- paste0("SimulationOutput/sim","_","center","_",center,"radius",radius,"_", "start",
                   "_",as.numeric(paste(timeperiod, collapse = "")),"_","rr","_",gsub("[.]","",risk.ratio),".RData")
save(res, file = filename)


#make maps
pdfname <- paste0("figures/simulations/sim","_","center","_",center,"radius",radius,"_", "start",
                  "_",as.numeric(paste(timeperiod, collapse = "")),"_","rr","_",gsub("[.]","",risk.ratio),".pdf")
easyplot.st(pdfname, res, mods)

#########################################################################################################
#########################################################################################################
#########################################################################################################
#
#Set Some Initial Conditions
x=dframe2$utmx/1000
y=dframe2$utmy/1000
rMax=30
Time=5
nsim=100
center=150
radius= 11
timeperiod=c(3:5)
risk.ratio=1.1


res <- clust.sim.all(x,y,rMax,dframe$period, dframe$expdeath, dframe$death, Time,
                     nsim,center, radius, risk.ratio, timeperiod, colors=TRUE, utm=TRUE, byrow=TRUE, threshold)

print(res$detect.out.qp.st)
print(res$detect.out.p.st)


#save results
filename <- paste0("SimulationOutput/sim","_","center","_",center,"radius",radius,"_", "start",
                   "_",as.numeric(paste(timeperiod, collapse = "")),"_","rr","_",gsub("[.]","",risk.ratio),".RData")
save(res, file = filename)


#make maps
pdfname <- paste0("figures/simulations/sim","_","center","_",center,"radius",radius,"_", "start",
                  "_",as.numeric(paste(timeperiod, collapse = "")),"_","rr","_",gsub("[.]","",risk.ratio),".pdf")
easyplot.st(pdfname, res, mods)


#########################################################################################################
#########################################################################################################
#########################################################################################################
#
#Set Some Initial Conditions
x=dframe2$utmx/1000
y=dframe2$utmy/1000
rMax=30
Time=5
nsim=100
center=150
radius= 9
timeperiod=c(3:5)
risk.ratio=2


res <- clust.sim.all(x,y,rMax,dframe$period, dframe$expdeath, dframe$death, Time,
                     nsim,center, radius, risk.ratio, timeperiod, colors=TRUE, utm=TRUE, byrow=TRUE, threshold)

print(res$detect.out.qp.st)
print(res$detect.out.p.st)


#save results
filename <- paste0("SimulationOutput/sim","_","center","_",center,"radius",radius,"_", "start",
                   "_",as.numeric(paste(timeperiod, collapse = "")),"_","rr","_",gsub("[.]","",risk.ratio),".RData")
save(res, file = filename)


#make maps
pdfname <- paste0("figures/simulations/sim","_","center","_",center,"radius",radius,"_", "start",
                  "_",as.numeric(paste(timeperiod, collapse = "")),"_","rr","_",gsub("[.]","",risk.ratio),".pdf")
easyplot.st(pdfname, res, mods)

#########################################################################################################
#########################################################################################################
#########################################################################################################
#
#Set Some Initial Conditions
x=dframe2$utmx/1000
y=dframe2$utmy/1000
rMax=30
Time=5
nsim=100
center=150
radius= 9
timeperiod=c(3:5)
risk.ratio=1.5


res <- clust.sim.all(x,y,rMax,dframe$period, dframe$expdeath, dframe$death, Time,
                     nsim,center, radius, risk.ratio, timeperiod, colors=TRUE, utm=TRUE, byrow=TRUE, threshold)

print(res$detect.out.qp.st)
print(res$detect.out.p.st)


#save results
filename <- paste0("SimulationOutput/sim","_","center","_",center,"radius",radius,"_", "start",
                   "_",as.numeric(paste(timeperiod, collapse = "")),"_","rr","_",gsub("[.]","",risk.ratio),".RData")
save(res, file = filename)


#make maps
pdfname <- paste0("figures/simulations/sim","_","center","_",center,"radius",radius,"_", "start",
                  "_",as.numeric(paste(timeperiod, collapse = "")),"_","rr","_",gsub("[.]","",risk.ratio),".pdf")
easyplot.st(pdfname, res, mods)


#########################################################################################################
#########################################################################################################
#########################################################################################################
#
#Set Some Initial Conditions
x=dframe2$utmx/1000
y=dframe2$utmy/1000
rMax=30
Time=5
#nsim=100
nsim=1
center=150
radius= 9
timeperiod=c(3:5)
#risk.ratio=1.1
risk.ratio=2


res <- clust.sim.all(x,y,rMax,dframe$period, dframe$expdeath, dframe$death, Time,
                     nsim,center, radius, risk.ratio, timeperiod, colors=TRUE, utm=TRUE, byrow=TRUE, threshold, space= "space")

#save results
sim.i <- paste0("sim","_","center","_",center,"radius",radius,"_", "start",
                   "_",as.numeric(paste(timeperiod, collapse = "")),"_","rr","_",gsub("[.]","",risk.ratio))
filename <- paste0("SimulationOutput/",sim.i,".RData")
save(res, file = filename)

#Print Detection for the Simulation
cbind(radius,risk.ratio,center,time=as.numeric(paste(timeperiod, collapse = "")), mod = "ST",
      rbind("QuasiPois",res$detect.out.qp.st), rbind("Pois",res$detect.out.qp.st))


#make maps
pdfname <- paste0("figures/simulations/sim","_","center","_",center,"radius",radius,"_", "start",
                  "_",as.numeric(paste(timeperiod, collapse = "")),"_","rr","_",gsub("[.]","",risk.ratio),".pdf")
easyplot.st(pdfname, res, mods)



sink()    
