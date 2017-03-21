########################################################################################################
########################################################################################################
########################################################################################################
#Space and Space-Time Cluster Detection Using the Lasso 
##Simulation Diagnostics

#Maria Kamenetsky
#3-6-2016
########################################################################################################
########################################################################################################
########################################################################################################
sink("diagnostics05.txt")
##############################################
##############################################
#Set-Up
##############################################
##############################################


#temporarily source my cluST.R file
source("scripts/cluST//R//cluST.R")

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

####################################################
#SET GLOBAL DETECTION THRESHOLD
####################################################
threshold <- 0.5

#####################################################
##NULL MODEL - Quasi-Poisson
#####################################################
#Set Some Initial Conditions
x=dframe2$utmx/1000
y=dframe2$utmy/1000
rMax=30
Time=5
nsim=2
#nsim=2
center=200
radius=18
timeperiod = c(1:5)
risk.ratio=1


(filename <- paste0("SimulationOutput/sim","_","center","_",center,"radius",radius,"_", "start",
                   "_",as.numeric(paste(timeperiod, collapse = "")),"_","rr","_",gsub("[.]","",risk.ratio),"NULL",".RData"))
load(filename)


#Detection
set <- detect.set(res$lassoresult, res$init.vec, res$rr.mat, Time=5, nb, x, y, rMax, 0, radius,nullmod=TRUE)
incluster <- detect.incluster(res$lassoresult, res$init.vec, res$rr.mat, set, timeperiod,Time=5, nsim,nb, x, y, rMax, 0, 
                               radius, IC = "ic", nullmod = TRUE)
(detect <- clust.diagnostics(incluster, threshold))


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


(filename <- paste0("SimulationOutput/sim","_","center","_",center,"radius",radius,"_", "start",
                   "_",as.numeric(paste(timeperiod, collapse = "")),"_","rr","_",gsub("[.]","",risk.ratio),"NULL",".RData"))
load(filename)

#Detection
set <- detect.set(res$lassoresult, res$init.vec, res$rr.mat, Time=5, nb, x, y, rMax, 0, radius, nullmod=TRUE)
incluster <- detect.incluster(res$lassoresult, res$init.vec, res$rr.mat, set, timeperiod,Time=5, nsim,nb, 
                               x, y, rMax, 0, radius, IC = "ic", nullmod=TRUE)
(detect <- clust.diagnostics(incluster, threshold))

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


(filename <- paste0("SimulationOutput/sim","_","center","_",center,"radius",radius,"_", "start",
                   "_",as.numeric(paste(timeperiod, collapse = "")),"_","rr","_",gsub("[.]","",risk.ratio),"NULLPOISSON",".RData"))
load(filename)

#Detection
set <- detect.set(res$lassoresult, res$init.vec, res$rr.mat, Time=5, nb, x, y, rMax, 0, radius, nullmod=TRUE)
incluster <- detect.incluster(res$lassoresult, res$init.vec, res$rr.mat, set, timeperiod,Time=5, 
                               nsim,nb, x, y, rMax, 0, radius, IC = "ic", nullmod=TRUE)
(detect <- clust.diagnostics(incluster, threshold))

########################################################################################################
########################################################################################################
########################################################################################################

###################################################
#NULL spaceonly - Quasi-Poisson
###################################################
#Set Some Initial Conditions
x=dframe2$utmx/1000
y=dframe2$utmy/1000
rMax=30
Time=1
nsim=100
center=200
radius=18
timeperiod=1
period_start = 1
risk.ratio=1


#create average dataframe for spaceonly

death <- with(dframe, tapply(death, id, function(x) round(mean(x))))
expdeath <- with(dframe, tapply(expdeath, id, function(x) mean(x)))
df <- cbind.data.frame(id = unique(dframe$id), period = rep("1", length(unique(dframe$id))), death = death, expdeath=expdeath)


(filename <- paste0("SimulationOutput/sim","_","center","_",center,"radius",radius,"_", "start",
                   "_",as.numeric(paste(timeperiod, collapse = "")),"_","rr","_",gsub("[.]","",risk.ratio),"spaceonlyQPOIS",".RData"))
load(filename)

#Detection
set <- detect.set(res$lassoresult, res$init.vec, res$rr.mat, Time=1, nb, x, y, rMax, 0, radius, nullmod=TRUE)
incluster <- detect.incluster(res$lassoresult, res$init.vec, res$rr.mat, set, timeperiod,Time=1, nsim,nb, x, y, rMax, 0, 
                               radius, IC = "ic", space=FALSE, nullmod=TRUE)
(detect <- clust.diagnostics(incluster, threshold))


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
period_start=1
risk.ratio=1


#create average dataframe for spaceonly

death <- with(dframe, tapply(death, id, function(x) round(mean(x))))
expdeath <- with(dframe, tapply(expdeath, id, function(x) mean(x)))
df <- cbind.data.frame(id = unique(dframe$id), period = rep("1", length(unique(dframe$id))), death = death, expdeath=expdeath)

(filename <- paste0("SimulationOutput/sim","_","center","_",center,"radius",radius,"_", "start",
                   "_",as.numeric(paste(timeperiod, collapse = "")),"_","rr","_",gsub("[.]","",risk.ratio),"spaceonlyPOISSONONLY",".RData"))
load(filename)

#Detection
set <- detect.set(res$lassoresult, res$init.vec, res$rr.mat, Time=1, nb, x, y, rMax, 0, radius, nullmod=TRUE)
incluster <- detect.incluster(res$lassoresult, res$init.vec, res$rr.mat, set, timeperiod,Time=1, nsim,nb,
                               x, y, rMax, 0, radius, IC = "ic", space=FALSE, nullmod=TRUE)
(detect <- clust.diagnostics(incluster, threshold))

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
#
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
#threshold = 0.9

(filename <- paste0("SimulationOutput/sim","_","center","_",center,"radius",radius,"_", "start",
                   "_",as.numeric(paste(timeperiod, collapse = "")),"_","rr","_",gsub("[.]","",risk.ratio),".RData"))
load(filename)

#Detection
set <- detect.set(res$lassoresult, res$init.vec, res$rr.mat, Time=5, nb, x, y, rMax, center, radius)
incluster <- detect.incluster(res$lassoresult, res$init.vec, res$rr.mat, set, timeperiod,Time=5, nsim,nb, x, y, rMax, center, radius, IC = "ic")
(detect <- clust.diagnostics(incluster, threshold))

#
#########################################################################################################
#########################################################################################################
#########################################################################################################
#
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

(filename <- paste0("SimulationOutput/sim","_","center","_",center,"radius",radius,"_", "start",
                   "_",as.numeric(paste(timeperiod, collapse = "")),"_","rr","_",gsub("[.]","",risk.ratio),".RData"))
load(filename)

#Detection
set <- detect.set(res$lassoresult, res$init.vec, res$rr.mat, Time=5, nb, x, y, rMax, center, radius)
incluster <- detect.incluster(res$lassoresult, res$init.vec, res$rr.mat, set, timeperiod,Time=5, nsim,nb, x, y, rMax, center, radius, IC = "ic")
(detect <- clust.diagnostics(incluster, threshold))


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
timeperiod=5
risk.ratio=2

(filename <- paste0("SimulationOutput/sim","_","center","_",center,"radius",radius,"_", "start",
                   "_",as.numeric(paste(timeperiod, collapse = "")),"_","rr","_",gsub("[.]","",risk.ratio),".RData"))
load(filename)

#Detection
set <- detect.set(res$lassoresult, res$init.vec, res$rr.mat, Time=5, nb, x, y, rMax, center, radius)
incluster <- detect.incluster(res$lassoresult, res$init.vec, res$rr.mat, set, timeperiod,Time=5, nsim,nb, x, y, rMax, center, radius, IC = "ic")
(detect <- clust.diagnostics(incluster, threshold))


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

(filename <- paste0("SimulationOutput/sim","_","center","_",center,"radius",radius,"_", "start",
                   "_",as.numeric(paste(timeperiod, collapse = "")),"_","rr","_",gsub("[.]","",risk.ratio),".RData"))
load(filename)

#Detection
set <- detect.set(res$lassoresult, res$init.vec, res$rr.mat, Time=5, nb, x, y, rMax, center, radius)
incluster <- detect.incluster(res$lassoresult, res$init.vec, res$rr.mat, set, timeperiod, Time=5, nsim,nb, x, y, rMax, center, radius, IC = "ic")
(detect <- clust.diagnostics(incluster, threshold))


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

(filename <- paste0("SimulationOutput/sim","_","center","_",center,"radius",radius,"_", "start",
                   "_",as.numeric(paste(timeperiod, collapse = "")),"_","rr","_",gsub("[.]","",risk.ratio),".RData"))
load(filename)

#Detection
set <- detect.set(res$lassoresult, res$init.vec, res$rr.mat, Time=5, nb, x, y, rMax, center, radius)
incluster <- detect.incluster(res$lassoresult, res$init.vec, res$rr.mat, set, timeperiod,Time=5, nsim,nb, x, y, rMax, center, radius, IC = "ic")
(detect <- clust.diagnostics(incluster, threshold))


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


(filename <- paste0("SimulationOutput/sim","_","center","_",center,"radius",radius,"_", "start",
                   "_",as.numeric(paste(timeperiod, collapse = "")),"_","rr","_",gsub("[.]","",risk.ratio),".RData"))
load(filename)

#Detection
set <- detect.set(res$lassoresult, res$init.vec, res$rr.mat, Time=5, nb, x, y, rMax, center, radius)
incluster <- detect.incluster(res$lassoresult, res$init.vec, res$rr.mat, set, timeperiod, Time=5, nsim,nb, x, y, rMax, center, radius, IC = "ic")
(detect <- clust.diagnostics(incluster, threshold))

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

(filename <- paste0("SimulationOutput/sim","_","center","_",center,"radius",radius,"_", "start",
                   "_",as.numeric(paste(timeperiod, collapse = "")),"_","rr","_",gsub("[.]","",risk.ratio),".RData"))
load(filename)

#Detection
set <- detect.set(res$lassoresult, res$init.vec, res$rr.mat, Time=5, nb, x, y, rMax, center, radius)
incluster <- detect.incluster(res$lassoresult, res$init.vec, res$rr.mat, set, timeperiod,Time=5, nsim,nb, x, y, rMax, center, radius, IC = "ic")
(detect <- clust.diagnostics(incluster, threshold))

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

(filename <- paste0("SimulationOutput/sim","_","center","_",center,"radius",radius,"_", "start",
                   "_",as.numeric(paste(timeperiod, collapse = "")),"_","rr","_",gsub("[.]","",risk.ratio),".RData"))
load(filename)

#Detection
set <- detect.set(res$lassoresult, res$init.vec, res$rr.mat, Time=5, nb, x, y, rMax, center, radius)
incluster <- detect.incluster(res$lassoresult, res$init.vec, res$rr.mat, set, timeperiod, multi_period = TRUE,Time=5, nsim,nb, x, y, rMax, center, radius, IC = "ic")
(detect <- clust.diagnostics(incluster, threshold))

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

(filename <- paste0("SimulationOutput/sim","_","center","_",center,"radius",radius,"_", "start",
                   "_",as.numeric(paste(timeperiod, collapse = "")),"_","rr","_",gsub("[.]","",risk.ratio),".RData"))
load(filename)

#Detection
set <- detect.set(res$lassoresult, res$init.vec, res$rr.mat, Time=5, nb, x, y, rMax, center, radius)
incluster <- detect.incluster(res$lassoresult, res$init.vec, res$rr.mat, set,timeperiod,Time=5, nsim,nb, x, y, rMax, center, radius, IC = "ic")
(detect <- clust.diagnostics(incluster, threshold))

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


#save results
(filename <- paste0("SimulationOutput/sim","_","center","_",center,"radius",radius,"_", "start",
                   "_",as.numeric(paste(timeperiod, collapse = "")),"_","rr","_",gsub("[.]","",risk.ratio),".RData"))
load(filename)

#Detection
set <- detect.set(res$lassoresult, res$init.vec, res$rr.mat, Time=5, nb, x, y, rMax, center, radius)
incluster <- detect.incluster(res$lassoresult, res$init.vec, res$rr.mat, set, timeperiod,Time=5, nsim,nb, x, y, rMax, center, radius, IC = "ic")
(detect <- clust.diagnostics(incluster, threshold))

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

#save results
(filename <- paste0("SimulationOutput/sim","_","center","_",center,"radius",radius,"_", "start",
                   "_",as.numeric(paste(timeperiod, collapse = "")),"_","rr","_",gsub("[.]","",risk.ratio),".RData"))
load(filename)

#Detection
set <- detect.set(res$lassoresult, res$init.vec, res$rr.mat, Time=5, nb, x, y, rMax, center, radius)
incluster <- detect.incluster(res$lassoresult, res$init.vec, res$rr.mat, set, timeperiod,Time=5, nsim,nb, x, y, rMax, center, radius, IC = "ic")
(detect <- clust.diagnostics(incluster, threshold))


sink()









