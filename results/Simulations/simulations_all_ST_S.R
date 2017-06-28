########################################################################################################
########################################################################################################
########################################################################################################
#Space and Space-Time Cluster Detection Using the Lasso 
##Running simulations and testing functions - running all models on same simulated values for comparison
##This should later be integrated into vignettes/examples

#Maria Kamenetsky
#3-21-2017
########################################################################################################
########################################################################################################
########################################################################################################

sink("multiclusters.txt")

##############################################
##############################################
#Set-Up
##############################################
##############################################
#sink("simulations_compareALL.txt")


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


#temporarily source my clustR files
file.sources = list.files(path="scripts/cluST/R/.",pattern="*.R", full.names = TRUE)
sapply(file.sources, source, .GlobalEnv)

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
##TESTER
#########################################################################################################
#########################################################################################################
#########################################################################################################
#####################################################
#TESTER - 1 SIM TO MAKE SURE IT WORKS
####################################################
#Set Some Initial Conditions
x=dframe2$utmx/1000
y=dframe2$utmy/1000
rMax=20
Time=5
#nsim=100
nsim=2
center=50
radius= 18
timeperiod=c(1:5)
risk.ratio=1.5

# center = 150
# radius =9
# timeperiod=c(3:5)
# risk.ratio = 1.1

#
# period <- dframe$period
# expected <- dframe$expdeath
# observed <- dframe$death

res <- clust.sim.all(x,y,rMax,dframe$period, dframe$expdeath, dframe$death, Time,
                    nsim,center, radius, risk.ratio, timeperiod, utm=TRUE, byrow=TRUE, threshold, space= "both", nullmod=TRUE)

#save results
(sim.i <- paste0("sim","_","center","_",center,"radius",radius,"_", "start",
                 "_",as.numeric(paste(timeperiod, collapse = "")),"_","rr","_",gsub("[.]","",risk.ratio)))
filename <- paste0("SimulationOutput/",sim.i,".RData")
save(res, file = filename)

#Print Detection for the Simulation
(tabn <- rbind("******",cbind(radius,risk.ratio,center,time=as.numeric(paste(timeperiod, collapse = "")), mod = "ST",
                              rbind("QuasiPois",res$detect.out.qp.st), rbind("Pois",res$detect.out.qp.st)), 
               cbind(radius,risk.ratio,center,time=as.numeric(paste(timeperiod, collapse = "")), mod = "Space",
                     rbind("QuasiPois",res$detect.out.qp.s), rbind("Pois",res$detect.out.qp.s))))

#make maps
pdfname <- paste0("figures/simulations/sim","_","center","_",center,"radius",radius,"_", "start",
                  "_",as.numeric(paste(timeperiod, collapse = "")),"_","rr","_",gsub("[.]","",risk.ratio),".pdf")
easyplot(pdfname, res, mods, space="both")


#########################################################################################################
#########################################################################################################
#########################################################################################################
##BEGIN NULL MODELS
#########################################################################################################
#########################################################################################################
#########################################################################################################


#Set Some Initial Conditions
x=dframe2$utmx/1000
y=dframe2$utmy/1000
rMax=20
Time=5
nsim=100

#need inputs for model to run
center=1
radius=18
timeperiod = c(1:5)
risk=1

table.detection.null <- NULL



system.time(res <- clust.sim.all(x,y,rMax,dframe$period, dframe$expdeath, dframe$death, Time,
                                 nsim,center, radius, risk, timeperiod,
                                 utm=TRUE, byrow=TRUE, threshold, space= "both", nullmod = TRUE))
#save results
(sim.i <- paste0("sim","_","center","_",center,"radius",radius,"_", "start",
                 "_",as.numeric(paste(timeperiod, collapse = "")),"_","rr","_",gsub("[.]","",risk), "NULL"))
filename <- paste0("SimulationOutput/",sim.i,".RData")
save(res, file = filename)

#Print Detection for the Simulation
(tabn <- rbind("******",cbind(radius,risk,center,time=as.numeric(paste(timeperiod, collapse = "")), mod = "ST",
                              rbind("QuasiPois",res$detect.out.qp.st), rbind("Pois",res$detect.out.p.st)),
               cbind(radius,risk,center,time=as.numeric(paste(timeperiod, collapse = "")), mod = "Space",
                     rbind("QuasiPois",res$detect.out.qp.s), rbind("Pois",res$detect.out.p.s))))
table.detection.null <- rbind(table.detection.null, tabn)

#make maps
pdfname <- paste0("figures/simulations/sim","_","center","_",center,"radius",radius,"_", "start",
                  "_",as.numeric(paste(timeperiod, collapse = "")),"_","rr","_",gsub("[.]","",risk), "NULL" ,".pdf")
easyplot(pdfname, res, mods, space="both")


#WRITE TO CSV
print(table.detection.null)
write.csv(table.detection.null, file="tabledetectionNULL.csv", row.names=TRUE)
save(table.detection.null, file="tabledetectionNULL.RData")





#########################################################################################################
#########################################################################################################
#########################################################################################################
##END NULL MODELS, BEGIN CLUSTER SIMULATIONS
#########################################################################################################
#########################################################################################################
#########################################################################################################
#########################################################################################################
#########################################################################################################
#########################################################################################################
#
#Set Some Initial Conditions
x=dframe2$utmx/1000
y=dframe2$utmy/1000
rMax=20
Time=5
nsim=100
#nsim=1

centers <- c(150, 35)
radii <- c(9, 11, 18)
timeperiods <- list(c(3:5), c(1:2), c(2:4))
risk.ratios <- c(1.1, 1.5, 2)
nullmod = NULL

table.detection <- NULL


for(cent in centers){
    for(rad in radii){
        for(tim in timeperiods){
            for(risk in risk.ratios){
                print(c(cent, rad, tim, risk))
                system.time(res <- clust.sim.all(x,y,rMax,dframe$period, dframe$expdeath, dframe$death, Time,
                                     nsim,cent, rad, risk, tim, utm=TRUE, byrow=TRUE, threshold, space= "both"))
                #save results
                (sim.i <- paste0("sim","_","center","_",cent,"radius",rad,"_", "start",
                                 "_",as.numeric(paste(tim, collapse = "")),"_","rr","_",gsub("[.]","",risk)))
                filename <- paste0("SimulationOutput/",sim.i,".RData")
                save(res, file = filename)

                #Print Detection for the Simulation
                (tabn <- rbind("******",cbind(rad,risk,cent,time=as.numeric(paste(tim, collapse = "")), mod = "ST",
                                               rbind("QuasiPois",res$detect.out.qp.st), rbind("Pois",res$detect.out.p.st)),
                                cbind(rad,risk,cent,time=as.numeric(paste(tim, collapse = "")), mod = "Space",
                                      rbind("QuasiPois",res$detect.out.qp.s), rbind("Pois",res$detect.out.p.s))))
                table.detection <- rbind(table.detection, tabn)
                
                #Print Descriptives

                #make maps
                pdfname <- paste0("figures/simulations/sim","_","center","_",cent,"radius",rad,"_", "start",
                                  "_",as.numeric(paste(tim, collapse = "")),"_","rr","_",gsub("[.]","",risk),".pdf")
                easyplot(pdfname, res, mods, space="both")

            }
        }
    }
}



#WRITE TO CSV
print(table.detection)
write.csv(table.detection, file="tabledetectionclustersim.csv", row.names=TRUE)
save(table.detection, file="tabledetectionclustersim.RData")

#Print Table Detection
table.detection
print(table.detection)

#########################################################################################################
#########################################################################################################
#########################################################################################################
##BEGIN CLUSTER SIMULATIONS - MULTICLUSTER MODELS
#########################################################################################################
#########################################################################################################
#########################################################################################################
#########################################################################################################
#########################################################################################################
#########################################################################################################

#Set Some Initial Conditions
x=dframe2$utmx/1000
y=dframe2$utmy/1000
rMax=30
Time=5
nsim=100
nullmod = NULL

centers <- c(list(c(150, 35)), list(c(50, 100)))
radii <- c(9, 18)
timeperiods <- list(c(2:4))
risk.ratios <- c(1.1, 1.5, 2)

table.detection <- NULL



for(cent in centers){
    for(rad in radii){
        for(tim in timeperiods){
            for(risk in risk.ratios){
                print(c(cent, rad, tim, risk))
                system.time(res <- clust.sim.all(x,y,rMax,dframe$period, dframe$expdeath, dframe$death, Time,
                                     nsim,cent, rad, risk, tim, colors=TRUE,
                                     utm=TRUE, byrow=TRUE, threshold, space= "both", nullmod=NULL))
                #save results
                (sim.i <- paste0("sim","_","center","_",as.numeric(paste(unlist(cent),collapse="")),"radius",rad,"_", "start",
                                 "_",as.numeric(paste(tim, collapse = "")),"_","rr","_",gsub("[.]","",risk), "multclust"))
                filename <- paste0("SimulationOutput/",sim.i,".RData")
                save(res, file = filename)

                #Print Detection for the Simulation
                (tabn <- rbind("******",cbind(rad,risk,cent=paste(cent, collapse=""),time=as.numeric(paste(tim, collapse = "")), mod = "ST",
                                               rbind("QuasiPois",res$detect.out.qp.st), rbind("Pois",res$detect.out.p.st)),
                                cbind(rad,risk,paste(unlist(cent), collapse=""),time=as.numeric(paste(tim, collapse = "")), mod = "Space",
                                      rbind("QuasiPois",res$detect.out.qp.s), rbind("Pois",res$detect.out.p.s))))
                table.detection <- rbind(table.detection, tabn)

                #make maps
                pdfname <- paste0("figures/simulations/sim","_","center","_",as.numeric(paste(unlist(cent),collapse="")),"radius",rad,"_", "start",
                                  "_",as.numeric(paste(tim, collapse = "")),"_","rr","_",gsub("[.]","",risk),"multclust",".pdf")
                easyplot(pdfname, res, mods, space="both")

            }
        }
    }
}

#WRITE TO CSV
print(table.detection)
write.csv(table.detection, file="tabledetectionmultclust.csv", row.names=TRUE)
save(table.detection, file="tabledetectionmultclust.RData")

#Print Table Detection
table.detection
print(table.detection)
#########################################################################################################
#########################################################################################################
#########################################################################################################

sink()

