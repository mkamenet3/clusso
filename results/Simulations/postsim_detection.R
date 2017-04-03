####################################################
#Post-Sim Detection: 
#M.Kamenetsky
#3-26-17
####################################################


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
library(stringr)
library(stringi)
library(plyr)

#Source .cpp files
sourceCpp("scripts/cluST/src/maxcol.cpp")
sourceCpp("scripts/cluST/src/st_matCpp.cpp")
sourceCpp("scripts/cluST/src/prod_yx.cpp")


#temporarily source my clustR files
file.sources = list.files(path="scripts/cluST/R/.",pattern="*.R", full.names = TRUE)
sapply(file.sources, source, .GlobalEnv)



####################################################
#LOAD JBC Data and Set Up
####################################################

#Load Data
dframe1 <- read.csv("data/JBC/jap.breast.F.9.10.11.csv")
dframe2 <- read.csv("data//JBC//utmJapan.csv")
dframe3 <- aggregate(dframe1, by=list(as.factor(rep(1:(nrow(dframe1)/4),each=4))), FUN="sum")
dframe=data.frame(id=as.factor(dframe3$id/4),period=as.factor(dframe3$year),death=dframe3$death,expdeath=dframe3$expdeath)
levels(dframe$period) <- c("1","2","3","4","5")

#Import Datasets with Map Polygons
dframe.poly2 <- read.csv("data/JBC/japan_poly2.csv")
japan.poly2 <- dframe.poly2[,2:3]
dframe.prefect2 <- read.csv("data/JBC/japan_prefect2.csv")
japan.prefect2 <- dframe.prefect2[,2:5]

#
#Set Some Initial Conditions
x=dframe2$utmx/1000
y=dframe2$utmy/1000
rMax=30
#Time=5
nsim=100
mods <- c("QuasiPoisson", "Poisson")
threshold <- c(0.9, 0.5)

tabn <- NULL

#Load each of the .RData files and recalculate

#files <- list.files(path="SimulationOutput/simoutputALL", pattern="*.RData", full.names=T, recursive=FALSE)
files <- list.files(path="SimulationOutput/", pattern="*.RData", full.names=T, recursive=FALSE)

table.detection <- lapply(1:length(files), function(x){
    load(files[[x]])
    fn <- files[[x]]
    print(fn)
    fn <- unlist(strsplit(fn, split="_"))
    
    rad <- as.numeric(str_extract(fn[[3]], "([^s]*)$"))
    risk <- as.numeric(gsub(".RData","", fn[[7]]))
    if(length(unlist(strsplit(as.character(risk), "")))>1){
        print("yes")
        risk <- as.numeric(gsub("(\\s\\d)","\\1.\\2",paste(" ",risk)))
    }
    cent <- as.integer(str_extract(fn[[3]], "\\d\\d"))
    if(cent==15){
        cent <- as.integer(150)
    }
    
    tim <- as.integer(unlist(strsplit(fn[[5]],split="")))
    
    print(c(rad, risk, cent, tim))
    
    #DETECTION
    # ##QP - Space
    set <- detect_set(res$lassoresult.qp.s, res$init.vec.s, as.matrix(res$rr.mat[,tim[1]]), Time=1, x, y, rMax, cent, rad)
    incluster.qp.s <- detect.incluster(res$lassoresult.qp.s, res$init.vec.s, as.matrix(res$rr.mat[,tim[1]]), set, 1, 1, nsim,
                                       x, y, rMax, cent,
                                       rad, IC = "ic")

    print("OKOKOK")
    detect.qp.s <- list(clust.diagnostics(incluster.qp.s, threshold[1]), clust.diagnostics(incluster.qp.s , threshold[2]))
    detect.out.qp.s <- (matrix(unlist(detect.qp.s),ncol=3, byrow=TRUE,
                               dimnames = list(c(
                                   paste0("incluster.any.", threshold[1]),
                                   paste0("outcluster.", threshold[1]),
                                   paste0("alldetect.",threshold[1]),
                                   paste0("potentialclusterdetect.",threshold[1]),
                                   paste0("trueclusterdetect.",threshold[1]),
                                   paste0("alldetect.summary.mean.", threshold[1]),
                                   paste0("alldetect.summary.median.", threshold[1]),
                                   paste0("alldetect.summary.sd", threshold[1]),
                                   paste0("potentialdetect.summary.mean.", threshold[1]),
                                   paste0("potentialdetect.summary.median.", threshold[1]),
                                   paste0("potentialdetect.summary.sd.", threshold[1]),
                                   paste0("truedetect.summary.mean.", threshold[1]),
                                   paste0("truedetect.summary.median.", threshold[1]),
                                   paste0("truedetect.summary.sd.", threshold[1]),

                                   paste0("incluster.any.", threshold[2]),
                                   paste0("outcluster.", threshold[2]),
                                   paste0("alldetect.",threshold[2]),
                                   paste0("potentialclusterdetect.",threshold[2]),
                                   paste0("trueclusterdetect.",threshold[2]),
                                   paste0("alldetect.summary.mean.", threshold[2]),
                                   paste0("alldetect.summary.median.", threshold[2]),
                                   paste0("alldetect.summary.sd", threshold[2]),
                                   paste0("potentialdetect.summary.mean.", threshold[2]),
                                   paste0("potentialdetect.summary.median.", threshold[2]),
                                   paste0("potentialdetect.summary.sd.", threshold[2]),
                                   paste0("truedetect.summary.mean.", threshold[2]),
                                   paste0("truedetect.summary.median.", threshold[2]),
                                   paste0("truedetect.summary.sd.", threshold[2])),
                                   c("aic","aicc","bic"))))
    ##P - Space
    set <- detect_set(res$lassoresult.p.s, res$init.vec.s, as.matrix(res$rr.mat[,tim[[1]]]), Time=1, x, y, rMax, cent, rad)
    incluster.p.s <- detect.incluster(res$lassoresult.p.s, res$init.vec.s, as.matrix(res$rr.mat[,tim[[1]]]), set, 1, 1, nsim, x, y, rMax, cent,
                                      rad, IC = "ic")
    detect.p.s <- list(clust.diagnostics(incluster.p.s, threshold[1]), clust.diagnostics(incluster.p.s, threshold[2]))
    detect.out.p.s <- (matrix(unlist(detect.p.s),ncol=3, byrow=TRUE,
                              dimnames = list(c(
                                  paste0("incluster.any.", threshold[1]),
                                  paste0("outcluster.", threshold[1]),
                                  paste0("alldetect.",threshold[1]),
                                  paste0("potentialclusterdetect.",threshold[1]),
                                  paste0("trueclusterdetect.",threshold[1]),
                                  paste0("alldetect.summary.mean.", threshold[1]),
                                  paste0("alldetect.summary.median.", threshold[1]),
                                  paste0("alldetect.summary.sd", threshold[1]),
                                  paste0("potentialdetect.summary.mean.", threshold[1]),
                                  paste0("potentialdetect.summary.median.", threshold[1]),
                                  paste0("potentialdetect.summary.sd.", threshold[1]),
                                  paste0("truedetect.summary.mean.", threshold[1]),
                                  paste0("truedetect.summary.median.", threshold[1]),
                                  paste0("truedetect.summary.sd.", threshold[1]),

                                  paste0("incluster.any.", threshold[2]),
                                  paste0("outcluster.", threshold[2]),
                                  paste0("alldetect.",threshold[2]),
                                  paste0("potentialclusterdetect.",threshold[2]),
                                  paste0("trueclusterdetect.",threshold[2]),
                                  paste0("alldetect.summary.mean.", threshold[2]),
                                  paste0("alldetect.summary.median.", threshold[2]),
                                  paste0("alldetect.summary.sd", threshold[2]),
                                  paste0("potentialdetect.summary.mean.", threshold[2]),
                                  paste0("potentialdetect.summary.median.", threshold[2]),
                                  paste0("potentialdetect.summary.sd.", threshold[2]),
                                  paste0("truedetect.summary.mean.", threshold[2]),
                                  paste0("truedetect.summary.median.", threshold[2]),
                                  paste0("truedetect.summary.sd.", threshold[2])),
                                  c("aic","aicc","bic"))))
    
    ##QP - SPACETIME
    set <- detect_set(res$lassoresult.qp.st, res$init.vec, res$rr.mat, Time, x, y, rMax, cent, rad)
    incluster.qp.st <- detect.incluster(res$lassoresult.qp.st, res$init.vec, res$rr.mat, set, tim, Time, nsim, x, y, rMax, cent, 
                                        rad, IC = "ic")
    detect.qp.st <- list(clust.diagnostics(incluster.qp.st , threshold[1]), clust.diagnostics(incluster.qp.st , threshold[2]))
    detect.out.qp.st <- (matrix(unlist(detect.qp.st),ncol=3, byrow=TRUE, 
                                dimnames = list(c(
                                    paste0("incluster.any.", threshold[1]),
                                    paste0("outcluster.", threshold[1]),
                                    paste0("alldetect.",threshold[1]), 
                                    paste0("potentialclusterdetect.",threshold[1]), 
                                    paste0("trueclusterdetect.",threshold[1]),
                                    paste0("alldetect.summary.mean.", threshold[1]),
                                    paste0("alldetect.summary.median.", threshold[1]),
                                    paste0("alldetect.summary.sd", threshold[1]),
                                    paste0("potentialdetect.summary.mean.", threshold[1]),
                                    paste0("potentialdetect.summary.median.", threshold[1]),
                                    paste0("potentialdetect.summary.sd.", threshold[1]),
                                    paste0("truedetect.summary.mean.", threshold[1]),
                                    paste0("truedetect.summary.median.", threshold[1]),
                                    paste0("truedetect.summary.sd.", threshold[1]),
                                    
                                    paste0("incluster.any.", threshold[2]),
                                    paste0("outcluster.", threshold[2]),
                                    paste0("alldetect.",threshold[2]), 
                                    paste0("potentialclusterdetect.",threshold[2]), 
                                    paste0("trueclusterdetect.",threshold[2]),
                                    paste0("alldetect.summary.mean.", threshold[2]),
                                    paste0("alldetect.summary.median.", threshold[2]),
                                    paste0("alldetect.summary.sd", threshold[2]),
                                    paste0("potentialdetect.summary.mean.", threshold[2]),
                                    paste0("potentialdetect.summary.median.", threshold[2]),
                                    paste0("potentialdetect.summary.sd.", threshold[2]),
                                    paste0("truedetect.summary.mean.", threshold[2]),
                                    paste0("truedetect.summary.median.", threshold[2]),
                                    paste0("truedetect.summary.sd.", threshold[2])),
                                    c("aic","aicc","bic"))))
    
    ##P - SPACETIME
    set <- detect_set(res$lassoresult.p.st, res$vectors.sim, res$rr.mat, Time, x, y, rMax, cent, rad)
    incluster.p.st <- detect.incluster(res$lassoresult.p.st, res$vectors.sim, res$rr.mat, set, tim, Time, nsim, x, y, rMax, cent, 
                                       rad, IC = "ic")
    detect.p.st <- list(clust.diagnostics(incluster.p.st, threshold[1]), clust.diagnostics(incluster.p.st, threshold[2]))
    detect.out.p.st <- (matrix(unlist(detect.p.st),ncol=3, byrow=TRUE, 
                               dimnames = list(c(
                                   paste0("incluster.any.", threshold[1]),
                                   paste0("outcluster.", threshold[1]),
                                   paste0("alldetect.",threshold[1]), 
                                   paste0("potentialclusterdetect.",threshold[1]), 
                                   paste0("trueclusterdetect.",threshold[1]),
                                   paste0("alldetect.summary.mean.", threshold[1]),
                                   paste0("alldetect.summary.median.", threshold[1]),
                                   paste0("alldetect.summary.sd", threshold[1]),
                                   paste0("potentialdetect.summary.mean.", threshold[1]),
                                   paste0("potentialdetect.summary.median.", threshold[1]),
                                   paste0("potentialdetect.summary.sd.", threshold[1]),
                                   paste0("truedetect.summary.mean.", threshold[1]),
                                   paste0("truedetect.summary.median.", threshold[1]),
                                   paste0("truedetect.summary.sd.", threshold[1]),
                                   
                                   paste0("incluster.any.", threshold[2]),
                                   paste0("outcluster.", threshold[2]),
                                   paste0("alldetect.",threshold[2]), 
                                   paste0("potentialclusterdetect.",threshold[2]), 
                                   paste0("trueclusterdetect.",threshold[2]),
                                   paste0("alldetect.summary.mean.", threshold[2]),
                                   paste0("alldetect.summary.median.", threshold[2]),
                                   paste0("alldetect.summary.sd", threshold[2]),
                                   paste0("potentialdetect.summary.mean.", threshold[2]),
                                   paste0("potentialdetect.summary.median.", threshold[2]),
                                   paste0("potentialdetect.summary.sd.", threshold[2]),
                                   paste0("truedetect.summary.mean.", threshold[2]),
                                   paste0("truedetect.summary.median.", threshold[2]),
                                   paste0("truedetect.summary.sd.", threshold[2])),
                                   c("aic","aicc","bic"))))
    
    tabn <- rbind("******",cbind(rad,risk,cent,time=as.numeric(paste(tim, collapse = "")), mod = "ST",
                                  rbind("QuasiPois",detect.out.qp.st), rbind("Pois",detect.out.p.st)), 
                   cbind(rad,risk,cent,time=as.numeric(paste(tim, collapse = "")), mod = "Space",
                         rbind("QuasiPois",detect.out.qp.s), rbind("Pois",detect.out.p.s)))
    tabn <- rbind(tabn)
    return(tabn)
})


    #load("SimulationOutput/sim_center_150radius11_start_234_rr_2.RData")
#fn <- "SimulationOutput/sim_center_150radius11_start_234_rr_2.RData"
#fn <- "SimulationOutput/simoutputALL/sim_center_150radius11_start_12_rr_11.RData"

#WRITE TO CSV
print(table.detection)
write.csv(table.detection, file="tabledetection1.csv", row.names=TRUE)
save(table.detection, file="tabledetection.RData")

mytable.df <- do.call("rbind", lapply(table.detection, as.data.frame))
write.csv(mytable.df, file="tabledetectiondf.csv", row.names=TRUE)




