########################################################################################################
########################################################################################################
########################################################################################################
#Space-time Analysis of Japanese Breast Cancer Data using cluST.R
#Maria Kamenetsky
#9-26-16
#updated: 6-27-17

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


mods <- c("QuasiPoisson", "Poisson")

########################################################################################################
########################################################################################################
########################################################################################################


########################################################################################################
#REAL DATA - Japanese Breast Cancer
########################################################################################################
#Initial inputs
x=dframe2$utmx/1000
y=dframe2$utmy/1000
rMax <- 20 
Time=5


#res.st.qp <- clust(x,y,rMax, dframe$period, dframe$expdeath, dframe$death, Time, spacetime=TRUE, pois=FALSE, utm=TRUE, byrow=TRUE)
system.time(res <- clust(x,y,rMax,dframe$period, dframe$expdeath, dframe$death, Time,
                                utm=TRUE, byrow=TRUE,space= "both"))

#save results
filename <- paste0("results/JBC/jbc_analysis",".RData")
save(res, file = filename)

#Create Empty PDF to Map Onto
pdfname <- paste0("figures/JBC/jbc_analysis",".pdf")

easyplot(pdfname, res, mods, space="both", obs=TRUE)


########################################################################################################
#Table of Number of clusters detected
########################################################################################################

model <- c(rep("Poisson",2), rep("Quasi-Poisson",2))
st <- rep(c("Space", "Space-Time"),2)
numclust.AIC <- c(res$lassoresult.p.s$numclust.qaic, res$lassoresult.p.st$numclust.qaic, 
                  res$lassoresult.qp.s$numclust.qaic, res$lassoresult.qp.st$numclust.qaic)
numclust.AICc <- c(res$lassoresult.p.s$numclust.qaicc, res$lassoresult.p.st$numclust.qaicc, 
                   res$lassoresult.qp.s$numclust.qaicc, res$lassoresult.qp.st$numclust.qaicc)
numclust.BIC <- c(res$lassoresult.p.s$numclust.qbic, res$lassoresult.p.st$numclust.qbic, 
                  res$lassoresult.qp.s$numclust.qbic, res$lassoresult.qp.st$numclust.qbic)

(table.clusters <- cbind(model, st, numclust.AIC, numclust.AICc, numclust.BIC))


#WRITE TO CSV
print(table.clusters)
write.csv(table.clusters, file="tables/tableclusters.csv", row.names=TRUE)

