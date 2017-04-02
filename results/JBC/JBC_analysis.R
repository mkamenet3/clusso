########################################################################################################
########################################################################################################
########################################################################################################
#Space-time Analysis of Japanese Breast Cancer Data using cluST.R
#Maria Kamenetsky
#9-26-16
#updated: 3-27-17

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
#SPACETIME ONLY
########################################################################################################

####################################################
#QUASI-POISSON ONLY
####################################################
#Initial inputs
x=dframe2$utmx/1000
y=dframe2$utmy/1000
rMax <- 30 
Time=5


res.st.qp <- clust(x,y,rMax, dframe$period, dframe$expdeath, dframe$death, Time, spacetime=TRUE, pois=FALSE, utm=TRUE, byrow=TRUE)

#save results
filename <- paste0("results/JBC/jbc_QP_ST",".RData")
save(res.st.qp, file = filename)

#Create Empty PDF to Map Onto
pdfname <- paste0("figures/JBC/jbc_QP_ST",".pdf")
easyplot(pdfname, res.st.qp, mods, space="spacetime")

####################################################
#POISSON ONLY
####################################################
res.st.p <- clust(x,y,rMax, dframe$period, dframe$expdeath, dframe$death, Time, spacetime=TRUE, pois=TRUE, utm=TRUE, byrow=TRUE)

#save results
filename <- paste0("results/JBC/jbc_P_ST",".RData")
save(res.st.p, file = filename)

#Create Empty PDF to Map Onto
pdfname <- paste0("figures/JBC/jbc_P_ST",".pdf")
easyplot(pdfname, res.st.p, mods, space="spacetime")


########################################################################################################
########################################################################################################
########################################################################################################

########################################################################################################
#SPACE ONLY
########################################################################################################


####################################################
#QUASI-POISSON
####################################################
#Initial inputs
x=dframe2$utmx/1000
y=dframe2$utmy/1000
rMax <- 30 
Time=1


#create average dataframe for spaceonly
death <- with(dframe, tapply(death, id, function(x) round(mean(x))))
expdeath <- with(dframe, tapply(expdeath, id, function(x) mean(x)))
df <- cbind.data.frame(id = unique(dframe$id), period = rep("1", length(unique(dframe$id))), death = death, expdeath=expdeath)


res.s.qp <- clust(x,y,rMax, df$period, df$expdeath, df$death, Time, spacetime=FALSE, pois=FALSE, utm=TRUE, byrow=TRUE)

#save results
filename <- paste0("results/JBC/jbc_QP_Space",".RData")
save(res.s.qp, file = filename)

#Create Empty PDF to Map Onto
pdfname <- paste0("figures/JBC/jbc_QP_Space",".pdf")
easyplot(pdfname, res.s.qp, mods, space="space")


####################################################
#POISSON
####################################################

res.s.p <- clust(x,y,rMax, df$period, df$expdeath, df$death, Time, spacetime=FALSE, pois=TRUE, utm=TRUE, byrow=TRUE)

#save results
filename <- paste0("results/JBC/jbc_P_Space",".RData")
save(res.s.qp, file = filename)

#Create Empty PDF to Map Onto
pdfname <- paste0("figures/JBC/jbc_P_Space",".pdf")
easyplot(pdfname, res.s.p, mods, space="space")





