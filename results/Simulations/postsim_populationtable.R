####################################################
#Post-Sim Population Table: 
#M.Kamenetsky
#4-24-17
####################################################

sink("postsimpopulation.txt")
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




####################################################
#Simulations Set-up
####################################################

#Set Some Initial Conditions
x=dframe2$utmx/1000
y=dframe2$utmy/1000
rMax=30
Time=5
nsim=100

centers <- c(150, 35)
radii <- c(9, 11, 18)
timeperiods <- list(c(1:2), c(2:4), c(1:5), c(3:5))

#create clusters
clusters <- clusters2df(x,y,rMax, utm=TRUE, length(x))
n <- length(x)

#create identifer
dframe$idx <- rep(seq(1:length(x)), each = 5)
which(unlist(table(dframe$idx, dframe$id))!=5 & unlist(table(dframe$idx, dframe$id)) !=0) #check

#create null
table.clusterpop <- NULL


###DO
#make one table
for(cent in centers){
    for(rad in radii){
        for(tim in timeperiods){
            print(c(cent, rad, tim))
            
            cells <- clusters[(clusters$center==cent & clusters$r <= rad ),]$last
            df <- subset(dframe, idx %in% cells & (as.integer(dframe$period) >=tim[[1]] & as.integer(dframe$period)<=tail(tim,n=1)))
            expected.agg <- aggregate(df$expdeath, list(df$idx), mean)
            observed.agg <- aggregate(df$death, list(df$idx), mean)
            #expected
            expected.mean <- mean(expected.agg$x)
            expected.sd <- sd(expected.agg$x)
            expected.sum <- sum(expected.agg$x)
            #observed
            observed.mean <- mean(observed.agg$x)
            observed.sd <- sd(observed.agg$x)
            observed.sum <- sum(observed.agg$x)
            
            tabn <- cbind(cent, time=as.numeric(paste(tim, collapse="")), rad,
                          expected.sum, expected.mean, expected.sd, observed.sum, observed.mean, observed.sd)
            
            table.clusterpop <- rbind(table.clusterpop, tabn)
        }
    }
}

#WRITE TO CSV
print(table.clusterpop)
write.csv(table.clusterpop, file="figures/tables/tableclustpopulation.csv", row.names=TRUE)


sink()





