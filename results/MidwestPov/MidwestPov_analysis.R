########################################################################################################
########################################################################################################
########################################################################################################
#Space-time Analysis of Midwest Poverty Data using cluST.R
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
#LOAD MPD Data and Set Up
####################################################

#Load Data
df <- read.csv("../../data/MidwestPov/povertydataNew.csv")
pov <- read.csv("../../data/MidwestPov//upMidWestpov_Iowa_cluster_names60_06_final_wide.csv")

#Clean the Dataset 
df$year <- factor(df$year, levels=c("1960","1970","1980","1990","2000"))
df$county <- as.factor(sub(".*,","", df$statecounty))
df$state <- as.factor(gsub("\\,.*","",df$statecounty))

#Set Coordinates Dataframe
coords <- pov[c("x","y","county","state")]

#Set Some Initial Conditions
x=coords$x
y=coords$y
rMax <- round((max(distm(cbind(x, y), fun=distHaversine))/10)/1000)
Time=5  

#Create Potential Clusters Dataframe
clusters <- clustersDF(x,y,rMax, utm=FALSE, length(x))

#Set initial expected and observed
pop <- as.vector(df$denom_poor)
Y.vec <- as.vector(df$poor)
period<- as.vector(df$year)

#Adjust for observed given expected counts as coming from negative binomial distribution
outinit <- glm.nb(Y.vec ~1)
out <- glm.nb(Y.vec ~ 1 + as.factor(period)  + offset(log(pop)), init.theta = outinit$theta, 
              link=log,control=glm.control(maxit=10000))


#Set initial expected to the fitted values
E0 <- out$fitted
MPDinit <- cbind.data.frame(period,E0, Y.vec)

####################################################
#RUN Model
####################################################

#set initial conditions for function (all centers can be potential origin of the cluster)
potentialClus <- max(clusters$center)
numberCenters <- max(clusters$center)

MPDresults <- spacetimeLasso(potentialClus, clusters, numberCenters, MPDinit, Time, spacetime=TRUE)


####################################################
#Set Risk Ratio Vectors Based on QIC
####################################################
rr <- setRR(MPDresults, MPDinit, Time=5)

####################################################
#Map RR to Colors
####################################################
rrcolors <- colormapping(rr, Time=5)


####################################################
#Map Colors to Maps
####################################################

#Import Map of the Midwest from 'maps' library
m <- map('county', region = c("Illinois", "Indiana", "Iowa", "Michigan", "Minnesota", "Wisconsin"))
m$names                                                

#Deal with issue where Shawano, Menominee, and Oconto counties were given the same initial values in the dataset
pcounty=as.character(df$county)
pcounty[533]<-"OCONTO"                                 
tmp=tolower(paste(df$state,',',pcounty,sep='')) ## State,County names (order in the poor dataset)

colSeq = c(1:504,533,505:519,533,520:532,533)   
tmp=tmp[colSeq] ## reorder the counties in the poor dataset

not=which(m$names!=tmp)              
m$names[not];tmp[not]    


#Create Empty PDF to Map Onto
pdf("../../figures/MidwestPov/MidwestPov_map.pdf", height=11, width=10)

#Set Plots
par(mfrow = c(4,5))

#Maps of Observed Counts
map('county', region = c("Illinois","Indiana","Iowa","Michigan","Minnesota","Wisconsin"), 
    fill=TRUE,col=rrcolors$colors.obs[,1][colSeq])
    title(main="Obs - 1960")

map('county', region = c("Illinois","Indiana","Iowa","Michigan","Minnesota","Wisconsin"), 
    fill=TRUE,col=rrcolors$colors.obs[,2][colSeq])
    title(main="Obs - 1970")

map('county', region = c("Illinois","Indiana","Iowa","Michigan","Minnesota","Wisconsin"), 
    fill=TRUE,col=rrcolors$colors.obs[,3][colSeq])
    title(main="Obs - 1980")

map('county', region = c("Illinois","Indiana","Iowa","Michigan","Minnesota","Wisconsin"), 
    fill=TRUE,col=rrcolors$colors.obs[,4][colSeq])
    title(main="Obs - 1990")

map('county', region = c("Illinois","Indiana","Iowa","Michigan","Minnesota","Wisconsin"), 
    fill=TRUE,col=rrcolors$colors.obs[,5][colSeq])
    title(main="Obs - 2000")


#Maps of AIC Path

map('county', region = c("Illinois","Indiana","Iowa","Michigan","Minnesota","Wisconsin"),
    fill=TRUE,col=rrcolors$color.qaic[,1][colSeq])
    title(main="AIC - 1960")

map('county', region = c("Illinois","Indiana","Iowa","Michigan","Minnesota","Wisconsin"), 
    fill=TRUE,col=rrcolors$color.qaic[,2][colSeq])
    title(main="AIC - 1970")

map('county', region = c("Illinois","Indiana","Iowa","Michigan","Minnesota","Wisconsin"), 
    fill=TRUE,col=rrcolors$color.qaic[,3][colSeq])
    title(main="AIC - 1980")

map('county', region = c("Illinois","Indiana","Iowa","Michigan","Minnesota","Wisconsin"),
    fill=TRUE,col=rrcolors$color.qaic[,4][colSeq])
    title(main="AIC - 1990")

map('county', region = c("Illinois","Indiana","Iowa","Michigan","Minnesota","Wisconsin"), 
    fill=TRUE,col=rrcolors$color.qaic[,5][colSeq])
    title(main="AIC - 2000")

#Maps of AICc Path

map('county', region = c("Illinois","Indiana","Iowa","Michigan","Minnesota","Wisconsin"), 
    fill=TRUE,col=rrcolors$color.qaicc[,1][colSeq])
    title(main="AICc - 1960")

map('county', region = c("Illinois","Indiana","Iowa","Michigan","Minnesota","Wisconsin"), 
    fill=TRUE,col=rrcolors$color.qaicc[,2][colSeq])
    title(main="AICc - 1970")

map('county', region = c("Illinois","Indiana","Iowa","Michigan","Minnesota","Wisconsin"), 
    fill=TRUE,col=rrcolors$color.qaicc[,3][colSeq])
    title(main="AICc - 1980")

map('county', region = c("Illinois","Indiana","Iowa","Michigan","Minnesota","Wisconsin"), 
    fill=TRUE,col=rrcolors$color.qaicc[,4][colSeq])
    title(main="AICc - 1990")

map('county', region = c("Illinois","Indiana","Iowa","Michigan","Minnesota","Wisconsin"), 
    fill=TRUE,col=rrcolors$color.qaicc[,5][colSeq])
    title(main="AICc - 2000")


#Maps of BIC Path

map('county', region = c("Illinois","Indiana","Iowa","Michigan","Minnesota","Wisconsin"), 
    fill=TRUE,col=rrcolors$color.qbic[,1][colSeq])
    title(main="BIC - 1960")

map('county', region = c("Illinois","Indiana","Iowa","Michigan","Minnesota","Wisconsin"), 
    fill=TRUE,col=rrcolors$color.qbic[,2][colSeq])
    title(main="BIC - 1970")

map('county', region = c("Illinois","Indiana","Iowa","Michigan","Minnesota","Wisconsin"), 
    fill=TRUE,col=rrcolors$color.qbic[,3][colSeq])
    title(main="BIC - 1980")

map('county', region = c("Illinois","Indiana","Iowa","Michigan","Minnesota","Wisconsin"), 
    fill=TRUE,col=rrcolors$color.qbic[,4][colSeq])
    title(main="BIC - 1990")

map('county', region = c("Illinois","Indiana","Iowa","Michigan","Minnesota","Wisconsin"), 
    fill=TRUE,col=rrcolors$color.qbic[,5][colSeq])
    title(main="BIC - 2000")


#Turn off pdf development
dev.off()







