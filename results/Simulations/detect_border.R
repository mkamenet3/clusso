
#load simulation



seg1 <- psp(dframe.prefect2$x1, dframe.prefect2$y1, dframe.prefect2$x2, dframe.prefect2$y2, window=owin(c(266, 436), c(3956,4109)))

x = c(0,1,1,0,0)
y = c(0,0,1,1,0)
bbox_matrix_sp = cbind(rep(x,13),rep(y,13))
sp_re_alle = SpatialPolygons(lapply(1:13, 
                                    function(x) Polygons(list(Polygon(bbox_matrix_sp[((x-1)*5+1):(x*5),])), paste0("reh",x))))
#SpatialPolygons(list(Polygons(list(Polygon(cbind(c(0, 0, 6, 6, 0), c(6, 0, 0, 6, 6)))), "peri")))

test = SpatialPolygons(lapply(1:nrow(dframe.prefect2),
                              function(x) Polygons(list(Polygon(dframe.prefect2[x,2:3], dframe.prefect2[x,1])))))
# 
# 
# #diagnostics
# set <- detect.set(res$lassoresult, res$init.vec, res$rr.mat, Time=5)
# 
# det <- detect(res$lassoresult, res$init.vec, res$rr.mat, period_start=1,period_end= 2, multi_period = FALSE, IC="aic", Time)
# 
# detin <- detect.incluster(res$lassoresult, res$init.vec, res$rr.mat, set,period_start=1, period_end=2, multi_period=FALSE, IC = "ic", Time=5) #this works
# detfalse <- detect.falsecluster(result$lassoresult, result$init.vec, result$rr.mat, set,period_start=1, period_end=2, multi_period=FALSE, IC = "aic", Time=5)
# detfalse <- detect.falsecluster(result$lassoresult, result$init.vec, result$rr.mat, set,period_start=1, period_end=2, multi_period=FALSE, IC = "ic", Time=5)
# 
# 
# detbk <- detect.inbackground(result$lassoresult, result$init.vec, result$rr.mat, set,period_start=1, period_end=2, multi_period=FALSE, IC = "bic", Time=5)


#TODO in background - clarify
#borders - best by plot?
#AIC, BIC should work
#Pois - space, quasi-pois - ST

###########################

set <- detect.set(result$lassoresult, result$init.vec, result$rr.mat, Time=5)
incluster <- detect.incluster.ic(result$lassoresult, result$init.vec, result$rr.mat, set, period=c(1,2),Time=5)

#############################
#borders

library(SpatialEpi)



dframe1 <- read.csv("data/JBC/jap.breast.F.9.10.11.csv")
dframe2 <- read.csv("data//JBC//utmJapan.csv")
dframe3 <- aggregate(dframe1, by=list(as.factor(rep(1:(nrow(dframe1)/4),each=4))), FUN="sum")
dframe=data.frame(id=as.factor(dframe3$id/4),period=as.factor(dframe3$year),death=dframe3$death,expdeath=dframe3$expdeath)
levels(dframe$period) <- c("1","2","3","4","5")

dframe.poly2 <- read.csv("data/JBC/japan_poly2.csv")
japan.poly2 <- dframe.poly2[,2:3]
dframe.prefect2 <- read.csv("data/JBC/japan_prefect2.csv")
japan.prefect2 <- dframe.prefect2[,2:5]


ind <- which(is.na(dframe.poly2[,2]))

names <- as.factor(seq(1,208,1))
coord.system <- '+proj=utm'
repeats <- rep(1, sum(is.na(dframe.poly2[, 2])) + 1)
poly <- as.matrix(dframe.poly2[, 2:3])
test <- polygon2spatial_polygon(poly, coord.system, names)


function (poly, coordinate.system, area.names = NULL, nrepeats = NULL) 
{
    if (missing(coordinate.system)) {
        stop("Coordinate system must be specified: '+proj=utm' or '+proj=longlat'.")
    }
    if (is.null(nrepeats)) {
        nrepeats <- rep(1, sum(is.na(poly[, 1])) + 1)
    }
    if (is.null(area.names)) {
        area.names <- as.character(1:length(nrepeats))
    }
    na.index <- which(is.na(poly[, 1]))
    n <- length(nrepeats)
    list.polygon <- NULL
    list.polygon <- list(Polygon(poly[1:(na.index[1] - 1), ], 
                                 hole = FALSE))
    for (i in 1:(length(na.index) - 1)) {
        list.polygon <- c(list.polygon, list(Polygon(poly[(na.index[i] + 
                                                               1):(na.index[i + 1] - 1), ], hole = FALSE)))
    }
    list.polygon <- c(list.polygon, list(Polygon(poly[(na.index[i + 1] + 1):length(poly[, 1]), ], hole = FALSE)))
    list.polygons <- NULL
    start <- 1
    for (i in 1:length(nrepeats)) {
        end <- start + nrepeats[i] - 1
        temp.polygon <- NULL
        for (j in start:end) {
            print(j)
            temp.polygon <- c(temp.polygon, list(list.polygon[[j]]))
        }
        list.polygons <- c(list.polygons, list(Polygons(temp.polygon, 
                                                        ID = area.names[i])))
        start <- end + 1
    }
    Spatial.Polygon <- SpatialPolygons(list.polygons, proj4string = CRS(coordinate.system))
    return(Spatial.Polygon)
}

library(sp)
library(spdep)
mypoly <- SpatialPolygons(list.polygons)
plot(mypoly)

nbtest <- poly2nb(mypoly)
plot(nbtest)

#create clusters
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


#Set parameters of your fake cluster
center <- 50
r_list <- 18
cluster_end <- 3
rr.ratio<- 2

#Create Potential Clusters Dataframe
clusters <- clusters.df(x1,y1,rMax, utm=TRUE, length(x1))
tmp <- clusters[clusters$center==center,]
cluster <- tmp[(tmp$r <= r_list),]

last <- sort(cluster$last)
myneigh <- last
for(i in last){
    myneigh <- c(myneigh, unlist(nbtest[[i]]))
}
neighs <- sort(unique(myneigh))

result <- res
res <- detect.set(result$lassoresult, result$init.vec, result$rr.mat, Time=5)
ix <- lapply(1: nsim, 
             function(j) sapply(1:Time, 
                                function(k) 
                                    which(round(matrix(res$rr.simAIC[[j]], ncol=Time)[,k],6) > round(as.numeric(attributes(res$alphaAIC[[j]][[k]])),6))))
in_cluster_background <- lapply(1:nsim, 
                     function(j) sapply(1:Time, 
                                        function(k) ix[[j]][[k]] %in% neighs))
#did it find all cells + some background?
clusterback_idx <- lapply(1:nsim,
                          function(j) unlist(sapply(1:Time,
                                                    function(k) which(in_cluster_background[[j]][[k]]==FALSE))))

clusterback <- (nsim - sum(unlist(lapply(1:nsim,
                                         function(j) any(clusterback_idx[[j]] >0)*1))))/100


#Clean attempt
ind <- which(is.na(dframe.poly2[,2]))

names <- as.factor(seq(1,208,1))
coord.system <- '+proj=utm'
nrepeats <- rep(1, sum(is.na(dframe.poly2[, 2])))
poly <- as.matrix(dframe.poly2[, 2:3])
#test <- polygon2spatial_polygon(poly, coord.system, names)


#function (poly, coordinate.system, area.names = NULL, nrepeats = NULL) 
#{
 #   if (missing(coordinate.system)) {
  #      stop("Coordinate system must be specified: '+proj=utm' or '+proj=longlat'.")
   # }
    #if (is.null(nrepeats)) {
        #nrepeats <- rep(1, sum(is.na(poly[, 1])) + 1)
#     }
#     if (is.null(area.names)) {
        area.names <- as.character(1:length(nrepeats))
    #}
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
period <- c(3:5)


set <- detect.set(res$lassoresult, res$init.vec, res$rr.mat, Time=5)
#incluster1 <- detect.incluster.aic(res$lassoresult, res$init.vec, res$rr.mat, set, period=c(3:5),Time=5, nsim,nb, x, y, rMax, center, radius)
incluster2 <- detect.incluster(res$lassoresult, res$init.vec, res$rr.mat, set, period_start = 3, period_end = 5, multi_period = TRUE,Time=5, nsim,nb, x, y, rMax, center, radius, IC = "ic")






