#' cluST.R
#' @name cluST
#' @title Spatial and Spatial-Temporal Clustering using Lasso
#' @description This package will host the functions and code for a spatial and spatial-temporal clustering model using Lasso regression. We assume a Poisson distribution for counts. The model includes a spatially-structured random effect to account for heterogeneity in a Poisson model. Model selection is based on QAIC/QAICc/QBIC.
#' @details This package contains C++ code which is compiled at run-time. These functions are located in the "src" file.
#' Dependencies and inputs are documented here.
#' Data used for this analysis comes from MidWest Poverty data.
#' @references Xu, Jiale, "Stagewise and Stepwise Methods for Space and Space-Time Cluster Detection"
#' 
#' 
#' 
#' 
#' Unlist2
#' 
#' This function is from https://github.com/Bioconductor-mirror/AnnotationDbi/blob/master/R/unlist2.R
#' It unlists the names without issues
#' 
make.name.tree <- function(x, recursive, what.names)
{
    if (!is.character(what.names) || length(what.names) != 1)
        stop("'what.names' must be a single string")
    what.names <- match.arg(what.names, c("inherited" , "full"))
    .make.name.tree.rec <- function(x, parent_name, depth)
    {
        if (length(x) == 0)
            return(character(0))
        x_names <- names(x)
        if (is.null(x_names))
            x_names <- rep.int(parent_name, length(x))
        else if (what.names == "full")
            x_names <- paste0(parent_name, x_names)
        else
            x_names[x_names == ""] <- parent_name
        if (!is.list(x) || (!recursive && depth >= 1L))
            return(x_names)
        if (what.names == "full")
            x_names <- paste0(x_names, ".")
        lapply(seq_len(length(x)),
               function(i) .make.name.tree.rec(x[[i]], x_names[i], depth + 1L))
    }
    .make.name.tree.rec(x, "", 0L)
}

unlist2 <- function(x, recursive=TRUE, use.names=TRUE, what.names="inherited")
{
    ans <- unlist(x, recursive, FALSE)
    if (!use.names)
        return(ans)
    if (!is.character(what.names) || length(what.names) != 1)
        stop("'what.names' must be a single string")
    what.names <- match.arg(what.names, c("inherited" , "full"))
    names(ans) <- unlist(make.name.tree(x, recursive, what.names), recursive, FALSE)
    ans
}

#' Create the clusters dataframe
#' 
#' @param vector xP x coordinates (easting/latitude); if utm coordinates, scale to km.
#' @param vector yP y coordinates (northing/longitude); if utm coordinates, scale to km.
#' @param r.max set max radius (in km)
#' @param utm TRUE/FALSE as to whether or not the x and y coordinates are in UTM (TRUE) or LAT/LONG(FALSE)
#' @param n Number of coordinate pairs/number of centers
#' @return This function returns a dataframe that contains 
#' @export
#' @examples
#' clustersDF(x1,y1,rMax, utm=TRUE, length(x1))
#' clustersDF(lat, long, utm=FALSE, length(lat))

clustersDF <- function(xP,yP, r.max, utm=FALSE,n){
    indR = (1:n)[!duplicated(cbind(xP,yP))] 
    if(utm==FALSE){
        tmpR <- (as.matrix(distm(cbind(xP, yP), fun=distHaversine))[indR,])/1000    
    } 
    else{
        tmpR = as.matrix(dist(cbind(xP,yP)))[indR,]
    }
    lastR = apply(tmpR, 1, function(x,r) order(x)[1:sum(x<=r)],r=r.max)
    ncR = unlist2(lapply(lastR, length))
    lastR = unlist2(lastR)
    rR=unlist2(apply(tmpR,1, function(x,r) { sort(x[x<=r]) },r=r.max))
    
    clustersR=data.frame(center=rep(indR,ncR),
                         x=xP[rep(indR,ncR)],y=yP[rep(indR,ncR)],
                         r=rR, 
                         n=unlist(lapply(ncR,seq)),
                         last=lastR)    
    return(clustersR)
}

#' Poisson Distribution Function
#' This is the main distriution function for our model. This assumes we have a Poisson fixed effect and Gamma random effect. In order to deal with constraints from the Lasso function, we use the Poisson distirbution function here and account for overdispersion in the QIC.
#' @param y observed values
#' @param lambda vector of expected outcomes * exp(each column of each potential path)
#' @param log whether or not the log-likelihood should be returned or the likelihood. Default is to be TRUE
#' @return returns a matrix of ________

dpoisson <- function(y, lambda, log = FALSE) {
    if(log == FALSE) 
        return(lambda^y * exp(-lambda)/factorial(y))
    else
        return(y*ifelse(lambda==0,1,log(lambda))-lambda)
}

###DO I need to document cpp functions?

#' sparseMe
#' This function creates a sparse matrix of 1's of all potential clusters for the Lasso algorithm to cycle over; this incorporates space
#' @param clusters clusters dataframe from (clustersDF function) that includes the center, x,y, r (radius), n (counter), and last (last observation in potential cluster)
#' @param numCenters the number of centers
#' @return returns sparse matrix of 1's
sparseMe <- function(clusters, numCenters){
    potClus <- numCenters
    mymat <- NULL
    for(i in 1:nrow(clusters)){
        myvec <- list(as(max_colCpp(numCenters, i, clusters$n, clusters$last), "sparseVector")) 
        mymat <- c(mymat,myvec) 
    }
    xx <- NULL
    jj <- NULL
    ii <- NULL
    
    for(k in 1:length(mymat)){
        xx<- c(xx, mymat[[k]]@x)
        jj <- c(jj, mymat[[k]]@i)
        ii <- c(ii, rep(k,length(mymat[[k]]@x)))
    }
    return(t(sparseMatrix(i = ii, j = jj, x =xx, dims = c(length(mymat), numCenters))))
}

###############################################################################
mymat <- spacesparse

ii = mymat@i +1
jj = findInterval(seq(mymat@x)-1, mymat@p[-1])+1
xx = mymat@x


matrixBlocks <- function(mymat, list, Time) {    
    block <- Matrix(data=0, nrow=1040, ncol=1, sparse=TRUE)  
    #for(m in 1:Time){
     #   newlist <- seq(1:(Time-(m-1)))
      #  print(newlist)
        for (k in newlist){
            blockNew <- sparseMatrix(i = (ii+((k-1)*208)), j = jj, x =xx, dims = c(1040, 8960))    
            block <- cBind(block, blockNew)
       #     print(m)
            print(k)
        } 
    }
    return(block[,-1])
}
test <- matrixBlocks(mymat, list, Time) 




###############################################################################
###############################################################################
#7-19 - attempt to create sparseme object st rbinds occur in sparsemet

sparseMe2 <- function(clusters, numCenters){
    potClus <- numCenters
    mymat <- NULL
    for(i in 1:nrow(clusters)){
        myvec <- list(as(max_colCpp(numCenters, i, clusters$n, clusters$last), "sparseVector")) 
        mymat <- c(mymat,myvec) 
    }
    xx <- NULL
    jj <- NULL
    ii <- NULL
    
    for(k in 1:length(mymat)){
        xx<- c(xx, mymat[[k]]@x)
        jj <- c(jj, mymat[[k]]@i)
        ii <- c(ii, rep(k,length(mymat[[k]]@x)))
    }
    return(cBind(t(sparseMatrix(i = ii, j = jj, x =xx, dims = c(length(mymat), numCenters))),
           t(sparseMatrix(i = ii, j = jj, x =xx, dims = c(length(mymat), numCenters)))))
}

matrixBlocks2 <- function(mymat, list, Time) {    
    block <- Matrix(data=0, nrow=(208*2), ncol=1, sparse=TRUE)  
    for(m in 1:Time){
       newlist <- seq(1:(Time-(m-1)))
      print(newlist)
    for (k in newlist){
        blockNew <- sparseMatrix(i = (ii+((k-1)*416)), j = jj, x =xx, dims = c(1040, 8960))    
        block <- cBind(block, blockNew)
        #     print(m)
        print(k)
    } 
}
return(block[,-1])
}


sparsetest2 <- sparseMe2(clusters, numCenters)
list <- c(1,2,3,4)
test <- matrixBlocks2(sparsetest2, list, 4)



###############################################################################

#' sparseTime
#' This function replicates the sparse matrix of potential clusters to account for time
#' @param clusters clusters dataframe from (clustersDF function) that includes the center, x,y, r (radius), n (counter), and last (last observation in potential cluster)
#' @param numCenters the number of centers
#' @param Time number of time periods in the space-time model

sparseTime <- function(clusters, numCenters, Time){
    mysparse <- sparseMe(clusters, numCenters, intercept=FALSE)
    
    myspacetime <- Matrix(data=0, nrow=1040, ncol=134400, sparse=TRUE)
    #myspacetime <- NULL
    #myspacetime <- Matrix(0, nrow=1040, sparse=TRUE)
    
    #myspacetime <- Matrix(data=NA, nrow=nrow(mysparse), ncol= ncol(mysparse)*Time)
    
truthjj <- findInterval(seq(spacetimemat@x)-1, spacetimemat@p[-1])+1
    
    ###store values
    ii = mymat@i +1
    jj = findInterval(seq(mymat@x)-1, mymat@p[-1])+1
    xx = mymat@x
    
    
   list <- c(1,2,3,4,5)

for(k in list){
    newlist<- list
}

###tester for correction
#Now, need to figure out how to replicate my sparse matrix n times going down by column and then have that cycle through  
iitruth <- spacetimemat[,44801:(44800+(8960*4))]@i +1
jjtruth <- findInterval(seq(spacetimemat[,44801:(44800+(8960*4))]@x)-1, spacetimemat[,44801:(44800+(8960*4))]@p[-1])+1
xxtruth <- spacetimemat[,44801:(44800+(8960*4))]@x

iiT = test@i +1
jjT = findInterval(seq(test@x)-1, test@p[-1])+1


as.vector(unlist(table(jjtruth)))
take <- as.vector(unlist(table(jjtruth)))/2
#mytest <- c(ii[1], ii[1]+208, ii[2:3], ii[2:3]+208, ii[4:6], ii[4:6]+208, ii[7:10], ii[7:10]+208) #this replicates iitruth
mytest2 <- c(ii[1], ii[1]+208, ii[2:(2+take[2]-1)], ii[2:(2+take[2]-1)]+208, ii[4:(4+take[3]-1)], ii[4:(4+take[3]-1)]+208, ii[7:(7+take[4]-1)], ii[7:(7+take[4]-1)]+208)

3+1
4+3
5+6

starter <- c(ii[1], ii[1]+208)
for(k in 2:10){
    #add <- c(ii[i:(i+take[i]-1)], ii[i:(i+take[i]-1)]+208)
    #bit <- c(ii[(k+take[k-1]+1):(k+take[i]-1)], ii[(k+take[k-1]+1):(k+take[k]-1)]+208)
    #starter <- c(starter, bit)
    print(k)
    print((k+take[k-1]+1))
}


str(test)

matrixBlocks <- function(mymat, list, Time) {    
    block <- Matrix(data=0, nrow=1040, ncol=1, sparse=TRUE)  
    for(m in 1:Time){
        newlist <- seq(1:(Time-(m-1)))
        #print(newlist)
        for (k in newlist){
            #blockNew <- sparseMatrix(i = (ii+((k-1)*208)), j = jj, x =xx, dims = c(208, 8960))    
            #blockRep <- sparseMatrix(i = (ii+((k)*208)), j = jj, x =xx, dims = c(208, 8960))    
            #test <- sparseMatrix(i = (ii+((k-m)*208)), j = jj, x =xx, dims = c(208, 8960))    
            #test <- sparseMatrix(i = (ii+(k-max(newlist))*208)), j = jj, x = (xx*(k)), dims = c(208, 8960)) 
            #test <- sparseMatrix(i = (ii+((m-1)*208)), j = jj, x = rep(xx,m), dims = c(1040, 8960)) #this works
            
            
            test <- sparseMatrix(i = (ii+((m-1)*208)), j = rep(jj, each=m), x = rep(xx,m), dims = c(1040, 8960)) #this works 7/5
            
            #interBlock <- rBind(blockNew, blockRep)
            #block <- cBind(block, test)
            print(newlist)
            print(k)
            print(m) #m is the number of column replicates I want of sparse matrix
            print(max(newlist))
        } 
    }
    return(block[,-1])
    #return(block)
}
metest <- matrixBlocks(mymat, list, Time) 

#truth
all.equal(rBind(mymat,mymat), spacetimemat[1:416, 8961:(8960*2)])

###this works
matrixBlocks <- function(mymat, list, Time) {    
    block <- Matrix(data=0, nrow=1040, ncol=1, sparse=TRUE)  
    for(m in 1:Time){
        newlist <- seq(1:(Time-(m-1)))
        print(newlist)
        for (k in newlist){
            blockNew <- sparseMatrix(i = (ii+((k-1)*208)), j = jj, x =xx, dims = c(1040, 8960))    
            block <- cBind(block, blockNew)
            print(m)
            print(k)
        } 
    }
    return(block[,-1])
}
test <- matrixBlocks(mymat, list, Time) 



    #cols: 
#     (8960*0)+1:(8960*1); 
#     (8960*1)+1:(8960*2);
#     (8960*2)+1: (8960*3)
#     (8960*3)+1: (8960*4)
#     (8960*4)+1: (8960*5)
    
    
    
##############################################
    test <- sparseMatrix(i = ii, j = jj, x =xx, dims = c(1040, 8960))
    test2 <- sparseMatrix(i = (ii+208), j =jj, x =xx ,dims =c(1040, 8960))
    test3 <- sparseMatrix(i = (ii+2*(208)), j =jj, x =xx ,dims =c(1040, 8960))
    test4 <- sparseMatrix(i = (ii+3*(208)), j =jj, x =xx ,dims =c(1040, 8960))
    test5 <- sparseMatrix(i = (ii+4*(208)), j =jj, x =xx ,dims =c(1040, 8960))
    
    mastertest<- cBind(test,test2,test3,test4,test5)
    
    myspacetime[1:dim(mastertest)[1], 1:dim(mastertest)[2]] <- mastertest
    #this works
############################################



    moop <- Matrix(data=0, nrow=20, ncol=10, sparse=TRUE)
    doop <- matrix(rnorm(n = 200,mean = 100, sd=50), nrow=20,ncol=10)
    dim(doop)
    moop[1:dim(doop)[1],1:dim(doop)[2]] <- doop
    
    
    
#     ii = head(ii)
#     jj = head(jj)
#     xx = head(xx)
    
    ###replace
    
    
    
    
    myspacetime@i = mymat@i
    #myspacetime@p = mymat@p
    myspacetime@x = mymat@x
    
    
    
    for(i in 1:Time){
        cBind(mysparse, Matrix(0, nrow=length(unique(clusters$center)), ncol=nrow(clusters)*(Time-1), sparse=TRUE))
        
    }
    
    
    
    
    
    
    
    
    #block1
    p1 <- cBind(mysparse, Matrix(0, nrow=length(unique(clusters$center)), ncol=nrow(clusters)*(Time-1), sparse=TRUE))
    p2 <- cBind(Matrix(0, nrow=length(unique(clusters$center)), ncol=nrow(clusters)*(Time-4), sparse=TRUE), mysparse, 
                Matrix(0, nrow=length(unique(clusters$center)), ncol=nrow(clusters)*(Time-2), sparse=TRUE))
    p3 <- cBind(Matrix(0, nrow=length(unique(clusters$center)), ncol=nrow(clusters)*(Time-3), sparse=TRUE), mysparse,
                Matrix(0, nrow=length(unique(clusters$center)), ncol=nrow(clusters)*(Time-3), sparse=TRUE))
    p4 <- cBind(Matrix(0, nrow=length(unique(clusters$center)), ncol=nrow(clusters)*(Time-2), sparse=TRUE), mysparse,
                Matrix(0, nrow=length(unique(clusters$center)), ncol=nrow(clusters)*(Time-4), sparse=TRUE))
    p5 <- cBind(Matrix(0, nrow=length(unique(clusters$center)), ncol=nrow(clusters)*(Time-1), sparse=TRUE), mysparse)
    block1 <- rBind(p1,p2,p3,p4,p5)
    #block2
    p1 <- cBind(mysparse, Matrix(0, nrow=length(unique(clusters$center)), ncol=nrow(clusters)*(Time-2), sparse=TRUE))
    p2 <- cBind(mysparse, mysparse, Matrix(0, nrow=length(unique(clusters$center)), ncol=nrow(clusters)*(Time-3), sparse=TRUE))
    p3 <- cBind(Matrix(0, nrow=length(unique(clusters$center)), ncol=nrow(clusters)*(Time-4), sparse=TRUE), mysparse, mysparse,
                Matrix(0, nrow=length(unique(clusters$center)), ncol=nrow(clusters)*(Time-4), sparse=TRUE))
    p4 <- cBind(Matrix(0, nrow=length(unique(clusters$center)), ncol=nrow(clusters)*(Time-3), sparse=TRUE), mysparse, mysparse)
    p5 <- cBind(Matrix(0, nrow=length(unique(clusters$center)), ncol=nrow(clusters)*(Time-2), sparse=TRUE), mysparse)
    block2 <- rBind(p1,p2, p3,p4,p5)
    #block3
    p1 <- cBind(mysparse, Matrix(0, nrow=length(unique(clusters$center)), ncol=nrow(clusters)*(Time-3), sparse=TRUE))
    p2 <- cBind(mysparse, mysparse, Matrix(0, nrow=length(unique(clusters$center)), ncol=nrow(clusters)*(Time-4), sparse=TRUE))
    p3 <- cBind(mysparse,mysparse,mysparse)
    p4 <- cBind(Matrix(0, nrow=length(unique(clusters$center)), ncol=nrow(clusters)*(Time-4), sparse=TRUE), mysparse,mysparse)
    p5 <- cBind(Matrix(0, nrow=length(unique(clusters$center)), ncol=nrow(clusters)*(Time-3), sparse=TRUE),mysparse)
    block3 <- rBind(p1,p2,p3,p4,p5)
    #block4
    p1 <- cBind(mysparse,Matrix(0, nrow=length(unique(clusters$center)), ncol=nrow(clusters)*(Time-4), sparse=TRUE))
    p2 <- cBind(mysparse, mysparse)
    p3 <- cBind(mysparse,mysparse)
    p4 <- cBind(mysparse, mysparse)
    p5 <- cBind(Matrix(0, nrow=length(unique(clusters$center)), ncol=nrow(clusters)*(Time-4), sparse=TRUE),mysparse)
    block4 <- rBind(p1,p2,p3,p4,p5)
    #block5
    p1 <- cBind(mysparse)
    p2 <- cBind(mysparse)
    p3 <- cBind(mysparse)
    p4 <- cBind(mysparse)
    p5 <- cBind(mysparse)
    block5 <- rBind(p1,p2,p3,p4,p5)
    
    sparseMAT <- cBind(block1,block2, block3,block4, block5)
    return(sparseMAT)
}


