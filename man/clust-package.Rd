\name{clust-package}
\alias{clust-package}
\alias{clust-package}
\docType{package}
\title{
Efficient spatial and spatio-temporal cluster detection using the LASSO
}
\description{
clust is an R package which can be used to efficiently detect spatial clusters
and spatio-temporal clusters with count data on a lattice. Examples include
rates of disease in a population broken down by a geographic unit. 
}
\details{
\tabular{ll}{
Package: \tab clust\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2017-08-20\cr
License: \tab GPL (>= 2)\cr
}
The main functions in this package are \code{clust} and \code{clust_sim}, which
both perform both spatial and spatio-temporal cluster detection using the
lasso. The former works with observed and expected data while the former is
used for the simulations of results. Diagnostic tools are also provided for
further analysis of background data.

Additionally, we provide the flexibility for the user to further explore
relative risk ratios from the output as well as provide colormapping. These can then be easily mapped onto a user-generated map.
}
\author{
Maintainer: Maria Kamenetsky <maria.kamenetsky@gmail.com>
}
\references{
Kamenetsky, M., Lee, J., Zhu, J., Gangnon, R.E. Space and Space-Time Cluster Detection Using the Lasso." (In preparation)
}
\keyword{ package }
\examples{
\donttest{
data(japanbreastcancer)
#Set Some Initial Conditions
x1=utmJapan$utmx/1000
y1=utmJapan$utmy/1000
#Set initial parameters
rMax <- 20 
Time=5
nsim=2 #number of simulations
center=1
radius=18
timeperiods = c(1:3)
risk.ratio=2
theta = 10000
overdispfloor = TRUE
nullmod <- TRUE
dframe <- dframe[,-1]
clst <- toclust(japanbreastcancer, expected = japanbreastcancer$expdeath, 
 observed = japanbreastcancer$death,timeperiod = japanbreastcancer$period, covars = FALSE)
res <- clust_sim(clst, x1,y1,rMax, Time, nsim,center, radius, risk.ratio, 
 timeperiod, utm=TRUE, byrow=TRUE, 
threshold, space= "both",theta = theta, nullmod = TRUE, overdispfloor)
}
\donttest{
data(japanbreastcancer)
#Set Some Initial Conditions
x1=utmJapan$utmx/1000
y1=utmJapan$utmy/1000
rMax <- 20 
Time=5
japanbreastcancer <- japanbreastcancer[,-1] #get rid of indicator column
clst <- toclust(japanbreastcancer, expected = japanbreastcancer$expdeath, 
 observed = japanbreastcancer$death,timeperiod = japanbreastcancer$period, covars = FALSE)
system.time(res <- clust(clst, x1,y1, rMax, Time, utm=TRUE, byrow=TRUE, 
 space="both", overdispfloor=TRUE, cv=NULL))
 }
%% ~~ simple examples of the most important functions ~~
}