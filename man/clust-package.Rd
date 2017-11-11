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
The main functions in this package are ```clust``` and ```clust_sim```, which
both perform both spatial and spatio-temporal cluster detection using the
lasso. The former works with observed and expected data while the former is
used for the simulations of results. Diagnostic tools are also provided for
further analysis of background data.

Additionally, we provide the flexibility for the user to further explore
relative risk ratios from the output as well as provide colormapping. These can then be easily mapped onto a user-generated map.
}
\author{
Maria Kamenetsky
Maintainer: Maria Kamenetsky <maria.kamenetsky@gmail.com>
}
\references{
~~ Literature or other references for background information ~~
}
\keyword{ package }
\seealso{
~~ Optional links to other man pages, e.g. ~~
~~ \code{\link[clust:clust-package]{clust}} ~~
}
\examples{
%% ~~ simple examples of the most important functions ~~
}