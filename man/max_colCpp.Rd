% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{max_colCpp}
\alias{max_colCpp}
\title{return maximum center based on a raidus supplied}
\usage{
max_colCpp(N, iMax, cn, clast)
}
\arguments{
\item{N}{total number of centroids}

\item{iMax}{each row of clusters dataframe (from clusters2df function)}

\item{cn}{identifier of centroid to which the count belongs to}

\item{clast}{last observation of that centroid group}
}
\value{
integer vector for the last observation in a sequence of observations by growing radius, given users rMax specifications.
}
\description{
Calculates maximum column for spatial radius
}
