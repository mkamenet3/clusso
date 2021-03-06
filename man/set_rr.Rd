% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/setRR.R
\name{set_rr}
\alias{set_rr}
\title{These functions will extract the risk ratios from the output}
\usage{
set_rr(lassoresult, vectors, Time, sim = FALSE)
}
\arguments{
\item{lassoresult}{List of QBIC, QAIC, QAICc estimates from the mylasso function}

\item{vectors}{dataframe of initial vectors of the observed and expected counts}

\item{Time}{number of time period}

\item{sim}{default is FALSE. If True, instead of returning observed over expected counts it will return the oracle.}
}
\value{
This returns a list of the risk ratios (observed over expected) as determined by 1) pure observed/expected counts,
2) observed based on QBIC path/expected; 3) observed based on QAIC path/expected; 4) observed based on QAICc path/expected.
}
\description{
set_rr
}
\details{
This function will create vectors of the risk ratios as determined by observed counts, QBIC, QAIC, and QAICc, respectively.
}
