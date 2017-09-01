#include <Rcpp.h>
using namespace Rcpp;

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar)

//' Calculates maximum column for spatial radius
//' @title return maximum center based on a raidus supplied
//' @param N total number of centroids
//' @param iMax each row of clusters dataframe (from clusters2df function)
//' @param cn identifier of centroid to which the count belongs to
//' @param clast last observation of that centroid group
//' @return integer vector for the last observation in a sequence of observations by growing radius, given users rMax specifications.
//' @export
// [[Rcpp::export]]
IntegerVector max_colCpp (const int N, int iMax, IntegerVector cn, IntegerVector clast){
  Rcpp::IntegerVector maxcol(N,0);    
  int num_ones = cn[iMax-1];
  for(int i = iMax; i > (iMax - num_ones); i--) {
      maxcol[clast[i-1]-1]=1;
      }
  return(maxcol);
}
