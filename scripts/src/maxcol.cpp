#include <Rcpp.h>
using namespace Rcpp;

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar)

// For more on using Rcpp click the Help button on the editor toolbar

// [[Rcpp::export]]
IntegerVector max_colCpp (const int N, int iMax, IntegerVector cn, IntegerVector clast){
  Rcpp::IntegerVector maxcol(N,0);    
  int num_ones = cn[iMax-1];
  for(int i = iMax; i > (iMax - num_ones); i--) {
      maxcol[clast[i-1]-1]=1;
      }
  return(maxcol);
}
