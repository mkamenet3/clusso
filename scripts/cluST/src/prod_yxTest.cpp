#include <Rcpp.h>
using namespace Rcpp;

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar)

// For more on using Rcpp click the Help button on the editor toolbar

// [[Rcpp::export]]
NumericVector prod_YxCpp (NumericVector vY, IntegerVector clast, IntegerVector ccenter){
  int ind;
  int M = clast.length();
  NumericVector cums(M);
  for(int i  = 0; i < M; i++){
    //  cums[i] = 0.0;
    double tmp;
    ind = clast[i] - 1;
	if((!i) || (ccenter[i]!=ccenter[i-1])){
			cums[i]=vY[ind];}
		else{
			cums[i]=tmp+vY[ind];
                  }
        tmp = cums[i];
  }
  return(cums);
}



