#include <Rcpp.h>

using namespace Rcpp;

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar)

// For more on using Rcpp click the Help button on the editor toolbar
//' @title
//' Calculates sequencce of values based on maximum radius
//' @name prod_YxCpp
//' @param vY numeric vector
//' @param clast last observation in sequence 
//' @param ccenter centroid value
//' @return numeric vector that feeds into st_matCpp function
//' @description helper function
//' @export
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

//' @title
//' Helper function
//' @name st_matCpp
//' @param X takes in numeric matrix
//' @param last takes in integer vector of last observation in a sequence given a radius
//' @param center centroid value
//' @param T integer
//' @return numeric vector
//' @description helper function
//' @export
// [[Rcpp::export]]
NumericVector st_matCpp (NumericMatrix X, IntegerVector last, IntegerVector center, int T){
    int n = X.nrow();
    int c = X.ncol();
    NumericMatrix U(n,c);
    NumericVector tempo(n);
    for(int j = 0; j < c; j++){
        tempo = prod_YxCpp(X(_,j), last, center);
        for(int i = 0; i < n; i++){
            U(i,j) = tempo(i);
        }
    }
    
    NumericMatrix W(n,T*(T+1)/2);
    int counter = 0;
    for(int i = 0; i < T; i++){
        NumericVector v(n);
        for(int j = 0; j < T-i; j++){
            v += U(_,i+j);
            for(int h = 0; h < n; h++){
                W(h,counter) = v(h);
            }
            counter++;
        }
    }
    return(W);
}

