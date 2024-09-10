#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix matrix_complement(NumericMatrix mat, double new_val, bool bin, bool dir, bool loops) {
  int nrow = mat.nrow();
  int ncol = mat.ncol();
  
  for (int i = 0; i < nrow; ++i) { //For each row
    int inner_limit = dir ? ncol : i + 1; //If symmetric, only one triangle
    for (int j = 0; j < inner_limit; ++j) { //For each column (up to `i`)
      if (mat(i, j) == 0) { // If there was no tie (value is 0)
        mat(i, j) = new_val; // Use new_val as predetermined in R
      } else { // If the tie already existed
        if (bin) { // For binary matrices
          mat(i, j) = 0; // Simply remove it
        } else { // For valued matrices
          mat(i, j) = 1 / mat(i, j); // Use the reciprocal
        }
      }
      if(!dir)mat(j, i) = mat(i, j); // If symmetric, fill the other triangle
    }
  }
  
  if (!loops) {
    for (int i = 0; i < std::min(nrow, ncol); ++i) {
      mat(i, i) = 0;
    }
  }
  
  return mat; // Return the complement matrix 
}
