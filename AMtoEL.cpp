#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
DataFrame adjMatrixToEdgeList(NumericMatrix mat) {
  // Get the dimensions of the matrix
  int n = mat.nrow();
  int m = mat.ncol();
  // First, count the number of non-zero elements
  int non_zero_count = 0; // Initialise counter
  for (int i = 0; i < n; i++) { // For each row
    for (int j = 0; j < m; j++) { // For each column
      if (mat(i, j) != 0) { // If there is an edge
        non_zero_count++; // Increment counter
      }
    }
  }
  
  // Pre-allocate vectors with the known size
  IntegerVector from(non_zero_count);
  IntegerVector to(non_zero_count);
  NumericVector weight(non_zero_count);
  
  // Fill vectors using indices
  int index = 0;
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < m; j++) {
      if (mat(i, j) != 0) {
        // Indexing in R starts from 1, in C++ from 0
        from[index] = i + 1; 
        to[index] = j + 1;
        // Extract the weight
        weight[index] = mat(i, j);
        // Advance the index
        index++;
      }
    }
  }
  
  // Create a DataFrame to hold the edge list
  return DataFrame::create(
    Named("from") = from,
    Named("to") = to,
    Named("weight") = weight
  );
}
