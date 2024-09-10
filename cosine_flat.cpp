#include <Rcpp.h>
#include <cmath>

using namespace Rcpp;

// [[Rcpp::export]]
double cosine_similarity(NumericMatrix mat1, NumericMatrix mat2) {
  // Check if dimensions match
  if (mat1.nrow() != mat2.nrow() || mat1.ncol() != mat2.ncol()) {
    stop("Matrices must have the same dimensions");
  }
  
  int n = mat1.nrow() * mat1.ncol();
  
  // Flatten matrices into vectors
  NumericVector vec1(n);
  NumericVector vec2(n);
  
  int index = 0;
  for (int i = 0; i < mat1.nrow(); i++) {
    for (int j = 0; j < mat1.ncol(); j++) {
      vec1[index] = mat1(i, j);
      vec2[index] = mat2(i, j);
      index++;
    }
  }
  
  // Compute dot product and magnitudes
  double dot_product = 0.0;
  double magnitude1 = 0.0;
  double magnitude2 = 0.0;
  
  for (int i = 0; i < n; i++) {
    dot_product += vec1[i] * vec2[i];
    magnitude1 += vec1[i] * vec1[i];
    magnitude2 += vec2[i] * vec2[i];
  }
  
  // Calculate cosine similarity
  if (magnitude1 == 0.0 || magnitude2 == 0.0) {
    return 0.0; // Return 0 if either vector is zero vector
  } else {
    return dot_product / (std::sqrt(magnitude1) * std::sqrt(magnitude2));
  }
}