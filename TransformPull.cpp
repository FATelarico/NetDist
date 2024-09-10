#include <Rcpp.h>
#include <algorithm>
#include <limits>

// [[Rcpp::export]]
// # Rcpp::NumericVector TransformPull(Rcpp::NumericVector input) {
Rcpp::NumericMatrix TransformPull(Rcpp::NumericVector input) {
  Rcpp::NumericVector result(input.size());
  double smallestPositive = std::numeric_limits<double>::max();
  
  // Find the smallest positive number in the vector
  for (double num : input) {
    if (num > 0 && num < smallestPositive) {
      smallestPositive = num;
    }
  }
  
  // If `smallestPositive` is zero, pick the largest non-zero number
  if(smallestPositive == 0) {
    smallestPositive = std::numeric_limits<double>::min();
    for (double num : input) {
      if (num != 0 && num > smallestPositive) {
        // Update if a smaller non-zero number is found
        smallestPositive = num; 
      }
    }
  }
  
  // Process each number in the input vector
  for (int i = 0; i < input.size(); ++i) {
    double num = input[i];
    if (num > 0) {
      // Raise positive numbers to the power of -1
      // <=! Zeroes would not work well here =>
      result[i] = 1.0 / num;
    } else {
      // Make negative numbers positive and add the inverse of
      // the absolute values of the smallest positive number
      // <=! Zeroes work fine here=>
      result[i] = -num + (1.0 / std::abs(smallestPositive));
    }
  }
  
  // <=! To return a vector=>
  
  // <= (1) Uncomment the line below=>
  // return result;
  
  // <= (2) Comment out lines 50-60=>
  int n = sqrt(result.size());
  Rcpp::NumericMatrix matrix(n, n);
  
  // Fill the matrix by columns
  for (int col = 0; col < n; ++col) {
    for (int row = 0; row < n; ++row) {
      matrix(row, col) = result[col * n + row];
    }
  }
  
  return matrix;
}
