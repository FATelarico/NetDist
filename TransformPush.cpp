#include <Rcpp.h>
#include <algorithm>
#include <limits>

// [[Rcpp::export]]
// # Rcpp::NumericVector TransformPush(Rcpp::NumericVector input) {
Rcpp::NumericMatrix TransformPush(Rcpp::NumericVector input) {
  Rcpp::NumericVector result(input.size());
  double LargestNegative = std::numeric_limits<double>::min();
  
  // Find the largest negative number in the vector
  for (double num : input) {
    if (num < 0 && num > LargestNegative) {
      LargestNegative = num;
    }
  }
  
  // If `Largest Negative` is zero, pick the smallest non-zero number
  if(LargestNegative == 0) {
    LargestNegative = std::numeric_limits<double>::max();
    for (double num : input) {
      if (num != 0 && num < LargestNegative) {
        // Update if a smaller non-zero number is found
        LargestNegative = num; 
      }
    }
  }
  
  // Process each number in the input vector
  for (int i = 0; i < input.size(); ++i) {
    double num = input[i];
    if (num < 0) {
      // Raise negative numbers to the power of -1
      // <=! Zeroes would not work well here =>
      result[i] = 1.0 / num;
    } else {
      // Add to positive numbers the inverse of the largest negative
      //  number's absolute value
      // <=! Zeroes work fine here=>
      result[i] = num + (1.0 / std::abs(LargestNegative));
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