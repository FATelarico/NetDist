#include <Rcpp.h>
#include <unordered_map>
using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::NumericVector tablecpp(NumericVector levels, NumericVector x) {
  // Create a map to store counts of levels elements
  std::unordered_map<double, int> counts;
  
  // Initialize the map with levels elements, setting initial counts to 0
  for (double num : levels) {
    counts[num] = 0;
  }
  
  // Count occurrences of levels elements in the x vector
  for (double num : x) {
    if (counts.find(num) != counts.end()) {
      counts[num]++;
    }
  }
  
  // Prepare the result as a named vector
  Rcpp::NumericVector result(levels.size());
  for (int i = 0; i < levels.size(); ++i) {
    result[i] = counts[levels[i]];
  }
  result.names() = levels;
  
  return result;
}