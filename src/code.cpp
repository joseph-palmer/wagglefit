/*
 * Example code for using Rcpp and C++
 * This is just a series of function which show how to get certain types in and
 * out of C++ from R. For example, altering an array in place.
 */
#include <Rcpp.h>
using namespace Rcpp;


double sum_array(double* x, int n)
{
  double res = 0;
  for (int i = 0; i < n; i++) {
    res += x[i];
  }
  return res;
}

// Multiply an array by 2 in place
// These functions (using pointers) cannot be linked directly with R
// Make sure you write tests for them in C++
void array_times_two(double* x, int n)
{
  for (int i = 0; i < n; i++) {
    x[i] = x[i]*2;
  }
}

//' Alter a NumericVector inplace
//'
//' @param x NumericArray
//' @export
// [[Rcpp::export]]
void alter_in_place(NumericVector x) {
  double* x_ = &x[0];

  // We can now pass the pointer around, e.g.:
  // pass the array pointer to function which adds them up
  int n = x.size();
  double res = sum_array(x_, n);
  Rcout << res << std::endl;

  array_times_two(x_, n);
}
