#include <iostream>
#include <cmath>
#include <Rcpp.h>
using namespace Rcpp;

//' Model CCDF for scout distribution
//'
//' @param x double Foraging distance
//' @param m double Minimum foraging distance
//' @param bs double Scout rate
//' @param as double scout alpha
//' @return result double the ccdf for that foraging distance (x)
//' @export
// [[Rcpp::export]]
double scout_ccdf(double x, double m,
                  double bs, double as)
{
  double result;
  if (x > 1/as) {
    result = 0;
  } else {
    result = ((1-as*x)*exp(-as*bs*x)-pow(bs, -1)*(exp(-as*bs*x)-exp(-bs))) /
    ((1-as*m)*exp(-as*bs*m)-pow(bs, -1)*(exp(-as*bs*m)-exp(-bs)));
  }
  return result;
}

//' Model CCDF for recruit distribution
//'
//' @param x double Foraging distance
//' @param m double Minimum foraging distance
//' @param br double Recruit rate
//' @param ar double recruit alpha
//' @return result double the ccdf for that foraging distance (x)
//' @export
// [[Rcpp::export]]
double recruit_ccdf(double x, double m,
                    double br, double ar)
{
  double result;
  if (x > 1/ar) {
    result = 0;
  } else {
    result = ((1-ar*x)*exp(-M_PI*br*pow(ar*x, 2)) +
      ((erf(ar*sqrt(M_PI*br)*x)-erf(sqrt(M_PI*br)))/(2*sqrt(br)))) /
      ((1-ar*m)*exp(-M_PI*br*pow(ar*m, 2)) +
      ((erf(ar*sqrt(M_PI*br)*m)-erf(sqrt(M_PI*br)))/(2*sqrt(br))));
  }
  return result;
}

//' Model ccdf function for scout and recruit superposition. Stores results in
//' given array (y)
//'
//' @param x NumericVector foraging distances
//' @param y NumericVector storage array for the results
//' @param p double Proportion of scouts (0<=p<=1)
//' @inheritParams scout_ccdf
//' @inheritParams recruit_ccdf
//' @export
// [[Rcpp::export]]
void ccdf_model_all(NumericVector x, NumericVector y,
                    double p, double bs, double br,
                    double as, double ar)
{
  const int x_size = x.size();
  const double m = min(x);
  for (int i = 0; i < x_size; i++) {
    y[i] = p * scout_ccdf(x[i], m, bs, as) +
            (1 - p) * recruit_ccdf(x[i], m, br, ar);
  }
}

//' Model ccdf function for scout model. Stores results in given array (y)
//'
//' @param x NumericVector foraging distances
//' @param y NumericVector storage array for the results
//' @inheritParams scout_ccdf
//' @export
// [[Rcpp::export]]
void ccdf_model_scout(NumericVector x, NumericVector y,
                      double bs, double as)
{
  const int x_size = x.size();
  const double m = min(x);
  for (int i = 0; i < x_size; i++) {
    y[i] = scout_ccdf(x[i], m, bs, as);
  }
}

//' Get model ccdf for a given model
//'
//' @param x NumericVector foraging distances
//' @param params NumericVector parameters to calculate the ccdfs for. E.g. if
//' the scout model their will only be three (ls, qn, a) but if all their will
//' be 5 (p, ls, ln, qn, a). The number of params determines which ccdf to make
//' @param model int The model to run. Must be 0 or 1 which means 'all',
//' or 'scout' respectively,  defaults to 0 ('all')
//' @return y NumericVector the ccdf
//' @export
// [[Rcpp::export]]
NumericVector model_ccdf(NumericVector x, NumericVector params,
                         int model = 0)
{
  NumericVector y(x.size());
  if (model == 0) {
    ccdf_model_all(x, y, params[0], params[1], params[2], params[3], params[4]);
  } else if (model == 1) {
    Rcout << params << std::endl;
    ccdf_model_scout(x, y, params[0], params[1]);
  } else {
        stop(
            "Model given is inconsistent with avaliable models.\n" \
            "Only 0 or 1 (meaning all or scout) is permited"
        );
    }
  return y;
}
