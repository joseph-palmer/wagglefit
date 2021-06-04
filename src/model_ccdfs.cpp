#include <iostream>
#include <cmath>
#include <Rcpp.h>
using namespace Rcpp;



//' Heaviside equation. Replace negatives with 0 and positives with 1.
//'
//' @param x double the parameter to convert with heaviside
//' @export
// [[Rcpp::export]]
double heaviside(double x)
{
  if (x < 0) return 0;
  else return 1;
}

//' Model CCDF function for scouts
//'
//' @param x double Foraging distance
//' @param p double Proportion of scouts (0<=p<=1)
//' @param ls double Scout rate
//' @param qn double Quality
//' @param a double alpha value
//' @param m double Minimum foraging distance
//' @return result double the ccdf for that foraging distance (x)
//' @export
// [[Rcpp::export]]
double scout_ccdf(double x, double p,
                  double ls, double qn,
                  double a, double m)
{
  double result;
  result = p*(
    (exp((m-x-(1/a))*ls)*(
      exp((x+(1/a))*ls)*a-exp((qn*ls)/a)*(a+ls-qn*ls+x*a*ls)
    )*heaviside((-m+x)*(-x+((-1+qn)/a)))
    ) / (
      exp(m*ls)*a-exp(((-1+qn)*ls)/a)*(a+ls-qn*ls+m*a*ls)
    )) + heaviside(m-x);
  return result;
}

//' Model CCDF function for recruits
//'
//' @param x double Foraging distance
//' @param m double Minimum foraging distance
//' @param p double Proportion of scouts (0<=p<=1)
//' @param ln double Recruit rate
//' @param qn double Quality
//' @param a double alpha value
//' @param m double Minimum foraging distance
//' @return result double the ccdf for that foraging distance (x)
//' @export
// [[Rcpp::export]]
double recruit_ccdf(double x, double p,
                    double ln, double qn,
                    double a, double m)
{
  double result;
  result = (1-p)*(
    (exp(M_PI*(m-x)*(m+x)*ln) * (
      2*(1+qn*(-1+x*a))*sqrt(ln)+exp(M_PI*pow(x, 2)*ln)*qn*a*(
        -erf(sqrt(M_PI)*x*sqrt(ln))+erf((sqrt(M_PI)*
          (-1+qn)*sqrt(ln))/(qn*a))
      )
    )) / (
    2*(1+qn*(-1+m*a))*sqrt(ln)+exp(pow(m, 2)*M_PI*ln)*qn*a*(
      -erf(m*sqrt(M_PI)*sqrt(ln))+erf((sqrt(M_PI)*
        (-1+qn)*sqrt(ln))/(qn*a))
    ))
    ) * heaviside((((-1+qn)/(qn*a))-x)*(x-m));
  return result;
}

//' Model ccdf function for scout and recruit superposition. Stores results in
//' given array (y)
//'
//' @param x NumericVector foraging distances
//' @param y NumericVector storage array for the results
//' @inheritParams scout_ccdf
//' @inheritParams recruit_ccdf
//' @export
// [[Rcpp::export]]
void ccdf_model_all(NumericVector x, NumericVector y,
                      double p, double ls,
                      double ln, double qn, double a)
{
  const int x_size = x.size();
  const double m = min(x);
  for (int i = 0; i < x_size; i++) {
    y[i] = scout_ccdf(x[i], p, ls, qn, a, m) +
            recruit_ccdf(x[i], p, ln, qn, a, m);
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
                      double ls, double qn, double a)
{
  const int x_size = x.size();
  const double m = min(x);
  for (int i = 0; i < x_size; i++) {
    y[i] = scout_ccdf(x[i], 1, ls, qn, a, m);
  }
}

//' Model ccdf function for recruit model. Stores results in given array (y)
//'
//' @param x NumericVector foraging distances
//' @param y NumericVector storage array for the results
//' @inheritParams recruit_ccdf
//' @export
// [[Rcpp::export]]
void ccdf_model_recruit(NumericVector x, NumericVector y,
                      double ln, double qn, double a)
{
  const int x_size = x.size();
  const double m = min(x);
  for (int i = 0; i < x_size; i++) {
    y[i] = recruit_ccdf(x[i], 0, ln, qn, a, m);
  }
}

//' Get model ccdf for a given model
//'
//' @param x NumericVector foraging distances
//' @param params NumericVector parameters to calculate the ccdfs for. E.g. if
//' the scout model their will only be three (ls, qn, a) but if all their will
//' be 5 (p, ls, ln, qn, a). The number of params determines which ccdf to make
//' @param model int The model to run. Must be 0, 1 or 2 which means 'all',
//' 'scout' or 'recruit' respectively,  defaults to 0 ('all')
//' @return y NumericVector the ccdf
//' @export
// [[Rcpp::export]]
NumericVector model_ccdf(NumericVector x, NumericVector params,
                         int model = 0)
{
  NumericVector y(x.size());
  if (model == 1) {
    ccdf_model_scout(x, y, params[0], params[1], params[2]);
  } else if (model == 2) {
    ccdf_model_recruit(x, y, params[0], params[1], params[2]);
  } else if (model == 0) {
    ccdf_model_all(x, y, params[0], params[1], params[2], params[3], params[4]);
  } else {
        stop(
            "Model given is inconsistent with avaliable models.\n" \
            "Only 0, 1 or 2 (meaning scout, recruit and all) is permited"
        );
    }
  return y;
}
