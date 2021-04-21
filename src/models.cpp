#include <iostream>
#include <cmath>
#include <Rcpp.h>
using namespace Rcpp;

//' Model function for scouts
//'
//' @param x double Foraging distance
//' @param m double Minimum foraging distance
//' @param p double Proportion of scouts (0<=p<=1)
//' @param ls double Scout rate
//' @param qn double Quality
//' @param a double alpha value
//' @export
// [[Rcpp::export]]
double scout_dist(double x, double m,
                  double p, double ls,
                  double qn, double a)
{
    double maxpart = (-1+qn-x*a)/a;
    if (maxpart < 0)
    {
        maxpart = 0;
    }
    double result =
    p * (
        (exp(((-1 + qn + m * a - x * a) * ls) / a) * a * pow(ls, 2) * maxpart) /
        (exp(m * ls) * a - exp(((-1 + qn) * ls) / a) *
            (a + ls - qn * ls + m * a * ls))
    );
    return result;
}

//' Model function for recruits
//'
//' @param ln double Recruit rate
//' @inheritParams scout_dist
//' @export
// [[Rcpp::export]]
double recruit_dist(double x, double m,
                    double p, double ln,
                    double qn, double a)
{
    double maxpart = -x+((-1+qn)/(qn*a));
    if (maxpart < 0)
    {
        maxpart = 0;
    }
    double result =
        (1 - p) * (
            (2 * exp(-M_PI * pow(x, 2) * ln) * M_PI * x * ln * maxpart) /
            (((exp(-pow(m, 2) * M_PI * ln) * (-1 + qn - m * qn * a)) /
                (qn * a)) + ((-1 + erf(m * sqrt(M_PI) * sqrt(ln)) +
                erfc((sqrt(M_PI) * (-1 + qn) * sqrt(ln)) /
                (qn * a))) / (2 * sqrt(ln))))
    );
    return result;
}

//' Log-likelihood function for scout and recruit superposition
//'
//' @param x double* Pointer to array of foraging distances
//' @param x_size int Number of foraging distances
//' @inheritParams scout_dist
//' @inheritParams recruit_dist
//' @export
// [[Rcpp::export]]
double loglike_model_all(NumericVector x, int x_size,
                         double p, double ls,
                         double ln, double qn, double a)
{
    double m = min(x);
    double ll = 0;
    for (int i = 0; i < x_size; i++)
    {
        ll += log(scout_dist(x[i], m, p, ls, qn, a) * +
            recruit_dist(x[i], m, p, ln, qn, a));
    }
    return ll;
}

//' Log-likelihood function for scouts
//'
//' @inheritParams loglike_model_all
//' @export
// [[Rcpp::export]]
double loglike_model_scout(NumericVector x, int x_size,
                           double ls, double qn, double a)
{
    double m = min(x);
    double ll = 0;
    for (int i = 0; i < x_size; i++)
    {
        ll += log(scout_dist(x[i], m, 1.0, ls, qn, a));
    }
    return ll;
}

//' Log-likelihood function for recruits
//'
//' @inheritParams loglike_model_all
//' @export
// [[Rcpp::export]]
double loglike_model_recruit(NumericVector x, int x_size,
                             double ln, double qn, double a)
{
    double m = min(x);
    double ll = 0;
    for (int i = 0; i < x_size; i++)
    {
        ll += log(recruit_dist(x[i], m, 0, ln, qn, a));
    }
    return ll;
}

//' Multiply a number by two
//'
//' @param x A single integer.
//' @export
// [[Rcpp::export]]
NumericVector times_two(NumericVector x) {
    return x * 2;
}
