#include <iostream>
#include <cmath>
#include <bits/stdc++.h>
#include <Rcpp.h>
using namespace Rcpp;

double scout_dist(double x, double m,
                  double p, double ls,
                  double qn, double a)
{
    // model function for scouts
    double maxpart = (1+qn-x*a)/a;
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

double recruit_dist(double x, double m,
                    double p, double ln,
                    double qn, double a)
{
    // model function for recruits
    double maxpart = x+((-1+qn)/(qn*a));
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

double loglike_model_all(double* x, int x_size,
                         double p, double ls,
                         double ln, double qn, double a)
{
    // loglikelihood for scout and recruit superposition
    double *m = std::min_element(x, x+x_size);
    double ll = 0;
    for (int i = 0; i < x_size; i++)
    {
        ll += log(scout_dist(x[i], *m, p, ls, qn, a) * +
            recruit_dist(x[i], *m, p, ln, qn, a));
    }
    return ll;
}

double loglike_model_scout(double* x, int x_size,
                           double ls, double qn, double a)
{
    // loglikelihood for scouts
    double *m = std::min_element(x, x+x_size);
    double ll = 0;
    for (int i = 0; i < x_size; i++)
    {
        ll += log(scout_dist(x[i], *m, 1.0, ls, qn, a));
    }
    return ll;
}

double loglike_model_recruit(double* x, int x_size,
                             double ln, double qn, double a)
{
    // loglikelihood for recruits
    double *m = std::min_element(x, x+x_size);
    double ll = 0;
    for (int i = 0; i < x_size; i++)
    {
        ll += log(recruit_dist(x[i], *m, 0, ln, qn, a));
    }
    return ll;
}

//' Multiply a number by two
//'
//' @param x A single integer.
//' @export
// [[Rcpp::export]]
int times_two(int x) {
    return x * 2;
}
