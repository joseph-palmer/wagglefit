#include <iostream>
#include <cmath>
#include <Rcpp.h>
#include <nloptrAPI.h>
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
//' @inheritParams scout_dist
//' @inheritParams recruit_dist
//' @export
// [[Rcpp::export]]
double loglike_model_all(NumericVector x,
                         double p, double ls,
                         double ln, double qn, double a)
{
    const int x_size = x.size();
    const double m = min(x);
    double ll = 0;
    for (int i = 0; i < x_size; i++)
    {
        ll += log(scout_dist(x[i], m, p, ls, qn, a) +
            recruit_dist(x[i], m, p, ln, qn, a));
    }
    return ll;
}

//' Log-likelihood function for scouts
//'
//' @inheritParams loglike_model_all
//' @export
// [[Rcpp::export]]
double loglike_model_scout(NumericVector x,
                           double ls, double qn, double a)
{
    const int x_size = x.size();
    const double m = min(x);
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
double loglike_model_recruit(NumericVector x,
                             double ln, double qn, double a)
{
    const int x_size = x.size();
    const double m = min(x);
    double ll = 0;
    for (int i = 0; i < x_size; i++)
    {
        ll += log(recruit_dist(x[i], m, 0.0, ln, qn, a));
    }
    return ll;
}

/*
------------------------------- optimizing functions ---------------------------
*/

static int fcount = 0;
static double bestres = -99999999;

typedef struct {
    NumericVector x;
} data_struct;

// Optimize function for scout and recruit superposition
// param: n unsigned, record parameter required for nlopt
// param: NumericVector, array of parameter estimates to run the model with
// param: grad double*, gradient value (unused but required as positional)
// param: x NumericVector, The data to load fit to
double objective_model_all(
    unsigned n, const double* params, double* grad, void* f_data
    )
{
    fcount++;
    data_struct data = *(data_struct *) (f_data);
    double result = loglike_model_all(
        data.x,
        params[0],
        params[1],
        params[2],
        params[3],
        params[4]
    );
    if (result > bestres) {bestres = result;}
    Rcout << "count: " << fcount << ", result: " << result <<  " Best res: " << bestres<< std::endl;
    Rcout << "p " << params[0] << ", ls " << params[1] << ", ln " << params[2] << ", qn " << params[3] << ", a " << params[4] << std::endl;
    Rcout << "-----" << std::endl;
    return result;
}

//' Optimise function for scout and recruit superposition
//'
//' @param x NumericVector Foraging distance
//' @param params NumericVector parameter estimates to run the model with
//' @param lb NumericVector, array of lower bounds for each paramater
//' @param ub NumericVector, array of upper bounds for each parameter
//' @export
// [[Rcpp::export]]
NumericVector optimise_all(
    NumericVector x, NumericVector params, NumericVector lb, NumericVector ub
)
{
    fcount = 0;
    bestres = -99999999;
    double maxf = 0;
    nlopt_opt opt;
    NumericVector results(params.size() + 1);
    data_struct ds;
    ds.x = x;

    opt = nlopt_create(NLOPT_LN_COBYLA, 5); //NLOPT_GN_CRS2_LM
    nlopt_set_max_objective(opt, objective_model_all, &ds);
    nlopt_set_lower_bounds(opt, lb.begin());
    nlopt_set_upper_bounds(opt, ub.begin());
    nlopt_set_xtol_rel(opt, 1e-06);
    if (nlopt_optimize(opt, params.begin(), &maxf) < 0) {
        Rcout << "nlopt fail" << std::endl;
    }

    nlopt_destroy(opt);
    results[1] = maxf;
    for (int i = 1; i < results.size()-1; i++) {
        results[i] = params[i];
    }
    return results;
}
