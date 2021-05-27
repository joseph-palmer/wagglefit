#include <iostream>
#include <cmath>
#include <Rcpp.h>
#include <nloptrAPI.h>
using namespace Rcpp;

/*
------------------------------- structs and data -- ----------------------------
*/


// struct to hold additional data for optimisation
//
// x: The distance data to fit to
// verbose: To display messages to console
// fcount: number of times the objective function has been called
// best_est: the best estimate up to that point
typedef struct {
    NumericVector x;
    bool verbose;
    int fcount;
    double best_est;
} data_struct;

/*
------------------------------- helper functions -------------------------------
*/

// Check if user has requested abort every 1000 iterations
//
// @param fcount double The current number of iterations of the objective fun
void check_user_input(int fcount) {
    if (fcount % 1000 == 0) {
        checkUserInterrupt();
    }
}

void check_optimise_result(int optimise_result) {
    if (optimise_result < 0) {
        Rcout << "NLOPT FAIL! Error code: " << optimise_result << std::endl;
        if (optimise_result == -4) {
            Rcout << "Fail is due to round off errors"
                "- result may still be usefull" << std::endl;
        }
    }
}

/*
------------------------------ distribution functions --------------------------
*/

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

/*
------------------------------- likelihood functions ---------------------------
*/

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
------------------------------- objective functions ----------------------------
*/

// Objective function for scout model
//
// param: n unsigned, record parameter required for nlopt
// param: NumericVector, array of parameter estimates to run the model with
// param: grad double*, gradient value (unused but required as positional)
// param: x NumericVector, The data to load fit to
double objective_model_all(
    unsigned n, const double* params, double* grad, void* f_data
    )
{
    data_struct* data = (data_struct*) (f_data);
    check_user_input(data->fcount);
    data->fcount++;
    double result = loglike_model_all(
        data->x,
        params[0],
        params[1],
        params[2],
        params[3],
        params[4]
    );
    if (result > data->best_est) {
        data->best_est = result;
    }
    if (data->verbose) {
        Rcout << "Iteration: "
            << data->fcount
            << ", Result: "
            << result
            << ", Best estimate: "
            << data->best_est
            << std::endl;
        Rcout << "ls = "
            << params[0]
            << ", qn = "
            << params[1]
            << ", a = "
            << params[2]
            << std::endl;
        Rcout << "-----" << std::endl;
    }
    return result;
}

// Objective function for scout and recruit superposition
//
// param: n unsigned, record parameter required for nlopt
// param: NumericVector, array of parameter estimates to run the model with
// param: grad double*, gradient value (unused but required as positional)
// param: x NumericVector, The data to load fit to
double objective_model_scout(
    unsigned n, const double* params, double* grad, void* f_data
    )
{
    data_struct* data = (data_struct*) (f_data);
    check_user_input(data->fcount);
    data->fcount++;
    double result = loglike_model_scout(
        data->x,
        params[0],
        params[1],
        params[2]
    );
    if (result > data->best_est) {
        data->best_est = result;
    }
    if (data->verbose) {
        Rcout << "Iteration: "
            << data->fcount
            << ", Result: "
            << result
            << ", Best estimate: "
            << data->best_est
            << std::endl;
        Rcout << "p = "
            << params[0]
            << ", ls = "
            << params[1]
            << ", ln = "
            << params[2]
            << ", qn = "
            << params[3]
            << ", a = "
            << params[4]
            << std::endl;
        Rcout << "-----" << std::endl;
    }
    return result;
}

/*
------------------------------- optimising functions ---------------------------
*/

//' Optimise function for fitting a model using NLOPT
//'
//' @param x NumericVector Foraging distance
//' @param params NumericVector parameter estimates to run the model with
//' @param lb NumericVector, array of lower bounds for each paramater
//' @param ub NumericVector, array of upper bounds for each parameter
//' @param verbose Bool, to display optimisation as it runs, defaults to FALSE
//' @param xtol double, The absolute tolerance on function value. If 0 (default)
//' then default to nlopt default value.
//' @export
// [[Rcpp::export]]
NumericVector optimise(
    NumericVector x, NumericVector params, NumericVector lb, NumericVector ub,
    bool verbose = false, double xtol = 0
)
{
    // put constant data into data_struct and initialise other data
    data_struct ds;
    ds.x = x;
    ds.fcount = 0;
    ds.best_est = -9999999;
    ds.verbose = verbose;

    // set up the model function to run
    double (*objective_fun)(unsigned, const double*, double*, void*);
    if (params.size() == 3) {
        objective_fun = &objective_model_scout;
    } else if (params.size() == 5) {
        objective_fun = &objective_model_all;
    } else {
        stop(
            "Number of parameters is inconsistent with avaliable models.\n" \
            "Only 3 and 5 parameters are permited, you provided %i",
            params.size()
        );
    }

    // set up NLOPT
    nlopt_opt opt;
    double maxf;
    opt = nlopt_create(NLOPT_LN_COBYLA, params.size());
    nlopt_set_max_objective(opt, objective_fun, &ds);
    nlopt_set_lower_bounds(opt, lb.begin());
    nlopt_set_upper_bounds(opt, ub.begin());
    if (xtol != 0) {
        nlopt_set_xtol_rel(opt, xtol);
    }

    // run optimisation
    nlopt_result optimise_result = nlopt_optimize(opt, params.begin(), &maxf);
    check_optimise_result(optimise_result);

    // store results
    NumericVector results(params.size() + 1);
    results[0] = maxf;
    for (int i = 0; i < results.size()-1; i++) {
        results[i + 1] = params[i];
    }
    nlopt_destroy(opt);
    return results;
}
