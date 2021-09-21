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
    double p;
} data_struct;

typedef struct {
    NumericVector x;
    bool verbose;
    int fcount;
    double best_est;
} data_struct_new;

typedef struct {
    NumericVector x;
    bool verbose;
    int fcount;
    double best_est;
    double ls;
    double ln;
    double qn;
    double a;
} data_struct_p;

/*
------------------------------- helper functions -------------------------------
*/

// Check if user has requested abort every 1000 iterations & if so exit
//
// @param fcount double The current number of iterations of the objective fun
void check_user_input(int fcount) {
    if (fcount % 1000 == 0) {
        checkUserInterrupt();
    }
}

// check optimise exit status and display important info about them to user
//
// @param optimise_result Integer The exit code from the nlopt optimisation
// @param verbose Bool To display messages or not (defaults to TRUE)
void check_optimise_result(int optimise_result, bool verbose = true) {
    if (verbose) {
        if (optimise_result < 0) {
            Rcout << "NLOPT FAIL! Error code: " << optimise_result << std::endl;
            if (optimise_result == -4) {
                Rcout << "Fail is due to round off errors"
                    "- result may still be usefull" << std::endl;
            }
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
//' @param ls double Scout rate
//' @param qn double Quality
//' @param a double alpha value
//' @export
// [[Rcpp::export]]
double scout_dist(double x, double m,
                  double ls, double qn, double a)
{
    double maxpart = (-1+qn-x*a)/a;
    if (maxpart < 0)
    {
        maxpart = 0;
    }
    double result =
        (
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
                    double ln, double qn, double a)
{
    double maxpart = -x+((-1+qn)/(qn*a));
    if (maxpart < 0)
    {
        maxpart = 0;
    }
    double result =
        (
            (2 * exp(-M_PI * pow(x, 2) * ln) * M_PI * x * ln * maxpart) /
            (((exp(-pow(m, 2) * M_PI * ln) * (-1 + qn - m * qn * a)) /
                (qn * a)) + ((-1 + erf(m * sqrt(M_PI) * sqrt(ln)) +
                erfc((sqrt(M_PI) * (-1 + qn) * sqrt(ln)) /
                (qn * a))) / (2 * sqrt(ln))))
    );
    return result;
}

//' New model function for scouts
//'
//' @param x double Foraging distance
//' @param m double Minimum foraging distance
//' @param bs double scout rate
//' @param as double scout alpha
//' @export
// [[Rcpp::export]]
double scout_dist_new(double x, double m, double bs, double as)
{
    double maxpart = 1-as*x;
    if (maxpart < 0)
    {
        maxpart = 0;
    }
    double result =
        as * (
            (bs * exp(-bs*as*(x-m))*maxpart) /
            ((1-as*m) - pow(bs, -1)*(1-exp(-bs*(1-as*m))))
        );
    return result;
}

//' New model function for recruits
//'
//' @param br double Recruit rate
//' @param ar double Recruit alpha
//' @inheritParams scout_dist_new
//' @export
// [[Rcpp::export]]
double recruit_dist_new(double x, double m, double br, double ar)
{
    double maxpart = 1-ar*x;
    if (maxpart < 0)
    {
        maxpart = 0;
    }
    double result = (maxpart*2*M_PI*pow(ar,2)*br*x*exp(-M_PI*br*pow(ar*x, 2))) /
        ((1-ar*m)*exp(-M_PI*br*pow(ar*m, 2)) +
            ((erf(ar*sqrt(M_PI*br)*m)-erf(sqrt(M_PI*br)))/(2*sqrt(br))));
    return result;
}

/*
------------------------------- likelihood functions ---------------------------
*/

//' Log-likelihood function for scout and recruit superposition using
//' re-parameterised model
//'
//' @param x double* Pointer to array of foraging distances
//' @param p double Proportion of scouts (0<=p<=1)
//' @inheritParams scout_dist_new
//' @inheritParams recruit_dist_new
//' @export
// [[Rcpp::export]]
double loglike_model_all_new(NumericVector x,
                         double p, double bs,
                         double br, double as, double ar)
{
    const int x_size = x.size();
    const double m = min(x);
    double ll = 0;
    double tmp_ll = 0;
    for (int i = 0; i < x_size; i++)
    {
        tmp_ll = (p * scout_dist_new(x[i], m, bs, as) +
            (1 - p) * recruit_dist_new(x[i], m, br, ar));
        if (tmp_ll < 1e-99) {
            ll += -1e99-ar;
        } else {
            ll += log(tmp_ll);
        };
    }
    return ll;
}

//' Log-likelihood function for scout and recruit superposition
//'
//' @param x double* Pointer to array of foraging distances
//' @param p double Proportion of scouts (0<=p<=1)
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
    double tmp_ll = 0;
    for (int i = 0; i < x_size; i++)
    {
        tmp_ll = (p * scout_dist(x[i], m, ls, qn, a) +
            (1 - p) * recruit_dist(x[i], m, ln, qn, a));
        if (tmp_ll < 1e-99) {
            ll += -1e8+(qn-1)/a;
        } else {
            ll += log(tmp_ll);
        };
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
        ll += log(scout_dist(x[i], m, ls, qn, a));
    }
    return ll;
}

/*
------------------------------- objective functions ----------------------------
*/

// Objective function for new collective model
//
// param: n unsigned, record parameter required for nlopt
// param: params NumericVector, array of parameter estimates to run the model with
// param: grad double*, gradient value (unused but required as positional)
// param: x NumericVector, The data to load fit to
double objective_model_all_new(
    unsigned n, const double* params, double* grad, void* f_data
    )
{
    data_struct_new* data = (data_struct_new*) (f_data);
    check_user_input(data->fcount);
    data->fcount++;
    double result = loglike_model_all_new(
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
        Rcout << "p = "
            << params[0]
            << ", bs = "
            << params[1]
            << ", br = "
            << params[2]
            << ", as = "
            << params[3]
            << ", ar = "
            << params[4]
            << std::endl;
        Rcout << "-----" << std::endl;
    }
    return result;
}

// Objective function for collective model for just p
//
// param: n unsigned, record parameter required for nlopt
// param: NumericVector, array of parameter estimates to run the model with
// param: grad double*, gradient value (unused but required as positional)
// param: x NumericVector, The data to load fit to
double objective_model_p(
    unsigned n, const double* params, double* grad, void* f_data
    )
{
    data_struct_p* data = (data_struct_p*) (f_data);
    check_user_input(data->fcount);
    data->fcount++;
    double result = loglike_model_all(
        data->x,
        params[0],
        data->ls,
        data->ln,
        data->qn,
        data->a
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
            << std::endl;
        Rcout << "-----" << std::endl;
    }
    return result;
}

// Objective function for collective model
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
        data->p,
        params[0],
        params[1],
        params[2],
        params[3]
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

// Objective function for scout model
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

//' Optimise function for fitting the new collective model using NLOPT
//'
//' @param x NumericVector Foraging distance
//' @param params NumericVector parameter estimates to run the model with
//' @param lb NumericVector, array of lower bounds for each paramater
//' @param ub NumericVector, array of upper bounds for each parameter
//' @param verbose Bool, to display optimisation as it runs, defaults to FALSE
//' @param xtol double, The absolute tolerance on function value. If 0 (default)
//' then default to nlopt default value.
//' @param model int The model to run. Must be 0 or 1 which means 'all' or
//' 'scout' respectively,  defaults to 0 ('all')
//' @export
// [[Rcpp::export]]
NumericVector optimise_model(
    NumericVector x, NumericVector params, NumericVector lb, NumericVector ub,
    bool verbose = false, double xtol = 0, int model = 0
)
{
    // put constant data into data_struct and initialise other data
    data_struct_new ds;
    ds.x = x;
    ds.fcount = 0;
    ds.best_est = -1e-90;
    ds.verbose = verbose;

    // set up NLOPT
    nlopt_opt opt;
    double maxf;

    // set up the model function to run


    // Add objective to NLOPT
    opt = nlopt_create(NLOPT_LN_NELDERMEAD, params.size());
    nlopt_set_max_objective(opt, objective_model_all_new, &ds);
    nlopt_set_lower_bounds(opt, lb.begin());
    nlopt_set_upper_bounds(opt, ub.begin());
    if (xtol != 0) {
        nlopt_set_xtol_rel(opt, xtol);
    }

    // run optimisation for all parameters but p
    nlopt_result optimise_result = nlopt_optimize(opt, params.begin(), &maxf);
    check_optimise_result(optimise_result, verbose);

    // store results
    NumericVector results(params.size()+1);
    results[0] = maxf;
    for (int i = 0; i < params.size(); i++) {
        results[i + 1] = params[i];
    }
    nlopt_destroy(opt);
    return results;
}

//' Optimise function for fitting a model using NLOPT
//'
//' @param x NumericVector Foraging distance
//' @param params NumericVector parameter estimates to run the model with
//' @param lb NumericVector, array of lower bounds for each paramater
//' @param ub NumericVector, array of upper bounds for each parameter
//' @param verbose Bool, to display optimisation as it runs, defaults to FALSE
//' @param xtol double, The absolute tolerance on function value. If 0 (default)
//' then default to nlopt default value.
//' @param model int The model to run. Must be 0 or 1 which means 'all' or
//' 'scout' respectively,  defaults to 0 ('all')
//' @export
// [[Rcpp::export]]
NumericVector optimise_model_old(
    NumericVector x, NumericVector params, NumericVector lb, NumericVector ub,
    bool verbose = false, double xtol = 0, int model = 0
)
{
    // put constant data into data_struct and initialise other data
    data_struct ds;
    ds.x = x;
    ds.fcount = 0;
    ds.best_est = -9999999999999;
    ds.verbose = verbose;

    // set up NLOPT
    NumericVector newparams(4);
    NumericVector newlb(4);
    NumericVector newub(4);
    nlopt_opt opt;
    double maxf;

    // set up the model function to run
    double (*objective_fun)(unsigned, const double*, double*, void*);
    if (model == 0) {
        ds.p = params[0];
        newparams = params[Range(1, params.size()-1)];
        newlb = lb[Range(1, lb.size()-1)];
        newub = ub[Range(1, ub.size()-1)];
        objective_fun = &objective_model_all;
    } else if (model == 1) {
        newparams = params[Range(0, params.size())];
        newlb = lb[Range(0, lb.size())];
        newub = ub[Range(0, ub.size())];
        objective_fun = &objective_model_scout;
    } else {
        stop(
            "Model given is inconsistent with avaliable models.\n" \
            "Only 0 or 1 (meaning all or scout) is permited"
        );
    }

    // Add objective to NLOPT
    opt = nlopt_create(NLOPT_LN_NELDERMEAD, newparams.size());
    nlopt_set_max_objective(opt, objective_fun, &ds);
    nlopt_set_lower_bounds(opt, newlb.begin());
    nlopt_set_upper_bounds(opt, newub.begin());
    if (xtol != 0) {
        nlopt_set_xtol_rel(opt, xtol);
    }

    // run optimisation for all parameters but p
    nlopt_result optimise_result = nlopt_optimize(opt, newparams.begin(), &maxf);
    check_optimise_result(optimise_result, verbose);

    // maximise p using Newtons method
    data_struct_p dsp;
    dsp.x = x;
    dsp.fcount = 0;
    dsp.best_est = -9999999;
    dsp.verbose = verbose;
    dsp.ls = newparams[0];
    dsp.ln = newparams[1];
    dsp.qn = newparams[2];
    dsp.a = newparams[3];

    nlopt_opt opt_p;
    NumericVector pest = {params[0]};
    NumericVector plb = {0};
    NumericVector pub = {1};
    double maxf_p;

    opt_p = nlopt_create(NLOPT_LN_NELDERMEAD, pest.size());
    nlopt_set_max_objective(opt_p, objective_model_p, &dsp);
    nlopt_set_lower_bounds(opt_p, plb.begin());
    nlopt_set_upper_bounds(opt_p, pub.begin());
    if (xtol != 0) {
        nlopt_set_xtol_rel(opt_p, xtol);
    }

    // run optimisation for just p
    nlopt_result optimise_result_p = nlopt_optimize(opt_p, pest.begin(), &maxf_p);
    check_optimise_result(optimise_result_p, verbose);

    // move optimised params back into params
    if (model == 0) {
        params[Range(1, params.size())] = newparams;
        params[0] = pest[0];
    } else {
        params = newparams;
    }

    // store results
    NumericVector results(params.size() + 1);
    results[0] = maxf_p;
    for (int i = 0; i < results.size()-1; i++) {
        results[i + 1] = params[i];
    }
    nlopt_destroy(opt);
    nlopt_destroy(opt_p);
    return results;
}
