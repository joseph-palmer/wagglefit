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
//' @param bs double scout rate
//' @param as double scout alpha
//' @export
// [[Rcpp::export]]
double scout_dist(double x, double m, double bs, double as)
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

//' Model function for recruits
//'
//' @param br double Recruit rate
//' @param ar double Recruit alpha
//' @inheritParams scout_dist
//' @export
// [[Rcpp::export]]
double recruit_dist(double x, double m, double br, double ar)
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

//' Log-likelihood function for scout and recruit superposition
//' (collective model)
//'
//' @param x double* Pointer to array of foraging distances
//' @param p double Proportion of scouts (0<=p<=1)
//' @inheritParams scout_dist
//' @inheritParams recruit_dist
//' @export
// [[Rcpp::export]]
double loglike_model_collective(NumericVector x,
                         double p, double bs,
                         double br, double as, double ar)
{
    const int x_size = x.size();
    const double m = min(x);
    double ll = 0;
    double tmp_ll = 0;
    for (int i = 0; i < x_size; i++)
    {
        tmp_ll = (p * scout_dist(x[i], m, bs, as) +
            (1 - p) * recruit_dist(x[i], m, br, ar));
        if (tmp_ll < 1e-99) {
            ll += -1e99-ar;
        } else {
            ll += log(tmp_ll);
        };
    }
    return ll;
}

//' Log-likelihood function for individual foraging
//'
//' @inheritParams loglike_model_collective
//' @export
// [[Rcpp::export]]
double loglike_model_individual(NumericVector x,
                           double bs, double as)
{
    const int x_size = x.size();
    const double m = min(x);
    double ll = 0;
    double tmp_ll = 0;
    for (int i = 0; i < x_size; i++)
    {
        tmp_ll = scout_dist(x[i], m, bs, as);
        if (tmp_ll < 1e-99) {
            ll += -1e99-as;
        } else {
            ll += log(tmp_ll);
        };
    }
    return ll;
}

/*
------------------------------- objective functions ----------------------------
*/

// Objective function for collective model
//
// param: n unsigned, record parameter required for nlopt
// param: params NumericVector, array of parameter estimates to run the model with
// param: grad double*, gradient value (unused but required as positional)
// param: x NumericVector, The data to load fit to
double objective_model_collective(
    unsigned n, const double* params, double* grad, void* f_data
    )
{
    data_struct* data = (data_struct*) (f_data);
    check_user_input(data->fcount);
    data->fcount++;
    double result = loglike_model_collective(
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

// Objective function for individual model
//
// param: n unsigned, record parameter required for nlopt
// param: NumericVector, array of parameter estimates to run the model with
// param: grad double*, gradient value (unused but required as positional)
// param: x NumericVector, The data to load fit to
double objective_model_individual(
    unsigned n, const double* params, double* grad, void* f_data
    )
{
    data_struct* data = (data_struct*) (f_data);
    check_user_input(data->fcount);
    data->fcount++;
    double result = loglike_model_individual(
        data->x,
        params[0],
        params[1]
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
        Rcout << "bs = "
            << params[0]
            << ", as = "
            << params[1]
            << std::endl;
        Rcout << "-----" << std::endl;
    }
    return result;
}

/*
------------------------------- optimising functions ---------------------------
*/

//' Optimise function for fitting the collective model using NLOPT
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
    data_struct ds;
    ds.x = x;
    ds.fcount = 0;
    ds.best_est = -1e-90;
    ds.verbose = verbose;

    // set up NLOPT
    nlopt_opt opt;
    double maxf;

    // set up the model function to run
    double (*objective_fun)(unsigned, const double*, double*, void*);
    if (model == 0) {
        objective_fun = &objective_model_collective;
    } else if (model == 1) {
        objective_fun = &objective_model_individual;
    } else {
        stop(
            "Model given is inconsistent with avaliable models.\n" \
            "Only 0 or 1 (meaning all or scout) is permited"
        );
    }

    // Add objective to NLOPT
    opt = nlopt_create(NLOPT_LN_NELDERMEAD, params.size());
    nlopt_set_max_objective(opt, objective_fun, &ds);
    nlopt_set_lower_bounds(opt, lb.begin());
    nlopt_set_upper_bounds(opt, ub.begin());
    if (xtol != 0) {
        nlopt_set_xtol_rel(opt, xtol);
    }

    // run optimisation for all parameters
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


/*
------------------------------- pdf plot functions -----------------------------
*/

//' Model pdf function for scout and recruit superposition. Stores results in
//' given array (y)
//'
//' @param x NumericVector foraging distances
//' @param y NumericVector storage array for the results
//' @param p double Proportion of scouts (0<=p<=1)
//' @inheritParams scout_dist
//' @inheritParams recruit_dist
//' @export
// [[Rcpp::export]]
void pdf_model_all(NumericVector x, NumericVector y,
                    double p, double bs, double br,
                    double as, double ar)
{
  const int x_size = x.size();
  const double m = min(x);
  for (int i = 0; i < x_size; i++) {
    y[i] = p * scout_dist(x[i], m, bs, as) +
            (1 - p) * recruit_dist(x[i], m, br, ar);
  }
}

//' Model pdf function for scout model. Stores results in given array (y)
//'
//' @param x NumericVector foraging distances
//' @param y NumericVector storage array for the results
//' @inheritParams scout_dist
//' @export
// [[Rcpp::export]]
void pdf_model_scout(NumericVector x, NumericVector y,
                      double bs, double as)
{
  const int x_size = x.size();
  const double m = min(x);
  for (int i = 0; i < x_size; i++) {
    y[i] = scout_dist(x[i], m, bs, as);
  }
}

//' Get model pdf for a given model
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
NumericVector model_pdf(NumericVector x, NumericVector params,
                         int model = 0)
{
  NumericVector y(x.size());
  if (model == 0) {
    pdf_model_all(x, y, params[0], params[1], params[2], params[3], params[4]);
  } else if (model == 1) {
    Rcout << params << std::endl;
    pdf_model_scout(x, y, params[0], params[1]);
  } else {
        stop(
            "Model given is inconsistent with avaliable models.\n" \
            "Only 0 or 1 (meaning all or scout) is permited"
        );
    }
  return y;
}