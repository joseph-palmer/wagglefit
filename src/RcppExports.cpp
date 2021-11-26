// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// alter_in_place
void alter_in_place(NumericVector x);
RcppExport SEXP _wagglefit_alter_in_place(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    alter_in_place(x);
    return R_NilValue;
END_RCPP
}
// scout_ccdf
double scout_ccdf(double x, double m, double bs, double as);
RcppExport SEXP _wagglefit_scout_ccdf(SEXP xSEXP, SEXP mSEXP, SEXP bsSEXP, SEXP asSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type m(mSEXP);
    Rcpp::traits::input_parameter< double >::type bs(bsSEXP);
    Rcpp::traits::input_parameter< double >::type as(asSEXP);
    rcpp_result_gen = Rcpp::wrap(scout_ccdf(x, m, bs, as));
    return rcpp_result_gen;
END_RCPP
}
// recruit_ccdf
double recruit_ccdf(double x, double m, double br, double ar);
RcppExport SEXP _wagglefit_recruit_ccdf(SEXP xSEXP, SEXP mSEXP, SEXP brSEXP, SEXP arSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type m(mSEXP);
    Rcpp::traits::input_parameter< double >::type br(brSEXP);
    Rcpp::traits::input_parameter< double >::type ar(arSEXP);
    rcpp_result_gen = Rcpp::wrap(recruit_ccdf(x, m, br, ar));
    return rcpp_result_gen;
END_RCPP
}
// ccdf_model_collective
void ccdf_model_collective(NumericVector x, NumericVector y, double p, double bs, double br, double as, double ar);
RcppExport SEXP _wagglefit_ccdf_model_collective(SEXP xSEXP, SEXP ySEXP, SEXP pSEXP, SEXP bsSEXP, SEXP brSEXP, SEXP asSEXP, SEXP arSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< double >::type p(pSEXP);
    Rcpp::traits::input_parameter< double >::type bs(bsSEXP);
    Rcpp::traits::input_parameter< double >::type br(brSEXP);
    Rcpp::traits::input_parameter< double >::type as(asSEXP);
    Rcpp::traits::input_parameter< double >::type ar(arSEXP);
    ccdf_model_collective(x, y, p, bs, br, as, ar);
    return R_NilValue;
END_RCPP
}
// ccdf_model_individual
void ccdf_model_individual(NumericVector x, NumericVector y, double bs, double as);
RcppExport SEXP _wagglefit_ccdf_model_individual(SEXP xSEXP, SEXP ySEXP, SEXP bsSEXP, SEXP asSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< double >::type bs(bsSEXP);
    Rcpp::traits::input_parameter< double >::type as(asSEXP);
    ccdf_model_individual(x, y, bs, as);
    return R_NilValue;
END_RCPP
}
// model_ccdf
NumericVector model_ccdf(NumericVector x, NumericVector params, int model);
RcppExport SEXP _wagglefit_model_ccdf(SEXP xSEXP, SEXP paramsSEXP, SEXP modelSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type params(paramsSEXP);
    Rcpp::traits::input_parameter< int >::type model(modelSEXP);
    rcpp_result_gen = Rcpp::wrap(model_ccdf(x, params, model));
    return rcpp_result_gen;
END_RCPP
}
// scout_dist
double scout_dist(double x, double m, double bs, double as);
RcppExport SEXP _wagglefit_scout_dist(SEXP xSEXP, SEXP mSEXP, SEXP bsSEXP, SEXP asSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type m(mSEXP);
    Rcpp::traits::input_parameter< double >::type bs(bsSEXP);
    Rcpp::traits::input_parameter< double >::type as(asSEXP);
    rcpp_result_gen = Rcpp::wrap(scout_dist(x, m, bs, as));
    return rcpp_result_gen;
END_RCPP
}
// recruit_dist
double recruit_dist(double x, double m, double br, double ar);
RcppExport SEXP _wagglefit_recruit_dist(SEXP xSEXP, SEXP mSEXP, SEXP brSEXP, SEXP arSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type m(mSEXP);
    Rcpp::traits::input_parameter< double >::type br(brSEXP);
    Rcpp::traits::input_parameter< double >::type ar(arSEXP);
    rcpp_result_gen = Rcpp::wrap(recruit_dist(x, m, br, ar));
    return rcpp_result_gen;
END_RCPP
}
// loglike_model_collective
double loglike_model_collective(NumericVector x, double p, double bs, double br, double as, double ar);
RcppExport SEXP _wagglefit_loglike_model_collective(SEXP xSEXP, SEXP pSEXP, SEXP bsSEXP, SEXP brSEXP, SEXP asSEXP, SEXP arSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type p(pSEXP);
    Rcpp::traits::input_parameter< double >::type bs(bsSEXP);
    Rcpp::traits::input_parameter< double >::type br(brSEXP);
    Rcpp::traits::input_parameter< double >::type as(asSEXP);
    Rcpp::traits::input_parameter< double >::type ar(arSEXP);
    rcpp_result_gen = Rcpp::wrap(loglike_model_collective(x, p, bs, br, as, ar));
    return rcpp_result_gen;
END_RCPP
}
// loglike_model_individual
double loglike_model_individual(NumericVector x, double bs, double as);
RcppExport SEXP _wagglefit_loglike_model_individual(SEXP xSEXP, SEXP bsSEXP, SEXP asSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type bs(bsSEXP);
    Rcpp::traits::input_parameter< double >::type as(asSEXP);
    rcpp_result_gen = Rcpp::wrap(loglike_model_individual(x, bs, as));
    return rcpp_result_gen;
END_RCPP
}
// optimise_model
NumericVector optimise_model(NumericVector x, NumericVector params, NumericVector lb, NumericVector ub, bool verbose, double xtol, int model);
RcppExport SEXP _wagglefit_optimise_model(SEXP xSEXP, SEXP paramsSEXP, SEXP lbSEXP, SEXP ubSEXP, SEXP verboseSEXP, SEXP xtolSEXP, SEXP modelSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type params(paramsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type lb(lbSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type ub(ubSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< double >::type xtol(xtolSEXP);
    Rcpp::traits::input_parameter< int >::type model(modelSEXP);
    rcpp_result_gen = Rcpp::wrap(optimise_model(x, params, lb, ub, verbose, xtol, model));
    return rcpp_result_gen;
END_RCPP
}
// pdf_model_all
void pdf_model_all(NumericVector x, NumericVector y, double p, double bs, double br, double as, double ar);
RcppExport SEXP _wagglefit_pdf_model_all(SEXP xSEXP, SEXP ySEXP, SEXP pSEXP, SEXP bsSEXP, SEXP brSEXP, SEXP asSEXP, SEXP arSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< double >::type p(pSEXP);
    Rcpp::traits::input_parameter< double >::type bs(bsSEXP);
    Rcpp::traits::input_parameter< double >::type br(brSEXP);
    Rcpp::traits::input_parameter< double >::type as(asSEXP);
    Rcpp::traits::input_parameter< double >::type ar(arSEXP);
    pdf_model_all(x, y, p, bs, br, as, ar);
    return R_NilValue;
END_RCPP
}
// pdf_model_scout
void pdf_model_scout(NumericVector x, NumericVector y, double bs, double as);
RcppExport SEXP _wagglefit_pdf_model_scout(SEXP xSEXP, SEXP ySEXP, SEXP bsSEXP, SEXP asSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< double >::type bs(bsSEXP);
    Rcpp::traits::input_parameter< double >::type as(asSEXP);
    pdf_model_scout(x, y, bs, as);
    return R_NilValue;
END_RCPP
}
// model_pdf
NumericVector model_pdf(NumericVector x, NumericVector params, int model);
RcppExport SEXP _wagglefit_model_pdf(SEXP xSEXP, SEXP paramsSEXP, SEXP modelSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type params(paramsSEXP);
    Rcpp::traits::input_parameter< int >::type model(modelSEXP);
    rcpp_result_gen = Rcpp::wrap(model_pdf(x, params, model));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_wagglefit_alter_in_place", (DL_FUNC) &_wagglefit_alter_in_place, 1},
    {"_wagglefit_scout_ccdf", (DL_FUNC) &_wagglefit_scout_ccdf, 4},
    {"_wagglefit_recruit_ccdf", (DL_FUNC) &_wagglefit_recruit_ccdf, 4},
    {"_wagglefit_ccdf_model_collective", (DL_FUNC) &_wagglefit_ccdf_model_collective, 7},
    {"_wagglefit_ccdf_model_individual", (DL_FUNC) &_wagglefit_ccdf_model_individual, 4},
    {"_wagglefit_model_ccdf", (DL_FUNC) &_wagglefit_model_ccdf, 3},
    {"_wagglefit_scout_dist", (DL_FUNC) &_wagglefit_scout_dist, 4},
    {"_wagglefit_recruit_dist", (DL_FUNC) &_wagglefit_recruit_dist, 4},
    {"_wagglefit_loglike_model_collective", (DL_FUNC) &_wagglefit_loglike_model_collective, 6},
    {"_wagglefit_loglike_model_individual", (DL_FUNC) &_wagglefit_loglike_model_individual, 3},
    {"_wagglefit_optimise_model", (DL_FUNC) &_wagglefit_optimise_model, 7},
    {"_wagglefit_pdf_model_all", (DL_FUNC) &_wagglefit_pdf_model_all, 7},
    {"_wagglefit_pdf_model_scout", (DL_FUNC) &_wagglefit_pdf_model_scout, 4},
    {"_wagglefit_model_pdf", (DL_FUNC) &_wagglefit_model_pdf, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_wagglefit(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
