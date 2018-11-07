// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// CrossCorr
double CrossCorr(arma::vec TestSignal, arma::vec TemplateSignal);
RcppExport SEXP _SpecMatch_CrossCorr(SEXP TestSignalSEXP, SEXP TemplateSignalSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type TestSignal(TestSignalSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type TemplateSignal(TemplateSignalSEXP);
    rcpp_result_gen = Rcpp::wrap(CrossCorr(TestSignal, TemplateSignal));
    return rcpp_result_gen;
END_RCPP
}
// test
double test(arma::vec hi);
RcppExport SEXP _SpecMatch_test(SEXP hiSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type hi(hiSEXP);
    rcpp_result_gen = Rcpp::wrap(test(hi));
    return rcpp_result_gen;
END_RCPP
}
// DeterStartEnd
arma::vec DeterStartEnd(arma::vec TestSignalWave, arma::vec TemplateSignalWave);
RcppExport SEXP _SpecMatch_DeterStartEnd(SEXP TestSignalWaveSEXP, SEXP TemplateSignalWaveSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type TestSignalWave(TestSignalWaveSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type TemplateSignalWave(TemplateSignalWaveSEXP);
    rcpp_result_gen = Rcpp::wrap(DeterStartEnd(TestSignalWave, TemplateSignalWave));
    return rcpp_result_gen;
END_RCPP
}
// StarletWT
arma::mat StarletWT(arma::vec X);
RcppExport SEXP _SpecMatch_StarletWT(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(StarletWT(X));
    return rcpp_result_gen;
END_RCPP
}
// StarletRC
arma::mat StarletRC(arma::mat X);
RcppExport SEXP _SpecMatch_StarletRC(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(StarletRC(X));
    return rcpp_result_gen;
END_RCPP
}
// PMT
Rcpp::List PMT(arma::vec X);
RcppExport SEXP _SpecMatch_PMT(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(PMT(X));
    return rcpp_result_gen;
END_RCPP
}
// PMTRC
arma::mat PMTRC(arma::mat C, arma::vec V);
RcppExport SEXP _SpecMatch_PMTRC(SEXP CSEXP, SEXP VSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type C(CSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type V(VSEXP);
    rcpp_result_gen = Rcpp::wrap(PMTRC(C, V));
    return rcpp_result_gen;
END_RCPP
}
// MAD
double MAD(arma::rowvec X);
RcppExport SEXP _SpecMatch_MAD(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::rowvec >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(MAD(X));
    return rcpp_result_gen;
END_RCPP
}
// MultiResSuppStarlet
arma::vec MultiResSuppStarlet(arma::vec X, int NoSimu);
RcppExport SEXP _SpecMatch_MultiResSuppStarlet(SEXP XSEXP, SEXP NoSimuSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type X(XSEXP);
    Rcpp::traits::input_parameter< int >::type NoSimu(NoSimuSEXP);
    rcpp_result_gen = Rcpp::wrap(MultiResSuppStarlet(X, NoSimu));
    return rcpp_result_gen;
END_RCPP
}
// HardThreshold
arma::mat HardThreshold(arma::mat C, arma::rowvec M);
RcppExport SEXP _SpecMatch_HardThreshold(SEXP CSEXP, SEXP MSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type C(CSEXP);
    Rcpp::traits::input_parameter< arma::rowvec >::type M(MSEXP);
    rcpp_result_gen = Rcpp::wrap(HardThreshold(C, M));
    return rcpp_result_gen;
END_RCPP
}
// sign
int sign(double x);
RcppExport SEXP _SpecMatch_sign(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(sign(x));
    return rcpp_result_gen;
END_RCPP
}
// SoftThreshold
arma::mat SoftThreshold(arma::mat C, arma::vec M);
RcppExport SEXP _SpecMatch_SoftThreshold(SEXP CSEXP, SEXP MSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type C(CSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type M(MSEXP);
    rcpp_result_gen = Rcpp::wrap(SoftThreshold(C, M));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_SpecMatch_CrossCorr", (DL_FUNC) &_SpecMatch_CrossCorr, 2},
    {"_SpecMatch_test", (DL_FUNC) &_SpecMatch_test, 1},
    {"_SpecMatch_DeterStartEnd", (DL_FUNC) &_SpecMatch_DeterStartEnd, 2},
    {"_SpecMatch_StarletWT", (DL_FUNC) &_SpecMatch_StarletWT, 1},
    {"_SpecMatch_StarletRC", (DL_FUNC) &_SpecMatch_StarletRC, 1},
    {"_SpecMatch_PMT", (DL_FUNC) &_SpecMatch_PMT, 1},
    {"_SpecMatch_PMTRC", (DL_FUNC) &_SpecMatch_PMTRC, 2},
    {"_SpecMatch_MAD", (DL_FUNC) &_SpecMatch_MAD, 1},
    {"_SpecMatch_MultiResSuppStarlet", (DL_FUNC) &_SpecMatch_MultiResSuppStarlet, 2},
    {"_SpecMatch_HardThreshold", (DL_FUNC) &_SpecMatch_HardThreshold, 2},
    {"_SpecMatch_sign", (DL_FUNC) &_SpecMatch_sign, 1},
    {"_SpecMatch_SoftThreshold", (DL_FUNC) &_SpecMatch_SoftThreshold, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_SpecMatch(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
