// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

// ComputeRmse
double ComputeRmse(SEXP Q1, SEXP Q2);
RcppExport SEXP tess3r_ComputeRmse(SEXP Q1SEXP, SEXP Q2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type Q1(Q1SEXP);
    Rcpp::traits::input_parameter< SEXP >::type Q2(Q2SEXP);
    rcpp_result_gen = Rcpp::wrap(ComputeRmse(Q1, Q2));
    return rcpp_result_gen;
END_RCPP
}
// ComputeAveragedCrossEntropy
double ComputeAveragedCrossEntropy(SEXP P, SEXP Q);
RcppExport SEXP tess3r_ComputeAveragedCrossEntropy(SEXP PSEXP, SEXP QSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type P(PSEXP);
    Rcpp::traits::input_parameter< SEXP >::type Q(QSEXP);
    rcpp_result_gen = Rcpp::wrap(ComputeAveragedCrossEntropy(P, Q));
    return rcpp_result_gen;
END_RCPP
}
// X2XBin
void X2XBin(const Rcpp::NumericMatrix& X, int ploidy, Eigen::Map<Eigen::MatrixXd>& XBin);
RcppExport SEXP tess3r_X2XBin(SEXP XSEXP, SEXP ploidySEXP, SEXP XBinSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type X(XSEXP);
    Rcpp::traits::input_parameter< int >::type ploidy(ploidySEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd>& >::type XBin(XBinSEXP);
    X2XBin(X, ploidy, XBin);
    return R_NilValue;
END_RCPP
}
// XBin2X
Eigen::MatrixXd XBin2X(const Eigen::Map<Eigen::MatrixXd> XBin, int ploidy);
RcppExport SEXP tess3r_XBin2X(SEXP XBinSEXP, SEXP ploidySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd> >::type XBin(XBinSEXP);
    Rcpp::traits::input_parameter< int >::type ploidy(ploidySEXP);
    rcpp_result_gen = Rcpp::wrap(XBin2X(XBin, ploidy));
    return rcpp_result_gen;
END_RCPP
}
// ComputeHeatKernelWeightSparse
Eigen::SparseMatrix<double> ComputeHeatKernelWeightSparse(const Eigen::Map<Eigen::MatrixXd> coord, double sigma);
RcppExport SEXP tess3r_ComputeHeatKernelWeightSparse(SEXP coordSEXP, SEXP sigmaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd> >::type coord(coordSEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    rcpp_result_gen = Rcpp::wrap(ComputeHeatKernelWeightSparse(coord, sigma));
    return rcpp_result_gen;
END_RCPP
}
// ComputeHeatKernelWeight
Eigen::MatrixXd ComputeHeatKernelWeight(const Eigen::Map<Eigen::MatrixXd> coord, double sigma);
RcppExport SEXP tess3r_ComputeHeatKernelWeight(SEXP coordSEXP, SEXP sigmaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd> >::type coord(coordSEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    rcpp_result_gen = Rcpp::wrap(ComputeHeatKernelWeight(coord, sigma));
    return rcpp_result_gen;
END_RCPP
}
// ComputeExponetialWeight
Eigen::MatrixXd ComputeExponetialWeight(const Eigen::Map<Eigen::MatrixXd> coord, double sigma);
RcppExport SEXP tess3r_ComputeExponetialWeight(SEXP coordSEXP, SEXP sigmaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd> >::type coord(coordSEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    rcpp_result_gen = Rcpp::wrap(ComputeExponetialWeight(coord, sigma));
    return rcpp_result_gen;
END_RCPP
}
// SampleGenoFromGenerativeModelTESS3
Rcpp::List SampleGenoFromGenerativeModelTESS3(const Rcpp::NumericMatrix& Q, const Rcpp::NumericMatrix& G, const Rcpp::NumericMatrix& coord, int ploidy, int openMP_core_num);
RcppExport SEXP tess3r_SampleGenoFromGenerativeModelTESS3(SEXP QSEXP, SEXP GSEXP, SEXP coordSEXP, SEXP ploidySEXP, SEXP openMP_core_numSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type Q(QSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type G(GSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type coord(coordSEXP);
    Rcpp::traits::input_parameter< int >::type ploidy(ploidySEXP);
    Rcpp::traits::input_parameter< int >::type openMP_core_num(openMP_core_numSEXP);
    rcpp_result_gen = Rcpp::wrap(SampleGenoFromGenerativeModelTESS3(Q, G, coord, ploidy, openMP_core_num));
    return rcpp_result_gen;
END_RCPP
}
// ComputeZHelper
Eigen::MatrixXi ComputeZHelper(const Eigen::Map<Eigen::MatrixXd> Q, int n, int L);
RcppExport SEXP tess3r_ComputeZHelper(SEXP QSEXP, SEXP nSEXP, SEXP LSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd> >::type Q(QSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type L(LSEXP);
    rcpp_result_gen = Rcpp::wrap(ComputeZHelper(Q, n, L));
    return rcpp_result_gen;
END_RCPP
}
// ComputeAdmixtedGeno
Eigen::MatrixXi ComputeAdmixtedGeno(const Rcpp::NumericVector& geno, const Eigen::Map<Eigen::MatrixXi> Z, int n, int L);
RcppExport SEXP tess3r_ComputeAdmixtedGeno(SEXP genoSEXP, SEXP ZSEXP, SEXP nSEXP, SEXP LSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type geno(genoSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXi> >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type L(LSEXP);
    rcpp_result_gen = Rcpp::wrap(ComputeAdmixtedGeno(geno, Z, n, L));
    return rcpp_result_gen;
END_RCPP
}
// ComputeFst
Eigen::MatrixXd ComputeFst(const Eigen::Map<Eigen::MatrixXd> Q, const Eigen::Map<Eigen::MatrixXd> G, int D);
RcppExport SEXP tess3r_ComputeFst(SEXP QSEXP, SEXP GSEXP, SEXP DSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd> >::type Q(QSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd> >::type G(GSEXP);
    Rcpp::traits::input_parameter< int >::type D(DSEXP);
    rcpp_result_gen = Rcpp::wrap(ComputeFst(Q, G, D));
    return rcpp_result_gen;
END_RCPP
}
// InitOpenMP
void InitOpenMP(int n);
RcppExport SEXP tess3r_InitOpenMP(SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    InitOpenMP(n);
    return R_NilValue;
END_RCPP
}
// ComputeMCPASolution
void ComputeMCPASolution(const Eigen::Map<Eigen::MatrixXd> X, int K, const Eigen::Map<Eigen::MatrixXd> Lapl, double lambdaPrim, int D, int maxIteration, double tolerance, Eigen::Map<Eigen::MatrixXd> Q, Eigen::Map<Eigen::MatrixXd> G, bool verbose);
RcppExport SEXP tess3r_ComputeMCPASolution(SEXP XSEXP, SEXP KSEXP, SEXP LaplSEXP, SEXP lambdaPrimSEXP, SEXP DSEXP, SEXP maxIterationSEXP, SEXP toleranceSEXP, SEXP QSEXP, SEXP GSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd> >::type X(XSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd> >::type Lapl(LaplSEXP);
    Rcpp::traits::input_parameter< double >::type lambdaPrim(lambdaPrimSEXP);
    Rcpp::traits::input_parameter< int >::type D(DSEXP);
    Rcpp::traits::input_parameter< int >::type maxIteration(maxIterationSEXP);
    Rcpp::traits::input_parameter< double >::type tolerance(toleranceSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type Q(QSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type G(GSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    ComputeMCPASolution(X, K, Lapl, lambdaPrim, D, maxIteration, tolerance, Q, G, verbose);
    return R_NilValue;
END_RCPP
}
// ComputeMCPASolutionNoCopyX
void ComputeMCPASolutionNoCopyX(const Eigen::Map<Eigen::MatrixXd> X, int K, const Eigen::Map<Eigen::MatrixXd> Lapl, double lambdaPrim, int D, int maxIteration, double tolerance, Eigen::Map<Eigen::MatrixXd> Q, Eigen::Map<Eigen::MatrixXd> G, bool verbose);
RcppExport SEXP tess3r_ComputeMCPASolutionNoCopyX(SEXP XSEXP, SEXP KSEXP, SEXP LaplSEXP, SEXP lambdaPrimSEXP, SEXP DSEXP, SEXP maxIterationSEXP, SEXP toleranceSEXP, SEXP QSEXP, SEXP GSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd> >::type X(XSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd> >::type Lapl(LaplSEXP);
    Rcpp::traits::input_parameter< double >::type lambdaPrim(lambdaPrimSEXP);
    Rcpp::traits::input_parameter< int >::type D(DSEXP);
    Rcpp::traits::input_parameter< int >::type maxIteration(maxIterationSEXP);
    Rcpp::traits::input_parameter< double >::type tolerance(toleranceSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type Q(QSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type G(GSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    ComputeMCPASolutionNoCopyX(X, K, Lapl, lambdaPrim, D, maxIteration, tolerance, Q, G, verbose);
    return R_NilValue;
END_RCPP
}