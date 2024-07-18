// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// DEx_Rcpp
NumericMatrix DEx_Rcpp(int PointNum, int FactorNum, int PosNum, int ParticleNum, int IterationNum, double PValue, double pMut, double pCR, double pGBest, double pSelf, int type, int seed, int maxTime, int maxNoImprove);
RcppExport SEXP _Meta4Design_DEx_Rcpp(SEXP PointNumSEXP, SEXP FactorNumSEXP, SEXP PosNumSEXP, SEXP ParticleNumSEXP, SEXP IterationNumSEXP, SEXP PValueSEXP, SEXP pMutSEXP, SEXP pCRSEXP, SEXP pGBestSEXP, SEXP pSelfSEXP, SEXP typeSEXP, SEXP seedSEXP, SEXP maxTimeSEXP, SEXP maxNoImproveSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type PointNum(PointNumSEXP);
    Rcpp::traits::input_parameter< int >::type FactorNum(FactorNumSEXP);
    Rcpp::traits::input_parameter< int >::type PosNum(PosNumSEXP);
    Rcpp::traits::input_parameter< int >::type ParticleNum(ParticleNumSEXP);
    Rcpp::traits::input_parameter< int >::type IterationNum(IterationNumSEXP);
    Rcpp::traits::input_parameter< double >::type PValue(PValueSEXP);
    Rcpp::traits::input_parameter< double >::type pMut(pMutSEXP);
    Rcpp::traits::input_parameter< double >::type pCR(pCRSEXP);
    Rcpp::traits::input_parameter< double >::type pGBest(pGBestSEXP);
    Rcpp::traits::input_parameter< double >::type pSelf(pSelfSEXP);
    Rcpp::traits::input_parameter< int >::type type(typeSEXP);
    Rcpp::traits::input_parameter< int >::type seed(seedSEXP);
    Rcpp::traits::input_parameter< int >::type maxTime(maxTimeSEXP);
    Rcpp::traits::input_parameter< int >::type maxNoImprove(maxNoImproveSEXP);
    rcpp_result_gen = Rcpp::wrap(DEx_Rcpp(PointNum, FactorNum, PosNum, ParticleNum, IterationNum, PValue, pMut, pCR, pGBest, pSelf, type, seed, maxTime, maxNoImprove));
    return rcpp_result_gen;
END_RCPP
}
// TAx_Rcpp
NumericMatrix TAx_Rcpp(int PointNum, int FactorNum, int PosNum, int ParticleNum, int IterationNum, double PValue, NumericVector vT, int SA, int type, int seed, int maxTime, int maxNoImprove);
RcppExport SEXP _Meta4Design_TAx_Rcpp(SEXP PointNumSEXP, SEXP FactorNumSEXP, SEXP PosNumSEXP, SEXP ParticleNumSEXP, SEXP IterationNumSEXP, SEXP PValueSEXP, SEXP vTSEXP, SEXP SASEXP, SEXP typeSEXP, SEXP seedSEXP, SEXP maxTimeSEXP, SEXP maxNoImproveSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type PointNum(PointNumSEXP);
    Rcpp::traits::input_parameter< int >::type FactorNum(FactorNumSEXP);
    Rcpp::traits::input_parameter< int >::type PosNum(PosNumSEXP);
    Rcpp::traits::input_parameter< int >::type ParticleNum(ParticleNumSEXP);
    Rcpp::traits::input_parameter< int >::type IterationNum(IterationNumSEXP);
    Rcpp::traits::input_parameter< double >::type PValue(PValueSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type vT(vTSEXP);
    Rcpp::traits::input_parameter< int >::type SA(SASEXP);
    Rcpp::traits::input_parameter< int >::type type(typeSEXP);
    Rcpp::traits::input_parameter< int >::type seed(seedSEXP);
    Rcpp::traits::input_parameter< int >::type maxTime(maxTimeSEXP);
    Rcpp::traits::input_parameter< int >::type maxNoImprove(maxNoImproveSEXP);
    rcpp_result_gen = Rcpp::wrap(TAx_Rcpp(PointNum, FactorNum, PosNum, ParticleNum, IterationNum, PValue, vT, SA, type, seed, maxTime, maxNoImprove));
    return rcpp_result_gen;
END_RCPP
}
// GAx_Rcpp
NumericMatrix GAx_Rcpp(int PointNum, int FactorNum, int PosNum, int PopulationNum, int IterationNum, int PValue, double pMut, int type, int seed, int maxTime, int maxNoImprove);
RcppExport SEXP _Meta4Design_GAx_Rcpp(SEXP PointNumSEXP, SEXP FactorNumSEXP, SEXP PosNumSEXP, SEXP PopulationNumSEXP, SEXP IterationNumSEXP, SEXP PValueSEXP, SEXP pMutSEXP, SEXP typeSEXP, SEXP seedSEXP, SEXP maxTimeSEXP, SEXP maxNoImproveSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type PointNum(PointNumSEXP);
    Rcpp::traits::input_parameter< int >::type FactorNum(FactorNumSEXP);
    Rcpp::traits::input_parameter< int >::type PosNum(PosNumSEXP);
    Rcpp::traits::input_parameter< int >::type PopulationNum(PopulationNumSEXP);
    Rcpp::traits::input_parameter< int >::type IterationNum(IterationNumSEXP);
    Rcpp::traits::input_parameter< int >::type PValue(PValueSEXP);
    Rcpp::traits::input_parameter< double >::type pMut(pMutSEXP);
    Rcpp::traits::input_parameter< int >::type type(typeSEXP);
    Rcpp::traits::input_parameter< int >::type seed(seedSEXP);
    Rcpp::traits::input_parameter< int >::type maxTime(maxTimeSEXP);
    Rcpp::traits::input_parameter< int >::type maxNoImprove(maxNoImproveSEXP);
    rcpp_result_gen = Rcpp::wrap(GAx_Rcpp(PointNum, FactorNum, PosNum, PopulationNum, IterationNum, PValue, pMut, type, seed, maxTime, maxNoImprove));
    return rcpp_result_gen;
END_RCPP
}
// SRSx_Rcpp
NumericMatrix SRSx_Rcpp(int PointNum, int FactorNum, int PosNum, int PopulationNum, int IterationNum, int PValue, int type, int seed, int maxTime, int maxNoImprove);
RcppExport SEXP _Meta4Design_SRSx_Rcpp(SEXP PointNumSEXP, SEXP FactorNumSEXP, SEXP PosNumSEXP, SEXP PopulationNumSEXP, SEXP IterationNumSEXP, SEXP PValueSEXP, SEXP typeSEXP, SEXP seedSEXP, SEXP maxTimeSEXP, SEXP maxNoImproveSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type PointNum(PointNumSEXP);
    Rcpp::traits::input_parameter< int >::type FactorNum(FactorNumSEXP);
    Rcpp::traits::input_parameter< int >::type PosNum(PosNumSEXP);
    Rcpp::traits::input_parameter< int >::type PopulationNum(PopulationNumSEXP);
    Rcpp::traits::input_parameter< int >::type IterationNum(IterationNumSEXP);
    Rcpp::traits::input_parameter< int >::type PValue(PValueSEXP);
    Rcpp::traits::input_parameter< int >::type type(typeSEXP);
    Rcpp::traits::input_parameter< int >::type seed(seedSEXP);
    Rcpp::traits::input_parameter< int >::type maxTime(maxTimeSEXP);
    Rcpp::traits::input_parameter< int >::type maxNoImprove(maxNoImproveSEXP);
    rcpp_result_gen = Rcpp::wrap(SRSx_Rcpp(PointNum, FactorNum, PosNum, PopulationNum, IterationNum, PValue, type, seed, maxTime, maxNoImprove));
    return rcpp_result_gen;
END_RCPP
}
// GeneticAlg_YY
NumericMatrix GeneticAlg_YY(int n, int m, int s, double gmp_beta, int popSize, double elitism_ratio, double crossover_ratio, double mutation_prob, int maxiter, int type, bool verbose);
RcppExport SEXP _Meta4Design_GeneticAlg_YY(SEXP nSEXP, SEXP mSEXP, SEXP sSEXP, SEXP gmp_betaSEXP, SEXP popSizeSEXP, SEXP elitism_ratioSEXP, SEXP crossover_ratioSEXP, SEXP mutation_probSEXP, SEXP maxiterSEXP, SEXP typeSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    Rcpp::traits::input_parameter< int >::type s(sSEXP);
    Rcpp::traits::input_parameter< double >::type gmp_beta(gmp_betaSEXP);
    Rcpp::traits::input_parameter< int >::type popSize(popSizeSEXP);
    Rcpp::traits::input_parameter< double >::type elitism_ratio(elitism_ratioSEXP);
    Rcpp::traits::input_parameter< double >::type crossover_ratio(crossover_ratioSEXP);
    Rcpp::traits::input_parameter< double >::type mutation_prob(mutation_probSEXP);
    Rcpp::traits::input_parameter< int >::type maxiter(maxiterSEXP);
    Rcpp::traits::input_parameter< int >::type type(typeSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(GeneticAlg_YY(n, m, s, gmp_beta, popSize, elitism_ratio, crossover_ratio, mutation_prob, maxiter, type, verbose));
    return rcpp_result_gen;
END_RCPP
}
// GeneticAlg_base
NumericMatrix GeneticAlg_base(int n, int m, int s, double beta, int popSize, double pMut, int maxiter, int type, bool verbose);
RcppExport SEXP _Meta4Design_GeneticAlg_base(SEXP nSEXP, SEXP mSEXP, SEXP sSEXP, SEXP betaSEXP, SEXP popSizeSEXP, SEXP pMutSEXP, SEXP maxiterSEXP, SEXP typeSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    Rcpp::traits::input_parameter< int >::type s(sSEXP);
    Rcpp::traits::input_parameter< double >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< int >::type popSize(popSizeSEXP);
    Rcpp::traits::input_parameter< double >::type pMut(pMutSEXP);
    Rcpp::traits::input_parameter< int >::type maxiter(maxiterSEXP);
    Rcpp::traits::input_parameter< int >::type type(typeSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(GeneticAlg_base(n, m, s, beta, popSize, pMut, maxiter, type, verbose));
    return rcpp_result_gen;
END_RCPP
}
// test
double test(int a);
RcppExport SEXP _Meta4Design_test(SEXP aSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type a(aSEXP);
    rcpp_result_gen = Rcpp::wrap(test(a));
    return rcpp_result_gen;
END_RCPP
}
// PSOx_Rcpp
NumericMatrix PSOx_Rcpp(int PointNum, int FactorNum, int PosNum, int ParticleNum, int IterationNum, double PValue, int SameNumP, int SameNumG, double pMut, int SwapNumR, double pGBest, double pPBest, int type, int seed, int maxTime, int maxNoImprove);
RcppExport SEXP _Meta4Design_PSOx_Rcpp(SEXP PointNumSEXP, SEXP FactorNumSEXP, SEXP PosNumSEXP, SEXP ParticleNumSEXP, SEXP IterationNumSEXP, SEXP PValueSEXP, SEXP SameNumPSEXP, SEXP SameNumGSEXP, SEXP pMutSEXP, SEXP SwapNumRSEXP, SEXP pGBestSEXP, SEXP pPBestSEXP, SEXP typeSEXP, SEXP seedSEXP, SEXP maxTimeSEXP, SEXP maxNoImproveSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type PointNum(PointNumSEXP);
    Rcpp::traits::input_parameter< int >::type FactorNum(FactorNumSEXP);
    Rcpp::traits::input_parameter< int >::type PosNum(PosNumSEXP);
    Rcpp::traits::input_parameter< int >::type ParticleNum(ParticleNumSEXP);
    Rcpp::traits::input_parameter< int >::type IterationNum(IterationNumSEXP);
    Rcpp::traits::input_parameter< double >::type PValue(PValueSEXP);
    Rcpp::traits::input_parameter< int >::type SameNumP(SameNumPSEXP);
    Rcpp::traits::input_parameter< int >::type SameNumG(SameNumGSEXP);
    Rcpp::traits::input_parameter< double >::type pMut(pMutSEXP);
    Rcpp::traits::input_parameter< int >::type SwapNumR(SwapNumRSEXP);
    Rcpp::traits::input_parameter< double >::type pGBest(pGBestSEXP);
    Rcpp::traits::input_parameter< double >::type pPBest(pPBestSEXP);
    Rcpp::traits::input_parameter< int >::type type(typeSEXP);
    Rcpp::traits::input_parameter< int >::type seed(seedSEXP);
    Rcpp::traits::input_parameter< int >::type maxTime(maxTimeSEXP);
    Rcpp::traits::input_parameter< int >::type maxNoImprove(maxNoImproveSEXP);
    rcpp_result_gen = Rcpp::wrap(PSOx_Rcpp(PointNum, FactorNum, PosNum, ParticleNum, IterationNum, PValue, SameNumP, SameNumG, pMut, SwapNumR, pGBest, pPBest, type, seed, maxTime, maxNoImprove));
    return rcpp_result_gen;
END_RCPP
}
// Move2Best_Rcpp
IntegerVector Move2Best_Rcpp(int PosNum, int PointNum, int SameNum, IntegerVector vParticle, IntegerVector vBest);
RcppExport SEXP _Meta4Design_Move2Best_Rcpp(SEXP PosNumSEXP, SEXP PointNumSEXP, SEXP SameNumSEXP, SEXP vParticleSEXP, SEXP vBestSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type PosNum(PosNumSEXP);
    Rcpp::traits::input_parameter< int >::type PointNum(PointNumSEXP);
    Rcpp::traits::input_parameter< int >::type SameNum(SameNumSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type vParticle(vParticleSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type vBest(vBestSEXP);
    rcpp_result_gen = Rcpp::wrap(Move2Best_Rcpp(PosNum, PointNum, SameNum, vParticle, vBest));
    return rcpp_result_gen;
END_RCPP
}
// gmp_log
double gmp_log(NumericMatrix D, double beta);
RcppExport SEXP _Meta4Design_gmp_log(SEXP DSEXP, SEXP betaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type D(DSEXP);
    Rcpp::traits::input_parameter< double >::type beta(betaSEXP);
    rcpp_result_gen = Rcpp::wrap(gmp_log(D, beta));
    return rcpp_result_gen;
END_RCPP
}
// gmp
double gmp(NumericMatrix D, double beta);
RcppExport SEXP _Meta4Design_gmp(SEXP DSEXP, SEXP betaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type D(DSEXP);
    Rcpp::traits::input_parameter< double >::type beta(betaSEXP);
    rcpp_result_gen = Rcpp::wrap(gmp(D, beta));
    return rcpp_result_gen;
END_RCPP
}
// gmp_diff
double gmp_diff(NumericMatrix D, int i0, int j0, int k0, double beta);
RcppExport SEXP _Meta4Design_gmp_diff(SEXP DSEXP, SEXP i0SEXP, SEXP j0SEXP, SEXP k0SEXP, SEXP betaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type D(DSEXP);
    Rcpp::traits::input_parameter< int >::type i0(i0SEXP);
    Rcpp::traits::input_parameter< int >::type j0(j0SEXP);
    Rcpp::traits::input_parameter< int >::type k0(k0SEXP);
    Rcpp::traits::input_parameter< double >::type beta(betaSEXP);
    rcpp_result_gen = Rcpp::wrap(gmp_diff(D, i0, j0, k0, beta));
    return rcpp_result_gen;
END_RCPP
}
// phip_log
double phip_log(NumericMatrix D, double p);
RcppExport SEXP _Meta4Design_phip_log(SEXP DSEXP, SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type D(DSEXP);
    Rcpp::traits::input_parameter< double >::type p(pSEXP);
    rcpp_result_gen = Rcpp::wrap(phip_log(D, p));
    return rcpp_result_gen;
END_RCPP
}
// phip
double phip(NumericMatrix D, double p);
RcppExport SEXP _Meta4Design_phip(SEXP DSEXP, SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type D(DSEXP);
    Rcpp::traits::input_parameter< double >::type p(pSEXP);
    rcpp_result_gen = Rcpp::wrap(phip(D, p));
    return rcpp_result_gen;
END_RCPP
}
// phip_diff
double phip_diff(NumericMatrix D, int i0, int j0, int k0, double p);
RcppExport SEXP _Meta4Design_phip_diff(SEXP DSEXP, SEXP i0SEXP, SEXP j0SEXP, SEXP k0SEXP, SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type D(DSEXP);
    Rcpp::traits::input_parameter< int >::type i0(i0SEXP);
    Rcpp::traits::input_parameter< int >::type j0(j0SEXP);
    Rcpp::traits::input_parameter< int >::type k0(k0SEXP);
    Rcpp::traits::input_parameter< double >::type p(pSEXP);
    rcpp_result_gen = Rcpp::wrap(phip_diff(D, i0, j0, k0, p));
    return rcpp_result_gen;
END_RCPP
}
// crt_log
double crt_log(NumericMatrix D, double beta, int type);
RcppExport SEXP _Meta4Design_crt_log(SEXP DSEXP, SEXP betaSEXP, SEXP typeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type D(DSEXP);
    Rcpp::traits::input_parameter< double >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< int >::type type(typeSEXP);
    rcpp_result_gen = Rcpp::wrap(crt_log(D, beta, type));
    return rcpp_result_gen;
END_RCPP
}
// crt_diff
double crt_diff(NumericMatrix D, int i0, int j0, int k0, double beta, int type);
RcppExport SEXP _Meta4Design_crt_diff(SEXP DSEXP, SEXP i0SEXP, SEXP j0SEXP, SEXP k0SEXP, SEXP betaSEXP, SEXP typeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type D(DSEXP);
    Rcpp::traits::input_parameter< int >::type i0(i0SEXP);
    Rcpp::traits::input_parameter< int >::type j0(j0SEXP);
    Rcpp::traits::input_parameter< int >::type k0(k0SEXP);
    Rcpp::traits::input_parameter< double >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< int >::type type(typeSEXP);
    rcpp_result_gen = Rcpp::wrap(crt_diff(D, i0, j0, k0, beta, type));
    return rcpp_result_gen;
END_RCPP
}
// generate_rbdesign
NumericMatrix generate_rbdesign(int n, int m, int s);
RcppExport SEXP _Meta4Design_generate_rbdesign(SEXP nSEXP, SEXP mSEXP, SEXP sSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    Rcpp::traits::input_parameter< int >::type s(sSEXP);
    rcpp_result_gen = Rcpp::wrap(generate_rbdesign(n, m, s));
    return rcpp_result_gen;
END_RCPP
}
// SA_XY
NumericMatrix SA_XY(int n, int m, int s, double beta, int Ntrials, int itermax, int maxTime, double T0, int Npairs, int nFail, int type, bool verbose);
RcppExport SEXP _Meta4Design_SA_XY(SEXP nSEXP, SEXP mSEXP, SEXP sSEXP, SEXP betaSEXP, SEXP NtrialsSEXP, SEXP itermaxSEXP, SEXP maxTimeSEXP, SEXP T0SEXP, SEXP NpairsSEXP, SEXP nFailSEXP, SEXP typeSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    Rcpp::traits::input_parameter< int >::type s(sSEXP);
    Rcpp::traits::input_parameter< double >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< int >::type Ntrials(NtrialsSEXP);
    Rcpp::traits::input_parameter< int >::type itermax(itermaxSEXP);
    Rcpp::traits::input_parameter< int >::type maxTime(maxTimeSEXP);
    Rcpp::traits::input_parameter< double >::type T0(T0SEXP);
    Rcpp::traits::input_parameter< int >::type Npairs(NpairsSEXP);
    Rcpp::traits::input_parameter< int >::type nFail(nFailSEXP);
    Rcpp::traits::input_parameter< int >::type type(typeSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(SA_XY(n, m, s, beta, Ntrials, itermax, maxTime, T0, Npairs, nFail, type, verbose));
    return rcpp_result_gen;
END_RCPP
}
// SA_MM
NumericMatrix SA_MM(int n, int m, int s, double beta, int Ntrials, int itermax, int maxTime, double T0, int Imax, int nFail, int type, bool verbose);
RcppExport SEXP _Meta4Design_SA_MM(SEXP nSEXP, SEXP mSEXP, SEXP sSEXP, SEXP betaSEXP, SEXP NtrialsSEXP, SEXP itermaxSEXP, SEXP maxTimeSEXP, SEXP T0SEXP, SEXP ImaxSEXP, SEXP nFailSEXP, SEXP typeSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    Rcpp::traits::input_parameter< int >::type s(sSEXP);
    Rcpp::traits::input_parameter< double >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< int >::type Ntrials(NtrialsSEXP);
    Rcpp::traits::input_parameter< int >::type itermax(itermaxSEXP);
    Rcpp::traits::input_parameter< int >::type maxTime(maxTimeSEXP);
    Rcpp::traits::input_parameter< double >::type T0(T0SEXP);
    Rcpp::traits::input_parameter< int >::type Imax(ImaxSEXP);
    Rcpp::traits::input_parameter< int >::type nFail(nFailSEXP);
    Rcpp::traits::input_parameter< int >::type type(typeSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(SA_MM(n, m, s, beta, Ntrials, itermax, maxTime, T0, Imax, nFail, type, verbose));
    return rcpp_result_gen;
END_RCPP
}
// bid_log
double bid_log(NumericMatrix D, double beta);
RcppExport SEXP _Meta4Design_bid_log(SEXP DSEXP, SEXP betaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type D(DSEXP);
    Rcpp::traits::input_parameter< double >::type beta(betaSEXP);
    rcpp_result_gen = Rcpp::wrap(bid_log(D, beta));
    return rcpp_result_gen;
END_RCPP
}
// bid
double bid(NumericMatrix D, double beta);
RcppExport SEXP _Meta4Design_bid(SEXP DSEXP, SEXP betaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type D(DSEXP);
    Rcpp::traits::input_parameter< double >::type beta(betaSEXP);
    rcpp_result_gen = Rcpp::wrap(bid(D, beta));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_Meta4Design_DEx_Rcpp", (DL_FUNC) &_Meta4Design_DEx_Rcpp, 14},
    {"_Meta4Design_TAx_Rcpp", (DL_FUNC) &_Meta4Design_TAx_Rcpp, 12},
    {"_Meta4Design_GAx_Rcpp", (DL_FUNC) &_Meta4Design_GAx_Rcpp, 11},
    {"_Meta4Design_SRSx_Rcpp", (DL_FUNC) &_Meta4Design_SRSx_Rcpp, 10},
    {"_Meta4Design_GeneticAlg_YY", (DL_FUNC) &_Meta4Design_GeneticAlg_YY, 11},
    {"_Meta4Design_GeneticAlg_base", (DL_FUNC) &_Meta4Design_GeneticAlg_base, 9},
    {"_Meta4Design_test", (DL_FUNC) &_Meta4Design_test, 1},
    {"_Meta4Design_PSOx_Rcpp", (DL_FUNC) &_Meta4Design_PSOx_Rcpp, 16},
    {"_Meta4Design_Move2Best_Rcpp", (DL_FUNC) &_Meta4Design_Move2Best_Rcpp, 5},
    {"_Meta4Design_gmp_log", (DL_FUNC) &_Meta4Design_gmp_log, 2},
    {"_Meta4Design_gmp", (DL_FUNC) &_Meta4Design_gmp, 2},
    {"_Meta4Design_gmp_diff", (DL_FUNC) &_Meta4Design_gmp_diff, 5},
    {"_Meta4Design_phip_log", (DL_FUNC) &_Meta4Design_phip_log, 2},
    {"_Meta4Design_phip", (DL_FUNC) &_Meta4Design_phip, 2},
    {"_Meta4Design_phip_diff", (DL_FUNC) &_Meta4Design_phip_diff, 5},
    {"_Meta4Design_crt_log", (DL_FUNC) &_Meta4Design_crt_log, 3},
    {"_Meta4Design_crt_diff", (DL_FUNC) &_Meta4Design_crt_diff, 6},
    {"_Meta4Design_generate_rbdesign", (DL_FUNC) &_Meta4Design_generate_rbdesign, 3},
    {"_Meta4Design_SA_XY", (DL_FUNC) &_Meta4Design_SA_XY, 12},
    {"_Meta4Design_SA_MM", (DL_FUNC) &_Meta4Design_SA_MM, 12},
    {"_Meta4Design_bid_log", (DL_FUNC) &_Meta4Design_bid_log, 2},
    {"_Meta4Design_bid", (DL_FUNC) &_Meta4Design_bid, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_Meta4Design(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
