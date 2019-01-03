
#define MaxTruncTrieds 100

#define Cq2Pi 0.3989423

#define CutCut  1.2534

#define MaxSliceLoop 16

//#ifndef RCPPH
//  #include <Rcpp.h>
//  #include <cmath>
//  #define RCPPH 0
//#endif

#ifndef RNeeds 
 #define RNeeds 0
 #include <R.h>
 #include <Rmath.h>
 #include <Rinternals.h>
 #include <R_ext/BLAS.h>
 #include <R_ext/Lapack.h>
 #include <Rdefines.h>
#endif


//SEXP rTNorm(SEXP sRet, SEXP smu, SEXP ssig, SEXP sL, SEXP sU);
inline double rTruncExpo(double lambda, double L, double U);
inline double rInsideHalfNorm(double U);
inline double rTailTruncHalfNorm(double L, double U);
double rTruncNorm2(double L, double U);
double rTruncNorm(double mu, double sig, double L, double U);
double rTruncT2(double nu, double L, double U);
inline double rInsideHalfT(double nu, double U);
inline double rTruncPareto(double nu, double lambda, double L, double U);
inline double rTailTruncHalfT(double nu, double L, double U);
inline double FindTLevel(double lfdens, double nu);


extern "C" {
SEXP rTT(SEXP sRet, SEXP snu, SEXP smu, SEXP sig, SEXP sL, SEXP sU);
}
