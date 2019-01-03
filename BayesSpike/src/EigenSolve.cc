/* ========================================================================== */
/*                                                                            */
/*   EigenSolve.c                                                             */
/*   (c) 2001 Author                                                          */
/*                                                                            */
/*   Test Eigenvalue LaPack capabilities                                      */
/*                                                                            */
/* ========================================================================== */
#define UpLo 'L'
#define TransWant 'N'
#ifndef RNeeds 
 #include <R.h>
 #include <Rmath.h>
 #include <Rinternals.h>
 #include <R_ext/BLAS.h>
 #include <R_ext/Lapack.h>
 #include <Rdefines.h>
 #define RNeeds 0
#endif

extern "C" {
SEXP LetsGoEigen(SEXP sMatA, SEXP sD, SEXP sE, SEXP sTau,
  SEXP sWORK) {
int N =  INTEGER(GET_DIM(sMatA))[0];
char MUpLo = UpLo;
int INFO = 0;
int llWork = length(sWORK);
F77_CALL(dsytrd)(&MUpLo, &N, REAL(sMatA), &N, REAL(sD), 
  REAL(sE), REAL(sTau), REAL(sWORK), &llWork, &INFO);

  return(sMatA);
}
}
