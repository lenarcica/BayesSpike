/* ========================================================================== */
/*                                                                            */
/*   EigenSolve.c                                                             */
/*   (c) 2011-2019 Alan Lenarcic                                              */
/*                                                                            */
/*   Test Eigenvalue LaPack capabilities                                      */
/*   This code is really only used for demonstration                          */
/* ========================================================================== */
/******************************************************************************/
//// LICENSE INFO: C CODE
//  This program is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 2 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  A copy of the GNU General Public License is available at
//  https://www.R-project.org/Licenses/
//
/******************************************************************************/
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
