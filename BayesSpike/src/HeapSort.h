/* ========================================================================== */
/*                                                                            */
/*   HeapSort.h                                                               */
/*   (c) 2013 Alan Lenarcic                                                   */
/*                                                                            */
/*  Based on Numerical Recipes in C++ Book                                    */
/*   Experiment in sorting XTResid in order                                   */
/*                                                                            */
/*   Adapted code based upon Numerical recipes Heap Sort                      */
/*  One version of these recipes found at:                                    */
// https://www.cec.uchile.cl/cinetica/pcordero/MC_libros/NumericalRecipesinC.pdf                                                                        */
/*                                                                            */
/*                                                                            */
/* ========================================================================== */


#ifndef RMATH
  #include <Rmath.h>
  #include <R.h>
  #include <Rdefines.h>
  #include <Rinternals.h>
  #include <stdio.h>
  #define RMATH 0
#endif

int *HeapSortAll(int N, double *Doubles);
int *HeapSortAll(int N, int *Doubles);
int HeapSortGot(int N, int *OutIntii, double *Doubles);
int ReverseString(int N, int *ReverseScores);

int FindInMe(int LengthList, double *MyDoubles, double FindMe);
