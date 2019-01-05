/******************************************************************************/
/*                                                                            */
/*  HeapSort.cc                                                               */
/*   Adapted 2013 by Alan Leanarcic                                           */
/*   Adapted code based upon Numerical recipes Heap Sort                      */
//  One version of these recipes found at: 
// https://www.cec.uchile.cl/cinetica/pcordero/MC_libros/NumericalRecipesinC.pdf
/*                                                                            */
/*                                                                            */
#ifndef HEAPSORTH
  #define HEAPSORTH 1
  #include "HeapSort.h"
#endif
extern "C" {
SEXP HeapSortShell(SEXP SDoubles, SEXP OutII) {
  int *Out = HeapSortAll(length(SDoubles), REAL(SDoubles));
  int ii;
  if (isReal(OutII)) {
   for (ii = 0; ii < length(OutII); ii++) {
     REAL(OutII)[ii] = (double) Out[ii];
   }
  } else if (isInteger(OutII)) {
   for (ii = 0; ii < length(OutII); ii++) {
     INTEGER(OutII)[ii] = (int) Out[ii];
   }  
  }
  Free(Out);
  return(OutII);  
}
}


int FindInMe(int LengthList, double *MyDoubles, double FindMe) {
  int iSmall = 0;  int iLarge = LengthList-1;
  int MoveSize  = (int) floor(LengthList / 2);

  if (LengthList  == 1) { return(0); }
  if (LengthList == 2)  {
    if (fabs(MyDoubles[1] - FindMe) < fabs(MyDoubles[0]-FindMe)) {
      return(1);
    } else { return(0); }
  }
  if (LengthList == 3) {
    if (MyDoubles[1] < FindMe) {
      if (fabs(MyDoubles[2] - FindMe) <
        fabs(MyDoubles[1] - FindMe)) { return(2); } else { return(1); }
    } else {
       if (fabs(MyDoubles[1] - FindMe) <
         fabs(MyDoubles[0] - FindMe) ) { return(1); } else {return(0); }
    }
  }
  if (MyDoubles[iLarge] <= FindMe) {
    return(iLarge);
  }
  if (MyDoubles[iSmall] >= FindMe) {
    return(iSmall);
  }
  int iProp;
  while(iLarge > iSmall +1) {
    iProp = (int) (floor( (iLarge - iSmall) / 2) + iSmall);
    if (MyDoubles[iProp] == FindMe) {
      return(iProp);
    } else if (MyDoubles[iProp] < FindMe) {
      iSmall = iProp;
    } else  {
      iLarge = iProp;
    }
  }
  if (fabs(MyDoubles[iLarge] - FindMe) < fabs(MyDoubles[iSmall] - FindMe)) {
    return(iLarge);
  } else {
    return(iSmall);
  }
  return(iLarge);
}

////////////////////////////////////////////////////////////////////////////////
// HeapSortAll
//
// This code is adapted from Numerical Recipes in C:
// https://www.cec.uchile.cl/cinetica/pcordero/MC_libros/NumericalRecipesinC.pdf
//
// Chapter 8. Sorting, page 335-336, hpsort algorithm
//
// It has been modified to return sort indexes instead of sort the actual Doubles array.
int *HeapSortAll(int N, double *Doubles) {
int *OutIntii = (int*)Calloc(N, int);
unsigned long i=0, ir, j, l;
double rra = 0.0;

int ii = 0;  int rraii = 0;
for (ii = 0; ii < N; ii++) { OutIntii[ii] = ii; }

if (N < 2) return( OutIntii );
l = (N >> 1) +1;  ir = N;
for (;;) {
  if (l > 1) {
    l--;
    rra = Doubles[ OutIntii[l-1] ];  rraii = OutIntii[l-1];
  } else {
    rra =  Doubles[ OutIntii[ir-1] ];  rraii = OutIntii[ir-1];
    OutIntii[ir-1] = OutIntii[0];
    if (--ir == 1) {
      OutIntii[0] = rraii;     // Least Competent worker of all?
      break;
    }
  }
  i=l; j = l + l;
  while( j <= ir) {
    if (j < ir && Doubles[OutIntii[j-1]] < Doubles[OutIntii[j]]) j++;
    if (rra < Doubles[OutIntii[j-1]]) {
      OutIntii[i-1] = OutIntii[j-1];
      i = j;
      j <<= 1;
    }  else {
      break;
    } 
  }
  OutIntii[i-1] = rraii;
}
return(OutIntii);
}

////////////////////////////////////////////////////////////////////////////////
// HeapSortAll
//
// This code is adapted from Numerical Recipes in C:
// https://www.cec.uchile.cl/cinetica/pcordero/MC_libros/NumericalRecipesinC.pdf
//
// Chapter 8. Sorting, page 335-336, hpsort algorithm
//
// It has been modified to return sort indexes instead of sort the actual Doubles array.
int *HeapSortAll(int N, int *Doubles) {
int *OutIntii = (int*)Calloc(N, int);
unsigned long i=0, ir, j, l;
double rra = 0.0;

int ii = 0;  int rraii = 0;
for (ii = 0; ii < N; ii++) { OutIntii[ii] = ii; }

if (N < 2) return( OutIntii );
l = (N >> 1) +1;  ir = N;
for (;;) {
  if (l > 1) {
    l--;
    rra = Doubles[ OutIntii[l-1] ];  rraii = OutIntii[l-1];
  } else {
    rra =  Doubles[ OutIntii[ir-1] ];  rraii = OutIntii[ir-1];
    OutIntii[ir-1] = OutIntii[0];
    if (--ir == 1) {
      OutIntii[0] = rraii;     // Least Competent worker of all?
      break;
    }
  }
  i=l; j = l + l;
  while( j <= ir) {
    if (j < ir && Doubles[OutIntii[j-1]] < Doubles[OutIntii[j]]) j++;
    if (rra < Doubles[OutIntii[j-1]]) {
      OutIntii[i-1] = OutIntii[j-1];
      i = j;
      j <<= 1;
    }  else {
      break;
    } 
  }
  OutIntii[i-1] = rraii;
}
return(OutIntii);
}

////////////////////////////////////////////////////////////////////////////////
// HeapSortGot
//
// This code is adapted from Numerical Recipes in C:
// https://www.cec.uchile.cl/cinetica/pcordero/MC_libros/NumericalRecipesinC.pdf
//
// Chapter 8. Sorting, page 335-336, hpsort algorithm
//
// It has been modified to return sort indexes instead of sort the actual Doubles array.
// In this case, the array for "OutIntii" is suplied by the user of the function.
int HeapSortGot(int N, int *OutIntii, double *Doubles) {
//int *OutIntii = (int*)Calloc(N, int);
unsigned long i=0, ir, j, l;
double rra = 0.0;

int ii = 0;  int rraii = 0;
for (ii = 0; ii < N; ii++) { OutIntii[ii] = ii; }

if (N < 2) return( 1 );
l = (N >> 1) +1;  ir = N;
for (;;) {
  if (l > 1) {
    l--;
    rra = Doubles[ OutIntii[l-1] ];  rraii = OutIntii[l-1];
  } else {
    rra =  Doubles[ OutIntii[ir-1] ];  rraii = OutIntii[ir-1];
    OutIntii[ir-1] = OutIntii[0];
    if (--ir == 1) {
      OutIntii[0] = rraii;     // Least Competent worker of all?
      break;
    }
  }
  i=l; j = l + l;
  while( j <= ir) {
    if (j < ir && Doubles[OutIntii[j-1]] < Doubles[OutIntii[j]]) j++;
    if (rra < Doubles[OutIntii[j-1]]) {
      OutIntii[i-1] = OutIntii[j-1];
      i = j;
      j <<= 1;
    }  else {
      break;
    } 
  }
  OutIntii[i-1] = rraii;
}
return(1);
}

// ReversString
//  Just a bad algorithm to revserse an integer string
int ReverseString(int N, int *ReverseScores) {
  int HalfN = ((int) floor(N/2.0));
  int ii;
  for (ii = 0; ii < HalfN; ii++) {
    if (N-ii-1 >= N) {
      Rprintf("ReverseString: Absoulte failure N-ii-1 = ", N-ii-1);
      R_FlushConsole();
    }
    ReverseScores[ii] += ReverseScores[N-ii-1];
    ReverseScores[N-ii-1] = ReverseScores[ii] - ReverseScores[N-ii-1];
    ReverseScores[ii] = ReverseScores[ii] - ReverseScores[N-ii-1];  
  }
  return(1);
}
