/* ========================================================================== */
/*                                                                            */
/*   BayesSpikeHelpers.c                                                               */
/*   (c) 2011 Alan Lenarcic                                                         */
/*                                                                            */
/*   Priors and Prior Derivatives for Bayes Spike                                                              */
/*                                                                            */
//  The function lFofX is the f_2(tau_k^2) density.  Container
//  MT which is of type "TauOfFContainer" contains all of the necessary
//  prior information to eigenrotate into the correct f_2(tau_k) form, and
//  then calculate lFofX and it's derivatives.
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

#ifndef BAYESSPIKECPPH
  #include "BayesSpikeCpp.h"
  #define BAYESSPIKECPPH 0
#endif

//#ifndef AllInclusiveSelectMCMCObject
//  #include "AllInclusiveSelectMCMCObject.h"
//  #define AllInclusiveSelectMCMCObject 0
//#endif

// This struct contains all of the prior information for the f_2(tau^2)
//  calculation, including what parameters will be used in f_3, f_4
//  bounding densities.
TauOfFContainer *NewTauContainer(SEXP sPriorOftau, 
  SEXP staubarnu, SEXP staupriordf, SEXP sPriorXScaling,
  SEXP PiA, SEXP sTemperature) {
  TauOfFContainer *MT = (TauOfFContainer*) malloc(sizeof(TauOfFContainer));
  if (sPriorOftau != NULL && !Rf_isNull(sPriorOftau) &&
    sPriorXScaling != NULL && !Rf_isNull(sPriorXScaling) &&
    Rf_length(sPriorOftau) > 1 && REAL(sPriorOftau)[0] >= 0 &&
    sPriorXScaling != NULL && !Rf_isNull(sPriorXScaling) && 
    Rf_length(sPriorXScaling) >= 3) {
    MT->sPriorOftau = sPriorOftau;  
    MT->sPriorXScaling = sPriorXScaling;
    //Rprintf("Including a Prior!\n");
  } else if (sPriorXScaling != NULL && !Rf_isNull(sPriorXScaling) 
    && Rf_length(sPriorXScaling) >= 3) {
    MT->sPriorXScaling = sPriorXScaling;
    MT->sPriorOftau = R_NilValue;
  } else {
    MT->sPriorOftau = R_NilValue;
    MT->sPriorXScaling = R_NilValue;
    //Rprintf("There Is no prior!\n");
  }
  MT->staubarnu = staubarnu;
  MT->staupriordf = staupriordf;
  MT->sPriorXScaling = sPriorXScaling;
  MT->D = NULL; MT->R = NULL; MT->NumJ = 0; 
  
  MT->q3X2 = 0.0; MT->f3X2 = 0.0; MT->f2X4 = 0.0; MT->q4X4 = 0.0;
  MT->X2 = 0.0; MT->X4 = 0.0; MT->q4X4=0.0; MT->xM3=0.0; MT->xM4=0.0;
  MT->lM3xM3 = 0.0; 
  MT->lM4xM4=0.0; MT->lM2xM2=0.0; MT->lA=0.0; MT->lLB42=0.0; MT->lUB32=0.0;
  MT->Df3=0.0; MT->Df4=0.0; MT->lLB42=0.0; MT->UsePrior=0.0; MT->lB = 0.0;
  
  MT->Onii = 0;
  if (Rf_length(PiA) == 2) {
    MT->logOddsPi = log(REAL(PiA)[1])  -  log(1.0-REAL(PiA)[1]);
  } else if (Rf_length(PiA) == 1) {
    MT->logOddsPi = log(REAL(PiA)[0])  -  log(1.0-REAL(PiA)[0]);
  }
  if (Rf_isNull(sTemperature) || !Rf_isReal(sTemperature) || REAL(sTemperature)[0] <= 0.0) {
    MT->Temperature = 1.0;  MT->invTemperature = 1.0;
  } else {
    MT->Temperature = REAL(sTemperature)[0]; MT->invTemperature = 1.0/MT->Temperature;
  }
  return(MT);
}


// Given necessary information, fills MT with the needed eigenvalue information.
int SetTauContainer(TauOfFContainer* MT, SEXP sPriorOftau, 
  SEXP staubarnu, SEXP staupriordf, SEXP sPriorXScaling,
  SEXP PiA, int Onii) {
  if (MT == NULL) {
    Rf_error("SetTauContainer: MT is NULL!\n");
  }
  //Rprintf("First Parts of the test!\n"); R_FlushConsole();
  MT->D = NULL; MT->R = NULL;
  MT->sPriorOftau = R_NilValue;
  MT-> sPriorXScaling = R_NilValue;  MT->staupriordf = R_NilValue;
  MT->staubarnu = R_NilValue;
  //Rprintf("Set Everything to NULL!\n"); R_FlushConsole();
  if (sPriorOftau == NULL || Rf_isNull(sPriorOftau)) {
    MT->sPriorOftau = R_NilValue;
    if (sPriorXScaling == NULL || Rf_isNull(sPriorXScaling)) {
      MT->sPriorXScaling = R_NilValue;
    } else if (Rf_length(sPriorXScaling) >= 3) {
      MT->sPriorXScaling = sPriorXScaling;
    }  else {
      MT->sPriorXScaling = R_NilValue;
    }
  } else if (!Rf_isReal(sPriorOftau)   || Rf_length(sPriorOftau) < 1 ||
    REAL(sPriorOftau)[0] <= 0) {
    MT->sPriorOftau = R_NilValue;
    if (sPriorXScaling == NULL || Rf_isNull(sPriorXScaling)) {
      MT->sPriorXScaling = R_NilValue;
    } else if (Rf_length(sPriorXScaling) >= 3) {
      MT->sPriorXScaling = sPriorXScaling;
    }  else {
      MT->sPriorXScaling = R_NilValue;
    }  
  } else if (sPriorOftau != NULL && !Rf_isNull(sPriorOftau) &&
    sPriorXScaling != NULL && !Rf_isNull(sPriorXScaling) &&
    Rf_length(sPriorOftau) > 1 && REAL(sPriorOftau)[0] >= 0 &&
    sPriorXScaling != NULL && !Rf_isNull(sPriorXScaling) && 
    Rf_length(sPriorXScaling) >= 3) {
    MT->sPriorOftau = sPriorOftau;  
    MT->sPriorXScaling = sPriorXScaling;
    //Rprintf("Including a Prior!\n");
  } else if (sPriorXScaling != NULL && !Rf_isNull(sPriorXScaling) 
    && Rf_length(sPriorXScaling) >= 3) {
    MT->sPriorXScaling = sPriorXScaling;
    MT->sPriorOftau = R_NilValue;
  } else {
    MT->sPriorOftau = R_NilValue;
    MT->sPriorXScaling = R_NilValue;
    //Rprintf("There Is no prior!\n");
  }
  MT->staubarnu = staubarnu;
  MT->staupriordf = staupriordf;
  //MT->sPriorXScaling = sPriorXScaling;
  MT->D = NULL; MT->R = NULL; MT->NumJ = 0;
  MT->Onii = Onii;
  if (Rf_length(PiA) == 2) {
    MT->logOddsPi = log(REAL(PiA)[1])  -  log(1.0-REAL(PiA)[1]);
  } else if (Rf_length(PiA) == 1) {
    MT->logOddsPi = log(REAL(PiA)[0])  -  log(1.0-REAL(PiA)[0]);
  }
  //Rprintf("Making it Out of MT Set!\n");  R_FlushConsole();
  return(1);
}

//////////////////////////////////////////////////////////////////////////////
// lFOfX is key function pointer for Splice sampler
//
//  This function, once set with information from container MT, will be the log density of f_2(tau_k).
//  
double lFOfX(double X, void *vMT) {
  TauOfFContainer*MT = (TauOfFContainer*) vMT;
  return((logF(MT->NumJ, X, -2,0,MT->D, MT->R) + 
    GivePriorOfX(X,MT) + MT->logOddsPi) * MT->invTemperature);
}

double GivePriorOfX(double X, TauOfFContainer* MT) {
  if (X == 0) {
    return(-999);
  }
  if (MT->sPriorOftau == NULL || Rf_isNull(MT->sPriorOftau) ||
    Rf_length(MT->sPriorOftau) <= 1) {
    if (REAL(MT->staupriordf)[0] <= 0) {
      return( -(.5 * REAL(MT->staupriordf)[0] + 1.0) * log(X));
    }
    return(  -(.5 * REAL(MT->staupriordf)[0] +1.0)*log(X) -
      .5 *REAL(MT->staubarnu)[0] * REAL(MT->staupriordf)[0] / X+
      -Rf_lgammafn(.5 *REAL(MT->staupriordf)[0]) + 
      .5*REAL(MT->staupriordf)[0] * log(
        REAL(MT->staupriordf)[0] * REAL(MT->staubarnu)[0] *.5));
  }
  Rprintf("PriorOfX:  We're Off using PriorOftau?"); R_FlushConsole();
  double firstX = REAL(MT->sPriorXScaling)[2];
  int ii = 0;
  double nX = X;
  if (REAL(MT->sPriorXScaling)[1] == 1) {
    nX = X;
  } else if (REAL(MT->sPriorXScaling)[1] == 2) {
    nX = sqrt(X);
  } else if (REAL(MT->sPriorXScaling)[1] == -1) {
    nX = log(X);
  }
  ii = (int) floor( (nX-firstX) / REAL(MT->sPriorXScaling)[0]);
  if (ii >= Rf_length(MT->sPriorOftau) || ii < 0) {
     return(NOPROB);
  }
  double RetProb = 0;
  if (Rf_length(MT->sPriorXScaling) > 3) {
    if (REAL(MT->staupriordf)[0] <= 0) {
      RetProb = (-.5 * REAL(MT->staupriordf)[0] -1)*log(X);     
    } else {
      RetProb = (-.5 * REAL(MT->staupriordf)[0] -1)*log(X) -
      .5 *REAL(MT->staubarnu)[0] * REAL(MT->staupriordf)[0] / X +
      -Rf_lgammafn( .5 *REAL(MT->staupriordf)[0]) + 
      .5*REAL(MT->staupriordf)[0] * log(
        REAL(MT->staupriordf)[0] * REAL(MT->staubarnu)[0] *.5); 
    }
  }
  RetProb += REAL(MT->sPriorOftau)[ii];
  return(RetProb);
}


double GiveDerivPriorOfX(double X, TauOfFContainer* MT) {
  if (X == 0) {return(-9999);}
  if (MT->sPriorOftau == NULL || Rf_isNull(MT->sPriorOftau) ||
    Rf_length(MT->sPriorOftau) < 1) {
    return(  (-.5 * REAL(MT->staupriordf)[0] -1) / X  +
      .5 *REAL(MT->staubarnu)[0] * REAL(MT->staupriordf)[0] / (X*X));
  }
  Rprintf("DerivPriorOfX:  We're Off using PriorOftau?"); R_FlushConsole();
  double firstX = REAL(MT->sPriorXScaling)[2];
  int ii = 0;
  double nX = X;
  double dX = .01;
  if (REAL(MT->sPriorXScaling)[1] == 1) {
    nX = X;
    dX = REAL(MT->sPriorXScaling)[0];
  } else if (REAL(MT->sPriorXScaling)[1] == 2) {
    nX = sqrt(X);
    dX = pow(nX + REAL(MT->sPriorXScaling)[0],2) - X;
  } else if (REAL(MT->sPriorXScaling)[1] == -1) {
    nX = log(X);
    dX = exp(nX + REAL(MT->sPriorXScaling)[0]) - X;
  }
  
  ii = (int) floor( (nX-firstX) / REAL(MT->sPriorXScaling)[0]);  
  if (ii >= Rf_length(MT->sPriorOftau) + 1 || ii < 0) {
     return(NOPROB);
  }
  double RetDeriv = 0;
  if (Rf_length(MT->sPriorXScaling) > 3) {
    RetDeriv = (-.5 * REAL(MT->staupriordf)[0] -1)/X +
      .5 *REAL(MT->staubarnu)[0] * REAL(MT->staupriordf)[0] / (X*X); 
  }
  RetDeriv += (REAL(MT->sPriorOftau)[ii+1] - 
    REAL(MT->sPriorOftau)[ii]) / dX;
  return(RetDeriv);
}


double GiveSDerivPriorOfX(double X, TauOfFContainer* MT) {
  if (MT->sPriorOftau == NULL || Rf_isNull(MT->sPriorOftau) ||
    Rf_length(MT->sPriorOftau) < 1) {
    return(  (.5 * REAL(MT->staupriordf)[0] +1) / (X*X)  -
      REAL(MT->staubarnu)[0] * REAL(MT->staupriordf)[0] / (X*X*X));
  }
  Rprintf("2ndDeriv PriorOfX:  We're Off using PriorOftau?"); R_FlushConsole();
  double firstX = REAL(MT->sPriorXScaling)[2];
  int ii = 0;
  double nX = X;
  double dX1 = .01; double dX2 = .01;
  if (REAL(MT->sPriorXScaling)[1] == 1) {
    nX = X;
    dX1 = REAL(MT->sPriorXScaling)[0];
    dX2 = dX1;
  } else if (REAL(MT->sPriorXScaling)[1] == 2) {
    nX = sqrt(X);
    dX1 = pow(nX + REAL(MT->sPriorXScaling)[0],2) - X;
    dX2 = X-pow(nX-REAL(MT->sPriorXScaling)[0],2);
  } else if (REAL(MT->sPriorXScaling)[1] == -1) {
    nX = log(X);
    dX1 = exp(nX + REAL(MT->sPriorXScaling)[0]) - X;
    dX2 = X-exp(nX - REAL(MT->sPriorXScaling)[0]);
  }
  
  ii = (int) floor( (nX-firstX) / REAL(MT->sPriorXScaling)[0]);  
  if (ii >= Rf_length(MT->sPriorOftau) + 1 || ii < 0) {
     return(NOPROB);
  }
  double RetSDeriv = 0;
  if (Rf_length(MT->sPriorXScaling) > 3) {
    RetSDeriv = (.5 * REAL(MT->staupriordf)[0] -1)/(X*X) -
      .5 *REAL(MT->staubarnu)[0] * REAL(MT->staupriordf)[0] / (X*X*X); 
  }
  if (ii > 0) {
    RetSDeriv += ( 
     (REAL(MT->sPriorOftau)[ii+1] - 
        REAL(MT->sPriorOftau)[ii]) / dX1 - 
     (REAL(MT->sPriorOftau)[ii] - 
       REAL(MT->sPriorOftau)[ii-1]) / dX2) / (.5 *( dX1+dX2));
  } else {
    RetSDeriv += ( 
     (REAL(MT->sPriorOftau)[ii+2] - 
        REAL(MT->sPriorOftau)[ii+1]) / dX1 - 
     (REAL(MT->sPriorOftau)[ii+1] - 
       REAL(MT->sPriorOftau)[ii]) / dX2) / (.5 *( dX1+dX2));  
  }
  return(RetSDeriv);
}


/////////////////////////////////////////////////////////////////////////////
//  logF()  :
//  logF density up to scale of f_2(tau_k)
//  
// As related in BayesSpike paper this is the parameteric form of likelihood
// in active density condtional for tau_k.
double logF(int J, double tau, double nu, double taubarnu, 
  double *D, double *RVals) {
  if (tau == 0) {return(-9999);}
  double A = -((J+nu+2))* log(tau); 
  int ii = 0;
  double itau = 1.0/tau;
  for (ii = 0; ii < J; ii++) {
    A -= log(D[ii] + itau);
  }
  A -= taubarnu * nu *(itau);
  for (ii = 0; ii < J; ii++) {
    A += RVals[ii] / ( D[ii] + itau);
  }
  return(.5 * A);
}

/////////////////////////////////////////////////////////////////////////////
// Identifying the maximum and the curvature at the maximum of f_2(tau_k)
// is critical in choosing best bounding densities.
//
// DerivF and SDerivF are first and second derivatives
//
double DerivF (int J, double tau, double nu, double taubarnu,
  double *D, double *RVals) {
  double itau = 1.0/tau;
  double A = -((J+nu+2)) * itau +  taubarnu * nu * itau*itau;
  int ii = 0;
  double tausq = tau * tau;
  for (ii = 0; ii < J; ii++) {
    A += pow( tausq * D[ii] + tau,-1.0);
  }
  for (ii = 0; ii < J; ii++) {
    A += RVals[ii] * pow(tau*D[ii] + 1.0,-2.0);
  }
  return(.5 * A);
}

double SDerivF(int J, double tau, double nu, double taubarnu,
  double *D, double *RVals) {
  double itau = 1.0/tau;
  double itausq = itau * itau;
  double itau3 = itau * itausq;
  double tausq = tau * tau;
  double A = ( J + nu+2) * itausq;
  int ii = 0;
  for (ii = 0; ii < J; ii++) {
    A -= (2.0 * tau * D[ii] + 1.0) *
      pow(tausq * D[ii] + tau,-2); 
  }
  A -= 2*taubarnu * nu * itau3;
  for (ii = 0; ii < J; ii++) {
    A+= RVals[ii] * (-2.0 * D[ii]) *
      pow(tau*D[ii] +1.0,-3.0); 
  }
  return(.5 * A);
}


int dgePcopy(double *SourceMat, double *NewMat, int L1, int L2,
  int LDA1, int LDA2) {
  int ii;
  for (ii = 0; ii < L1; ii++) {
    F77_CALL(dcopy)(&L2, SourceMat + ii, &LDA1, NewMat+ii, &LDA2);
  }
  return(1);
}

// Alternate Dependencies Tree
// Sometime, we only wish tau_k to be possible if other fixed or random effects
// are already non-zero.  This is, we don't normally want higher order or interactive
// terms to be non zero if there was no observed simpler effects.
int CheckDependencies(SEXP DependenceList, int iFirstRandom,
  short *BFixed, SEXP sOnTau, int iOnii) {
  int ii;
  //if (Rf_length(DependenceList) < length(BFixed)) {
  //  Rf_error("Dependence List fails for %d, not long enough here!", 
  //    GetFirstInteger(Onii));
  //}
  if (Rf_length(DependenceList) == 0) { return(1);}
  for (ii = 0; ii < Rf_length(DependenceList); ii++) {
    if ( (ANINT(DependenceList, ii)) < 0 &&
      (abs(ANINT(DependenceList,ii)) < iFirstRandom)) {
      if ((int)BFixed[(abs(ANINT(DependenceList, ii)))]  == 0 ) {
        return(0);
      }
    }  else if (ANINT(DependenceList,ii) < Rf_length(sOnTau)) {
      if (REAL(sOnTau)[ANINT(DependenceList, ii)] <= 0) {
        return(0);
      }
    }
  }
  return(1);
}


void TestMTMaterials(SEXP sPriorOftau, SEXP sPriorXScaling,
  SEXP snu, SEXP staubarnu) {
    Rprintf("SampleNewTaus: Inspecting MT Materials\n"); R_FlushConsole();
    if (sPriorOftau == NULL) {
      Rprintf("MTTest: sPriorOftau is Really NULL\n"); R_FlushConsole();
      Rprintf("  This is a big Error, we're going to return an error now!\n");
      R_FlushConsole();
      Rf_error(" We want at least a NIL Object Doh!!\n");
    } else if (Rf_isNull(sPriorOftau)) {
      Rprintf("MTTest: sPriorOftau is NIL Object\n"); R_FlushConsole();
    } else {
      Rprintf("MTTest: sPriorOftau has length(sPriorOftau) =%d \n",
        Rf_length(sPriorOftau)); R_FlushConsole();
    }
    if (sPriorXScaling == NULL) {
      Rprintf("MTTest:  sPriorXScaling is Really NULL\n"); R_FlushConsole();
      Rprintf("  This is a big Error, we're going to return an error now!\n");
      R_FlushConsole();
      Rf_error(" We want at least a NIL Object Doh!!\n");
    } else if (Rf_isNull(sPriorXScaling)) {
      Rprintf("MTTest: sPriorXScaling is NIL Object\n"); R_FlushConsole();
    } else {
      Rprintf("MTTest: sPriorXScaling has length(sPriorXScaling) = %d\n",
        Rf_length(sPriorXScaling)); R_FlushConsole();
    }
    if (staubarnu == NULL) {
      Rprintf("MTTest: staubarnu is Really NULL\n"); R_FlushConsole();
      Rprintf("  This is a big Error, we're going to return an error now!\n");
      R_FlushConsole();      
      Rf_error(" We want at least a NIL Object Doh!!\n");
    } else if (Rf_isNull(staubarnu)) {
      Rprintf("MTTest: staubarnu is NIL Object\n"); R_FlushConsole();
    } else {
      Rprintf("MTTest: staubarnu has length(staubarnu) =%d, = %f \n",
        Rf_length(staubarnu), REAL(staubarnu)[0]); R_FlushConsole();
    }
    if (snu == NULL) {
      Rprintf("MTTest: snu is Really NULL\n"); R_FlushConsole();
      Rf_error(" We want at least a NIL Object Doh!!\n");
    } else if (Rf_isNull(snu)) {
      Rprintf("MTTest: snu is NIL Object\n"); R_FlushConsole();
    } else {
      Rprintf("MTTest:  snu has length(snu) = %d\n",
        Rf_length(snu)); R_FlushConsole();
      Rprintf("MTTest: snu[0] = %f\n", REAL(snu)[0]);
    }    
  return;
}

///////////////////////////////////////////////////////////////////////////
// Expand in Reverse
//
//  Not clear this function ever gets called.  Mostly testing quality of
//  memory in the PBeta vector.
int ExpandInReverse(SEXP PBeta, SEXP FirstRandom, SEXP tauEndList,
  SEXP CurrentTau, SEXP BFixed) {
int IPT = Rf_length(PBeta)-1;
int ii; int jj;
int St; int k;
int CKill = 0;
if (Rf_length(CurrentTau) >= 1) {
  for (jj = Rf_length(CurrentTau)-1; jj >= 0; jj--) {
    if (jj == 0) {
        St = GetFirstInteger(FirstRandom);
    } else {
        St = (ANINT(tauEndList, (jj-1))) +1;
    }
    k = (ANINT(tauEndList, jj)) - St + 1;
    if (REAL(CurrentTau)[jj] > 0) {
      CKill +=k;
    }
  }
}
if (Rf_length(BFixed) >=1) {
  for (jj = Rf_length(BFixed)-1; jj>=0; jj--) {
    if ((Rf_isReal(BFixed) && REAL(BFixed)[jj] > .2) ||
        (Rf_isInteger(BFixed) && INTEGER(BFixed)[jj] > 0)) {
      CKill++;    
    }
  }
}
IPT = CKill-1;
if (CKill > Rf_length(PBeta)) {
  Rprintf("PrintBack: CKill Error, PBeta of wrong length");
}
if (Rf_length(CurrentTau) >= 1) {
  for (jj = Rf_length(CurrentTau)-1; jj >= 0; jj--) {
    if (jj == 0) {
        St = GetFirstInteger(FirstRandom);
    } else {
        St = (ANINT(tauEndList, (jj-1))) +1;
    }
    k = (ANINT(tauEndList, jj)) - St + 1;
    if (REAL(CurrentTau)[jj] > 0) {
      for (ii = k-1; ii >= 0; ii--) {
        if (IPT < 0) {
         Rf_error("PrintBack, IPT is less than zero!");
        }
        REAL(PBeta)[St + ii] = REAL(PBeta)[IPT]; IPT--;
      }
    }
  }
}
if (Rf_length(BFixed) >=1) {
  for (jj = Rf_length(BFixed)-1; jj>=0; jj--) {
    if ((Rf_isReal(BFixed) && REAL(BFixed)[jj] > .2) ||
        (Rf_isInteger(BFixed) && INTEGER(BFixed)[jj] > 0)) {
      if (IPT < 0) {
        Rf_error("PrintBack, IPT is less than zero!");
      }
      REAL(PBeta)[jj] = REAL(PBeta)[IPT]; IPT--;   
    }
  }
}
if (IPT > 0) {
  Rprintf("PrintBack, error IPT > 0!");
}
return(1); 
} 


int PrintDetailOfRsprec(SEXP MySEXP, char* MyC) {
  Rprintf("  Testing now:  %s  \n", MyC);  R_FlushConsole();
  if (MySEXP == NULL) {
    Rprintf("  %s == True NULL \n", MyC); R_FlushConsole(); return(-1);
  }
  if (Rf_isNull(MySEXP)) {
    Rprintf("  %s is NilValue \n", MyC); R_FlushConsole(); return(-1);
  }
  if (Rf_isInteger(MySEXP)) {
    Rprintf("  %s is Integer, Rf_length %d :", MyC, Rf_length(MySEXP));
    R_FlushConsole();
    if (Rf_length(Rf_getAttrib(MySEXP, R_DimSymbol)) >= 2) {
      Rprintf("  %s is a matrix of high dim=(%d,%d)\n", MyC,
        INTEGER(Rf_getAttrib(MySEXP, R_DimSymbol))[0],INTEGER(Rf_getAttrib(MySEXP, R_DimSymbol))[1]);
      PrintRMatrix(INTEGER(MySEXP),
        INTEGER(Rf_getAttrib(MySEXP, R_DimSymbol))[0],INTEGER(Rf_getAttrib(MySEXP, R_DimSymbol))[1]);
      R_FlushConsole(); return(1);
    }
    PrintVector(INTEGER(MySEXP), Rf_length(MySEXP));
    Rprintf("\n"); R_FlushConsole(); return(2);
  }
  
  if (Rf_isReal(MySEXP)) {
    Rprintf("  %s is Real, Rf_length %d :", MyC, Rf_length(MySEXP));
    R_FlushConsole();
    if (Rf_length(Rf_getAttrib(MySEXP, R_DimSymbol)) >= 2) {
      Rprintf("  %s is a matrix of high dim=(%d,%d)\n", MyC,
        INTEGER(Rf_getAttrib(MySEXP, R_DimSymbol))[0],INTEGER(Rf_getAttrib(MySEXP, R_DimSymbol))[1]);
      PrintRMatrix(REAL(MySEXP),
        INTEGER(Rf_getAttrib(MySEXP, R_DimSymbol))[0],INTEGER(Rf_getAttrib(MySEXP, R_DimSymbol))[1]);
      R_FlushConsole(); return(1);
    }
    PrintVector(REAL(MySEXP), Rf_length(MySEXP));
    Rprintf("\n"); R_FlushConsole(); return(2);
  }
  Rprintf("  Oops:  I don't know what %s, is!!\n", MyC);
  R_FlushConsole();
  return(-4);
}


/*
double CheckStuff(SEXP sOnBeta, SEXP XtResid, SEXP MCMCSexp) {
AllMCMCSexp *AMS = (AllMCMCSexp*) R_ExternalPtrAddr(MCMCSexp);
  int ii = 0; int NotBeta= 0;
  if (AMS == NULL) {
    Rf_error("CheckStuff: AMS is Null!");
  }
  Rprintf(" Checking Stuff! "); R_FlushConsole();
  if (Rf_length(AMS->sOnBeta) != Rf_length(sOnBeta)) {
    Rf_error("AMS->sOnBeta not same length as OnBeta!");
  }
  for (ii = 0; ii < Rf_length(sOnBeta); ii++) {
    if (REAL(AMS->sOnBeta)[ii] != REAL(sOnBeta)[ii]) {
       NotBeta = 1; break;
    }
  }
  //int NotXTResid = 0;
  //for (ii = 0; ii < Rf_length(XtResid); ii++) {
  //  if (MCMCSexp->CDO->XTResid[ii] != REAL(XtResid)[ii]) {
  //     NotXTResid = 1; break;
  //  }
  //}
  if (NotBeta == 1) {
    Rprintf("Rf_error Checkstuff, sOnBeta = ");  PrintVector(REAL(sOnBeta), Rf_length(sOnBeta));
    Rprintf("\n  But in MCMCSexp = ");  PrintVector(REAL(AMS->sOnBeta), Rf_length(AMS->sOnBeta));
  }
  //if (NotXTResid == 1) {
  //  Rprintf("Rf_error Checkstuff, XTResid = ");  PrintVector(REAL(XtResid), Rf_length(XtResid));
  //  Rprintf("\n  But in MCMCSexp = ");  PrintVector(AMS->CDO->XtResid, Rf_length(AMS->sOnBeta));
  //}
  if (NotBeta == 0) {
    Rprintf("Clean Bill of Health!\n"); R_FlushConsole();
  }
  return(2.0);
}
*/
