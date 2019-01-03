/* ========================================================================== */
/*                                                                            */
/*   BayesSpikeCpp.h                                                          */                              
/*   (c) 2012 Alan Lenarcic                                                   */
/*                                                                            */
/*   Description                                                              */
/*   Header File For BayesSpike Project, Defines BayesSpikeCL RCpp class      */
/*                                                                            */
/*  This file defines the BayesSpikeCL class which is externed as a "Rcpp"    */
/*    module                                                                  */
/*   through "BayesSpike.cpp"                                                 */
/*                                                                            */
/*  Most of the public class functions are accessible at the R prompt.        */
/*  This file contains the constructor and destructor functions.              */
/*     as well as getter/setters                                              */
/*                                                                            */
/*  Most Algorithm functions are in "BayesSpikeGibbs.cpp"                     */
/*  The Group Selector slice sampler is in "BayesSpikeTau.cpp"                */
/*  Nonclass Helper functions are in "BayesSpikeHelpers.cpp"                  */
/*                                                                            */
/*                                                                            */
/*                                                                            */
/*                                                                            */
/* ========================================================================== */

// Uncomment to activate timing systems
//#define DOTIMEH 1

#ifdef DOTIMEH
#ifndef TIMEH
 #include <time.h>
 #define TIMEH 1
#endif
#endif

#ifndef RCPPH
  #include <Rcpp.h>
  #include <cmath>
  #define RCPPH 0
#endif



#define UpLo 'L'
#define TransWant 'N'
#ifndef RNeeds 
 #define RNeeds 0
 #include <R.h>
 #include <Rmath.h>
 #include <Rinternals.h>
 #include <R_ext/BLAS.h>
 #include <R_ext/Lapack.h>
 #include <Rdefines.h>
#endif

//  MyMarixOp2009.h contains mainly printing functions for diagnosing matrix issues.
#ifndef MyMatrixOp2009DD
  #define MyMatrixOp2009DD 0
  #include "MyMatrixOp2009.h"
#endif

// A HeapSort function
#ifndef HEAPSORTH
  #define HEAPSORTH 1
  #include "HeapSort.h"
#endif

//  Basic Square Macro
#ifndef Sq
  #define Sq(X)  (X) * (X)
#endif

//  Maximum integration value for TAU
//  Note, if you have to integrate up to tau = 10000,
//    then it is likely that  Tau != 0!
#ifndef MAXTAU
  #define MAXTAU 10000
#endif


#ifndef ACPPH
  #include "Acpp.h"
#endif 

//  We might have to include an alternative RInternals if Rcpp does not compile correctly.
//
//
#ifndef SEXPREC_HEADER2
  #include "RInternals2.h"
#endif
    
#define TooMuchBuffer 100000
#define MaxSquaresAllocation 20480
//#define DefaultMaximumAllocation 40960
#define DefaultMaximumAllocation 4096
#define DefaultStartRemoval   1024
#define DefaultCheckRemovable 100
#define DefaultKillPct .2
#define DefaultMinimumProbKeeper .01
#define DefaultMaxInvert 500

typedef  long int bufferILocType;
//  All of this macro is only loaded if BAYESSPIKECPPH is not defined.
//
#ifndef BAYESSPIKECPPH
  #define BAYESSPIKECPPH 0
  
  

#define logSqrt2PI 0.9189385

#define NOPROB 999

#ifndef MRF
   #ifdef DOTIMEH
 #define MRF(FF, charS, step, Loc) if (Verbose > 1) { \
   Rprintf("  %d,  Num %d: Running Function %s \n", tt, step, charS); R_FlushConsole(); \
   } \
   if (Loc >= 0) { Clockt1 = clock(); } \
   SF = FF(); \
   if ( ((step) == EarlyEndStep) && tt == EarlyEndtt)  { \
     Rf_error("  Quitting Based upon Early End after %d, %d, %s \n",  \
       tt, (step), charS ); \
   } \
   if (SF < 0) { \
     Rprintf("BayesSpikeGibbs.cpp:::RunAlgorithm():: ");                      \
     Rprintf(" We get a fail=%d in trying to run function %s \n", SF, charS); \
     R_FlushConsole(); \
     Rf_error("--- ERROR RunAlgoirthm Error returned in function %s \n", charS);  \
   } \
   if (Loc >= 0) { Clockt2 = clock();     \
      if (RsTimePartsList != NULL) { \
        REAL(RsTimePartsList->asSexp())[Loc] += ((double)(Clockt2-Clockt1) / CLOCKS_PER_SEC);  \
      } \
   }  \
   SF = 1;
   #else 
   #define MRF(FF, charS, step, Loc) if (Verbose > 1) { \
   Rprintf("  %d,  Num %d: Running Function %s \n", tt, step, charS); R_FlushConsole(); \
   } \
   SF = FF(); \
   if ( ((step) == EarlyEndStep) && tt == EarlyEndtt)  { \
     Rf_error("  Quitting Based upon Early End after %d, %d, %s \n",  \
       tt, (step), charS ); \
   } \
   if (SF < 0) { \
     Rprintf("  We get a fail=%d in trying to run function %s \n", SF, charS); \
     R_FlushConsole(); \
     Rprintf(" We get a fail=%d in trying to run function %s \n", SF, charS); \
     R_FlushConsole(); \
     Rf_error("--- ERROR RunAlgoirthm Error returned in function %s \n", charS);  \
   } \
   SF = 1
   #endif
#endif

#ifndef MRFRefreshOrderedActive
 #ifdef DOTIMEH 
 #define MRFRefreshOrderedActive(FF, charS, step, Loc) if (Verbose > 1) { \
   Rprintf("  %d,  Num %d: Running Function %s \n", tt, step, charS); R_FlushConsole(); \
   } \
   if (Loc >= 0) { Clockt1 = clock(); } \
   SF = RefreshOrderedActive(FF); \
   if ( ((step) == EarlyEndStep) && tt == EarlyEndtt)  { \
     Rf_error("  Quitting Based upon Early End after %d, %d, %s \n",  \
       tt, (step), charS ); \
   } \
   if (SF < 0) { \
     Rprintf("  We get a fail=%d in trying to run function %s \n", SF, charS); \
     R_FlushConsole(); \
   } \
   if (Loc >= 0) { Clockt2 = clock();     \
      if (RsTimePartsList != NULL) { \
        REAL(RsTimePartsList->asSexp())[Loc] += ((double)(Clockt2-Clockt1) / CLOCKS_PER_SEC);  \
      } \
   }  \
   SF = 1
  #else 
     #define MRFRefreshOrderedActive(FF, charS, step, Loc) if (Verbose > 1) { \
   Rprintf("  %d,  Num %d: Running Function %s \n", tt, step, charS); R_FlushConsole(); \
   } \
   SF = RefreshOrderedActive(FF); \
   if ( ((step) == EarlyEndStep) && tt == EarlyEndtt)  { \
     Rf_error("  Quitting Based upon Early End after %d, %d, %s \n",  \
       tt, (step), charS ); \
   } \
   if (SF < 0) { \
     Rprintf("  We get a fail=%d in trying to run function %s \n", SF, charS); \
     R_FlushConsole(); \
   } \
   SF = 1
  #endif
#endif
   
#define DefaultMaxCODA 10000
// These functions are self defined deletion/free macros for memory management
//
#ifndef DDelete
  #define DDelete(X, charS)  if (Verbose > 2) {  \
    Rprintf("deleting  %s \n ", charS); R_FlushConsole(); } \
    if ( (X) != NULL ) {  delete((X));  (X) = NULL; }
    
  #define DDeleteA(X, charS)  if (Verbose > 2) {  \
    Rprintf("deleting  %s location %d\n ", charS, X->getInVectorLoc()); R_FlushConsole(); } \
    if ( (X) != NULL ) {  delete((X));  (X) = NULL; }
#endif
#ifndef FFree
  #define FFree(X, charS)  if (Verbose > 2) {  \
    Rprintf("Freeing  %s \n ", charS); R_FlushConsole(); } \
    if ( (X) != NULL ) {  Free((X));  (X) = NULL; }
#endif
#ifndef ffree
  #define ffree(X, charS)  if (Verbose > 2) {  \
    Rprintf("Freeing  %s \n ", charS); R_FlushConsole(); } \
    if ( (X) != NULL ) {  free((X));  (X) = NULL; }
#endif
#ifndef GetFirstInteger
#define GetFirstInteger(x) (( int)( Rf_isInteger((x)) ? INTEGER((x))[0] : (int) REAL((x))[0] ) )
#endif
#ifndef ANINT
  #define ANINT(x,ii) ( (int) Rf_isInteger((x)) ? INTEGER((x))[(ii)] : (int) REAL((x))[(ii)] )
#endif
#ifndef MemError
  #define MemError(X, charS) if ( (X) == NULL ) {  \
    Rf_error("%s could not be allocated \n", (char*) charS); }


  #define RMemGetD( X , charS , Count) \
    if ( Verbose > 3 ) { \
      Rprintf("Allocating %s \n", (char*) (charS) ); \
      R_FlushConsole(); \
    } \
    if ((X) != NULL ) { \
      Rf_error("RMemGetD: Error, trying to allocate to %s, but is not null!\n", charS); \
    }\
    (X) = (double*) Calloc( (Count), double ); \
    if ( (X) == NULL) {     \
      Rf_error("%s, could not be allocated \n", (char*) charS); \
    }

  #define RMemGetI( X , charS , Count) \
    if ( Verbose > 3 ) { \
      Rprintf("Allocating %s \n", (char*) (charS) ); \
      R_FlushConsole(); \
    } \
    (X) = (int*) Calloc( (Count), int ); \
    if ( (X) == NULL) {     \
      Rf_error("%s, could not be allocated \n", (char*) charS); \
    }
    
  #define RMemGetsuI( X , charS , Count) \
    if ( Verbose > 3 ) { \
      Rprintf("Allocating %s \n", (char*) (charS) ); \
      R_FlushConsole(); \
    } \
    (X) = (unsigned short int*) Calloc( (Count), unsigned short int ); \
    if ( (X) == NULL) {     \
      Rf_error("%s, could not be allocated \n", (char*) charS); \
    }          
  #define RMemGetS( X , charS , Count) \
    if ( Verbose > 3 ) { \
      Rprintf("Allocating %s \n", (char*) (charS) ); \
      R_FlushConsole(); \
    } \
    (X) = (short int*) Calloc( (Count), unsigned int ); \
    if ( (X) == NULL) {     \
      Rf_error("%s, could not be allocated \n", (char*) charS); \
    } 
  #define RMemGet( X , charS , Count , type ) \
    if ( Verbose > 3 ) { \
      Rprintf("Allocating %s \n", (char*) (charS) ); \
      R_FlushConsole(); \
    } \
    (X) = ((type)*) Calloc( (Count), (type) ); \
    if ( (X) == NULL) {     \
      Rf_error("%s, could not be allocated \n", (char*) charS); \
    }
    
    
  #define RMemRealloc(X, charS, Count, type) if ( Verbose > 3 ) { \
    Rprintf("ReAllocating %s \n", (char*) (charS)); R_FlushConsole(); } \
    X = ((type)*) Realloc( (X), (Count), (type) ); \
    if ((X) == NULL) {     \
      Rf_error("%s, could not be reallocated \n", (char*) charS); \
    }

  #define RMemReallocD(X, charS, Count) if ( Verbose > 3 ) { \
    Rprintf("ReAllocating %s \n", (char*) (charS)); R_FlushConsole(); } \
    X = (double*) Realloc( (X), (Count), double ); \
    if ((X) == NULL) {     \
      Rf_error("%s, could not be reallocated \n", (char*) charS); \
    }

  #define RMemReallocI(X, charS, Count) if ( Verbose > 3 ) { \
    Rprintf("ReAllocating %s \n", (char*) (charS)); R_FlushConsole(); } \
    X = ( int *) Realloc( (X), (Count), int ); \
    if ((X) == NULL) {     \
      Rf_error("%s, could not be reallocated \n", (char*) charS); \
    }
  #define RMemReallocS(X, charS, Count) if ( Verbose > 3 ) { \
    Rprintf("ReAllocating %s \n", (char*) (charS)); R_FlushConsole(); } \
    X = ( short *) Realloc( (X), (Count), short ); \
    if ((X) == NULL) {     \
      Rf_error("%s, could not be reallocated \n", (char*) charS); \
    }  
#endif

int dgePcopy(double *SourceMat, double *NewMat, int L1, int L2,
  int LDA1, int LDA2);

/***************************************************************************
 *  TauOfFContainer
 *
 *    This structure contains information on P(tau_k|tau_k != 0) group prior
 *    
 *    Typically this should be inverse gamma with 
 *       staubarnu (mean), staupriordf (degrees of freedom) 
 *    As prior, however, it is possible for the user to define a complex probability
 *    Distribution which goes into sPriorXScaling as a double vector.
 *    
 *    D and R are two diagonal vectors  For the entire group of parameters
 *    The algorithm decomposes XJkResid =  X_Jk^T(Y-X_{not Jk} Beta_{not Jk)) 
 *      this is the residual correlation of group J(k) with all other
 *      Beta coefficients not in J(k).  Then an Eigen decomposition is made of
 *                 X_JkJk =   VDV^T 
 *      where D is diagonal of eigenvalues of X_JkJk.   We then take
 *              R = V^T %*%  XJkResid
 *      this gives the correlation of the eigenvectors of X_JkJk with XJkResid.
 *      Both of these vectors combined 
 *         make it O|Jk| to calculate the posterior P(tau^2 | Beta_{not Jk}, Y))         
 *
 *    This container is used in "BayesSpikeTau.cpp"
 *      There is a getter function "get_MT"
 *      The current TauOfFContainer for the BayesSpikeCL is element MT  
 *
 *****************************************************************************/       
typedef struct TauOfFContainer_t {
  int NumJ;
  SEXP sPriorOftau;
  SEXP staubarnu;
  SEXP staupriordf;
  SEXP sPriorXScaling;   
  double *D; 
  double *R;
  int Onii;  double logOddsPi;
  double Temperature;  double invTemperature;
  double q3X2, f3X2, f2X2,  f2X4, X2, X4, q4X4, xM2, xM3, xM4, lM3xM3, lB;
  double lM4xM4, lM2xM2, lA, lLB42, lUB32, Df3, Df4, UsePrior, rUnif1;
} TauOfFContainer;

//double FunInt(double (*GiveMelFOfX)(double, void*), 
//  TauOfFContainer* MT, int MaxN, double Tweak);
 
//  Given a tau value, R, D values what is posterior F (up to a constant)
//    What is derivative of that log posterior?
//    What is second derivative of that log Posterior?
//                        
double logF(int J, double tau, double nu, double taubarnu,
  double *D, double *RVals);
double DerivF(int J, double tau, double nu, double taubarnu,
  double *D, double *RVals);
double SDerivF (int J, double tau, double nu, double taubarnu,
  double *D, double *RVals);

//  Given a tau value, R, D values what is prior F (up to a constant)
//    What is derivative of that log Prior?
//    What is second derivative of that log Prior?
//
double GivePriorOfX(double X, TauOfFContainer* MT);
double lFOfX(double X, void *vMT);
double GiveDerivPriorOfX(double X, TauOfFContainer* MT);
double GiveSDerivPriorOfX(double X, TauOfFContainer* MT);
TauOfFContainer *NewTauContainer(SEXP sPriorOftau, 
  SEXP staubarnu, SEXP staupriordf, SEXP sPriorXScaling, SEXP PiA,
  SEXP sTemperature);
int PrintDetailOfRsprec(SEXP MySEXP, char* MyC);


/////////////////////////////////////////////
//  DerivesA suggested degrees of freedom for dominating density; 
inline double DeriveAlpha(double D, double xm2); 

////////////////////////////////////////////////////////////////////////////
//
// Class:: BayesSpikeCL
//  New Rcpp linked class which contains practical info for BayesSpike
//   and TwoLasso Data
//
//  This class is adaptable to do BayesianGroupSelection for large covariate sets X
//   It uses dynamic memory management similar to Coordinate Descent Lasso 
//   (Friedman et. al 2007) to keep a relevant set of X^TX or wX^TX for
//   calculation of the key XTResid vector = X^T(Y-XBeta)
//
//
//
//using namespace Rcpp;
//using namespace Rcpp;
class BayesSpikeCL {
private :

  SEXP BayesSpikeNameSpace;  AObject *RBayesSpikeNameSpace;
  
  //  Permanent, private constants: Dataset in sX, sY
  int p, n;   // X is a p x n  matrix, Y is n length vector
  SEXP sX, sY;  AObject *RsX, *RsY;
  
  //  For non gaussian Robit data
  //    intY is the Y (0,1) logit vector, dfRobit is Robit degrees of freedom
  int *intY; double dfRobit;
  SEXP sRandomInfoList;  RRObject *RsRandomInfoList;
  
  ///////////////////////////////////////////
  // SmallXtX, rIChol, RICholBack
  //   Since BayesSpike Selects a subset of vectors, rSmallXtX is the small
  //   matrix XtX of relevant selected vectors for each gibbs sample
  //
  int MaxLargeSizeSquares;
  double *rSmallXtX, *rI, *rIChol, *rICholBack; 
  double *PropBeta; double *XtXRecorded;  int LengthXtXRecorded;
  SEXP InitFlags; AObject *RInitFlags;
  
  ///////////////////////////////////////////
  //  XtX is dynamically allocated according to Coordinate Descent Lasso 
  //   allocation methods.  If non-Gaussian (Robit or t-noise) data is used
  //   XtX will really be sum_i w_i X_ij X_ij'   weighted sum elements
  //
  double **pXtX; double *cXtY;   double *smXtY;
  double *XtY;
  
  /////////////////////////////////////////////////////////
  //  XtResid = XtY - XtX Beta
  //
  //  This is the key selection vector for regression in BayesSpike
  double *XtResid;
  int MaxInvert;

  
  // Indicates which vectors if XtX are properly weighted (for non-Gauss)
  //   and which are not
  short int *iWeightedXtX;
  
  //  TBSRoo is a TBSR5 object, aka R5 Reference Class
  //   This is a storage container that is non Cpp of relevant features of R data
  SEXP TBSRoo; AObject *RTBSRoo;
  AObject *RTBSR5;  SEXP TBSR5;
  
  //  Fixed Effects coefficients are not shrunk by group
  //     These are priors for the fixed effects coefficients.
  int LengthNoShrinkFixed, LengthNoShrinkRandom;
  int *NoShrinkFixed;  double *NoShrinkFixedPrior;
  int *NoShrinkRandom;  double *NoShrinkRandomPrior;
  
  //////////////////////////////////////////////////////////////////////////
  //  If Non-Gauss data is used, iiWeight is vector of all weights of 
  //    each datapoint Y_i
  //  dfTNoise is the degrees of freedom of noise.
  //
  double dfTNoise;
  double *iiWeight;
  int StartRecordingiiWeight;  double *SumiiWeight; int CountiiWeight;
  

  
  // kXFinder is vector of length OnKammaMem (OnKappaS is filled)
  //  It maps the columns of XtX to the variables they belong to.
  //
  // XLC is a vector of length p.  For variables that do not have existance
  //  in XtX, the XLC[ii] > 0.  For variables that do have existance in XtX
  //  then XLC is a map up to that value.
  int *kXFinder; int *XLC;
  int OnKappaS;  int OnKappaMem;  
  int InitKKs;
  int Verbose;

  SEXP sBeta; AObject *RsBeta;


  int DoTheKills(int *RankMyUseless, int FindPivot, int LeftRight);
  
  double *TemperatureList;  int LengthTemperatureList; int Tempii;
  double TemperatureDecreasingRate;
  double Temperature, invTemperature;
  int *LengthTempDraws;
  int *NumMerges;
  int FileCodaChainIter;

  double *WW;
  int TrueOrder;
  int TotalGibbsIters;
  
  int *OrderAttack;

  int *Codajj; int *CodaTjj; bufferILocType *CodaTLocjj;   int OnCodaTLocjj;
  int LengthCodaTjj;
  int LengthCodaTLoc;  int TotCodaTLocjj;
  int NewWriteTLoc; 
  double *CodaBeta;  double *CodaTau;  double *CodaP;
  int CodaTisSetup;
  int OnCoda; int OnCodaT;  int TotOnCodaT;
  

  
  
  SEXP sOnPiA, PiAPrior;   AObject *RsOnPiA, *RPiAPrior;
  SEXP sOnSigma, SigmaPrior; AObject *RsOnSigma, *RSigmaPrior;
  SEXP tauFixed; SEXP sOnTau;  AObject *RtauFixed, *RsOnTau;
  SEXP tauFixedPrior;  AObject *RtauFixedPrior;
  int iFirstRandom; SEXP tauEndList;  AObject *RtauEndList;
  AObject *RsCenteredColumns;
  AObject *RsDeCenteredCodaList;
  short int *BFixed; double *ProbFixed; double *ProbTau; double CurrentProb;
  int NumActive; int *OrderedActive;
  
  int LengthMergeActive; int MaxLengthMergeActive; int *MergeActive; double *MergePropBeta;

  SEXP staubarnu; SEXP staupriordf;
  SEXP sPriorXScaling; SEXP sPriorOftau;
  AObject *Rstaubarnu, *Rstaupriordf, *RsPriorXScaling, *RsPriorOftau;    
  TauOfFContainer *MT; 
  
  double YResidSq;  double YSq;  
  SEXP DependenciesList;
  AObject *InitiateTime, *SecondInitiateTime, *CompleteTime;
  
  int MaxTauList;  double *SmallXtResid, *SmallRVec;
  
  SEXP AllEigenValues;  SEXP AllEigenVectors;
  AObject *RAllEigenValues, *RAllEigenVectors;
  AObject*RDependenciesFixed; AObject*RDependenciesTau;
  SEXP DependenciesFixed;  SEXP DependenciesTau;
  
  // Coda Lists attempt to save the chains in objects suitable for the CODA package
  SEXP CodaTable; AObject *RCodaTable;
  AObject *ROldColNames, *ROldDeCenterColNames;
  AObject *ROtherNameCodaList, *ROtherNameDeCenterCodaList;
  AObject *RSubCodaList;  AObject *RSubCodaSubSetCoords;
  AObject *RTauCodaList;  AObject *RSubSetTau;
  AObject *RsTimePartsList;
  #ifdef DOTIMEH
    time_t Rect1, Rect2, AllRect1, AllRect2, TauRect1, TauRect2;
    unsigned long TauClock1, TauClock2, Clockt1, Clockt2;
    //double MyAns;
  #endif
  SEXP DoRecord;  AObject *RDoRecord;
  int NumActiveFixed;  int NumActiveTau;
  double *XjSq;
  int OnRNGState;

///////////////////////////////////////////////////////////////////
//  These are just info passed to SetupTauSamplers  
  int k, St;  SEXP MyEigenValueList;  SEXP MyEigenVectorList;
  double TweakCons;
  
  SEXP CodaList;  // List of Coda Output MCMC tables if saved 
  AObject *RAllTempCodaLists;
  SEXP AFD;    // AFD is a "AFullDiallelObject" - in case one wants it attached to MBS
    AObject *RCodaList, *RAFD;
 //   AObject **RCodaChainsList;  int RCodaChainsList;
  
  //double MaxX;

  //CodaIFile: Integer Draws of coordinates j for nonzero Beta_j per iter-t
  //CodaJFile: Draws of Double values of non-zero Beta_j
  //CodaiTFile: Draws of groups of coordinates k for nonzero tau^2_k per iter-t
  //CodadTFile: Double draws of groups of nonzero- tau^2_k
  //CodaPiAFile: Draws of PiA (either for both or one Fixed one Random)
  //CodaSigFile: Draws of sigma (for noise) 
  SEXP sCodaIFile, sCodaJFile; 
  SEXP sCodadTFile, sCodaiTFile;
  
 
  SEXP  sCodaPFile;
  SEXP sSaveDir;  AObject *RsSaveDir;
  RRObject *RsCodaIFile, *RsCodaJFile;
  RRObject *RsCodadTFile, *RsCodaiTFile, 
    *RsCodaPFile;
  RRObject *RsCodaiTLocFile;    

  RRObject *RsCodaILocFile;  RRObject *RsCodaProbFile;
  RRObject *RsPostProbBufferFile;
  RRObject *RsOldCodaILocFile;  RRObject *RsOldCodaProbFile;
  RRObject *RsCodaOldIFile;  RRObject *RsCodaOldJFile;
  RRObject *RsPiACodaFile; 
  
  double *PiACodaBuffer;   int LengthPiACodaBuffer;
  int LengthWrittenPiACodaBuffer;  int NewWritePiACodaBuffer;
  int LengthTotalWrittenPiACodaBuffer;

  RRObject *RsCodaBetaAllDrawBuffer; 
  RRObject *SubCodaLongList;   
  RRObject *RsSigCodaFile;  double *SigCodaBuffer;   int LengthSigCodaBuffer;
  int LengthWrittenSigCodaBuffer;  int NewWriteSigCodaBuffer;
  int LengthTotalWrittenSigCodaBuffer;
  

  double *BetaAllDrawBuffer;
  int LengthWrittenBetaAllDrawBuffer;
  int LengthTotalWrittenBetaAllDrawBuffer;
  int LengthBetaAllDrawBuffer;
  int WriteToBetaAllDrawBuffer();
  int ClearToBetaAllDrawBuffer();
       
  
  AObject *RsWeightFile, *RsYFile;     
  double *WeightBuffer, *YBuffer;
  int LengthWeightBuffer, LengthYBuffer, LengthWrittenWeightBuffer, LengthWrittenYBuffer;   
  int LengthTotalWrittenWeightBuffer;  int LengthTotalWrittenYBuffer;
  int NewWeightBufferWrite, NewYBufferWrite;
  AObject *RsYCodaList, *RsWeightCodaList, *RsPiACodaList, *RsSigCodaList;
  
  AObject *RsAlterWeightFile;
  AObject *RsAlterWeightCodaList;
  double *AlterWeightBuffer;  int LengthAlterWeightBuffer; int DoLogitNonePostPreProb;
  int LengthWrittenAlterWeightBuffer;
  int LengthTotalWrittenAlterWeightBuffer, NewAlterWeightBufferWrite;  
  double AlterWeightdfRobit, AlterWeightdfTNoise;
  double AlterWeightTemperature;   int AlterInitRun;
    
  ////////////////////////////////////////////////////////////////////////////
  // Files to write running algorithm steps to drive.
  //
  //  Many of these files are compressed with ILoc, IFile, JFile
  //  ILoc File is the integer location of each sample in IFile, JFile
  //  IFile is integer identifier of each activated coefficient
  //  JFile is value of these activated coefficients.
  //
  //  OldFiles are being drawn from for EquiEnergy buffer and were written
  //   in a previous chain at a higher temperature
  //
  //  ProbCoda contains the posterior model inclusion draws.

  double *ProbOldCoda;  double *ProbCoda; 
  bufferILocType *ICoda; bufferILocType *IOldCoda;  int LengthOldCoda;
  int LengthProbCoda;   int LengthPostProb;  int LengthWrittenPostProb;
  
  int LengthPostProbCodaBuffer; double *PostProbCodaBuffer;
  int NewPostProbWrite; int LengthTotalWrittenPostProb;
  int iiOnILocProb;  long int LocationOfBetaWrite;  long int MaxLengthOldIDFiles;
  double OldTemperature;  double invOldTemperature;   int IterOldCoda;
  int SampleTausOnly;
    
  //  A running total of Posterior Model Inclusion probabilities, which we
  //   Average for RMIP
  double *RunProbVector; int StartRunProbVector;  int TotInProbVector; 
  double *RunProbRegionVector;
  double *RunSumBeta;
  unsigned int *RunMoreZero;  unsigned int *RunLessZero;
  //double *MedianUp;  double *MedianDown;
  int *TotEveryProbVector; int DoShortMIP;
  int RegionWidth;  // A Length for +- X coordinates in case one wishes to combine
                    //   associated regions.
   
  int OnCodaTable;         // Which CodaTable of CodaList are we writing too?  
  int ChainIter;    
  int TotalOpenedFiles, TotalR5OpenedFiles;
  int TotalClosedFiles, TotalR5ClosedFiles;

  //  Prior for ProbFixed overruns sOnPiA if it is filled
  SEXP sPriorProbFixed;
  SEXP sPriorProbTau; 
  double gPrior; 
  AObject *rPriorProbFixed;  int HadSetPriorProbFixed;
  AObject *rPriorProbTau;   int HadSetPriorProbTau;
  AObject *RPostProbCodaList; // Posterior Probabilities of Random Effects
  AObject *RMIP;              // Model Inclusion Probabilities uses RunProbVector
  AObject *RegionMIP;          // Model Inclusion Probabilities using RunProbVector for regions using RegionWidth
    
public :

  int BeingDestroyed;//  Flag whether this object is in the process of destruction
  int IamBayesSpike; // Simple flag noting this is a valid BayesSpikeCL Object
  int PrintIter;     // If Verbose >= 0, prints Report on this Iter;
  int ReorderIter;   // Reorders Indices on this iteration
  int NewWrite;      // Flag whether to Write Beta vector draws to new file on harddrive or append
  int RevertTemperatureEvery; // Flag how often to swich low temperature to high temperature.
  
  //  EEProbSort:  Equi-Energy sampling  How many energy levels do we need.
  int DoingEEProbSort; double EEProbSortWidth;
  
  //  fd, fd1, fd2: These are file location markers for writing to file
  int fd, fd1, fd2;
  int RanRecordCoda, RanFillTau, WroteTotalTau, WroteAllTau;
  
  int MaxCODABeta;  int DontPartBetaNoise;

  int NInstance; // Number of Object as locked to Global Environment

  // Increase Memory Size by Factor
  int Resize(double Factor); 
    
  ////////////////////////////////////////////////////////////////////////
  //  Functions of Algorithm: 
  //    
  //   Most of these are executed in RunAlgorithm() which is defined
  //    in BayesSpikeGibbs.cpp
  // 
  int SetupMTForOnii(int iOnii); // Setup MT so that a Tau of current Group Onk can be drawn  
  int SampleNewTaus();  // SampleNewTaus() is in "BayesSpikeTaus.cpp"                        
  //  UpdatePiA, UpdateSigma:  New draws for PiA, Sigma
  int AddAllNewCoords();  int WillAddNewTau; int AllNewCoords; int NewFixedCoords;
  int TestAllNewCoords();

  int DoSampleBFixedii;
  int UpdatePiA();  
  int LenRandom, LenFixed, CountFreeFixed, CountFreeRandom, OnShrinkFixed, OnShrinkRandom;  
  
  int UpdateSigma();
  int SampleFixedB();     // Samples B_j indicators for fixed effects.
  int SamplePropBeta();   // Samples Beta for current active set.
  int PrepareForRegression();  // Prep for SampleProbBeta()
  
  int OnKSqPBack;  int NumActiveBack;  int JustDidPartBeta;
  int JustDidPartBFixed;   double BackSigma;
  
  
  //  RunAlgorithm in "BayesSpikeGibbs.cpp": is chief algorithm
  int RunAlgorithm();
  
  ///////////////////////////////////////////////////////////////////
  //  Timing of algorithm using clock ticks.
  //  Must define DOTIMEH at top of this page to activate this code.
  SEXP get_RsTimePartsList();   void set_RsTimePartsList(SEXP iIn);
  int SetupRsTimePartsList();
  
  //  Kills off coefficients with low MIP
  int GiveItAChanceAndKillUselessCoeficients(int PleaseKill);
  double KillPct; // How many to kill per iteration.
  double MinimumProbKeeper;
  int MaximumAllocation;
  int StartRemoval; int ItersSinceRemoval;
  int CheckRemovable;  int DoAddCoordsOnSetBeta; 
  
  //  TweakCons is set for SampleNewTaus
  double get_TweakCons() { return(TweakCons); }

  // Resamples the active Beta
  int SamplePartBeta(int StF, int MaxDo); int SetupRICholBack();
  int SamplePartBeta2(int StF, int MaxDo);  int LastStF, LastGoFor;
  int CopyIn(int StF, int GoFor);   int CopyIn2(int StF, int GoFor);
  int FillPropBetaFromBeta();  // Fills all Beta with active Beta
  int UpdateTNoise();  // Update weights for t-noise

  int RecordHistory(); // Save current state to buffers/files
  int UpdateXtY();     // Update XtY vector for new info
  int RefillXtX();     // Refill XtX SubMatrix for new weights
  int ReorderXtX();    // Reorder XtX SubMatrix if out of order by kXFinder;
  int AddCoord(int NewCoord);  //  Add a new coordinate to XtX SubMatrix
  int UpdateFreshXtResid();    // Refresh XtResid
  int UpdateOnlyDeleteXtResid(); // Delete members of XtResid where BFixed = 0 but sBeta != 0
  int SetXtResidNoXtX(); // Create XtResid without using XtX, found in 
  int LengthWrittenTauCodaBuffer, LengthTotalWrittenTauCodaBuffer;
  
  int TestMemAllocateFailure();  int TestXtResid();
 
  // Flag whether to start with zeroes (for big p)
  //  Code for this is actually done in R n ABayesSpike.r
  int ZeroOutBeforeSample; 
  int AssignTau(SEXP sAssTau);   // Assigns Tau to Zeroed out Betas without changing Beta.
  
  int SetCodaFileNames(SEXP sCodaIFile_, SEXP sCodaJFile_, SEXP siFileCodaChainIter) {
    if (Verbose >= -1) {
      Rprintf("BayesSpikeCpp.h:SetCodaFileNames(): Will sset sCodaIFile, JFile, siFileCodaChainIter set\n");
      R_FlushConsole();
    }
    if (Rf_isNull(sCodaIFile_) || Rf_isNull(sCodaJFile_)) {
      if (Verbose >= 2) {
        Rprintf("SetCodaFileNames: I/J/Coda is to cause NULL effect\n");
        R_FlushConsole();
      }
      sCodaIFile = R_NilValue;  sCodaJFile = R_NilValue;  OnCoda = 0;
      FFree(CodaBeta, "CodaBeta");  FFree(Codajj, "Codajj");
      return(1);
    }
    if (RsCodaIFile != NULL && RsCodaOldIFile == NULL) {
      if (Verbose >= 2) {
        Rprintf("SetCodaFileNames: OldI and set to old, IFile to NULL!\n");
      }
      RsCodaOldIFile = RsCodaIFile;  RsCodaIFile = NULL;
    } else if (RsCodaIFile != NULL) {
      DDelete(RsCodaOldIFile, "RsCodaOldIFile"); 
      RsCodaOldIFile = RsCodaIFile;  RsCodaIFile = NULL;
    }
    if (RsCodaJFile != NULL && RsCodaOldJFile == NULL) {
      RsCodaOldJFile = RsCodaJFile;  RsCodaJFile = NULL;
    }  else if (RsCodaJFile != NULL) {
      DDelete(RsCodaOldJFile, "RsCodaOldJFile");
      RsCodaOldJFile = RsCodaJFile; RsCodaJFile = NULL;
    }
    if (Rf_isInteger(sCodaIFile_)) {
      Rf_error("SetCodaFileNames: error, sCodaIFile is an Integer?\n");
    } else if (Rf_isReal(sCodaIFile_)) {
      Rf_error("SetCodaFileNames: error, sCodaIFile is an Real?\n");
    } else if (!Rf_isString(sCodaIFile_)) {
      Rf_error("SetCodaFileNames: error, sCodaIFile to be set but not string?\n");
    }
    if (Rf_isInteger(sCodaJFile_)) {
      Rf_error("SetCodaFileNames: error, sCodaJFile is an Integer?\n");
    } else if (Rf_isReal(sCodaJFile_)) {
      Rf_error("SetCodaFileNames: error, sCodaJFile is an Real?\n");
    } else if (!Rf_isString(sCodaJFile_)) {
      Rf_error("SetCodaFileNames: error, sCodaJFile to be set but not string?\n");
    }
    DDelete(RsCodaIFile, "RsCodaIFile");
    DDelete(RsCodaJFile, "RsCodaJFile");
    if (Verbose >= -1) {
      Rprintf("BayesSpikeCpp.h:SetCodaFileName(): Begin new objects for sCodaIFile, sCodaJFile. \n");
      R_FlushConsole();
    }
    RsCodaIFile = new RRObject(sCodaIFile_);
    sCodaIFile = RsCodaIFile->asSexp();  
    RsCodaJFile = new RRObject(sCodaJFile_);
    sCodaJFile = RsCodaJFile->asSexp(); OnCoda = 0;
    FileCodaChainIter = GetFirstInteger(siFileCodaChainIter);
    LocationOfBetaWrite = (long int) 0; NewWrite = 1; 
    FFree(CodaBeta, "CodaBeta");  FFree(Codajj, "Codajj");
    if (MaxCODABeta <= 0) {
      MaxCODABeta = DefaultMaxCODA;
    }
    if (MaxCODABeta < Rf_length(sBeta)+50) {
      MaxCODABeta = Rf_length(sBeta) + 50;
    }
    if (MaxCODABeta >= 3*TooMuchBuffer) {
      Rprintf("Error: Hey: MaxCODABeta=%d, but TooMuchBuffer=%d, maybe better space work?\n",
        MaxCODABeta, TooMuchBuffer);
      MaxCODABeta = TooMuchBuffer * 3;
    }
    //if (GoMaxCODA > 3*TooMuchBuffer) {
    //  Rprintf("Error: Hey: MaxCODA=%d, but TooMuchBuffer=%d, maybe better space work?\n",
    //    GoMaxCODA, TooMuchBuffer);
    //  Rf_error("Error: MaxCODA=%d, TooMuchBuffer=%d\n", MaxCODA, TooMuchBuffer);
    //}
    RMemGetD(CodaBeta, "CodaBeta", MaxCODABeta);
    RMemGetI(Codajj, "Codajj", MaxCODABeta);
    if (CodaBeta == NULL || Codajj == NULL) {
      Rprintf("BayesSpikeCpp.h: Error issue, Codajj or CodaBeta did not allocate!\n");
      R_FlushConsole();
    }
    if (Verbose >= 3) {
      Rprintf("------------------------------------------------------\n");
      Rprintf("SetupSaveFile: We have that IFile/JFile were set up.\n");
      R_FlushConsole();
      Rprintf("  We got IFile is %s\n",
        CHAR(STRING_ELT(sCodaIFile,0))); R_FlushConsole();
      Rprintf("  We got JFile is %s\n",
        CHAR(STRING_ELT(sCodaJFile,0))); R_FlushConsole();      
    }
    return(1);
  }
  int get_LocationOfBetaWrite() {return((int) LocationOfBetaWrite); }
  int SetupBFixedAndFriends(); // Setup BFixed vector as well as OrderedActive Vector
  int RefreshOrderedActive(int ForceFill);  // Fill OrderedActive as set of BFixed>0, Ontau>0 coordinates in order
  int FillSmallXtX();  //  Fill SmallXtX compact triangular matrix using XtX and OrderedActive
  int FillSmallXtResid();  //  Fill SmallXtX compact triangular matrix using XtX and OrderedActive
  int FillsBetaFromPropBetaAndCompute(); // Restore sBeta from PropBeta, compute residuals
  
  double ThisMaxTau;  int burnin;
  
  int  OnTauIndex;  int get_OnTauIndex() { return(OnTauIndex); }

  SEXP ReadFileCodas(SEXP StartIter, SEXP EndIter,
    SEXP SubSetCoords);
    
  SEXP ReadSetPctBetas(SEXP StartIter, SEXP EndIter);
  

  
  double ProbLDraw();
  
  SEXP sLOfX(SEXP SampleStart);
  
  //////////////////////////////////////////////
  // Running Percentages
  double *PctBetas;  double *PctTaus;
  
  
  int EEMergeEvery; // If temperature list is set, asn ProbOldCoda, then merge from old distributions every the iters
  int EESamplerMerge(); // Run that Merge
  int TestFindInMe(SEXP FindMe);
  
  // Reweighting of XtX,  If t Statistics we need to reweight XtX
  int ReweightedOneXtX(int iOneX);
  int ReweightOnlyActiveXtX();   int QuickCheckWeightXtX();
  int BackReweightOnlyActiveXtX();  
  
  void set_Tempii(int iTempii) {
    if (TemperatureList == NULL) { Rf_error(" No SetTempii with no TemperatureList! "); }
    if (LengthTemperatureList < iTempii || iTempii < 0) {
      Rf_error("set_Tempii: Sorry, iTempii not valid; ");
    }
    Tempii = iTempii;
    Temperature = TemperatureList[Tempii]; invTemperature = 1.0 / Temperature;
    if (Tempii == 0) {
      OldTemperature = TemperatureList[Tempii];  invOldTemperature = 1.0 / OldTemperature;
    } else if (TemperatureList[Tempii-1] != 1.0) {
      OldTemperature = TemperatureList[Tempii-1]; invOldTemperature = 1.0 / OldTemperature;
    }
    if (MT != NULL) {
      MT->Temperature = Temperature;  MT->invTemperature = invTemperature;
    }
    if (TemperatureDecreasingRate < 1.0) {
      CurrentMaxGibbsIters = (int) round(MaxGibbsIters * 
      powl(TemperatureDecreasingRate, LengthTemperatureList-Tempii-1));
      if (LengthTempDraws != NULL) {
        LengthTempDraws[Tempii] = CurrentMaxGibbsIters;
      }
    }
    return;
  }
  int get_Tempii(){ return(Tempii); }
  SEXP get_LengthTempDraws() {
    if (LengthTempDraws == NULL) { return(R_NilValue); }
    if (LengthTemperatureList  <= 0) { return(R_NilValue); }
    SEXP sOut = R_NilValue;
    Rf_protect(sOut = Rf_allocVector(INTSXP, LengthTemperatureList));
    if (sOut == NULL || Rf_isNull(sOut)) {
      Rf_error("Error: get_LengthTempDraws, could not allocate\n");
    }
    int ii = 0;
    if (Tempii >= LengthTemperatureList) { Tempii = LengthTemperatureList -1; }
    if (Tempii < 0) { Tempii = 0; }
    for (ii = 0; ii < Tempii+1; ii++) {
      if (LengthTempDraws[Tempii] < 0) {
        if (TemperatureDecreasingRate < 1.0) {
          LengthTempDraws[Tempii] = (int) round(MaxGibbsIters * 
            powl(TemperatureDecreasingRate, LengthTemperatureList-ii-1));
        }
      }
    }
    for (ii = 0; ii < LengthTemperatureList; ii++) {
      INTEGER(sOut)[ii] = LengthTempDraws[ii];
    }
    Rf_unprotect(1);
    return(sOut);
  }
  SEXP get_Temperature() { SEXP sOn = R_NilValue;
    Rf_protect(sOn = Rf_allocVector(REALSXP, 1));
    REAL(sOn)[0] = Temperature;  Rf_unprotect(1); return(sOn); }
  double get_invTemperature() { return(invTemperature); }
  double get_OldTemperature() { return(OldTemperature); }
  void set_OldTemperature(double iOldTemperature) {
    if (iOldTemperature <= 0) {
      Rf_error("Sorry, we can't set old temperature to less than zero");
    }
    OldTemperature = iOldTemperature;  invOldTemperature = 1.0  / OldTemperature;
  }
  double get_invOldTemperature() { return(invOldTemperature); }
  SEXP get_TemperatureList() {
    if (TemperatureList == NULL) { return(R_NilValue); }
    SEXP sOut = R_NilValue;
    Rf_protect(sOut = Rf_allocVector(REALSXP, LengthTemperatureList));
    for (int jj = 0; jj < LengthTemperatureList; jj++) {
      REAL(sOut)[jj] = TemperatureList[jj];
    }
    Rf_unprotect(1); return(sOut); 
  }
  void set_TemperatureList(SEXP iTemperatureList) {
    if (Rf_isNull(iTemperatureList)) {
      if (Verbose > 2) {
        Rprintf("Note: TemperatureList input is NULL! \n");
        R_FlushConsole();
      }
      FFree(TemperatureList, "TemperatureList");
      FFree(LengthTempDraws, "LengthTempDraws"); FFree(NumMerges, "NumMerges");
      TemperatureList = NULL;  return;
    }
    if (Verbose >= -4) {
      Rprintf("set_TemperatureList: Starting \n"); R_FlushConsole();
    } 
    if (!Rf_isReal(iTemperatureList)) {
      Rf_error("Error set_TemperatureList, iTemperatureList is not real list!");
    }
    if (Rf_length(iTemperatureList) < 1) {
      Rf_error("Error: set_TemperatureList, give something of length greater than zero");
    }
    int One = 1;
    if (TemperatureList != NULL) {
      if (Verbose >= 2) {
        Rprintf("set_TemperatureList: Weird, TemperatureList is Not NULL!\n");
      }
    }
    FFree(TemperatureList, "TemperatureList");
    FFree(LengthTempDraws, "LengthTempDraws"); FFree(NumMerges, "NumMerges");
    if (LengthTempDraws != NULL) {
    FFree(LengthTempDraws, "LengthTempDraws");   }
    if (REAL(iTemperatureList)[Rf_length(iTemperatureList) - 1]  != 1.0) {
      LengthTemperatureList = Rf_length(iTemperatureList) +1;
      if (Verbose > 1) {
        Rprintf("Last member is not one, allocating LengTemperatureList = %d",
          LengthTemperatureList); R_FlushConsole();
      }
      RMemGetD(TemperatureList, "TemperatureList", LengthTemperatureList);
      int pT = Rf_length(iTemperatureList);
      F77_CALL(dcopy)(&pT, REAL(iTemperatureList), &One, TemperatureList, &One);
      TemperatureList[LengthTemperatureList -1] = 1.0;
    } else {
       LengthTemperatureList = Rf_length(iTemperatureList);
      if (Verbose > 2) {
        Rprintf("Allocating regular TemperatureList, length %d \n",
          LengthTemperatureList); R_FlushConsole();
      }
      RMemGetD(TemperatureList, "TemperatureList", LengthTemperatureList);      
      F77_CALL(dcopy)(&LengthTemperatureList, REAL(iTemperatureList), &One, TemperatureList, &One);
    }
    if(Verbose >= -2) {
      Rprintf("SetTemperatureList, here we are at the end , Alloc NumMerges, LengthtempDraws\n"); R_FlushConsole();
    }
    Tempii = 0;  Temperature = TemperatureList[0];  invTemperature = 1.0 / Temperature;
    FFree(NumMerges, "NumMerges");
    NumMerges = (int*) Calloc( LengthTemperatureList, int);
    for (int iti = 0; iti < LengthTemperatureList; iti++) { NumMerges[iti] = 0; }
    LengthTempDraws = (int* )Calloc(LengthTemperatureList, int);
    for (int iti2 = 0; iti2 < LengthTemperatureList; iti2++) { LengthTempDraws[iti2] = -1; }
    if (Verbose >=2) {
      Rprintf("SetTempreatureList: Reached End \n"); R_FlushConsole();
    }
    return;
  }
  SEXP get_NumMerges() {
    if (TemperatureList == NULL || NumMerges == NULL) {
      return(R_NilValue);
    }
    SEXP sOn = R_NilValue;  Rf_protect(sOn = Rf_allocVector(INTSXP, LengthTemperatureList));
    for (int iti = 0; iti < LengthTemperatureList; iti++) {
      INTEGER(sOn)[iti] = NumMerges[iti];
    }
    Rf_unprotect(1); return(sOn);
  }
  int LoadPCodaICoda(); //PCodaICoda, record info from previous chain (old Temperature),
     //  Useful for ladder draws
     
    
  double ProbPosterior(); // In BayesSpikecpp.cpp Calculate current Posterior probability of eleemnts
  int RobitReplace();  // In BayesSpikeGibbs.cpp : For Robit regression, replace Y with truncated t draws
  int RobitNotReplaceFlag;
  double MaximizeOn(); // Attempts to find tau^2 that maximizes P(tau^2_j | y, \beta_{\backslashJ(k)}
  double ValAtMaximum(); double get_Df3(); double get_Df4();
  double get_lUB32(); double get_lLB42();
  double get_lM4xM4(); double get_lM3xM3();
  double get_xM4(); double get_xM3(); //double get_xM2();
  
  double get_UsePrior();  double get_InitlA();  double get_lA();  double InitlA;
  double get_lMyHLB();  double get_NewlA();  double lA; double q4X4; double f2X4;
  double lB;  double get_lB() { return(lB); }
  int ItsNewX4;  double ScaleOfSmallXtResid;  double ScaleOfEigenVectors;
  double X4;  double get_NewX4(); double get_OldX4() { return(X4); }
  double get_f2X4();  double get_q4X4();
  double get_ScaleOfSmallXtResid() { return(ScaleOfSmallXtResid); }
  double get_ScaleOfEigenVectors() { return(ScaleOfEigenVectors); }
    
    
  int NukeBetaTau, NukeBetaNotTau;
  int TestOrderedActive();
  double GiveLOfX(double SampleStart);
  SEXP BPriorLOfX(SEXP BVector);
  double BDerivlF(double Outtau);
  double BDerivPrior(double Outtau);
  double BSDerivlF(double Outtau); 
  double BSDerivPrior(double Outtau);

  int RecordPandT(); int WriteP(int sP, int eP); 
  int WriteTauFromTable(int sT, int eT);
  int WriteTau();  int FillTau();
  int RecordCoda();  
  int WriteCoda(); 
  
  int NoSave; int DoSave;
  
  
  double BinFindLOfX(double StartX, double GoalVal);
  double BinFindROfX(double StartX, double GoalVal);
  double get_MaxX();
  SEXP GRunif(int Counts);
  SEXP BAllDeriv(SEXP sIn);
  SEXP BAllSDeriv(SEXP sIn);
  
  int UpdateNonFreshXtResid();   // Found in 
  
  int HowSample; int NumIntegrations; double FakeOuttau;
  double IntegratedDensity; 
  int NumFails;  int MaxAlsErrors; int GoBackAfter;
  int TotalFails;
  int NumlAZero;  int TotalLAZero;  int GoBackAfterlA;  int MaxlAErrors;
  
  int tt; int NumSteps;
  
  int MaxIters;  double CauchyEpsilon; int MaximizeMeCauchyTotalIters;
  double SpliceMove;  int NumSpliceSampleReps;  int DoMax;
  double MaxTauEst;  int MaxGibbsIters; int CurrentMaxGibbsIters;
  
  int TypePrior; int NumberNew;
  
  int EarlyEndStep;  int EarlyEndtt;
  double nWeight;
  int RecrodPandT();  int NewTWrite;  int NewPWrite;
  int StartIter, EndIter;
  void Finalizer();
  void DeleteMe();
  int AmDestroying;  
  int get_FileCodaChainIter() { 
    if (FileCodaChainIter < 0 ||  FileCodaChainIter >= 1000) {
      Rprintf("get_FileCodaChainIter, we have it is %d which is probably not what you want.\n", FileCodaChainIter); R_FlushConsole();
      return(FileCodaChainIter);
    }
    return(FileCodaChainIter); 
  }  
      
  ///////////////////////////////////////////////////////////////////
  // Destructor
  ~BayesSpikeCL() {
     Verbose = 5;
     
     int PTECTED = 0;
     Rprintf("----------------------------------------------------------------\n");
     Rprintf("-  Still a Big Deal: We're deleting the BayesSpikeCL, Enjoy !   \n");
     Rprintf("-     Oh what fun it will be. \n"); 
     if (RsSaveDir == NULL || Rf_isNull(sSaveDir)) {
       Rprintf("-  But there is only null SaveDir! \n"); R_FlushConsole();
     } else if (!Rf_isString(sSaveDir)) {
       Rprintf("- But SaveDir is not a satisfactory string! \n"); 
     } else {
       Rprintf("-  Our Save Dir is %s.\n",  CHAR(STRING_ELT(sSaveDir, 0)) ); 
     }
     if (TotalOpenedFiles != TotalClosedFiles) {
Rprintf("***************************************************************************************** \n");
Rprintf("**  BAD ERROR ON THE Destruction of BayesSpike Object!!!!!!!! \n");
Rprintf("**\n");
Rprintf("**\n");
Rprintf("**  We have that TotalOpenedFiles = %d, TotalClosedFiles = %d, quit on error!\n"); R_FlushConsole();     
     Rf_error("Error On a BayesSpike Destructor. \n", 
     TotalOpenedFiles, TotalClosedFiles);    
     }
     if (TotalR5OpenedFiles != TotalR5ClosedFiles) {
Rprintf("***************************************************************************************** \n");
Rprintf("**  BAD ERROR ON THE Destruction of BayesSpike Object!!!!!!!! \n");
Rprintf("**\n");
Rprintf("**\n");
Rprintf("**  We have that TotalR5OpenedFiles = %d, TotalR5ClosedFiles = %d, quit on error!\n",
  TotalR5OpenedFiles, TotalR5ClosedFiles); R_FlushConsole();     
     Rf_error("Error On a BayesSpike Destructor. \n");    
     }
     Rprintf("-  Note, ABayesSpikeCL, We're getting Destroyed now and I think thats a big deal\n");
     Rprintf("---------------------------------------------------------------------------------\n");
     DeleteMe();
     R_FlushConsole();
     if (Verbose > 0) {
       Rprintf("~BayesSpikeCL:  Destructing.\n"); R_FlushConsole();
     }
     if (Verbose > 0 && BeingDestroyed ==1) {
       Rprintf("~BayesSpikeCL: Won't destroy, because,we're being destroyed"); 
       R_FlushConsole(); return;
     }
     AmDestroying = 1;
     BeingDestroyed = 1;
     
     DDelete(RsBeta, "RsBeta");   DDelete(RsRandomInfoList, "RsRandomInfoList");
     DDelete(RsTimePartsList, "RsTimePartsList");
    if (RBayesSpikeNameSpace != NULL && !Rf_isNull(BayesSpikeNameSpace)) {
      ////////////////////////////////////////////////////////////
      //  Re-enter R and Kill TBSRoo
      if (TBSRoo != NULL && !Rf_isNull(TBSRoo) && TBSRoo != NULL) {
      SEXP FreeTBSRoo = Rf_protect(Rf_findVarInFrame( BayesSpikeNameSpace, Rf_install("FreeTBSRoo")));
      PTECTED++;
      if (Rf_isNull(FreeTBSRoo)) {
        Rf_unprotect(PTECTED);
        Rf_error("BayesSpikeCL: Error: no FreeTBSRoo");
      }
      SEXP call = Rf_protect(Rf_lcons( FreeTBSRoo, TBSRoo));
      if (!Rf_isNull(call)) {
        Rf_unprotect(PTECTED);
        Rf_error("BayesSpikeCL Destroy: Failed to delete TBSRoo");
      }
      SEXP RTM = R_NilValue;
      Rf_protect(RTM = Rf_eval(call, BayesSpikeNameSpace));
      Rf_unprotect(3);    PTECTED = 0;
      } 
      SEXP FreeTBSR5=R_NilValue;  SEXP call = R_NilValue;  SEXP RTM = R_NilValue;
      SEXP sVerbose = R_NilValue;
      if (RTBSR5 != NULL && !Rf_isNull(RTBSR5->asSexp())) {
         Rf_protect(FreeTBSR5 = Rf_findVarInFrame( BayesSpikeNameSpace, Rf_install("KillFreeTBSR5")));
         PTECTED++;
         Rf_protect(sVerbose = Rf_allocVector(INTSXP,1));  PTECTED++;
         INTEGER(sVerbose)[0] = Verbose;
         if (Rf_isNull(FreeTBSR5)) {
           Rprintf("~ABayesSpikeCL(): Error FreeTBSR5 cannot find function. \n");
         } else {
           call = R_NilValue;
           if (Verbose > 0) {
             Rprintf("DeleteMe: We still managed to find TBSR5 to try to kill!\n");
           }
           Rf_protect(call = Rf_lcons( FreeTBSR5, 
           Rf_cons(TBSR5, 
             Rf_cons(sVerbose, R_NilValue)))); PTECTED++;
           if (Rf_isNull(call)) {
             Rf_unprotect(PTECTED); PTECTED = 0;
             Rf_error("BayesSpikeCL Destroy: Call is a null to delete TBSR5\n");
           }
           RTM = R_NilValue;
           if (Verbose > 0) {
             Rprintf("DeleteMe: We've got the call to delete TBSR5, now getting to work\n "); R_FlushConsole();
           }
           Rf_protect(RTM = Rf_eval(call, get_BayesSpikeNameSpace())); PTECTED++;
           if (Verbose > 0) {
             Rprintf("DeleteMe: Just conducted FIRST RTBSR5 KillFree on current TBSR5. \n");
           }
           DDelete(RTBSR5, "RTBSR5");  TBSR5 = R_NilValue;
         }
      }
      SEXP tTBSR5 = get_TBSR5();
      if (!Rf_isNull(tTBSR5)) {
         FreeTBSR5 = Rf_protect(Rf_findVarInFrame( BayesSpikeNameSpace, Rf_install("FreeTBSR5")));  PTECTED++;
         if (Rf_isNull(FreeTBSR5)) {
           Rprintf("~ABayesSpikeCL(): Deleter, cannot go further. Error we tried to find FreeTBSR5 in both spaces but failed.");
           Rf_error("BayesSpikeCL: Error: no FreeTBSR5");
         }
         Rf_protect(sVerbose = Rf_allocVector(INTSXP,1));  PTECTED++;
         INTEGER(sVerbose)[0] = Verbose;
         call = Rf_protect(Rf_lcons( FreeTBSR5, TBSR5));   PTECTED++;
         if (!Rf_isNull(call)) {
          Rf_error("BayesSpikeCL Destroy: Failed to delete TBSR5");
         }
         RTM = R_NilValue;
         if (Verbose >= 2) {
           Rprintf("~ABayesSpikeCL(): Deleter: About to run FreeTBSR5. \n");
           R_FlushConsole();
         }
         Rf_protect(RTM = Rf_eval(call, BayesSpikeNameSpace));  PTECTED++;
         Rf_unprotect(3);      
         if (Verbose >= 2) {
           Rprintf("~ABayesSpikeCL(): Deleter: We have ran FreeTBSR5 and it should be deleted. \n");
           R_FlushConsole();
         }
         Rf_unprotect(PTECTED); PTECTED = 0;
       }
       tTBSR5 = get_TBSR5();
       if (!Rf_isNull(tTBSR5)) {
         if (Verbose >= 0) {
           Rprintf("~ABayesSpikeCL(): We have to play clear off for tTBSR5.\n");
         }
           Rf_protect(FreeTBSR5 = Rf_findVarInFrame( BayesSpikeNameSpace, Rf_install("KillFreeTBSR5")));
           PTECTED++;
           if (Rf_isNull(FreeTBSR5)) {
             Rprintf("~ABayesSpikeCL(): Error FreeTBSR5 cannot find function. \n");
           } else {
             Rf_protect(sVerbose = Rf_allocVector(INTSXP,1));  PTECTED++;
             INTEGER(sVerbose)[0] = Verbose;
             call = R_NilValue;
             if (Verbose > 0) {
               Rprintf("DeleteMe: We still managed to find TBSR5 to try to kill!\n");
             }
             Rf_protect(call = Rf_lcons( FreeTBSR5, 
               Rf_cons(TBSR5, 
               Rf_cons(sVerbose, R_NilValue)))); PTECTED++;
             if (Rf_isNull(call)) {
               Rf_unprotect(PTECTED); PTECTED = 0;
               Rf_error("BayesSpikeCL Destroy: Call is a null to delete TBSR5\n");
             }
             RTM = R_NilValue;
            if (Verbose > 0) {
               Rprintf("~ABayesSpikeCL(): We've got the call to delete TBSR5, now getting to work\n "); R_FlushConsole();
             }
             Rf_protect(RTM = Rf_eval(call, get_BayesSpikeNameSpace())); PTECTED++;
             if (Verbose > 0) {
               Rprintf("~ABayesSpikeCL(): Just conducted LAST RTBSR5 KillFree on current TBSR5. \n");
             }
             Rf_unprotect(PTECTED); PTECTED = 0;
          }
       }
      
       DDelete(RBayesSpikeNameSpace, "RBayesSpikeNameSpace");
       BayesSpikeNameSpace = R_NilValue;
     } else if (RBayesSpikeNameSpace != NULL) {
       DDelete(RBayesSpikeNameSpace, "RBayesSpikeNameSpace");
     }
     DDeleteA(RTBSRoo, "RTBSRoo"); DDeleteA(RTBSR5, "RTBSR5");
     TBSRoo = NULL;  TBSR5 = NULL;
     
     DDelete(RtauEndList, "RtauEndList"); tauEndList= R_NilValue; DDeleteA(RsY, "RsY"); sY = R_NilValue;
     DDeleteA(RsX, "RsX"); sX = R_NilValue;
     DDelete(RtauFixed, "RtauFixed");  DDeleteA(RsOnTau, "RsOnTau");  sOnTau = R_NilValue;
     DDelete(RtauFixedPrior, "RtauFixedPrior");
     DDelete(RsCenteredColumns, "RsCenteredColumns");
     DDelete(RsDeCenteredCodaList, "RsDeCenteredCodaList");
     DDelete(RsOnSigma, "RsOnSigma"); sOnSigma = R_NilValue; DDelete(RSigmaPrior, "RSigmaPrior");  SigmaPrior = R_NilValue;
     DDelete(RInitFlags, "RInitFlags"); DDelete(RsOnPiA, "RsOnPiA");  sOnPiA = R_NilValue;
     
     DDelete(RPiAPrior, "RPiAPrior");
     FFree(TemperatureList, "TemperatureList");
     FFree(NumMerges, "NumMerges");  FFree(LengthTempDraws, "LengthTempDraws");
     
     DDelete(RDependenciesFixed, "RDependenciesFixed"); DependenciesFixed = R_NilValue;
     DDelete(RDependenciesTau, "RDependenciesTau"); DependenciesTau = R_NilValue;
     DDelete(RBayesSpikeNameSpace, "RBayesSpikeNameSpace");
     DDelete(RAllEigenVectors, "RAllEigenVectors");
     DDelete(RAllEigenValues, "RAllEigenValues");  
     DDelete(Rstaupriordf, "Rstaupriordf");
     DDelete(Rstaubarnu, "Rstaubarnu");  
     DDelete(RsPriorXScaling, "RsPriorXScaling");
     DDelete(RsPriorOftau, "RsPriorOftau");     
     
     if (MergeActive != NULL) {
       LengthMergeActive = 0; FFree(MergeActive, "MergeActive"); 
       FFree(MergePropBeta, MergePropBeta);   MaxLengthMergeActive = 0;
     }
     
     DDelete(InitiateTime, "InitiateTime");
     DDelete(SecondInitiateTime, "SecondInitiateTime");
     DDelete(CompleteTime, "CompleteTime");
  
  
     DDelete(RCodaTable, "RCodaTable");  
     DDelete(ROldColNames, "ROldColNames"); 
     DDelete(ROldDeCenterColNames, "ROldDeCenterColNames"); 
     DDelete(RSubCodaList, "RSubCodaList");
     DDelete(RSubCodaSubSetCoords, "RSubCodaSubSetCoords");
     DDelete(ROtherNameCodaList, "ROtherNameCodaList"); 
     DDelete(RTauCodaList, "RTauCodaList");
     DDelete(RSubSetTau, "RSubSetTau");    
     DDelete(RDoRecord, "RDoRecord");  
     DDelete(RCodaList, "RCodaList");
     DDelete(RAllTempCodaLists, "AllTempCodaLists");
     DDelete(RAFD, "RAFD");     
     FFree(intY, "intY");
     
     DDelete(RsCodaIFile, "RsCodaIFile"); 
     DDelete(RsCodaJFile, "RsCodaJFile"); 
     DDelete(RsCodaOldIFile, "RsCodaOldIFile");
     DDelete(RsCodaOldJFile, "RsCodaOldJFile");
     DDelete(RsCodadTFile, "RsCodadTFile"); 
     DDelete(RsCodaiTLocFile, "RsCodaiTLocFile");
     DDelete(RsCodaiTFile, "RsCodaiTFile"); 
     DDelete(RsCodaPFile, "RsCodaPFile");
     
     DDelete(RsCodaILocFile, "RsCodaILocFile");
     DDelete(RsCodaProbFile, "RsCodaProbFile");
     DDelete(RsOldCodaILocFile, "RsOldCodaILocFile");
     DDelete(RsOldCodaProbFile, "RsOldCodaProbFile");   
     
     DDelete(RsPiACodaFile, "RsPiACodaFile");  
     FFree(PiACodaBuffer, "CoaPiABuffer"); 
     
     DDelete(RsSigCodaFile, "RsSigCodaFile");
     FFree(SigCodaBuffer, "SigCodaBuffer");  LengthSigCodaBuffer = 0;
  
     
     FFree(PostProbCodaBuffer, "PostProbCodaBuffer"); 
     DDelete(RsPostProbBufferFile, "RsPostProbBufferFile");
     DDelete(RPostProbCodaList, "RPostProbCodaList");
     DDelete(RMIP, "RMIP");    DDelete(RegionMIP, "RegionMIP");
     
     DDelete(rPriorProbFixed, "rPriorProbFixed"); 
     DDelete(rPriorProbTau, "rPriorProbTau");
     
     DDelete(RsWeightFile, "RsWeightFile");
     DDelete(RsYFile, "RsYFile");
     FFree(YBuffer, "YBuffer");  FFree(WeightBuffer, "WeightBuffer");
     DDelete(RsWeightCodaList, "RsWeightCodaList");
     DDelete(RsYCodaList, "RsYCodaList");
     DDelete(RsSigCodaList, "RsSigCodaList");
     DDelete(RsPiACodaList, "RsPiACodaList");
     
     FFree(AlterWeightBuffer, "AlterWeightBuffer");
     DDelete(RsAlterWeightCodaList, "RsAlterWeightCodaList");
     DDelete(RsAlterWeightFile, "RsAlterWeightFile");
              
     //FFree(MedianUp, "MedianUp");  FFree(MedianDown, "MedianDown");
     if (RunMoreZero != NULL) {  Free(RunMoreZero);  } 
     if (RunLessZero != NULL) { Free(RunLessZero);  RunLessZero = NULL; }
     FFree(RunSumBeta, "RunSumBeta");  FFree(RunSumBeta, "RunSumBeta");
     FFree(RunProbVector, "RunProbVector"); FFree(RunProbRegionVector, "RunProbRegionVector");
     FFree(TotEveryProbVector, "TotEveryProbVector");
     FFree(PropBeta, ("PropBeta" ));  FFree(rI, "rI");  FFree(rIChol, "rIChol");
     FFree(rICholBack, "rICholBack");
     MaxLargeSizeSquares= 0;
     FFree(rSmallXtX, "SmallXtX" );  FFree(XLC, "XLC"); FFree(cXtY, "cXtY"); 
     FFree(smXtY, "smXtY");
     FFree(kXFinder, "kXFinder" ); FFree(iiWeight, "iiWeight");
     FFree(SumiiWeight, "SumiiWeight");
     FFree(CodaBeta, "CodaBeta"); FFree(Codajj, "Codajj");
     FFree(XtResid, "XtResid");  
     if (Verbose >= 3) {
       Rprintf("Clearing pXtX!\n"); R_FlushConsole();
     }
     if (pXtX != NULL) {
       for (int ii = 0; ii < OnKappaMem; ii++) {
         if (pXtX[ii] != NULL) {
           if (Verbose >= 4) {
             Rprintf("Clearing pXtX[ii=%d]\n", ii); R_FlushConsole();
           }
           Free(pXtX[ii]);  pXtX[ii] = NULL;
         }
       }
       if (Verbose >= 3) {
         Rprintf("pXtX cleared internally, clearing at end.\n"); R_FlushConsole();
       }
       Free(pXtX); pXtX = NULL;
     }
     FFree(XtY, "XtY"); FFree(WW, "WW"); FFree(ProbFixed, "ProbFixed");  
     FFree(ProbTau, "ProbTau"); ffree(MT, "MT");  FFree(OrderAttack, "OrderAttack");
     FFree(SmallXtResid, "SmallXtResid");  FFree(SmallRVec, "SmallRVec");
     FFree(cXtY, "cXtY");  FFree(XjSq, "XjSq"); FFree(smXtY, "smXtY"); 
     FFree(CodaP, "CodaP");  FFree(CodaTau, "CodaTau"); CodaTau = NULL; CodaTisSetup = 0;
     FFree(CodaTjj, "CodaTjj"); 
     FFree(CodaTLocjj, "CodaTLocjj");
     FFree(NoShrinkFixed, "NoShrinkFixed");
     FFree(NoShrinkFixedPrior, "NoShrinkFixedPrior");
     FFree(NoShrinkRandom, "NoShrinkRandom");  
     FFree(NoShrinkRandomPrior, "NoShrinkRandomPrior");
     FFree(XtXRecorded, "XtXRecorded");
     FFree(XjSq, "XjSq"); FFree(BFixed, "BFixed");  FFree(iWeightedXtX, "iWeightedXtX");
     FFree(NumMerges, "NumMerges");   FFree(AlterWeightBuffer, "AlterWeightBuffer");
     
     DDelete(RsCodaBetaAllDrawBuffer, "RsCodaBetaAllDrawBuffer");
     DDelete(SubCodaLongList, "SubCodaLongList");
     FFree(BetaAllDrawBuffer, "BetaAllDrawBuffer");  
     //if (RCodaChainsList != NULL && LengthRCodaChainsList > 0) {
     //  for (tt  = 0; tt < LengthRCodaChainsList; tt++) {
    //     DDelete(RCodaChainsList[tt], "RCodaChainsList[tt]");
    //   }
    //   LenghRCodaChainsList = 0;
    //   DDelete(RCodaChainsList, "RCodaChainsList");
    //   RCodaChainsList = NULL;
    // }
     
     FFree(PctBetas, "PctBetas");  FFree(PctTaus, "PctTaus");
     

    
     if (Verbose > 0) {
       Rprintf("~BayesSpikeCL:  Finish Destruction\n\n "); R_FlushConsole();
     }
  }

  BayesSpikeCL(SEXP X_, SEXP Y_, SEXP Beta_, SEXP RandomInfoList_,
    SEXP InitFlags_, SEXP _BayesSpikeNameSpace) {
    
    AmDestroying = 0;
    RsX = NULL; RsY = NULL; RsBeta = NULL;
    sX = R_NilValue;  sY = R_NilValue;  sBeta = R_NilValue;
    XjSq = NULL;  CodaTau = NULL;   LengthWrittenTauCodaBuffer = 0;  
    CodaTisSetup = 0;    iiWeight = NULL; iWeightedXtX = NULL; pXtX = NULL;
    LengthTotalWrittenTauCodaBuffer = 0;   RsDeCenteredCodaList = NULL;
    ZeroOutBeforeSample=0;  KillPct = DefaultKillPct;
    MinimumProbKeeper = DefaultMinimumProbKeeper;  ItersSinceRemoval = 0;
    DoAddCoordsOnSetBeta = 0;  RunProbRegionVector = NULL;  RegionWidth = -1;
    MaximumAllocation = DefaultMaximumAllocation;
    StartRemoval = DefaultStartRemoval;
    CheckRemovable =DefaultCheckRemovable;   DontPartBetaNoise = 0;
    MaxCODABeta = DefaultMaxCODA;  LengthCodaTjj = 0;
    LastStF = 0; LastGoFor = 0;
    ScaleOfSmallXtResid=1.0;   ScaleOfEigenVectors = 1.0;
    RsTimePartsList = NULL;  AlterWeightBuffer = NULL;
    RanFillTau =0; RanRecordCoda = 0;   WroteTotalTau = 0;  WroteAllTau = 0;
    MaxInvert = DefaultMaxInvert;   RevertTemperatureEvery = -1;
    LengthMergeActive = 0; MaxLengthMergeActive = 0;  MergeActive = NULL; MergePropBeta = NULL;
    RsSaveDir = NULL; sSaveDir = R_NilValue;  SampleTausOnly = 0;
    if (Rf_length(InitFlags_) >= 2 ) {
      if (Rf_isInteger(InitFlags_)) {  Verbose = INTEGER(InitFlags_)[0]; 
      } else if (Rf_isReal(InitFlags_)) { 
        Verbose = (int) REAL(InitFlags_)[0]; 
      }
    }
    //Verbose = 5;
    if (Verbose > 1) {
      Rprintf("BayesSpikeCL: Welcome\n"); R_FlushConsole();
    }
    if (Rf_isNull(X_)) {
      Rprintf("BayesSpikeCL: Hey, bad result X_ is NULL!\n"); R_FlushConsole();
      return;
    }
    if (!Rf_isReal(X_)) {
      Rprintf("BayesSpikeCL: Hey, bad result X_ is not real!\n");  R_FlushConsole();
      return;
    }
    if (Rf_isNull(Y_)) {
      Rprintf("BayesSpikeCL: Hey, bad result Y_ is null!\n"); R_FlushConsole();
    }
    
    if (Verbose > 1) {
      Rprintf("BayesSpikeCL:: Checking RBayesSpikeNameSpace\n"); R_FlushConsole();
    }
    RBayesSpikeNameSpace = new AObject(_BayesSpikeNameSpace);
    BayesSpikeNameSpace = RBayesSpikeNameSpace->asSexp();
        
    if (Verbose > 1) {
      Rprintf("Starting by getting Beta \n"); R_FlushConsole();
    }
    RsBeta = NULL;  
    if (Rf_isInteger(Beta_)) {
      RsBeta = new AObject(Rf_allocVector(REALSXP, Rf_length(Beta_)));
      for (int iit = 0; iit < Rf_length(Beta_); iit++) {
        REAL(RsBeta->asSexp())[iit] = ((double) INTEGER(Beta_)[iit]);
      }
    } else {
      RsBeta = new AObject(Beta_);
    }
    sBeta = RsBeta->asSexp();
    
    if (Verbose > 1) {
      Rprintf("Starting by getting X,Y \n"); R_FlushConsole();
    }
    RsX = new AObject(X_); sX = RsX->asSexp();
    if (Rf_isInteger(Y_)) {
      RsY = new AObject(Rf_allocVector(REALSXP, Rf_length(Y_)));
      for (int iit = 0; iit < Rf_length(Y_); iit++) {
        REAL(RsY->asSexp())[iit] = ((double) INTEGER(Y_)[iit]);
      }
    } else {
      RsY = new AObject(Y_); sY = RsY->asSexp();
    }
    
    RInitFlags = new AObject(Rf_allocVector(INTSXP, Rf_length(InitFlags_)));
    InitFlags = RInitFlags->asSexp();
    SampleTausOnly = 0; DoSampleBFixedii = 1;
    for (int ii = 0; ii < Rf_length(InitFlags_); ii++) {
      if (Rf_isInteger(InitFlags_)) {
        INTEGER(InitFlags)[ii] = INTEGER(InitFlags_)[ii];
      } else if (Rf_isReal(InitFlags_)) {
        INTEGER(InitFlags)[ii] = ((int) REAL(InitFlags_)[ii]);
      } else {
      
      }
    }
    if (!Rf_isNull(Rf_getAttrib(InitFlags_, R_NamesSymbol))) {
      char *OnName;
      SEXP sStrings = R_NilValue;
      Rf_protect(sStrings = Rf_allocVector(STRSXP, Rf_length(InitFlags_)));
      for (int iti = 0; iti < Rf_length(InitFlags_); iti++) {
        OnName = (char*) CHAR(STRING_ELT(Rf_getAttrib(InitFlags_, R_NamesSymbol), iti));
        SET_STRING_ELT(sStrings, iti, Rf_mkChar(OnName));
      }
      Rf_setAttrib(InitFlags, R_NamesSymbol, sStrings);
      Rf_unprotect(1);
    }

    sBeta = RsBeta->asSexp(); 
    if (Verbose > 1) {
      Rprintf("BayesSpikeCL:: Now RandomInfoList\n"); R_FlushConsole();
    }
    RsRandomInfoList = new RRObject(RandomInfoList_);
    sRandomInfoList = RsRandomInfoList->asSexp();
    if (Rf_isNull(sX)) {
      Rf_error("sX is NULL!\n");
    }
    if (Rf_isNull(sY)) {
      Rf_error("sY is NULL!\n");
    }
    if (Rf_length(sY) <= 0) {
      Rf_error("sY has zero length!\n");
    }
    if (!R_finite(Rf_length(sY))) {
      Rf_error("sY has infinite length!\n");
    }    
    if (Rf_isNull(Rf_getAttrib(sX, R_DimSymbol))) {
      Rf_error("No Dimension to sX!");
    }
    if (Verbose > 1) {
      if (Rf_isReal(Rf_getAttrib(sX, R_DimSymbol))) {
        Rprintf("Dimension of X is, %d, %d \n",
         REAL(Rf_getAttrib(sX, R_DimSymbol))[0], 
         REAL(Rf_getAttrib(sX, R_DimSymbol))[1]
        );      
        Rprintf("Also, crazy dim of sX is real!\n"); R_FlushConsole();
      } else if (Rf_isInteger(Rf_getAttrib(sX, R_DimSymbol))) {
        Rprintf("Dimension of X is, %d, %d \n",
          INTEGER(Rf_getAttrib(sX, R_DimSymbol))[0], 
         INTEGER(Rf_getAttrib(sX, R_DimSymbol))[1]
        ); 
      }
       Rprintf("Length of Y is %d \n", Rf_length(sY));
       Rprintf("Length of sBeta is %d\n", Rf_length(sBeta));  R_FlushConsole();
       Rprintf("Length of RandomInfoList is is %d\n", Rf_length(sRandomInfoList));  R_FlushConsole();       
       Rprintf("Length of InitFlags is is %d\n", Rf_length(InitFlags_));  R_FlushConsole();  
    }


    if (Verbose > 1) {
      Rprintf("BayesSpikeCL: Assigned SEXPs\n"); R_FlushConsole();
    }
    if (!Rf_isInteger(Rf_getAttrib(X_, R_DimSymbol))) {
      Rprintf("BayesSpikeCL: Weird, do not have INTEGER attirbute to X_ dimension!\n");
      R_FlushConsole();
    }
    p = INTEGER(Rf_getAttrib(X_, R_DimSymbol))[1];
    n = Rf_length(Y_);
    if (Rf_length(sBeta) != p) {
      Rf_error("Beta Error has length %d, but p %d\n", 
        Rf_length(sBeta), p);
    }
    if (Verbose > 1) {
      Rprintf("BayesSpikeCL: p = %d, n = %d, Putting in all results\n", p,n); R_FlushConsole();
    }

     DoShortMIP = 1;  RunProbRegionVector = NULL;
     RunProbVector = NULL; StartRunProbVector = -100;  RunProbRegionVector = NULL;
     TotInProbVector = 0;  RunProbRegionVector = NULL;  RegionWidth = -5;
     TotEveryProbVector = NULL;  RsTimePartsList = NULL;
     //MedianUp = NULL;  MedianDown=NULL;
     RunMoreZero = NULL;  RunLessZero = NULL; RunSumBeta = NULL;
     RtauEndList = NULL; tauEndList = R_NilValue; RtauFixed = NULL; RtauFixedPrior = NULL; RTBSRoo = NULL; RTBSR5 = NULL;
     RsOnTau = NULL; RsOnSigma = NULL;  sOnSigma = R_NilValue; RSigmaPrior = NULL; SigmaPrior = R_NilValue;
      WillAddNewTau = 0;
     RsCenteredColumns = NULL; RsDeCenteredCodaList = NULL;
     RsOnPiA = NULL; sOnPiA = R_NilValue; RPiAPrior = NULL;  PiAPrior = R_NilValue; 
     OnTauIndex = 0;
     Rstaubarnu = NULL;  Rstaupriordf = NULL;   staubarnu = R_NilValue; staupriordf = R_NilValue;
     sPriorXScaling = R_NilValue; sPriorOftau = R_NilValue;
     RsPriorXScaling = NULL;  RsPriorOftau = NULL;   
     MyEigenVectorList = R_NilValue;  MyEigenValueList = R_NilValue;
     RAllEigenVectors = NULL;  RAllEigenValues = NULL;  AllEigenVectors=R_NilValue;
     AllEigenValues = R_NilValue;
     RCodaTable = NULL; RDoRecord = NULL;  
     ROldColNames = NULL; ROldDeCenterColNames = NULL;
     ROtherNameCodaList = NULL;  ROtherNameDeCenterCodaList = NULL;
     RCodaList = NULL; RAFD = NULL; RSubCodaList = NULL;  RSubCodaSubSetCoords=NULL;
     RSubSetTau = NULL;  RTauCodaList = NULL;
     RAllTempCodaLists = NULL;
     
     RsCodaILocFile = NULL;  RsCodaProbFile = NULL;
     RsOldCodaILocFile = NULL; RsOldCodaProbFile = NULL;
     RsCodaIFile = NULL; RsCodaJFile = NULL; 
     RsCodadTFile = NULL; RsCodaiTFile = NULL;  RsCodaiTLocFile=NULL;
     RsCodaPFile = NULL;  XtXRecorded = NULL; LengthXtXRecorded = 0;

    InitiateTime = NULL; SecondInitiateTime = NULL; CompleteTime = NULL;
    RsTimePartsList = NULL;
  
    sOnTau = R_NilValue; CodaTable = R_NilValue;
    BFixed = NULL; OrderedActive = NULL;     tauEndList = R_NilValue;
    rIChol = NULL; rICholBack = NULL;   OnKSqPBack = 0;
    rSmallXtX = NULL; rI = NULL;  MaxLargeSizeSquares = 0;
    PropBeta = NULL; iiWeight = NULL;   SumiiWeight = NULL; XtY = NULL;
    Codajj = NULL;  CodaBeta = NULL; WW = NULL;
    sCodaIFile = R_NilValue;  sCodaJFile = R_NilValue;
    OnCoda = 0; OnCodaT = 0;  XjSq = NULL;  TotOnCodaT= 0;
    ProbFixed = NULL;  ProbTau = NULL;
    iFirstRandom = -1; MT = NULL;
    staubarnu = R_NilValue; staupriordf = R_NilValue;
    sPriorXScaling = R_NilValue; sPriorOftau = R_NilValue;
    pXtX = NULL;  XtY = NULL;  NoSave = 0; DoSave = 1;
    XtResid = NULL; iiWeight = NULL;  SumiiWeight = NULL; DependenciesList = R_NilValue;
    RDependenciesFixed=NULL;  RDependenciesTau=NULL;
    DependenciesFixed = R_NilValue;  DependenciesTau = R_NilValue;
    SmallXtResid = NULL; SmallRVec = NULL;  AllNewCoords = 0; NewFixedCoords=0;
    BFixed = NULL; OrderedActive = NULL;  NumActive = 0;
    OrderAttack = NULL;  HowSample = 3; NumIntegrations = 10000; FakeOuttau = 0.0;
    IntegratedDensity = 0.0;  iWeightedXtX = NULL; 
    MaxIters = 100;  CauchyEpsilon = .00001;
    SpliceMove = 1;  NumSpliceSampleReps = 5;  DoMax = 1;
    MyEigenVectorList = R_NilValue;  MyEigenValueList = R_NilValue;
    RAllEigenVectors = NULL;  RAllEigenValues = NULL;
    AllEigenVectors = R_NilValue;  AllEigenValues = R_NilValue;
    tauFixed = R_NilValue; tauFixedPrior = R_NilValue;
    NInstance = 0; MaxGibbsIters = 100; CurrentMaxGibbsIters=MaxGibbsIters;CodaTable = R_NilValue;
    NumActiveFixed = 0;  NumActiveTau = 0; cXtY = NULL; smXtY = NULL; NumActiveBack = 0;
    JustDidPartBeta = -1;  JustDidPartBFixed = -10;
    PiAPrior = R_NilValue; sOnSigma = R_NilValue;
    sOnPiA = R_NilValue;  RSigmaPrior = NULL; SigmaPrior = R_NilValue;  RDoRecord = NULL; DoRecord = R_NilValue;
    SmallRVec = NULL; MaxTauList = 0;  gPrior = 0.0; 
    NewWrite = 1;    EarlyEndStep = -1;  EarlyEndtt = -1;   OnCodaTable = 0;
    PrintIter = 200; ReorderIter = 200; nWeight = Rf_length(sY);
    TBSRoo = R_NilValue; IamBayesSpike = 1; CodaList = R_NilValue; RSubCodaList = NULL;
    RSubCodaSubSetCoords = NULL;   RTauCodaList = NULL; RSubSetTau = NULL;
    iWeightedXtX = NULL; dfTNoise = -1; RsSaveDir =NULL; sSaveDir = R_NilValue;
    sCodadTFile = R_NilValue; sCodaiTFile = R_NilValue; sCodaPFile = R_NilValue;
    RsCodaiTLocFile = NULL; CodaTisSetup = 0;
    NewTWrite = 1;  NewPWrite = 1; CodaTjj = NULL;  CodaTau = NULL; CodaTLocjj = NULL;
    LengthCodaTjj = 0;  LengthWrittenTauCodaBuffer = 0;  LengthTotalWrittenTauCodaBuffer = 0;
    LengthCodaTLoc = 100;  DoLogitNonePostPreProb = -1;
    CodaP = NULL; tt = 0;  AFD = R_NilValue; BeingDestroyed = 0; NumFails=0; TotalFails= 0;
    MaxAlsErrors=10; GoBackAfter=10;  RPostProbCodaList = NULL; RMIP = NULL;  RegionMIP = NULL;
    LengthNoShrinkFixed = 0; LengthNoShrinkRandom = 0; sOnSigma = R_NilValue;
    NoShrinkFixed = NULL;  NoShrinkFixedPrior = NULL;
    NoShrinkRandom = NULL; NoShrinkRandomPrior = NULL;
    PctBetas = NULL;  PctTaus = NULL; StartIter = 0; EndIter = 0;
    //RCodaChainsList = NULL;  LengthRCodasList = 0;
    NumlAZero = 0;   TotalLAZero = 0;  GoBackAfterlA = 20;  MaxlAErrors = 30;
    sPriorProbFixed = R_NilValue; sPriorProbTau = R_NilValue; rPriorProbFixed = NULL; rPriorProbTau = NULL;
    HadSetPriorProbFixed = 0;
    HadSetPriorProbFixed = 0; HadSetPriorProbTau = 0;
    ThisMaxTau = MAXTAU;   TemperatureList = NULL; Temperature = 1.0;   TemperatureDecreasingRate = 1.0;
    invTemperature = 1.0;  LengthTempDraws = NULL;   NumMerges = NULL;
    TotalOpenedFiles = 0;  TotalClosedFiles = 0;  TotalR5OpenedFiles = 0;
    TotalR5ClosedFiles = 0;
    Tempii = 0; LengthTemperatureList = 0; CurrentProb = 0.0; LengthProbCoda = 0;  
    LengthTempDraws = NULL;
    LengthPostProb = 0; LengthPostProbCodaBuffer = 0; LengthWrittenPostProb = 0;
    NewPostProbWrite = 1;
    PostProbCodaBuffer = NULL;   RsPostProbBufferFile=NULL;
    ICoda = NULL; ProbCoda = NULL;  ProbOldCoda = NULL; IOldCoda = NULL;
    RsCodaOldIFile = NULL;  RsCodaOldJFile = NULL;  IterOldCoda = -1;
    intY = NULL;   dfRobit = -1.0;   EEMergeEvery = -1; burnin = 100;
    OnRNGState = 0; DoingEEProbSort = 0; EEProbSortWidth = .5;
    MaxLengthOldIDFiles = 0;  NumMerges = NULL;   fd = -1; fd1 = -1; fd2 = -1;
    RsWeightFile = NULL;  RsYFile = NULL; RsYCodaList = NULL;  RsWeightCodaList=NULL;
    WeightBuffer = NULL;  YBuffer = NULL; RsPiACodaList=NULL;  RsSigCodaList =NULL;
    LengthYBuffer = 100;  LengthWrittenYBuffer = 0; 
    LengthTotalWrittenYBuffer = 0;  LengthTotalWrittenWeightBuffer = 0;
    LengthWeightBuffer = 100;  LengthWrittenWeightBuffer = 0;   
    NewYBufferWrite = 1;  NewWeightBufferWrite = 1; 
    RobitNotReplaceFlag = 0;
    
    RsCodaBetaAllDrawBuffer = NULL;  BetaAllDrawBuffer = NULL;
    SubCodaLongList = NULL;  LengthBetaAllDrawBuffer = 0;
    LengthWrittenBetaAllDrawBuffer = 0;
    LengthTotalWrittenBetaAllDrawBuffer = 0;  rSmallXtX = NULL;

    AlterInitRun = 0;
    AlterWeightBuffer = NULL;  RsAlterWeightFile = NULL;  LengthWeightBuffer = 0;
    LengthWrittenAlterWeightBuffer = 0;  LengthTotalWrittenAlterWeightBuffer = 0;
    NewAlterWeightBufferWrite = 0;  RsAlterWeightCodaList = NULL;
    AlterWeightdfRobit = -1.0;  AlterWeightdfTNoise = -1.0;
    AlterWeightTemperature = 1.0;
    
    RsPiACodaFile=NULL;  PiACodaBuffer=NULL; LengthPiACodaBuffer=0;
    LengthWrittenPiACodaBuffer=0; NewWritePiACodaBuffer = 0;
    LengthTotalWrittenPiACodaBuffer = 0;
  
    RsSigCodaFile=NULL;  SigCodaBuffer=NULL; LengthSigCodaBuffer=0;
    LengthWrittenSigCodaBuffer=0;  NewWriteSigCodaBuffer = 0;
    LengthTotalWrittenSigCodaBuffer = 0;
    
    int ii = 0;
    if (InitFlags == NULL || Rf_isNull(InitFlags) || Rf_length(InitFlags) < 1) {
      if (InitFlags == NULL) {
        Rprintf("BayesSpikeCL:: InitFlags is NULL\n");
      } else if (Rf_isNull(InitFlags)) { 
        Rprintf("BayesSpikeCL:: InitFlags is NilValue\n");
      } else {
        Rprintf("BayesSpikeCL:: Length InitFlags is %d \n", Rf_length(InitFlags));
      }
      Rf_error("BayesSpikeCL: InitFlags has no length!");
    }
    if (Verbose > 1) {
      Rprintf("BayesSpikeCL: Setting Items from InitFlags\n", p,n); R_FlushConsole();
    }
    if (Verbose > 1) {
      Rprintf("  BayesSpikeCL:: Note that InitFlags NInstance from first");
      Rprintf(" InitFlags = %d, setting.\n", ANINT(InitFlags,0));
      R_FlushConsole();
    }
    set_NInstance(ANINT(InitFlags, 0));
    if (Rf_length(InitFlags) < 2) {
      Rf_error("BayesSpikeCL: InitFlags is not long enough for Verbose, NInstance, InitKKs\n");
    }
    if (Verbose > 2) {
      Rprintf("  BayesSpikeCL:: Note that InitFlags Verbose = %d \n", ANINT(InitFlags,1));
      R_FlushConsole();
    }
    if (Verbose > 2) {
      Rprintf("  BayesSpikeCL:: Length InitFlags = %d, Now setting Verbose \n", Rf_length(InitFlags));
      R_FlushConsole();
    }
    Verbose = ANINT( InitFlags, 1 );
    if( Rf_length(InitFlags) < 3 ) {
      Rf_error("  BayesSpikeCL:: Error, InitKKs less than length 3 \n"); R_FlushConsole();
    }
    if (Verbose > 1) {
      Rprintf("  BayesSpikeCL:: Length InitFlags = %d, Setting InitKKs \n");
      R_FlushConsole();
    }
    InitKKs = ANINT( InitFlags, 2);
    if (Verbose > 1) {
      Rprintf("BayesSpikeCL: Got Stuff from InitKKs\n", p,n); R_FlushConsole();
    }
    if (Rf_length(InitFlags) >= 4 && ANINT(InitFlags, 3) > 0 &&
      ANINT(InitFlags,3) <= p) {
      iFirstRandom = ANINT(InitFlags,3);  
    } 
    if (InitKKs <= 0) {
      InitKKs = 5;
    }
    OnKappaMem = InitKKs; OnKappaS = 0;
    InitKKs = InitKKs;
  
    int OnKSq = OnKappaMem * OnKappaMem;
    int OnKSqP = (int) ceil( (OnKappaMem)*(OnKappaMem+1) / 2.0) + 2;
    if (OnKappaMem >= MaxSquaresAllocation) {
      Rprintf("BayesSpikeCL: MaxSquaresAllocation = %d, OnKappaMem =%d is already too large!\n",
        MaxSquaresAllocation, OnKappaMem);
      R_FlushConsole();
      MaxLargeSizeSquares = MaxSquaresAllocation;
      OnKSq = MaxLargeSizeSquares*MaxLargeSizeSquares;
      OnKSqP = (int) ceil( (MaxSquaresAllocation)*(MaxSquaresAllocation+1) / 2.0) + 2;
    } else {
      MaxLargeSizeSquares = OnKappaMem;
      OnKSqP = MaxLargeSizeSquares * (MaxLargeSizeSquares+1)/2+2;
    }
    if (Verbose > 1) {
      Rprintf("BayesSpikeCL: Allocate Mem, OnKappaMem = %d, OnKSq = %d, OnKSqP = %d\n", 
        OnKappaMem, OnKSq, OnKSqP); R_FlushConsole();
    }  
    
    if (OnKSqP <= 0) {
      Rf_error("BayesSpikeCL: AllocateMem, this won't work if OnKSqP=%d!\n", OnKSqP);
      R_FlushConsole();
    }
    RMemGetD( rIChol, ("rIChol") , OnKSqP );
    if (Verbose > 1) {
      Rprintf("BayesSpikeCL: Allocate first time inside start rSmallXtX\n"); R_FlushConsole();
    }
    if (rSmallXtX != NULL) {
      RMemGetD( rSmallXtX, "rSmallXtX", OnKSqP );
    }
    if (Verbose > 1) {
      Rprintf("BayesSpikeCL: Allocate first time inside start rI\n"); R_FlushConsole();
    }
    RMemGetD( rI, "rI",  OnKSqP);
      if (Verbose > 1) {
      Rprintf("BayesSpikeCL: Allocate first time inside start PropBeta\n"); R_FlushConsole();
    }
    RMemGetD( PropBeta, "PropBeta", OnKappaMem+1);

    if (Verbose > 1) {
      Rprintf("BayesSpikeCL: Allocate first time inside start XLC\n"); R_FlushConsole();
    }    
    RMemGetI(XLC, "XLC", p);
    for (ii = 0; ii < p; ii++) { XLC[ii] = -1;}
    RMemGetI(kXFinder, "kXFinder" , OnKappaMem);
    
    RMemGetD( XtY, ("XtY"),  p);
    RMemGetD( XtResid, ("XtResid"),  p);
    if (OnKappaMem >= MaxSquaresAllocation) {
      Rprintf("Error going to happen in XtX Allocation MaxSquaresAllocation=%d, OnKappaMem = %d\n",
        MaxSquaresAllocation, OnKappaMem);   R_FlushConsole();
    }
    //RMemGetD( XtX, ("XtX"), (p * OnKappaMem) );
    
    if (pXtX != NULL) {
      Rprintf("Warning: No, pXtX is not NULL but we don't know its length to kill it!\n");
      R_FlushConsole();
    }
    pXtX = (double**) Calloc(OnKappaMem, double *);
    for (int ii = 0; ii < OnKappaMem; ii++) {
      pXtX[ii] = NULL;
    }
    //double *Block;
    if (Verbose>=3) {
      Rprintf("BayesSpike.h:: Allocation of Block of Memory \n"); R_FlushConsole();
    }
    if (p < 10000) {
      //Block = NULL;
      //Block = (double*) Calloc( p * OnKappaMem, double);
      //RMemGetD(Block, "Block", (p*OnKappaMem));
      for (int ii = 0; ii < OnKappaMem; ii++) {
        pXtX[ii] = (double *) Calloc(p, double);
      }
      //Block = NULL;
    } else {
      for (int ii = 0; ii < OnKappaMem; ii++) {
        pXtX[ii] = (double *) Calloc(p, double);
      }
    }
    RMemGetD( cXtY, ("cXtY"), OnKappaMem);
    TrueOrder = 1;

    if (Verbose > 1) {
      Rprintf("BayesSpikeCL: Update XtY\n", 
        OnKappaMem, OnKSq, OnKSqP); R_FlushConsole();
    }
    UpdateXtY();
    UpdateNonFreshXtResid();
    RMemGetD(XjSq, ("XjSq"), p);
    SetupXjSq();    


    if (Verbose > 1) {
      Rprintf("BayesSpikeCL: Check sBeta\n", 
        OnKappaMem, OnKSq, OnKSqP); R_FlushConsole();
    }   
    if (Rf_length(sBeta) != p) {
      Rf_error("BayesSpikeCL:  p = %d, but length sBeta = %d\n",
        p, Rf_length(sBeta));
    }   
    //int One = 1; double ZeroD = 0.0; char Trans = 'T'; double OneD = 1.0;   
    int One = 1; 
    double CountSBeta = 0.0;
    SEXP nsBeta = R_NilValue;  AObject *RnsBeta;
    if (Verbose >= 0 && !Rf_isNull(nsBeta)) {
      Rprintf("BayesSpikeCL: weird for nsBeta to not be null!\n"); R_FlushConsole();
    }
    if (Rf_isNull(sBeta)) {
      Rprintf("BayesSpikeCL: sBeta is Null, weird warning! \n"); R_FlushConsole();
    }
    if (Rf_length(sBeta) != p) {
      Rprintf("BayesSpikeCL: sBeta is weird, warning! length %d \n",
        Rf_length(sBeta)); R_FlushConsole();
    }
    if (Rf_isInteger(sBeta)) {
      RnsBeta = new AObject(Rf_allocVector(REALSXP,p));
      for (int itt = 0; itt < p; itt++) {
        REAL(RnsBeta->asSexp())[itt] = INTEGER(sBeta)[itt];
      }
      if (RsBeta != NULL) {
        DDelete(RsBeta, "RsBeta");
      }
      RsBeta = RnsBeta;  sBeta = RsBeta->asSexp();
    }   
    CountSBeta = F77_CALL(dasum)(&p, REAL(sBeta), &One);
    
    if (Verbose > 2 && CountSBeta > 0.0) {
      Rprintf("BayesSpikeCL: Beta is Not Zero, going to calc XtResid\n"); 
      R_FlushConsole();
    }
    if (CountSBeta > 0.0) {
      SetXtResidNoXtX();
    }
    if (Verbose > 2 && CountSBeta > 0.0) {
      Rprintf("BayesSpikeCL: Beta was Not Zero, calculated XtResid\n");  
      R_FlushConsole();
    }

      //F77_CALL(dgemv)((char*) &Trans, &p, &LenCopy, &NegOneD,
  //  REAL(XtX), &p, REAL(PropBeta), &One,
  //  &OneD, REAL(XTResid), &One);
    if (Verbose > 2 && iFirstRandom < 0) {
      Rprintf("BayesSpikeCL: Setting UpBFixedAndFriends\n"); 
      R_FlushConsole();
    }
    if (iFirstRandom > 0 && iFirstRandom <= p) {
      SetupBFixedAndFriends();
    } 
    if (Verbose > 2) {
      Rprintf("BayesSpikeCL:: BFixedAndFriends Are Setup, now YSq\n");
       R_FlushConsole();
    }
    
    SetupYSq(); 
    if (Verbose > 2) {
      Rprintf("BayesSpikeCL:: BFixedAndFriends Are Setup, now YResidSq\n");
       R_FlushConsole();
    } 
    //SetupYResidSq();
    if (Verbose > 2 && iFirstRandom < 0) {
      Rprintf("BayesSpikeCL: Returning Home\n"); 
      R_FlushConsole();
    }  
    
    NumSteps = 100; tt = 0; TypePrior = 1; NumberNew = 0; MaxTauEst = 0;
    MaximizeMeCauchyTotalIters = 0;
  }
  BayesSpikeCL() {
    Rprintf("BayesSpikeCL: You Initialized me Blank!"); R_FlushConsole();
    AmDestroying = 0;      RsY = NULL; sY = R_NilValue; RsX = NULL; sX = R_NilValue;
         RtauEndList = NULL; tauEndList = R_NilValue; RtauFixed = NULL; RtauFixedPrior = NULL;
     RsOnTau = NULL; RsOnSigma = NULL;  sOnSigma = R_NilValue; RSigmaPrior = NULL; SigmaPrior = R_NilValue; 
     RsCenteredColumns = NULL;  RsDeCenteredCodaList = NULL; 
     WillAddNewTau = 0; AllNewCoords = 0;  NewFixedCoords = 0;
     RsOnPiA = NULL; RPiAPrior = NULL;  KillPct = DefaultKillPct;
     MaximumAllocation = DefaultMaximumAllocation;  gPrior = 0.0;
     StartRemoval = DefaultStartRemoval; DoAddCoordsOnSetBeta = 0;
     CheckRemovable =DefaultCheckRemovable; iiWeight = NULL; iWeightedXtX = NULL;
     ItersSinceRemoval = 0;  MaxCODABeta= DefaultMaxCODA;
     MinimumProbKeeper = DefaultMinimumProbKeeper;
     Rstaubarnu = NULL;  Rstaupriordf = NULL;
     RsPriorXScaling = NULL;  RsPriorOftau = NULL;   
     RAllEigenVectors = NULL;  RAllEigenValues = NULL;  AllEigenVectors = R_NilValue;
     AllEigenValues = R_NilValue;  RsTimePartsList = NULL;
     MyEigenVectorList = R_NilValue;  MyEigenValueList = R_NilValue;
     RCodaTable = NULL; RDoRecord = NULL; DoRecord = R_NilValue;
     ROldColNames = NULL; ROldDeCenterColNames=NULL;
     ROtherNameCodaList = NULL;  ROtherNameDeCenterCodaList=NULL;
     CodaTau = NULL; CodaTisSetup = 0;
     LengthWrittenTauCodaBuffer = 0;  LengthTotalWrittenTauCodaBuffer = 0;
     RCodaList = NULL; RAFD = NULL;  ScaleOfEigenVectors = 1.0;
     RAllTempCodaLists = NULL;  RunProbRegionVector = NULL;  RegionWidth = -1;
     RsCodaIFile = NULL; RsCodaJFile = NULL;  RsCodaOldIFile = NULL;  
     RsCodaOldJFile = NULL;  RsTimePartsList = NULL;
     RsCodadTFile = NULL; RsCodaiTFile = NULL;  RsCodaiTLocFile=NULL;
     RsCodaPFile = NULL;  XjSq = NULL; RsCodaILocFile = NULL; RsCodaProbFile = NULL;
     RsOldCodaILocFile = NULL;  RsOldCodaProbFile = NULL;
     RsSaveDir = NULL;  sSaveDir = R_NilValue;
     SampleTausOnly = 0;  DoSampleBFixedii = 1;
     LastStF=0, LastGoFor = 0;  AlterWeightBuffer = NULL;
     MaxInvert = DefaultMaxInvert;
     RInitFlags = NULL;  InitFlags = R_NilValue;
    BFixed = NULL; OrderedActive = NULL; NoSave = 0;  DoSave = 1;
    sX = NULL;  sY = NULL; sBeta = NULL; sRandomInfoList = NULL;
    InitFlags = NULL; p = 0; n = 0;
    rIChol = NULL; rICholBack = NULL; rSmallXtX = NULL; rI = NULL; OnKSqPBack = 0;
    PropBeta = NULL; pXtX = NULL;  XtY = NULL; iiWeight = NULL; SumiiWeight = NULL;
    XtResid = NULL; CodaBeta = NULL; Codajj = NULL; WW = NULL;
    sCodaIFile = R_NilValue;  sCodaJFile = R_NilValue;  
    sCodadTFile = R_NilValue; sCodaiTFile = R_NilValue; sCodaPFile = R_NilValue;
    OnCoda = 0; OnCodaT = 0; TotOnCodaT = 0;  sOnSigma = NULL; RTBSRoo = NULL; RTBSR5 = NULL;
    OrderedActive = NULL;  NumActive = 0; OnKappaMem = 0; OnKappaS = 0; 
    ProbFixed = NULL;  ProbTau = NULL;  iFirstRandom = -1;
    staubarnu = R_NilValue; staupriordf = R_NilValue;
    sPriorXScaling = R_NilValue; sPriorOftau = R_NilValue;   MT = NULL;
    DependenciesList = R_NilValue;  tauFixed = R_NilValue;  tauFixedPrior = R_NilValue;
    OrderAttack = NULL;  XjSq = NULL;   XtXRecorded = NULL;   LengthXtXRecorded = 0;
    YSq = 0; YResidSq = 0;  tauEndList = R_NilValue; DontPartBetaNoise = 0;
    SmallXtResid = NULL; SmallRVec = NULL; HowSample = 3;  NumIntegrations = 10000;
    NumSteps = 100; tt = 0; TypePrior = 0; NumberNew = 0;
    MaximizeMeCauchyTotalIters = 100; NInstance = 0;
    RDependenciesFixed=NULL;  RDependenciesTau=NULL;
    DependenciesFixed = R_NilValue;  DependenciesTau = R_NilValue;
    MaxGibbsIters = 0;  CurrentMaxGibbsIters = MaxGibbsIters; CodaTable = R_NilValue;
    NumActiveFixed = 0; NumActiveTau = 0; cXtY = NULL; smXtY = NULL; NumActiveBack = 0;
    JustDidPartBeta = -1;   JustDidPartBFixed = -10;
    NewWrite = 1; TBSRoo = R_NilValue; IamBayesSpike = 1; CodaList = R_NilValue;
    iWeightedXtX = NULL; dfTNoise = 0;  BayesSpikeNameSpace = R_NilValue;  tt = 0;
    CodaTjj = NULL; CodaTLocjj=NULL; CodaTau = NULL; CodaP = NULL;  AFD = R_NilValue;
    CodaTisSetup = 0;    DoLogitNonePostPreProb = -1;
    LengthCodaTLoc = 100;   LengthWrittenTauCodaBuffer = 0;
    LengthTotalWrittenTauCodaBuffer = 0;
    NumlAZero = 0;   TotalLAZero = 0;  GoBackAfterlA = 20;  MaxlAErrors = 30;
    BeingDestroyed = 0;  NumFails= 0; TotalFails = 0; GoBackAfter = 10;  MaxAlsErrors = 10;
    LengthNoShrinkFixed = 0; LengthNoShrinkRandom = 0; PctBetas = NULL; PctTaus = NULL;
   NoShrinkFixed = NULL; NoShrinkFixedPrior = NULL;
   NoShrinkRandom = NULL; NoShrinkRandomPrior = NULL; StartIter = 0; EndIter = 0;
    RsPiACodaFile=NULL;  PiACodaBuffer=NULL; LengthPiACodaBuffer=0;
    LengthWrittenPiACodaBuffer=0;    LengthTotalWrittenPiACodaBuffer = 0;
   InitiateTime = NULL; SecondInitiateTime = NULL; CompleteTime = NULL;
   ZeroOutBeforeSample = 0;  RevertTemperatureEvery = -1;
   
   LengthMergeActive = 0; MaxLengthMergeActive = 0;  MergeActive = NULL; MergePropBeta = NULL;
   
   RsCodaBetaAllDrawBuffer = NULL;  BetaAllDrawBuffer = NULL;
   SubCodaLongList = NULL;
   LengthWrittenBetaAllDrawBuffer = 0;
   LengthTotalWrittenBetaAllDrawBuffer = 0;
  
    RsSigCodaFile=NULL;  SigCodaBuffer=NULL; LengthSigCodaBuffer=0;
    LengthWrittenSigCodaBuffer=0;   LengthTotalWrittenSigCodaBuffer = 0;
   rPriorProbFixed = NULL; rPriorProbTau = NULL; HadSetPriorProbFixed = 0;
   sPriorProbFixed=R_NilValue;  sPriorProbTau = R_NilValue;
   HadSetPriorProbTau = 0;
   TotalOpenedFiles = 0; TotalClosedFiles = 0;  
   TotalR5OpenedFiles = 0; TotalR5ClosedFiles = 0;fd = -1;  fd1 = -1; fd = -2;
   NumMerges = NULL;
   OnRNGState = 0;
   IterOldCoda = -1;
   burnin = 200;
   ThisMaxTau = MAXTAU;
   intY = NULL;
   TemperatureList = NULL; Temperature = 1.0;  invTemperature = 1.0;
   OldTemperature = 1.0;  invOldTemperature = 1.0;  NumMerges=NULL;
    Tempii = 0; LengthTemperatureList = 0;
    CurrentProb = 0.0; LengthProbCoda = 0;  LengthTempDraws = NULL; 
    LengthPostProb = 0; LengthPostProbCodaBuffer = 0; LengthWrittenPostProb = 0;
    NewPostProbWrite = 1;
    PostProbCodaBuffer = NULL;   RsPostProbBufferFile=NULL;
   OnCodaTable = 0;  RobitNotReplaceFlag = 0;
   ICoda = NULL; ProbCoda = NULL;  ProbOldCoda = NULL;  IOldCoda = NULL;
   dfRobit = -1;   EEMergeEvery = -1;  DoingEEProbSort = 0; EEProbSortWidth = .5;
   MaxLengthOldIDFiles = 0;
   
   YBuffer = NULL;  RsYFile = NULL;  WeightBuffer = NULL;  RsWeightFile = NULL;
   RsYCodaList = NULL;  RsWeightCodaList=NULL; RsSigCodaList=NULL; RsPiACodaList=NULL;
   LengthYBuffer = 100;  LengthWrittenYBuffer = 0;  LengthWeightBuffer = 100;
   LengthTotalWrittenYBuffer = 0;  LengthTotalWrittenWeightBuffer = 0;
   LengthWrittenWeightBuffer = 0; NewWeightBufferWrite = 1;  NewYBufferWrite = 1;
    AlterInitRun = 0;
    AlterWeightBuffer = NULL;  RsAlterWeightFile = NULL;  LengthWeightBuffer = 0;
    LengthWrittenAlterWeightBuffer = 0;  LengthTotalWrittenAlterWeightBuffer = 0;
    NewAlterWeightBufferWrite = 0;  RsAlterWeightCodaList = NULL;
    AlterWeightdfRobit = -1.0;  AlterWeightdfTNoise = -1.0;
    AlterWeightTemperature = 1.0;
      
    Rprintf("You Know that I'm completely dead in the water, right?"); R_FlushConsole();
  }
  
  int get_NumActive() {return(NumActive);}
  int get_NumActiveFixed() {return(NumActiveFixed);}
  int get_NumActiveTau() {return(NumActiveTau);}
  int get_p() {return(p);}  int get_n() {return(n);}
  void RenameCodaList(SEXP dimames); 
  
  int get_MaxLengthOldIDFiles() { return(MaxLengthOldIDFiles); }

  void set_DoRecord(SEXP DoRecord);  int CountRecord();
  SEXP get_DoRecord(){return(DoRecord);}
  double SampleANewTau(int iOnii);
  SEXP get_MT();
  int AlsSwitchingSampler(double(*lFOfX)(double, void*), int iOnii);
  int SetupProbFixed() {
    int ii;
    if (Verbose >= 4) {
      Rprintf("SetupProbFixed, Try to setup with p=%d, iFirstRandom = %d. \n",
        p, iFirstRandom);  R_FlushConsole();
    }
    if (sOnTau == NULL || Rf_isNull(sOnTau) || Rf_length(sOnTau) <= 0 ||
      tauEndList == NULL || Rf_isNull(tauEndList) || 
      Rf_length(tauEndList) <= 0) {
      FFree(ProbFixed, "ProbFixed");
      RMemGetD(ProbFixed, "ProbFixed", p);
      if (Verbose >= 2) {
        Rprintf("Setup ProbFixed to length p = %d\n!", p); R_FlushConsole();
      }
      if (Rf_isReal(sOnPiA)) {
        for (ii = 0; ii < p; ii++) {
          ProbFixed[ii] = REAL(sOnPiA)[0];
        }
      } else if (Rf_isInteger(sOnPiA)) {
        for (ii = 0; ii < p; ii++) {
          ProbFixed[ii] = ((double) INTEGER(sOnPiA)[0]);
        }
      }
      return(1);
    } else  if (iFirstRandom < 0) {
      Rprintf("SetupProbFixed: Please Set FirstRandom First\n");
      R_FlushConsole();
      return(-1);
    }
    if (Rf_isNull(sOnPiA)|| Rf_length(sOnPiA) < 0) {
      Rf_error("SetupProbFixed: Please Setup OnPiA First\n");
    }
    if (Verbose >= 4) {
      Rprintf("SetupProbFixed: Try to allocate for sOnTau != NULL, %d \n", 
        iFirstRandom); R_FlushConsole();
    }
    if (iFirstRandom <= 0) {
      Rprintf("SetupProbFixed: it doesn't even make any ");
      Rprintf(" sense, iFirstRandom = %d!\n", iFirstRandom);  R_FlushConsole();
      return(-1);
    }
    FFree(ProbFixed, "ProbFixed");
    RMemGetD(ProbFixed, "ProbFixed", iFirstRandom);
    if (Rf_isReal(sOnPiA)) {
    for (ii = 0; ii < iFirstRandom; ii++)  {
      ProbFixed[ii] = REAL(sOnPiA)[0];
    }
    } else if (Rf_isInteger(sOnPiA)) {
      for (ii = 0; ii < iFirstRandom; ii++)  {
        ProbFixed[ii] = ((double) INTEGER(sOnPiA)[0]);
      }    
    }
    if (Verbose >= 4) {
      Rprintf("SetupPRobFixed: Successful Setup of ProbFixed of length %d\n",
        iFirstRandom); R_FlushConsole();
    }
    return(1);
  }
  void set_MaxGibbsIters(int MaxGibbsIters_) {
    if (MaxGibbsIters_ < 0) {
      Rf_error("MaxGibbsIters cannot be set less than Zero!");
    }
    MaxGibbsIters = MaxGibbsIters_;  CurrentMaxGibbsIters = MaxGibbsIters;
    if (EndIter <= 0) { EndIter = MaxGibbsIters; }
  }  
  int get_MaxGibbsIters(){return(MaxGibbsIters);}
  void set_tt(SEXP tt_) {
    if (Rf_isNull(tt_) || Rf_length(tt_) <= 0) {
      Rf_error("Error tt_ cannot be set to that value. \n");
    }
    if (Rf_isInteger(tt_)) {
     tt = INTEGER(tt_)[0];
    } else if (Rf_isReal(tt_)) {
      tt = REAL(tt_)[0];
    }
    if (tt <= -2) { tt = 0; return; }
    if (tt == -1) {
      tt = -1;
      Rprintf("set_tt you have set to -1 to get debug mode. \n"); R_FlushConsole();
    }
    if (tt >= MaxGibbsIters) { tt = MaxGibbsIters-1; return; }
  }
  SEXP get_tt() { SEXP sOut = Rf_allocVector(INTSXP,1);
    if (tt == -1) {
      INTEGER(sOut)[0] = -1;
      Rprintf("get_tt: MBS is in debug mode. \n");
      return(sOut);
    }
    INTEGER(sOut)[0] = tt+1;
    return(sOut); 
  }
  SEXP get_DoAddCoordsOnSetBeta() {
    SEXP sOut = Rf_allocVector(INTSXP,1);
    INTEGER(sOut)[0] = DoAddCoordsOnSetBeta;
    return(sOut); 
  }
  void set_DoAddCoordsOnSetBeta(SEXP DoAddCoordsOnSetBeta_) {
    if (Rf_isNull(DoAddCoordsOnSetBeta_) || Rf_length(DoAddCoordsOnSetBeta_) <= 0) {
      Rprintf("set_DoAddCoords: Invalid input. \n"); R_FlushConsole();
      return;
    }
    int In = -666;
    if (Rf_isInteger(DoAddCoordsOnSetBeta_) && Rf_length(DoAddCoordsOnSetBeta_) >= 1) {
       In = INTEGER(DoAddCoordsOnSetBeta_)[0];
    } else if (Rf_isReal(DoAddCoordsOnSetBeta_) && Rf_length(DoAddCoordsOnSetBeta_) >= 1) {
      In = REAL(DoAddCoordsOnSetBeta_)[0];
    }
    if (In == 0 || In == 1) {
      DoAddCoordsOnSetBeta = In;
    } else {
      Rprintf("set_DoAddCoordsOnSetBeta: please either supply 1 or 0"); R_FlushConsole();
    }
  }
  void set_AllTempCodaLists(SEXP AllTempCodaLists_) {
    DDelete(RAllTempCodaLists, "AllTempCodaLists");
    RAllTempCodaLists = new AObject(AllTempCodaLists_);
  }
  SEXP get_AllTempCodaLists() {
    if (RAllTempCodaLists == NULL) { return(R_NilValue); }
    return(RAllTempCodaLists->asSexp());
  }
  int CheckResizeIntegrity();
  void set_OnCodaTable(int OnT) {
    if (Rf_isNull(CodaList)) {
      if (CodaList == NULL || Rf_isNull(CodaList) || Rf_length(CodaList) <= 0) {
         int TI = 0;
         for (int ii = 0; ii < Rf_length(DoRecord); ii++) {
           if (Rf_isInteger(DoRecord)) { TI += INTEGER(DoRecord)[ii]; }
           if (Rf_isReal(DoRecord)) { TI += (int) REAL(DoRecord)[ii]; }
         }
         if (TI >= 1) {
           Rprintf("ERROR ERROR ERROR ERROR ERROR ERROR \n"); R_FlushConsole();
           Rprintf("BayesSpike: Really Weird!!!! \n");
           Rprintf("DoRecord = "); 
           if (Rf_isReal(DoRecord)) {
             PrintVector(REAL(DoRecord), Rf_length(DoRecord));
           } else if (Rf_isInteger(DoRecord)) { 
             PrintVector(INTEGER(DoRecord), Rf_length(DoRecord));
           }
           Rprintf("But CodaList is  definitely NULL!\n"); R_FlushConsole();
         }
      }
      Rf_error("BayesSpike: Cannot set OnCodaTable with Null OnCodaList\n");
    }
    if (Rf_length(CodaList) < OnT) {
      Rf_error("BayesSpike:  Error: CodaList not long enough for this info \n");
    }
    if (OnT <= 0) {
      Rf_error("BayesSpike: Set OnCodaTable  anywhere from 1 to CodaList length");
    }
    DDelete(RCodaTable, "RCodaTable");
    //DDelete(ROldColNames, "ROldColNames");
    OnCodaTable = OnT -1;
    RCodaTable = new AObject(VECTOR_ELT(CodaList, OnCodaTable));
    CodaTable = RCodaTable->asSexp();
  }
  void set_OtherNameCodaList(SEXP setMe) {
    if (Rf_isNull(setMe)) {
      DDelete(ROtherNameCodaList, "ROtherNameCodaList");
      return;
    }
    DDelete(ROtherNameCodaList, "ROtherNameCodaList");
    ROtherNameCodaList = new AObject(setMe);
    return;
  }
  int get_TooMuchBuffer() {
    int IWant = TooMuchBuffer;
    return(IWant);
  }
  SEXP get_OtherNameDeCenterCodaList() {
    if (ROtherNameDeCenterCodaList == NULL) { return(R_NilValue); }
    return(ROtherNameDeCenterCodaList->asSexp());
  }
    void set_OtherNameDeCenterCodaList(SEXP setMe) {
    if (Rf_isNull(setMe)) {
      DDelete(ROtherNameDeCenterCodaList, "ROtherNameDeCenterCodaList");
      return;
    }
    DDelete(ROtherNameDeCenterCodaList, "ROtherNameDeCenterCodaList");
    ROtherNameDeCenterCodaList = new AObject(setMe);
    return;
  }
  SEXP get_OtherNameCodaList() {
    if (ROtherNameCodaList == NULL) { return(R_NilValue); }
    return(ROtherNameCodaList->asSexp());
  }
  void set_OldCodaNames(SEXP setMe) {
    if (Rf_isNull(setMe)) {
      DDelete(ROldColNames, "ROldColNames");
      return;
    }
    DDelete(ROldColNames, "ROldColNames");
    ROldColNames = new AObject(setMe);
    return;
  }
  void set_OldDeCenterCodaNames(SEXP setMe) {
    if (Rf_isNull(setMe)) {
      DDelete(ROldDeCenterColNames, "ROldDeCenterColNames");
      return;
    }
    DDelete(ROldDeCenterColNames, "ROldDeCenterColNames");
    ROldDeCenterColNames = new AObject(setMe);
    return;
  }
  SEXP get_InitiateTime() {
    if (InitiateTime == NULL) { return(R_NilValue); }
    return(InitiateTime->asSexp());
  }
  void set_InitiateTime(SEXP iInitiateTime) {
    InitiateTime = new AObject(iInitiateTime);
  }
  SEXP get_SecondInitiateTime() {
    if (SecondInitiateTime == NULL) { return(R_NilValue); }
    return(SecondInitiateTime->asSexp());
  }
  void set_SecondInitiateTime(SEXP iInitiateTime) {
    SecondInitiateTime = new AObject(iInitiateTime);
  }
  SEXP get_CompleteTime() {
    if (CompleteTime == NULL) { return(R_NilValue); }
    return(CompleteTime->asSexp());
  }
  void set_CompleteTime(SEXP iCompleteTime) {
    CompleteTime = new AObject(iCompleteTime);
  }
   
  void set_SampleTausOnly(SEXP iInsert) {
    if (Rf_isNull(iInsert) || Rf_length(iInsert) <= 0) {
    } else if (Rf_isInteger(iInsert) && Rf_length(iInsert) == 1) {
      SampleTausOnly = INTEGER(iInsert)[0];
    } else if (Rf_isReal(iInsert) && Rf_length(iInsert) == 1) {
      SampleTausOnly = REAL(iInsert)[0];
    } else {
       Rprintf("set_SampleTausOnly: invalid SEXP given. \n"); R_FlushConsole();
    }
  }
  SEXP get_SampleTausOnly() {
    SEXP sOut = Rf_allocVector(INTSXP, 1);
    INTEGER(sOut)[0] = SampleTausOnly; return(sOut);
  }
  
  SEXP get_OldCodaNames() {
    if (ROldColNames == NULL) { return(R_NilValue); }
    return(ROldColNames->asSexp());
  }
  SEXP get_OldDeCenterCodaNames() {
    if (ROldDeCenterColNames == NULL) { return(R_NilValue); }
    return(ROldDeCenterColNames->asSexp());
  }
  int get_OnCodaTable() { return( OnCodaTable +1 );  }
  void set_CodaTable(SEXP CodaTable_) {
    if (CodaTable_ == NULL || Rf_isNull(CodaTable_)) {
      Rf_error("set_CodaTable: Error: cannot set Null Table!");
    }
    if (DoRecord == NULL || Rf_isNull(DoRecord) || 
     Rf_length(DoRecord) != 7) {
      Rf_error("set_CodaTable:  Error, DoRecord must be set first!"); 
    }                                            
    if (CodaList == NULL || Rf_isNull(CodaList)) {
      Rf_error("Error, cannot refresh CodaTable with Null CodaList \n");
    }
    if (RCodaTable == NULL || Rf_isNull(CodaTable)) {
      Rf_error("Error, we never set up CodaTable!");
    }
    SEXP MyODim = R_NilValue;
    MyODim = Rf_getAttrib(CodaTable, R_DimSymbol);
    if (Rf_length(MyODim) != 2) {
      Rf_error("set_CodaTable, CurrentCodaTable does not have 2 dimensions!");
    }
    if (!Rf_isInteger(MyODim)) {
      Rf_error("set_CodaTable: Why is MyODim not an integer!\n");
    }
    SEXP MyTDim = Rf_getAttrib(CodaTable_, R_DimSymbol);


    if (Rf_isNull(MyTDim)) {
      Rf_error("set_CodaTable: InsertTable must have dimensions!");
    }
    if (Rf_length(MyTDim)!= 2) {
      Rf_error("set_CodaTable: Dimension of Insert CodaTable is not 2");
    }
    if (!Rf_isInteger(MyTDim)) {
      Rf_error("set_CodaTable: Why is MyTDim not an integer!\n"); R_FlushConsole();
    }
    if (Verbose > 2) {
      Rprintf("set_CodaTable: Doing CountRecord()\n"); R_FlushConsole();
    }
    int CR = CountRecord();
    if (INTEGER(MyTDim)[1] != CR) { //Count Record is in BayesSpikeGibbs
      Rprintf("set_CodaTable: InsertTable has bad dims \n");
      Rprintf("Length(sBeta) = %d, Length(sOntau) = %d \n",
       Rf_length(sBeta), Rf_length(sOnTau)); R_FlushConsole();
      Rprintf("Length(sOnPiA) = %d, Length(sOnSigma) = %d \n",
        Rf_length(sOnPiA), Rf_length(sOnSigma)); R_FlushConsole();
      Rprintf("Length(tauFixed) = %d, Length(tauFixedPrior) = %d, iFirstRandom = %d \n",
        Rf_length(tauFixed), Rf_length(tauFixedPrior), iFirstRandom); 
      Rf_error("set_CodaTable: InsertTable has dim, %d, %d, but Do Record insist on = %d\n",
        INTEGER(MyTDim)[0], INTEGER(MyTDim)[1], CountRecord());
    }
    int RecordLength = CountRecord();
    int Mini = INTEGER(MyODim)[0];  if (INTEGER(MyTDim)[0] < Mini) {
      Mini = INTEGER(MyTDim)[0];
    }
    if (!Rf_isReal(CodaTable_)) {
      Rprintf("set_CodaTable: Error, CodaTable_ is not real!\n");
      R_FlushConsole();
    }
    for (int ii = 0; ii < Mini; ii++) {
      F77_CALL(dcopy)(&RecordLength, REAL(CodaTable_), INTEGER(MyTDim),
        REAL(CodaTable), INTEGER(MyODim));
    } 
  }
  SEXP get_CodaTable() { return(CodaTable);}
  void SetupProbTau() {
    if (Rf_isNull(tauEndList) || Rf_length(tauEndList) < 0) {
      Rf_error("SetupProbTau: Please Set tauEndList First\n");
    }
    if (Rf_isNull(sOnPiA) || Rf_length(sOnPiA) < 0) {
      Rf_error("SetupProbTau: Please Setup OnPiA First\n");
    }
    FFree(ProbTau, "ProbTau");
    RMemGetD(ProbTau, "ProbTau", Rf_length(tauEndList));
    double PP;
    if (Rf_isReal(sOnPiA)) {
      PP = REAL(sOnPiA)[0];
      if (Rf_length(sOnPiA) >= 2) {
        PP = REAL(sOnPiA)[1];
      }
    } else {
      PP = .5;
    }
    int ii = 0; for (ii = 0; ii < Rf_length(tauEndList); ii++)  {
      ProbTau[ii] = PP;
    }
  }
  double SetupYSq() {
    //int One;  
    if (Verbose > 2) {
      Rprintf("SetupYSq: Starting, n = %d\n", n); R_FlushConsole();
      if (Rf_isNull(sY)) {
        Rprintf("Hey: sY is NULL on SetupYSq!\n"); R_FlushConsole();
      }
      if (Rf_isReal(sY)) {
        Rprintf("Y = "); PrintVector(REAL(sY), n); Rprintf("\n"); R_FlushConsole();
      } else {
        Rprintf("Hey, Y is not real!! on SetupYSq!"); R_FlushConsole();
      }
    }
    if (Rf_length(sY) != n) {
      Rf_error("SetupYSq: huge error, we've a bad length %d for n = %d\n", 
        Rf_length(sY), n);
    }
    int ii;  YSq = 0;     //int C = 0;
    if (iiWeight == NULL) {
      for (ii = 0; ii < n; ii++) {
        //Rprintf(" %d ", ii);  C++; if (C >= 10) {Rprintf("\n"); C = 0;} R_FlushConsole();
        YSq += REAL(sY)[ii] * REAL(sY)[ii];
      }
    } else {
      for (ii = 0; ii < n; ii++) {
        //Rprintf(" %d ", ii);  C++; if (C >= 10) {Rprintf("\n"); C = 0;} R_FlushConsole();
        YSq += REAL(sY)[ii] * REAL(sY)[ii] * iiWeight[ii];
      }    
    }
    //Rprintf(" Okay we got to the end, %f, now to try it with a ddot, \n", YSq); R_FlushConsole(); 
    //YSq = F77_CALL(ddot)(&n, REAL(sY), &One, REAL(sY), &One);
    if (Verbose > 2) {
      Rprintf("SetupYSq: Finishing\n"); R_FlushConsole();
    }
    return(YSq);
  }
  double get_TemperatureDecreasingRate() {
    return(TemperatureDecreasingRate);
  }
  void set_TemperatureDecreasingRate(double TemperatureDecreasingRate_) {
    if (TemperatureDecreasingRate_ <= 0) {
      Rf_error("set_TemperatureDecreasingRate: can't do that for %f < 0!\n",
        TemperatureDecreasingRate_);
    }
    if (TemperatureDecreasingRate_ > 1.0) {
      Rf_error("set_TemperatureDecrasingRate: should't do that for %f > 1.0!\n",
        TemperatureDecreasingRate_);
    }
    TemperatureDecreasingRate = TemperatureDecreasingRate_;
  }
  int KillNanBeta() {
    for (int jj = 0; jj < p; jj++) {
      if (R_isnancpp(REAL(sBeta)[jj])) {
        REAL(sBeta)[jj] = 0.0;
      }
    }
    UpdateFreshXtResid();
    if (AllNewCoords > 0) {
      AddAllNewCoords();
    }
    RefreshOrderedActive(1);
    return(1);
  }
  int get_CurrentMaxGibbsIters() {
    return(CurrentMaxGibbsIters);
  }
  void set_CurrentMaxGibbsIters(int aCMGI) {
    if (aCMGI  > 0) {
      CurrentMaxGibbsIters = aCMGI;
    }
  }
  int get_NoSave() { return(NoSave); }
  void set_NoSave(int NoSave_) {
    if (NoSave_ == 1) {
      NoSave = 1;
    } else {
      NoSave = 0;
    }
  }
  int get_DoSave() { return(DoSave); }
  void set_DoSave(int DoSave_) {
    if (DoSave_ == 1) {
      DoSave = 1;
    } else {
      DoSave = 0;
    }
  }
  ////////////////////////////////////////////////////////////
  //  WW is behind the scenes junk weight vector of length n
  SEXP get_WW() {
    if (WW == NULL) { return(R_NilValue); }
    SEXP sOut;
    Rf_protect(sOut = Rf_allocVector(REALSXP, n));
    int One = 1;
    F77_CALL(dcopy)(&n, WW, &One, REAL(sOut), &One); 
    Rf_unprotect(1); return(sOut);
  }
  
  int ReweightEigenvalues();
  SEXP get_DoShortMIP() {
    SEXP sOut = R_NilValue;
    Rf_protect(sOut = Rf_allocVector(INTSXP, 1));
    INTEGER(sOut)[0] = DoShortMIP;
    Rf_unprotect(1);
    return(sOut);
  }
  void set_DoShortMIP(SEXP iIn) {
    int AIn = -1;
    if (Rf_isNull(iIn) || (!Rf_isReal(iIn) && !Rf_isInteger(iIn)) || Rf_length(iIn) <= 0) {
      Rf_error("set_DoShortMIP: Invalid input");
    }
    AIn = GetFirstInteger(iIn);
    if (AIn == 0) { DoShortMIP = 0; }
    if (AIn == 1) { DoShortMIP = 1; }
  }
  int get_StartRunProbVector() {return(StartRunProbVector);}
  int FillRunProbVector(SEXP AV, SEXP BV) {
    if (RunProbVector == NULL) {
      Rf_error("FillRunProbVector: We do not run this without setting ProbVector up!\n");
    }
    if (RunProbVector == NULL) {
      Rf_error("FillRunProbVector: We do not run this without setting TotEveryProbVector up!\n");
    }    
    if (Rf_length(AV) != Rf_length(BV)) {
      Rf_error("FillRunProbVector: can't run unless AV, BV same length!\n");
    }
    int lAV = Rf_length(AV);
    int dLen = 0;
    if (Rf_isNull(sOnTau) || Rf_length(sOnTau) <= 0) {
      dLen = p;
    } else {
      dLen = iFirstRandom  + Rf_length(sOnTau);
    }
    if (Rf_length(AV) !=  dLen) {
      Rf_error("Error, we need AV[len=%d] to have length %d \n", lAV, dLen);
    }
    if (Rf_isInteger(BV)) {
      for (int ii = 0; ii < dLen;ii++) {
        TotEveryProbVector[ii] = INTEGER(BV)[ii]; 
      }
      TotInProbVector = INTEGER(BV)[0];
    } else if (Rf_isReal(BV)) {
      for (int ii = 0; ii < dLen;ii++) {
        TotEveryProbVector[ii] = (int) REAL(BV)[ii]; 
      }
      TotInProbVector = (int) REAL(BV)[0];    
    }
    if (Rf_isInteger(AV)) {
     for (int ii = 0; ii < dLen; ii++) {
       RunProbVector[ii] = (double) INTEGER(AV)[ii];
     }
    } else if (Rf_isReal(AV)) {
      for (int ii = 0; ii < dLen; ii++) {
        RunProbVector[ii] = (double) REAL(AV)[ii];
      }
    }
    return(1);
  }
  void SetupRunProbVector(int StartRunProbVector_) {
    if (StartRunProbVector_ < -1) {
      FFree(RunProbVector, "RunProbVector");  TotInProbVector = 0;  StartRunProbVector = -1;
      FFree(TotEveryProbVector, "TotEveryProbVector");
      return;
    }
    FFree(RunProbVector, "RunProbVector");  FFree(TotEveryProbVector, "TotEveryProbVector");
    if (RunLessZero != NULL) { 
      Free(RunLessZero); RunLessZero = NULL; 
    }
    if (RunMoreZero != NULL) {
      Free(RunMoreZero);  RunMoreZero = NULL;
    }
    //FFree(MedianUp, "MedianUp");  FFree(MedianDown, "MedianDown");
    if (Rf_isNull(sOnTau) || Rf_length(sOnTau) <= 0) {
      RMemGetD(RunProbVector, "RunProbVector", p); 
      RMemGetI(TotEveryProbVector, "TotEveryProbVector", p);
      for (int iti = 0; iti < p; iti++) {
        TotEveryProbVector[iti] = 0;
      }
      for (int iti = 0; iti < p; iti++) {
        RunProbVector[iti] = 0.0;
      }
      if (RegionWidth > 1 && ProbFixed != NULL) {
        if (RunProbRegionVector != NULL) { FFree(RunProbRegionVector, "RunProbRegionVector"); }
        RMemGetD(RunProbRegionVector, "RunProbRegionVector", p);
        for (int iti = 0; iti < p; iti++) {
          RunProbRegionVector[iti] = 0.0;
        }
      }
    } else {
      RMemGetD(RunProbVector, "RunProbVector", iFirstRandom + Rf_length(sOnTau)); 
      RMemGetI(TotEveryProbVector, "TotEveryProbVector", iFirstRandom+Rf_length(sOnTau));
      for (int iti = 0; iti < iFirstRandom+Rf_length(sOnTau); iti++) {
        TotEveryProbVector[iti] = 0;
      }
      for (int iti = 0; iti < iFirstRandom+Rf_length(sOnTau); iti++) {
        RunProbVector[iti] = 0;
      }
      if (RegionWidth > 1 && ProbTau != NULL) {
        if (RunProbRegionVector != NULL) { FFree(RunProbRegionVector, "RunProbRegionVector"); }
        RMemGetD(RunProbRegionVector, "RunProbRegionVector", (iFirstRandom+Rf_length(sOnTau)));
        for (int iti = 0; iti < iFirstRandom+Rf_length(sOnTau); iti++) {
          RunProbRegionVector[iti] = 0.0;
        }
      }
    } 
    StartRunProbVector = StartRunProbVector_;  TotInProbVector = 0;  
    RunLessZero = (unsigned int *) Calloc(p, unsigned int);
    RunMoreZero = (unsigned int *) Calloc(p, unsigned int);
    RunSumBeta = (double *)  Calloc(p, double);
    //MedianUp = (double *)  Calloc(p, double);
    //MedianDown = (double *) Calloc(p, double);
    for (int iti = 0; iti < p; iti++) {
      RunLessZero[iti] = (unsigned int) 0;
    }
    for (int iti = 0; iti < p; iti++) {
      RunMoreZero[iti] = (unsigned int) 0;
    }
    for (int iti = 0; iti < p; iti++) {
      RunSumBeta[iti] = 0.0;
    }
    //for (int iti = 0; iti < p; iti++) {
    //  MedianUp[iti] = -666.6;
    //}
    //for (int iti = 0; iti < p; iti++) {
    //  MedianDown[iti] = 666.6;
    //}
  }
  int SetupRunLessZeroMoreZeroSumBeta(SEXP sRunLessZero, SEXP sRunMoreZero, SEXP sRunSumBeta) {
    if (Rf_isNull(sRunLessZero) || Rf_length(sRunLessZero) != p) {
      Rf_error("SetupRunLessZeroMoreZeroSumBeta: No RunLessZero is NULL");
    }
    if (Rf_isNull(sRunMoreZero) || Rf_length(sRunMoreZero) != p) {
      Rf_error("SetupRunLessZeroMoreZeroSumBeta: No sRunMoreZero is NULL");
    }
    if (Rf_isNull(sRunSumBeta) || Rf_length(sRunSumBeta) != p) {
      Rf_error("SetupRunLessZeroMoreZeroSumBeta: No sRunSumBeta is NULL");
    }
    if (Verbose >= 2) {
      Rprintf("SetupRunLessZeroMoreZeroSumBeta: Setting up\n"); R_FlushConsole();
    }
    if (RunLessZero != NULL) { Free(RunLessZero); RunLessZero = NULL; }
    if (RunMoreZero != NULL) { Free(RunMoreZero); RunMoreZero = NULL; }
    if (RunSumBeta != NULL) { Free(RunSumBeta); RunSumBeta = NULL; }
    RunLessZero = (unsigned int *) Calloc(p, unsigned int);
    RunMoreZero = (unsigned int *) Calloc(p, unsigned int);
    RunSumBeta = (double *) Calloc(p, double);
    int iti;
    for (iti = 0; iti < p; iti++) {
      RunLessZero[iti] = Rf_isInteger(sRunLessZero) ? (unsigned int) INTEGER(sRunLessZero)[iti] : (unsigned int) REAL(sRunLessZero)[iti];
    }
    for (iti = 0; iti < p; iti++) {
      RunMoreZero[iti] = Rf_isInteger(sRunMoreZero) ? (unsigned int)  INTEGER(sRunMoreZero)[iti] : (unsigned int)  REAL(sRunMoreZero)[iti];
    }
    for (iti = 0; iti < p; iti++) {
      RunSumBeta[iti] = Rf_isInteger(sRunSumBeta) ? (unsigned int)  INTEGER(sRunSumBeta)[iti] : (unsigned int)  REAL(sRunSumBeta)[iti];
    }
    
    return(1);
  }
  SEXP get_RunLessZero() {
    if (RunLessZero == NULL) { return(R_NilValue); }
    SEXP sOn = R_NilValue;
    Rf_protect(sOn = Rf_allocVector(INTSXP, p));
    for (int iti = 0; iti < p; iti++) {
      INTEGER(sOn)[iti] = RunLessZero[iti];
    }
    Rf_unprotect(1);
    return(sOn);
  }
  SEXP get_RunSumBeta() {
    if (RunSumBeta == NULL) { return(R_NilValue); }
    SEXP sOn = R_NilValue;
    Rf_protect(sOn = Rf_allocVector(REALSXP, p));
    for (int iti = 0; iti < p; iti++) {
      REAL(sOn)[iti] = RunSumBeta[iti];
    }
    Rf_unprotect(1);
    return(sOn);  
  }
  SEXP get_QMedianBeta() {
     if (RunSumBeta == NULL || RunMoreZero == NULL || RunLessZero == NULL  || TotEveryProbVector == NULL) {
       return(R_NilValue);
     }
     SEXP sOn = R_NilValue;
     Rf_protect(sOn = Rf_allocVector(REALSXP, p));
     for (int iti = 0; iti < p; iti++) {
       if ((double) RunMoreZero[iti] > ((double) TotEveryProbVector[iti])/ 2.0 ||
         (double) RunLessZero[iti] > ((double) TotEveryProbVector[iti]) /2.0 ) {
         REAL(sOn)[iti] = RunSumBeta[iti] / TotEveryProbVector[iti];
       } else {
         REAL(sOn)[iti] = 0.0;
       }
     }
     Rf_unprotect(1);
     return(sOn);
  }
  SEXP BSDeDouble(SEXP saIn);
  SEXP get_RunMoreZero() {
    if (RunMoreZero == NULL) { return(R_NilValue); }
    SEXP sOn = R_NilValue;
    Rf_protect(sOn = Rf_allocVector(INTSXP, p));
    for (int iti = 0; iti < p; iti++) {
      INTEGER(sOn)[iti] = RunMoreZero[iti];
    }
    Rf_unprotect(1);
    return(sOn);
  }
  int UpdateRunMoreZero() {
    if (RunMoreZero == NULL || RunLessZero == NULL || RunSumBeta == NULL) {
      return(-2);
    }
    if (sBeta == NULL || Rf_isNull(sBeta) || !Rf_isReal(sBeta)  || Rf_length(sBeta) != p) {
      Rprintf("UpdateRunMoreZero(): no way, sBeta is NULL!\n");
      return(-666);
    }
    if (tt < StartRunProbVector) {
      return(-3);
    }
    for (int iti = 0; iti < p; iti++) {
      if (REAL(sBeta)[iti] > 0.0) {
        RunMoreZero[iti]++;
      } else if (REAL(sBeta)[iti] < 0.0) {
        RunLessZero[iti]++;
      }
      if (REAL(sBeta)[iti] != 0.0) {
        RunSumBeta[iti] += REAL(sBeta)[iti];
      }
    }
    return(1);
  }
  int get_TotInProbVector() { return(TotInProbVector); }
  SEXP GetTotEveryProbVector() {
    if (TotEveryProbVector == NULL) { return(R_NilValue); }
    SEXP sOut = R_NilValue;        int CC = 0;
    if (Rf_isNull(sOnTau) || Rf_length(sOnTau) <= 0) {
     CC = p;
    } else {
      CC = iFirstRandom + Rf_length(sOnTau);
    }
    Rf_protect(sOut = Rf_allocVector(INTSXP, CC));
    for (int ii = 0; ii < CC; ii++) {
      INTEGER(sOut)[ii] = TotEveryProbVector[ii];
    }
    Rf_unprotect(1);
    return(sOut);
  }
  SEXP GetRunProbVector() {
    if (RunProbVector == NULL) { return(R_NilValue); }
    int One = 1;  int CC = 0;
    SEXP sOut = R_NilValue;
    if (Rf_isNull(sOnTau) || Rf_length(sOnTau) <= 0) {
      CC = p;
    } else {
      CC = iFirstRandom + Rf_length(sOnTau);
    }
    Rf_protect(sOut = Rf_allocVector(REALSXP, CC));
    F77_CALL(dcopy)(&CC, RunProbVector, &One, REAL(sOut), &One);
    Rf_unprotect(1);
    return(sOut);
  }
  SEXP GetProbVector() {
    if (RunProbVector == NULL) {
      Rprintf("Can't give ProbVector because it wasn't set up \n"); R_FlushConsole();
      return(R_NilValue);
    }
    SEXP sOut = R_NilValue;
    SEXP sStrings = R_NilValue;  int One = 1;
    double Scale = 1.0 / TotInProbVector;
    char VectorName[100];
    if (Verbose >= 8) {
       Rprintf("GetProbVector: Scale is %f.\n", Scale); R_FlushConsole();
    }
    int CC = 0;
    if (Rf_isNull(sOnTau) || Rf_length(sOnTau) <= 0) {
      Rf_protect(sOut = Rf_allocVector(REALSXP, p));
      F77_CALL(dcopy)(&p, RunProbVector, &One, REAL(sOut), &One);
      for (int iti = 0; iti < p; iti++) {
         //if (TotInProbVector != TotEveryProbVector[iti]) {
           //Rprintf("GetProbVector:  Note TotEveryProbVector[%d]=%d, TotInProbVector = %d \n",
           //  iti, TotEveryProbVector[iti], TotInProbVector); R_FlushConsole();
         //}
         REAL(sOut)[iti]  = REAL(sOut)[iti] / TotEveryProbVector[iti];
      }
      Rf_protect(sStrings=Rf_allocVector(STRSXP,p));
      for (int iti = 0; iti < p; iti++) {
        for (int jtj2 = 0; jtj2 < 100-1; jtj2++) {
          VectorName[jtj2] = '\0';
        }
        sprintf(VectorName, "ProbFixed:%d", iti+1);
        SET_STRING_ELT(sStrings, iti, Rf_mkChar(VectorName));
      } 
      Rf_setAttrib(sOut, R_NamesSymbol, sStrings);
      Rf_unprotect(2); 
      return(sOut);
    }
    Rf_protect(sOut = Rf_allocVector(REALSXP, iFirstRandom + Rf_length(sOnTau)));
    
    CC  = iFirstRandom + Rf_length(sOnTau);
    F77_CALL(dcopy)(&CC, RunProbVector, &One, REAL(sOut), &One);
    for (int iti = 0; iti < CC; iti++) {
         //if (TotInProbVector != TotEveryProbVector[iti]) {
         //  Rprintf("GetProbVector:  Note TotEveryProbVector[%d]=%d, TotInProbVector = %d \n",
         //    iti, TotEveryProbVector[iti], TotInProbVector); R_FlushConsole();
         //}
         REAL(sOut)[iti]  = REAL(sOut)[iti] / TotEveryProbVector[iti];
    }
    Rf_protect(sStrings=Rf_allocVector(STRSXP,CC));
    if (iFirstRandom > 0) {
      for (int iti = 0; iti < iFirstRandom; iti++) {
        for (int jtj2 = 0; jtj2 < 100-1; jtj2++) {
          VectorName[jtj2] = '\0';
        }
        sprintf(VectorName, "ProbFixed:%d", iti+1);
        SET_STRING_ELT(sStrings, iti, Rf_mkChar(VectorName));
      } 
    } 
    for (int iti = 0; iti < Rf_length(sOnTau); iti++) {
        for (int jtj2 = 0; jtj2 < 100-1; jtj2++) {
          VectorName[jtj2] = '\0';
        }
        sprintf(VectorName, "ProbTau:%d", iti+1);
        SET_STRING_ELT(sStrings, iFirstRandom+iti, Rf_mkChar(VectorName));
    } 
    
    Rf_setAttrib(sOut, R_NamesSymbol, sStrings);
    Rf_unprotect(2); return(sOut);
  }
  ////////////////////////////////////////////////////////////
  //  YResidSq is the square sum residuals, or sum( (Y-X %*% Beta)^2 )
  void SetupYResidSq() {
    // sum( (Y-Xb)^2 ) = sum(Y(Y-Xb)) -  sum(bX^T(Y-Xb))
    if (Verbose > 2) {
      Rprintf("SetupYResidSq: Starting\n"); R_FlushConsole();
    }
    if (iiWeight != NULL) {SetupYSq();}
    int One = 1;
    
    if (XtY == NULL) {
      Rf_error("SetupYResidSq:  Error, XtY is NULL!\n");
    }
    if (XtResid == NULL) { Rf_error("SetupYResidSq: SEtup XtResid First!\n"); }
    YResidSq = YSq - 
      F77_CALL(ddot)(&p, XtY, &One, REAL(sBeta), &One) -
      F77_CALL(ddot)(&p, XtResid, &One, REAL(sBeta), &One);
    if (Verbose > 2) {
      Rprintf("SetupYResidSq: finished\n"); R_FlushConsole();
    }
    if (YResidSq < 0 || R_isnancpp(YResidSq)) {
      Rprintf("Error in SetupYResidSq: sum YSq = %f, XtY*Beta=%f, XtResid*sBeta=%f\n",
        YSq, F77_CALL(ddot)(&p, XtY, &One, REAL(sBeta), &One), 
        F77_CALL(ddot)(&p, XtResid, &One, REAL(sBeta), &One)); R_FlushConsole();
      Rf_error("SetupYResidSq: Uh Oh, sXtResid is wrong and YResidSq = %f.\n, check the problem!", YResidSq);
    }
  }
  SEXP get_PriorProbs()  {
    //SEXP sOn = R_NilValue;
    if (rPriorProbFixed == NULL) {
      Rprintf("PriorProbFixed is NULL!");
    } else if (Rf_isNull(rPriorProbFixed->asSexp())) {
      Rprintf("rPriorProbFixed->asSexp() is NilValue!");
    } else {
      if (rPriorProbTau == NULL) {
        Rprintf("rPriorProbfixed is not null, but rPriorProbTau is Null!");
        R_FlushConsole();
        return(rPriorProbFixed->asSexp());
      } else if (Rf_isNull(rPriorProbTau->asSexp()) ) {
        Rprintf("rPriorProbTau->asSexp is a Nilvalue. but rPriorProbFixed is good");
        return(rPriorProbFixed->asSexp());
      } else {
        Rprintf("Returning All PriorProbs.");
        return(Rf_cons(rPriorProbFixed->asSexp(), rPriorProbTau->asSexp()));
      }
    }
    if (rPriorProbTau == NULL) {
       Rprintf("rPriorProbfixed is null, and rPriorProbTau is Null!");
       R_FlushConsole();
       return(R_NilValue);
    } else if (Rf_isNull(rPriorProbTau->asSexp()) ) {
      Rprintf("rPriorProbTau->asSexp is a Nilvalue. but rPriorProbFixed is also null");
      return(R_NilValue);
    }
    Rprintf("returning just rPriorProbTau!");
    return(rPriorProbTau->asSexp());
  }
  void set_PriorProbFixed(SEXP ssPriorProbFixed) {
    if (Rf_isNull(ssPriorProbFixed)) {
     DDelete(rPriorProbFixed, "rPriorProbFixed");
     rPriorProbFixed = NULL; sPriorProbFixed = R_NilValue;
     return;
    }
    if (!Rf_isReal(ssPriorProbFixed)) {
      Rprintf("set_PriorProbFixed: ERROR Supplied ssPriorProbFixed is not REAL!\n");   R_FlushConsole();
      if (Rf_isInteger(ssPriorProbFixed)) {
        Rf_error("set_PriorProbFixed: Apparently ssPriorProbFixed is an Integer!\n");
      } else {
        Rf_error("I don't know what ssPriorProbFixed is !\n");
      }
    }
    if (sOnTau == NULL || Rf_isNull(sOnTau) || Rf_length(sOnTau) <= 0 ||
      tauEndList == NULL || Rf_isNull(tauEndList) || Rf_length(tauEndList) <= 0) {
      if (p != Rf_length(ssPriorProbFixed)) {
        Rprintf("set_PriorProbFixed: Hey, no sOnTau but p=%d, ", p);
        Rprintf("and length Sub = %d!\n",
           Rf_length(ssPriorProbFixed));   R_FlushConsole();
        Rf_error("ERROR: set_PriorProbFixed But that's that.\n");
      }  
    } else {
      if (iFirstRandom <= 0) {
        Rf_error("set_PriorProbFixed: Sorry, There are no fixed paramters");
      }
      if (Rf_length(ssPriorProbFixed) != iFirstRandom) {
        if (ssPriorProbFixed == NULL || Rf_isNull(ssPriorProbFixed)) {
          Rprintf("set_rsPriorProbFixed: you supplied NULL ssPriorProbFixed!\n");
          R_FlushConsole();
        }
        Rprintf("set_rsPriorProbFixed: Sorry, length of ssPriorProbFixed %d\n", 
          Rf_length(ssPriorProbFixed));
        Rprintf("But iFirst Random is : %d. \n", iFirstRandom);
        Rprintf("We Just can't set things that way."); R_FlushConsole();
        rPriorProbFixed = (AObject*) new AObject(
          Rf_allocVector(REALSXP, iFirstRandom));
        if (Rf_length(ssPriorProbFixed) >= 1) {
          for (int tjt = 0; tjt < iFirstRandom; tjt++) {
            if (tjt < Rf_length(ssPriorProbFixed)) {
              REAL(rPriorProbFixed->asSexp())[tjt] = REAL(ssPriorProbFixed)[tjt];
            } else {
              REAL(rPriorProbFixed->asSexp())[tjt] = REAL(ssPriorProbFixed)[
                Rf_length(ssPriorProbFixed)-1];
            }
          }
        }
        sPriorProbFixed=rPriorProbFixed->asSexp();
        Rprintf("PriorProbFixed Well we did a hack to save it!\n"); R_FlushConsole();
      }
    } 
    DDelete(rPriorProbFixed, "rPriorProbFixed");
    rPriorProbFixed = (AObject*) new AObject(ssPriorProbFixed);
    sPriorProbFixed = (SEXP) rPriorProbFixed->asSexp();
    HadSetPriorProbFixed = 1;
    if (Verbose >= 1) {
      Rprintf("set_PriorProbFixed: Successful!\n"); R_FlushConsole();
    }
  }
  void set_gPrior(double In) {
    if (In <= 0.0) {
      gPrior = 0.0; return;
    }
    if (Verbose >= 1) {
      Rprintf("set_gPrior: We will set Zellner g prior to %f \n", In); R_FlushConsole();
    }
    gPrior = In; return;
  }
  double get_gPrior() { return(gPrior); }
  SEXP get_PriorProbFixed() {
    if( rPriorProbFixed == NULL ) { return(R_NilValue); }
    if (Rf_length(rPriorProbFixed->asSexp()) != iFirstRandom) {
      Rprintf("get_PriorProbFixedError: the length is only %d, but iFirsRandom = %d\n",
        Rf_length(rPriorProbFixed->asSexp()), iFirstRandom);
      Rprintf(" Set to null now!");
      DDelete(rPriorProbFixed, "rPriorProbFixed");
      sPriorProbFixed = R_NilValue;
    }
    return(rPriorProbFixed->asSexp());
  }
  int AddToRunProbRegionVector();
  SEXP get_RunProbRegionVector() {
    if (RegionWidth <= 1 && RunProbRegionVector != NULL) {
      Rf_error("get_RunProbRegionVector: Error, this should not be operated upon!, RegionWidth = %d!\n", RegionWidth);
    }
    if (RunProbRegionVector == NULL) { return(R_NilValue);}
    SEXP sOn = R_NilValue;
    int TotFixed = p; int TotalAll = p;
    if (sOnTau == NULL || Rf_isNull(sOnTau) || Rf_length(sOnTau) <= 0) {
    
    } else if (iFirstRandom >= 0 && iFirstRandom < p) {
      TotFixed = iFirstRandom; TotalAll = TotFixed + Rf_length(sOnTau);
    } else if (iFirstRandom < 0 && sOnTau != NULL && !Rf_isNull(sOnTau) && Rf_length(sOnTau) >= 1) {
      TotFixed = 0;  TotalAll = Rf_length(sOnTau);
    }
    if (Verbose >= 1) {
      Rprintf("get_RunProbRegionVector(), Allocate TotalAll = %d \n", TotalAll);
      R_FlushConsole();
    }
    Rf_protect(sOn = Rf_allocVector(REALSXP, TotalAll));
    for (int ii = 0; ii < TotalAll; ii++) {
      REAL(sOn)[ii] = (double) RunProbRegionVector[ii];
    }
    Rf_unprotect(1); return(sOn);
  }
  int SetRegionWidthAndRunProbRegionVector(SEXP sRegionWidth, SEXP sRunProbRegionVector) {
    if (sRegionWidth == NULL  || Rf_isNull(sRegionWidth) || Rf_length(sRegionWidth) <= 0) {
      Rf_error("SetRegionWidthAndRunProbRegionVector, no this will not go!\n");
    }
    int TotFixed = p; int TotalAll = p;
    if (sOnTau == NULL || Rf_isNull(sOnTau) || Rf_length(sOnTau) <= 0) {
    
    } else if (iFirstRandom >= 0 && iFirstRandom < p) {
      TotFixed = iFirstRandom; TotalAll = TotFixed + Rf_length(sOnTau);
    } else if ((iFirstRandom < 0 || iFirstRandom > p) &&
      sOnTau != NULL && !Rf_isNull(sOnTau) && Rf_length(sOnTau) >= 1) {
      TotFixed = 0;  TotalAll = Rf_length(sOnTau);  
    }
    if (sRunProbRegionVector == NULL || Rf_isNull(sRunProbRegionVector) ||
      !Rf_isReal(sRunProbRegionVector)  || !(Rf_length(sRunProbRegionVector) == TotalAll)) {
      Rf_error("Set RunProbRegionVecto, no you supplied bad sRunProbRegionVector\n"); 
    }
    if (TotEveryProbVector == NULL) {
      Rf_error("SetRegionWidthAndRunProbRegionVector: TotEveryProbVector Is Not Set!\n");
    }
    int iRegionWidth = 0;
    iRegionWidth = Rf_isInteger(sRegionWidth) ? INTEGER(sRegionWidth)[0] : (int) REAL(sRegionWidth)[0];
    if (iRegionWidth <= 1) {
      Rprintf("SetRegionWidthAndRunProbRegionVector: No, iRegionWidth is %d\n", iRegionWidth);
    }
    RegionWidth = iRegionWidth;
    if (RunProbRegionVector != NULL) { Free(RunProbRegionVector);  RunProbRegionVector = NULL; }
    RunProbRegionVector = (double *) Calloc(TotalAll, double);
    for (int ii = 0; ii < TotalAll; ii++) {
      RunProbRegionVector[ii] = REAL(sRunProbRegionVector)[ii];
    }
    return(1);
  }
  SEXP get_PriorProbTau() {
    if (rPriorProbTau == NULL) { return(R_NilValue); }
    if (Rf_isNull(tauEndList)) { return(R_NilValue); }
    if (Rf_length(rPriorProbTau->asSexp()) != Rf_length(tauEndList) ) {
      Rprintf("get_PriorProbFixedTau: the length is only %d, but length(tau) = %d\n",
        Rf_length(rPriorProbTau->asSexp()), Rf_length(tauEndList) );
      Rprintf(" Set to null now!");
      DDelete(rPriorProbTau, "rPriorProbTau");
      sPriorProbTau = R_NilValue;
    }
    if (rPriorProbTau == NULL) { return(R_NilValue); }
    return(rPriorProbTau->asSexp());
  }
  void set_PriorProbTau(SEXP rsPriorProbTau) {
    if (Rf_isNull(rsPriorProbTau)) {
     DDelete(rPriorProbTau, "rPriorProbTau");
     sPriorProbTau = R_NilValue;
     return;
     rPriorProbTau = NULL; sPriorProbTau = R_NilValue;
    }
    if (!Rf_isReal(rsPriorProbTau)) {
      Rprintf("set_PriorProbTau: Error PriorProbTau is not a real!");
      R_FlushConsole();
      if (Rf_isInteger(rsPriorProbTau)) {
        Rf_error("set_PriorProbTau: It is an integer!\n");
      } else {
        Rf_error("set_PriorProbTau: I don't know what it is!\n");
      }
    }
    if (Rf_isNull(tauEndList) || tauEndList == NULL) {
      DDelete(rPriorProbTau, "rPriorProbTau");
      sPriorProbTau = R_NilValue;
      return;
      DDelete(rPriorProbTau, "rPriorProbTau");
      Rf_error("set_PriorProbTau: Sorry, There are no random paramters");
    }
    if (Rf_length(rsPriorProbTau) != Rf_length(tauEndList)) {
      Rf_error("set_rsPriorProbTau: Sorry, given PriorProbTau = Length %d, use length %d\n",
        Rf_length(rsPriorProbTau), Rf_length(tauEndList));
    }
    DDelete(rPriorProbTau, "rPriorProbTau");
    rPriorProbTau = (AObject*) new AObject(rsPriorProbTau);
    sPriorProbTau = rPriorProbTau->asSexp();
    HadSetPriorProbTau = 1;
  }
  ////////////////////////////////////////////////////////////////
  // tauFixed():  Role is to represent standard deviation of 
  //    ``Fixed'' coefficients, so that Bayesian inference can be performed
  //    comparing Select-out and Select-in hypotheses
  SEXP get_tauFixed() {  return(tauFixed); }
  void set_tauFixed(SEXP tauFixed_) { 
    DDelete(RtauFixed, "RTauFixed"); 
    if (Rf_isNull(tauFixed_) || Rf_length(tauFixed_) <= 0) {
      if (Verbose >= 1) {
        Rprintf("set_tauFixed: NOTE Setting to NULL!\n");    R_FlushConsole();
      }
      RtauFixed = NULL; tauFixed = R_NilValue;
      return;
    } else if (!Rf_isReal(tauFixed_)) {
      if (Verbose >= 1) {
         Rprintf("BayesSpikeCpp.h:set_tauFixed: note that tauFixed_ is not real!\n");
      }
      if (!Rf_isInteger(tauFixed_)) {
        Rf_error("BayesSpikeCpp.h:set_tauFixed: tau_fixed_ isn't even integer!\n");
      }
      RtauFixed = new AObject(Rf_allocVector(REALSXP, Rf_length(tauFixed_)));
      for (int jj = 0; jj < Rf_length(tauFixed_); jj++) {
        REAL(RtauFixed->asSexp())[jj] = INTEGER(tauFixed_)[jj];
      }
      tauFixed = RtauFixed->asSexp();
    } else {
      RtauFixed = new AObject(tauFixed_);
      tauFixed = RtauFixed->asSexp();
    } 
    if (Verbose >= 1) {
       printf("set_tauFixed: tauFixed Just set, to length %d and first element %f \n",
         Rf_length(tauFixed), REAL(tauFixed)[0]);   R_FlushConsole();
    }
  }
  SEXP get_tauFixedPrior() {  return(tauFixedPrior); }
  void set_tauFixedPrior(SEXP tauFixedPrior_) { 
    DDelete(RtauFixedPrior, "RTauFixedPrior");  tauFixedPrior = R_NilValue;
    if (Rf_isNull(tauFixedPrior_) || Rf_length(tauFixedPrior_) <= 0) {
      if (Verbose >= 1) {
        Rprintf("set_tauFixedPrior: NOTE Setting to NULL!\n");    R_FlushConsole();
      }
      RtauFixedPrior = NULL; tauFixedPrior = R_NilValue;
      return;
    } else if (!Rf_isReal(tauFixedPrior_)) {
      if (Verbose >= 1) {
        Rprintf("set_tauFixedPrior: Issue, RtauFixedPrior given is not actually real!\n"); R_FlushConsole();
      }
      if (!Rf_isInteger(tauFixedPrior_)) {
        Rf_error("BayesSpikeCpp.h:set_tauFixeedPrior: RtauFixedPrior isn't Integer either, big error!\n");
      }
      SEXP NewtauFixedPrior = Rf_allocVector(REALSXP, Rf_length(tauFixedPrior_));
      RtauFixedPrior = new AObject(NewtauFixedPrior);
      for (int jj = 0; jj < Rf_length(tauFixedPrior_); jj++) {
        REAL(RtauFixedPrior->asSexp())[jj] = INTEGER(tauFixedPrior_)[jj];
      }
      tauFixedPrior = RtauFixedPrior->asSexp();      
    } else {
      RtauFixedPrior = new AObject(tauFixedPrior_);
      tauFixedPrior = RtauFixedPrior->asSexp();
    } 
    if (Verbose >= 1) {
       printf("set_tauFixedPrior: tauFixedPrior Just set, to length %d and first element %f \n",
         Rf_length(tauFixedPrior), REAL(tauFixedPrior)[0]);   R_FlushConsole();
    }
  }
  double get_YSq() { return(YSq);}
  double get_YResidSq() {return(YResidSq); }
  
  int get_HadSetPriorProbFixed() { return(HadSetPriorProbFixed); }
  int get_HadSetPriorProbTau() {return(HadSetPriorProbTau); }
  //////////////////////////////////////////////////////////////
  //  XjSq
  //     XjSq is the square of diagonals of XtX, which changes based upon weights
  SEXP get_XjSq() {
    if (XjSq == NULL ||  p <= 0) { return(R_NilValue); }
    SEXP sOut;
    Rf_protect(sOut = Rf_allocVector(REALSXP, p));
    int One = 1; 
    F77_CALL(dcopy)(&p, XjSq, &One, REAL(sOut), &One);
    Rf_unprotect(1);
    return(sOut);
  }  
  int SetupXjSq();
  SEXP get_SmallXtX() {
    if (rSmallXtX == NULL || NumActive == 0) { return(R_NilValue); }
    SEXP sOut;
    Rf_protect(sOut = Rf_allocMatrix(REALSXP, NumActive, NumActive));
    int One = 1;  double ZeroD = 0.0;  int NASq = NumActive * NumActive;
    F77_CALL(dscal)(&NASq, &ZeroD, REAL(sOut), &One);
    int ii = 0; int CI; int NASqPQ = NumActive * (NumActive + 1) / 2;
    int dd = NumActive;  int kk = 1; int ii2 = 0;
    for (CI = 0; CI < NASqPQ; CI++) {
      REAL(sOut)[ii] = rSmallXtX[CI];
      REAL(sOut)[ii2] = rSmallXtX[CI];
      ii++;  ii2 += NumActive;
      if (ii >= dd) {  ii = kk * NumActive + kk;  dd+=NumActive; 
        ii2 = ii; kk++; }
    }
    Rf_unprotect(1);
    return(sOut);
  }
  SEXP RICholSubDim(int SubDim) {
    if (rIChol == NULL || NumActive == 0) { return(R_NilValue); }
    if (SubDim > NumActive) { return(R_NilValue); }
    SEXP sOut;
    Rf_protect(sOut = Rf_allocMatrix(REALSXP, SubDim, SubDim));
    int One = 1;  double ZeroD = 0.0;  int NASq = SubDim * SubDim;
    F77_CALL(dscal)(&NASq, &ZeroD, REAL(sOut), &One);
    int ii = 0;
    for (ii = 0; ii < NASq; ii++) {
      REAL(sOut)[ii] = 0.0;
    }
    int CI; int NASqPQ = SubDim * (SubDim + 1) / 2;
    int dd = SubDim;  int kk = 1;
    ii = 0;
    for (CI = 0; CI < NASqPQ; CI++) {
      REAL(sOut)[ii] = rIChol[CI];
      ii++;  if (ii >= dd) {  ii = kk * SubDim + kk;  dd+=SubDim; kk++; }
    }
    Rf_unprotect(1);
    return(sOut);  
  }
  SEXP get_rICholSym() {
    if (rIChol == NULL || NumActive == 0) { return(R_NilValue); }
    SEXP sOut;
    Rf_protect(sOut = Rf_allocMatrix(REALSXP, NumActive, NumActive));
    int One = 1;  double ZeroD = 0.0;  int NASq = NumActive * NumActive;
    F77_CALL(dscal)(&NASq, &ZeroD, REAL(sOut), &One);
    int ii = 0; int CI; int NASqPQ = NumActive * (NumActive + 1) / 2;
    int dd = NumActive;  int kk = 1;
    int gg = 0;  int goG = 0;
    for (CI = 0; CI < NASqPQ; CI++) {
      REAL(sOut)[ii] = rIChol[CI];
      REAL(sOut)[gg] = rIChol[CI];
      gg+=NumActive;
      ii++;  if (ii >= dd) {  ii = kk * NumActive + kk;  dd+=NumActive; kk++; }
      if (gg >= NumActive*NumActive) {
        goG++;  gg=goG;
      }
    }
    Rf_unprotect(1);
    return(sOut);  
  }
  SEXP get_rIChol() {
    if (rIChol == NULL || NumActive == 0) { return(R_NilValue); }
    SEXP sOut;
    Rf_protect(sOut = Rf_allocMatrix(REALSXP, NumActive, NumActive));
    int One = 1;  double ZeroD = 0.0;  int NASq = NumActive * NumActive;
    F77_CALL(dscal)(&NASq, &ZeroD, REAL(sOut), &One);
    int ii = 0; int CI; int NASqPQ = NumActive * (NumActive + 1) / 2;
    int dd = NumActive;  int kk = 1;
    for (CI = 0; CI < NASqPQ; CI++) {
      REAL(sOut)[ii] = rIChol[CI];
      ii++;  if (ii >= dd) {  ii = kk * NumActive + kk;  dd+=NumActive; kk++; }
    }
    Rf_unprotect(1);
    return(sOut);
  }
  SEXP get_rICholPart2() {
      if (rIChol == NULL || NumActive == 0) { return(R_NilValue); }
    SEXP sOut;
    Rf_protect(sOut = Rf_allocMatrix(REALSXP, NumActiveBack, NumActiveBack));
    int One = 1;  double ZeroD = 0.0;  int NASq = NumActiveBack * NumActiveBack;
    F77_CALL(dscal)(&NASq, &ZeroD, REAL(sOut), &One);
    int ii = 0; int CI; int NASqPQ = NumActiveBack * (NumActiveBack + 1) / 2;
    int dd = NumActiveBack;  int kk = 1;
    for (CI = 0; CI < NASqPQ; CI++) {
      REAL(sOut)[ii] = rIChol[CI];
      ii++;  if (ii >= dd) {  ii = kk * NumActiveBack + kk;  dd+=NumActiveBack; kk++; }
    }
    Rf_unprotect(1);
    return(sOut);
  }
  SEXP get_rICholBack() {
    if (rICholBack == NULL || NumActive == 0) { return(R_NilValue); }
    SEXP sOut;
    Rf_protect(sOut = Rf_allocMatrix(REALSXP, NumActive, NumActive));
    int One = 1;  double ZeroD = 0.0;  int NASq = NumActive * NumActive;
    F77_CALL(dscal)(&NASq, &ZeroD, REAL(sOut), &One);
    int ii = 0; int CI; int NASqPQ = NumActive * (NumActive + 1) / 2;
    int dd = NumActive;  int kk = 1;
    for (CI = 0; CI < NASqPQ; CI++) {
      REAL(sOut)[ii] = rICholBack[CI];
      ii++;  if (ii >= dd) {  ii = kk * NumActive + kk;  dd+=NumActive; kk++; }
    }
    Rf_unprotect(1);
    return(sOut);
  }
  void set_BayesSpikeNameSpace(SEXP BSNS) {
    DDelete(RBayesSpikeNameSpace, "RBayesSpikeNameSpace"); 
    RBayesSpikeNameSpace = new AObject(BSNS); 
    BayesSpikeNameSpace = RBayesSpikeNameSpace->asSexp();
  }
  SEXP get_BayesSpikeNameSpace()  { return(BayesSpikeNameSpace);}
  SEXP get_D() {
    if (MT->NumJ <= 0) {  return(R_NilValue); }
    SEXP sOut;
    Rf_protect(sOut = Rf_allocVector(REALSXP, MT->NumJ));
    int One = 1;  
    if (MT->D == NULL) {
      Rf_error("get_D: sorry, we're NULL \n");
    }
    F77_CALL(dcopy)(&MT->NumJ, MT->D, &One, REAL(sOut), &One);
    double iSigma = 1.0 / REAL(sOnSigma)[0];
    F77_CALL(dscal)(&MT->NumJ, &iSigma, REAL(sOut), &One);
    Rf_unprotect(1);
    return(sOut);
  }
  double SolveAlpha();
  double get_d4();
  double get_d3() {
    return(MT->NumJ + REAL(MT->staupriordf)[0]+1);
  }
  SEXP r2o(); SEXP r2(int num);  // Random Draws from proposal distribution from SampleTauSpliced
  SEXP get_R() {
    if (MT->NumJ <= 0) {  return(R_NilValue); }
    SEXP sOut;
    Rf_protect(sOut = Rf_allocVector(REALSXP, MT->NumJ));
    int One = 1;  
    F77_CALL(dcopy)(&MT->NumJ, MT->R, &One, REAL(sOut), &One);
    Rf_unprotect(1);
    return(sOut);
  }
  SEXP RISubDim(int SubDim) {
   if (rI == NULL || NumActive == 0) { return(R_NilValue); }
   if (SubDim > NumActive) { return(R_NilValue); }
    SEXP sOut;
    Rf_protect(sOut = Rf_allocMatrix(REALSXP, SubDim, SubDim));
    int One = 1;  double ZeroD = 0.0;  int NASq = SubDim * SubDim;
    F77_CALL(dscal)(&NASq, &ZeroD, REAL(sOut), &One);
    int ii = 0; int CI; int NASqPQ = SubDim * (SubDim + 1) / 2;
    int dd = SubDim;  int kk = 1; int ii2 = 0;
    for (CI = 0; CI < NASqPQ; CI++) {
      REAL(sOut)[ii] = rI[CI];
      REAL(sOut)[ii2] = rI[CI];
      ii++;  ii2 += SubDim;
      if (ii >= dd) {  ii = kk * SubDim + kk;  dd+=SubDim; ii2 = ii; kk++;}
    }
    Rf_unprotect(1);
    return(sOut);
  }
  SEXP get_rI() {
    if (rI == NULL || NumActive == 0) { return(R_NilValue); }
    SEXP sOut;
    Rf_protect(sOut = Rf_allocMatrix(REALSXP, NumActive, NumActive));
    int One = 1;  double ZeroD = 0.0;  int NASq = NumActive * NumActive;
    F77_CALL(dscal)(&NASq, &ZeroD, REAL(sOut), &One);
    int ii = 0; int CI; int NASqPQ = NumActive * (NumActive + 1) / 2;
    int dd = NumActive;  int kk = 1; int ii2 = 0;
    for (CI = 0; CI < NASqPQ; CI++) {
      REAL(sOut)[ii] = rI[CI];
      REAL(sOut)[ii2] = rI[CI];
      ii++;  ii2 += NumActive;
      if (ii >= dd) {  ii = kk * NumActive + kk;  dd+=NumActive; ii2 = ii; kk++;}
    }
    Rf_unprotect(1);
    return(sOut);
  }
  void set_AFD(SEXP _AFD) {
    if (Rf_isNull(_AFD)) {  Rprintf("AFD set, error null given \n"); }
    DDelete(RAFD, "RAFD");
    RAFD = new AObject(_AFD);
    AFD = RAFD->asSexp();
  }
  SEXP get_AFD() { return(AFD); }
  SEXP get_SmallXtResid() {
    if (SmallXtResid == NULL) { return(R_NilValue); }
    SEXP sOut;
    if (MT == NULL || MT->Onii < 0 || MT->Onii >= MT->NumJ || MT->NumJ <= 0) {
    Rf_protect(sOut = Rf_allocVector(REALSXP, MaxTauList));
    int One = 1;
    F77_CALL(dcopy)(&MaxTauList, SmallXtResid, &One, REAL(sOut), &One);
    } else {
      Rf_protect(sOut = Rf_allocVector(REALSXP, MT->NumJ));
      int One = 1;
      F77_CALL(dcopy)(&MT->NumJ, SmallXtResid, &One, REAL(sOut), &One);    
    }
    Rf_unprotect(1);
    return(sOut);
  }
  int get_JustDidPartBeta() {
    return(JustDidPartBeta);
  }
  int get_JustDidPartBFixed() {
    return(JustDidPartBFixed);
  }
  SEXP get_smXtY() {
    if (NumActiveBack < 0 || NumActiveBack > OnKappaMem) {return(R_NilValue); }
    SEXP sOut;
    Rf_protect(sOut = Rf_allocVector(REALSXP, NumActiveBack));
    int One = 1; F77_CALL(dcopy)(&NumActiveBack, smXtY, &One, REAL(sOut), &One);
    Rf_unprotect(1);
    return(sOut);
  }
  SEXP get_cXtY() {
    if (NumActive < 0 || NumActive > OnKappaMem) {return(R_NilValue); }
    SEXP sOut;
    Rf_protect(sOut = Rf_allocVector(REALSXP, NumActive));
    int One = 1; F77_CALL(dcopy)(&NumActive, cXtY, &One, REAL(sOut), &One);
    Rf_unprotect(1);
    return(sOut);
  }
  SEXP get_SmallRVec() {
    if (SmallRVec == NULL) { return(R_NilValue); }
    SEXP sOut;
    Rf_protect(sOut = Rf_allocVector(REALSXP, MaxTauList));
    int One = 1;
    F77_CALL(dcopy)(&MaxTauList, SmallRVec, &One, REAL(sOut), &One);
    Rf_unprotect(1);
    return(sOut);
  }
  void SetupMT() {
    if (MT != NULL) {
      if (Verbose > 2) { Rprintf(" Freeing MT !\n"); R_FlushConsole(); }
      ffree(MT, "MT");
    }
    if (Rf_isNull(staubarnu)) { Rf_error("SetupMT: Please set taubarnu!\n"); }
    if (Rf_isNull(staupriordf)) { Rf_error("SetupMT: Please set staupriordf!\n"); }
    SEXP sTemp = R_NilValue;
    Rf_protect(sTemp = Rf_allocVector(REALSXP, 1));
    REAL(sTemp)[0] = 1.0;
    MT = NewTauContainer(sPriorOftau, staubarnu, staupriordf, sPriorXScaling,
      sOnPiA, sTemp);
    Rf_unprotect(1);
    if (MT == NULL) { Rf_error(" SetupMT: Failed! \n"); }
  }
  int SetupOrderAttack();
  SEXP get_OrderAttack() {
   int ii;
   if (OrderAttack == NULL) { return(R_NilValue); }
    SEXP sOut;  
    int AA = iFirstRandom;
    if (tauEndList == NULL || Rf_isNull(tauEndList) || Rf_length(tauEndList) <= 0 ||
      sOnTau == NULL || Rf_isNull(sOnTau) || Rf_length(sOnTau) <= 0) {
      AA = p;  
    } 
    Rf_protect(sOut = Rf_allocVector(INTSXP,AA));   
    for (ii = 0; ii < AA; ii++) { INTEGER(sOut)[ii] = OrderAttack[ii]; }
    Rf_unprotect(1);
    return(sOut);   
  }   
  SEXP get_WeightedXtX() {
    if (iWeightedXtX == NULL) {
      if (iiWeight != NULL) {
        Rf_error("WeigthedXtX, uh, oh iWeightedXtX is NULL!\n");
      }
      return(R_NilValue);
    }
    int ii;
    SEXP sOut;
    if (OnKappaS == 0) {
      Rprintf("WeightedXtX: OnKappaS = %d, so there are no weighted members!\n"); 
      return(R_NilValue);
    }
    Rf_protect(sOut = Rf_allocMatrix(INTSXP, 2, OnKappaS));
    if (OnKappaS > 0) {
      for (ii=0; ii < OnKappaS; ii++) {
         INTEGER(sOut)[ii*2] = kXFinder[ii];
         INTEGER(sOut)[ii*2 + 1] = (int) iWeightedXtX[ii];
      }
    }
    Rf_unprotect(1);
    return(sOut);
  }
  
  ////////////////////////////////////////////////////
  //  EEProbSortWidth
  //    If doing EEProbability merging, and one wants
  //    only to draw EE draws within closeness of the draw,
  //    this is up down ratio width.
  void set_EEProbSortWidth(SEXP sEEProbSortWidth) {
    if (Rf_isInteger(sEEProbSortWidth)) {
      if (INTEGER(sEEProbSortWidth)[0] > 0) {
         EEProbSortWidth = INTEGER(sEEProbSortWidth)[0];
      }
    } else if (Rf_isReal(sEEProbSortWidth)) {
      if (REAL(sEEProbSortWidth)[0] > 0 ) {
       EEProbSortWidth = REAL(sEEProbSortWidth)[0];
      }
    }
  }
  SEXP get_EEProbSortWidth() {
    SEXP sOn = R_NilValue;
    Rf_protect( sOn = Rf_allocVector(REALSXP, 1));
    REAL(sOn)[0] = EEProbSortWidth;
    Rf_unprotect(1); return(sOn);
  }
  
  int SetupIntY(SEXP fsIntY, SEXP sdfRobit) {
  
    if (Rf_length(fsIntY) != Rf_length(sY)) {
      Rf_error("Hey: Y doesn't have IntY length!");
    }
    if (Rf_isNull(sdfRobit)) {
      Rf_error("SetupIntY:  Error: sdfRobit supplied is NULL!");
    }
    if (Rf_isInteger(sdfRobit)) {
      dfRobit = (double) INTEGER(sdfRobit)[0];
      if (dfRobit < 0) {
        dfRobit = -1;
        Rf_error("SetupIntY, Error, sdfRobit supplied is < 0 is %d \n",
          INTEGER(sdfRobit)[0]);
      }
    } else if (Rf_isReal(sdfRobit)) {
      dfRobit = (double) REAL(sdfRobit)[0];
      if (dfRobit < 0) {
        dfRobit = -1;
        Rf_error("SetupIntY, Error, sdfRobit supplied is < 0 is %f \n",
          REAL(sdfRobit)[0]);
      }    
    } else {
      Rf_error("SetupIntY: I don't know what you supplied for dfRobit");
    }
    int ii;
    FFree( intY, "intY");
    intY = (int*) Calloc(n, int);
    if (Rf_isInteger(fsIntY)) {
      for (ii = 0; ii < Rf_length(sY); ii++) {
        intY[ii] = (int) INTEGER(fsIntY)[ii];
      } 
    } else {
     for (ii = 0; ii < Rf_length(sY); ii++) {
        intY[ii] = (int) REAL(fsIntY)[ii];
     }
    }
    if (dfRobit > 0.0) {
      if (iiWeight == NULL) {
        RMemGetD(iiWeight, "iiWeight", n);
        for (int iti = 0; iti < n; iti++) {
          iiWeight[iti] = 1.0;
        }
        if (OnKappaMem <= 0) {
          Rprintf("SetupIntY: UhOh, can't do if OnKappaMem = %d.\n", OnKappaMem);
          R_FlushConsole();
          Rf_error("SetupIntY, want positive OnKappaMem. \n");
        }
        if (iWeightedXtX == NULL) {
          RMemReallocS(iWeightedXtX, "iWeightedXtX", OnKappaMem); 
          for (int iti = 0; iti < OnKappaMem; iti++) {
            iWeightedXtX[iti] = 0;
          }
        }
      }
    }
    return(1);
  }
  SEXP get_dfRobit() {
    SEXP sOn = R_NilValue;
    Rf_protect(sOn = Rf_allocVector(REALSXP, 1));
    REAL(sOn)[0] = dfRobit;  Rf_unprotect(1); return(sOn); }
  void set_dfRobit(SEXP adfRobit) {
    if (adfRobit == NULL || Rf_isNull(adfRobit)) {
      dfRobit = -1;
    }
    if (!Rf_isReal(adfRobit) || !Rf_isInteger(adfRobit)) {
      dfRobit = -1;
    }
    if (Rf_isReal(adfRobit)) {
      dfRobit = REAL(adfRobit)[0];
    } else if (Rf_isInteger(adfRobit)) { dfRobit = INTEGER(adfRobit)[0]; }
    if (dfRobit >= 0.0 && iiWeight == NULL) {
      RMemGetD(iiWeight, "iiWeight", n);
      for (int iti = 0.0; iti < n; iti++) {
        iiWeight[iti] = 1.0;
      }
    }
  }
  int SetupSumiiWeight(SEXP iStartRecordingiiWeight) {
    CountiiWeight = 0;  
    if (Rf_isNull(iStartRecordingiiWeight) || !(Rf_isReal(iStartRecordingiiWeight) ||
      Rf_isInteger(iStartRecordingiiWeight))) {
      Rprintf("SetupSumiiWeight: iStartRecordingiiWeight is not good. \n");
      R_FlushConsole(); return(-1);
    }
    int AG = GetFirstInteger(iStartRecordingiiWeight);
    if (AG < 0) {
      StartRecordingiiWeight = -1;  CountiiWeight = 0;
      if (SumiiWeight != NULL) { FFree(SumiiWeight, "SumiiWeight"); }
      return(0);
    }
    StartRecordingiiWeight = AG;
    CountiiWeight = 0;
    if (SumiiWeight != NULL) { FFree(SumiiWeight, "SumiiWeight"); }
    RMemGetD(SumiiWeight, "SumiiWeight", n);
    for (int iti = 0; iti < n; iti++) { SumiiWeight[iti] = 0.0; }
    return(StartRecordingiiWeight);

  }
  int get_StartRecordingiiWeight() {
    return(StartRecordingiiWeight);
  }
  int get_CountiiWeight() { return(CountiiWeight); }
  SEXP GetAverageiiWeight() {
    if (StartRecordingiiWeight < 0) {
      Rprintf("GetAverageiiWeight: Sorry, StartRecordingiiWeight = %d",
        StartRecordingiiWeight); R_FlushConsole();  return(R_NilValue);
    }
    if (SumiiWeight == NULL) {
      Rprintf("GetAverageiiWeight: Sorry, SumiiWeight is not set up.\n");
      R_FlushConsole(); return(R_NilValue);
    }
    if (CountiiWeight <= 0) {
      Rprintf("GetAverageiiWeight: Sorry, CountiiWeight = %d\n", CountiiWeight);
      R_FlushConsole(); return(R_NilValue);
    }
    SEXP sOn = R_NilValue; 
    Rf_protect(sOn = Rf_allocVector(REALSXP, n));
    for (int iti = 0; iti < n; iti++) {
      REAL(sOn)[iti] = SumiiWeight[iti] / CountiiWeight;
    }
    Rf_unprotect(1);
    return(sOn);
  }
  SEXP get_intY() {
    SEXP sOn = R_NilValue;
    if (intY == NULL) { return(R_NilValue); }
    Rf_protect(sOn = Rf_allocVector(INTSXP, n));
    for (int iti = 0; iti < n; iti ++) {
      INTEGER(sOn)[iti] = intY[iti];
    }
    Rf_unprotect(1); return(sOn);
  }
  SEXP get_ProbTau() {
   if (ProbTau == NULL) { return(R_NilValue); }
    SEXP sOut;
    Rf_protect(sOut = Rf_allocVector(REALSXP, Rf_length(tauEndList)));
    int One = 1;   int LL = Rf_length(tauEndList);
    F77_CALL(dcopy)(&LL, ProbTau, &One, REAL(sOut), &One);
    Rf_unprotect(1);
    return(sOut);   
  }
  SEXP get_ProbFixed() {
   if (ProbFixed == NULL) { return(R_NilValue); }
    SEXP sOut;     int One = 1; 
    if (sOnTau == NULL || Rf_isNull(sOnTau) || Rf_length(sOnTau) <= 0 ||
      tauEndList == NULL || Rf_isNull(tauEndList) || Rf_length(tauEndList) <= 0) {
      if (Rf_isNull(tauEndList) && !Rf_isNull(sOnTau)) {
        Rf_error("BayesSpikeCpp.h::get_ProbFixed(), why is tauEndList NULL but not sOnTau!\n");
      }
      if (Rf_isNull(sOnTau) && !Rf_isNull(tauEndList)) {
        Rf_error("BayesSpikeCpp.h::get_ProbFixed(), why is sOnTau NULL but not tauEndList!\n");
      }
      Rf_protect(sOut = Rf_allocVector(REALSXP, p));
      if (Verbose >= -1) {
        Rprintf("BayesSpikeCpp.h():get_ProbFixed: Made the manual fill.\n"); R_FlushConsole();
      }
      for (int jj = 0; jj < p; jj++) {
        REAL(sOut)[jj] = ProbFixed[jj];
      }
      if (Verbose >= -1) {
        Rprintf("BayesSpikeCpp.h():get_ProbFixed: Made the manual fill.\n"); R_FlushConsole();
      }
      REAL(sOut)[p-1] = REAL(sOut)[p-1]+1.0-1.0;  
      // F77_CALL(dcopy)(&p, ProbFixed, &One, REAL(sOut), &One);
    } else if (iFirstRandom == 0) {
      return(R_NilValue);
    } else if (iFirstRandom > p) {
      Rf_error("BayesSpikeCpp.h():get_ProbFixed: why is p=%d but iFirstRandom = %d!?\n",
        p, iFirstRandom);
    }  else {
      Rf_protect(sOut = Rf_allocVector(REALSXP, iFirstRandom));
      F77_CALL(dcopy)(&iFirstRandom, ProbFixed, &One, REAL(sOut), &One);
    } 
    Rf_unprotect(1);
    return(sOut);   
  }
  
  
  double TestPt(SEXP sP, SEXP sDf, SEXP sLTF, SEXP sLog);
  double TestPnorm(SEXP sP, SEXP sMu, SEXP sSig, SEXP sLTF, SEXP sLog);
  ///////////////////////////////////
  //  Dependencies Fixed.
  //   Unimplemented, but this would a linked list of dependencies that govern
  //   What fixed effects/random effects can be turned on in order.
  SEXP get_DependenciesFixed() { if (RDependenciesFixed==NULL) { return(R_NilValue); }
    return(RDependenciesFixed->asSexp()); }
  SEXP get_DependenciesTau() { if (RDependenciesTau==NULL) { return(R_NilValue); }
    return(RDependenciesTau->asSexp()); } 
  
  void set_DependenciesFixed(SEXP DependenciesFixed_) {
    if (DependenciesFixed_ == NULL || Rf_isNull(DependenciesFixed_)) {
      DependenciesFixed = R_NilValue;
    }
    RDependenciesFixed = new AObject(DependenciesFixed_);
    DependenciesFixed = RDependenciesFixed->asSexp();
  }
  void set_DependenciesTau(SEXP DependenciesTau_) {
    if (DependenciesTau_ == NULL || Rf_isNull(DependenciesTau_)) {
      DependenciesTau = R_NilValue;  return;
    }
    if (Rf_length(DependenciesTau_)  != Rf_length(tauEndList)) {
      Rf_error("DependenciesTau: %d, Not same length as tauEndList = %d",
        Rf_length(DependenciesTau_), Rf_length(tauEndList));
    }
    RDependenciesTau = new AObject(DependenciesTau_);
    DependenciesTau = RDependenciesTau->asSexp();
  }
  void set_sCodaILocAndProbFile(SEXP _sCodaILocFile, SEXP _sCodaProbFile,
    SEXP _sIters) {
   if (!Rf_isString(_sCodaILocFile)) {
      Rf_error("set_sCodaILocFile: Set to a string please!");
    }
    if (!Rf_isString(_sCodaProbFile)) {
      Rf_error("set_sCodaProbFile: Set to a string please!");
    }
    int iIters = GetFirstInteger(_sIters);
    
    if (iIters <= 0) {
      Rf_error("set_sCodaProbFile: can't set iIters to %d \n", iIters);
    }
    DDelete(RsCodaProbFile, "RsCodaProbFile");  
    RsCodaProbFile = new RRObject(_sCodaProbFile);
    DDelete(RsCodaILocFile, "RsCodaILocFile");  
    RsCodaILocFile = new RRObject(_sCodaILocFile);
    if (ProbCoda == NULL || ICoda == NULL) {
      SetupProbCoda(iIters);
    }
    iiOnILocProb = 0;  
  }
  int LoadProbOldCodaICoda(int DoSort, int StartIter, int LastIter);
  void set_sCodaOldLocAndProbFile(SEXP _sCodaOldILocFile, SEXP _sCodaOldProbFile,
   SEXP StartLoad, SEXP ReadLength, SEXP iOldIterCoda, SEXP DoSort) {
    if (Rf_isNull(_sCodaOldILocFile) || !Rf_isString(_sCodaOldILocFile)) {
      Rf_error("set_sCodaOldILocFile: Set to a string please!");
    }
    if (Rf_isNull(_sCodaOldProbFile) || !Rf_isString(_sCodaOldProbFile)) {
      Rf_error("set_sCodaOldProbFile: Set to a string please!");
    }
    DDelete(RsOldCodaProbFile, "RsOldCodaProbFile");  
    RsOldCodaProbFile = new RRObject(_sCodaOldProbFile);
    DDelete(RsOldCodaILocFile, "RsOldCodaILocFile");  
    RsOldCodaILocFile = new RRObject(_sCodaOldILocFile);
    IterOldCoda = GetFirstInteger(iOldIterCoda);
    LoadProbOldCodaICoda( GetFirstInteger(DoSort), 
      GetFirstInteger(StartLoad), GetFirstInteger(ReadLength));
  }
  int WritePICoda();
  
  SEXP get_CodaILocFile() {
    if (RsCodaILocFile == NULL) { return(R_NilValue); }
    return(RsCodaILocFile->asSexp());
  }
  SEXP get_CodaProbFile() {
   if (RsCodaProbFile == NULL) { return(R_NilValue); }
    return(RsCodaProbFile->asSexp());
  }
  void set_CodaProbFile(SEXP RsCodaProbFile_) {
    if (RsCodaProbFile_ == NULL || Rf_isNull(RsCodaProbFile_) ||
      !Rf_isString(RsCodaProbFile_)) {
      DDelete(RsCodaProbFile, "RsCodaProbFile");  
    }
    RsCodaProbFile = new RRObject(RsCodaProbFile_);
    return;
  }
  SEXP get_OldCodaILocFile() {
    if (RsOldCodaILocFile == NULL) { return(R_NilValue); }
    return(RsOldCodaILocFile->asSexp());
  }
  SEXP get_OldCodaProbFile() {
    if (RsOldCodaProbFile == NULL) { return(R_NilValue); }
    return(RsOldCodaProbFile->asSexp());
  }
  void set_sCodaPFile(SEXP _sCodaPFile) {
    if (Rf_isNull(_sCodaPFile)) {
      if (Verbose >= 2) {
        Rprintf("Setting CodaPFile to NULL."); R_FlushConsole();
      }
      DDelete(RsCodaPFile, "RsCodaPFile");  
      RsCodaPFile = NULL;  sCodaPFile = R_NilValue;  NewPWrite = 1;
      //RsCodaPFile = new AObject(R_NilValue);
      //sCodaPFile = RsCodaPFile->asSexp();  NewPWrite = 1;
    }
    if (!Rf_isString(_sCodaPFile)) {
      Rf_error("set_sCodaPFile: Set to a string please!");
    }
    if (DoRecord == NULL || Rf_isNull(DoRecord)) {
      Rprintf("set_CodaPFile, weird, DoRecord is NULL!\n");
      R_FlushConsole();
    } else if (Rf_isNull(sOnTau) || Rf_length(sOnTau) <= 0 ||
      Rf_isNull(tauEndList) || Rf_length(tauEndList) <= 0) {
      if (Verbose >= 2) {
        Rprintf("set_sCodaPFile, set CodaPFile to non null, though No Tau is present.");
        R_FlushConsole();
      }
    } else if (!Rf_isNull(DoRecord) && 
      Rf_length(DoRecord) >= 4 && ANINT(DoRecord, 2) == 0 && ANINT(DoRecord,3) == 0) {
      Rf_error("set_sCodaPFile: Do not set sCodaPFile, unless recording in CodaTakes place");
      R_FlushConsole();
    }
    DDelete(RsCodaPFile, "RsCodaPFile");  
    RsCodaPFile = new RRObject(_sCodaPFile);
    sCodaPFile = RsCodaPFile->asSexp();  NewPWrite = 1;
  }
  SEXP get_sCodaPFile() {
    return(RsCodaPFile->asSexp());
  }
  int set_sCodaTFile(SEXP _sCodaiTFile, SEXP _sCodadTFile, SEXP _sCodaiTLocFile) {
    if (sOnTau == NULL || Rf_isNull(sOnTau) || Rf_length(sOnTau) <= 0 ||
      tauEndList == NULL || Rf_isNull(tauEndList) || Rf_length(tauEndList) <= 0) {
      if (CodaTisSetup == 1) {
        Rprintf("set_sCodaTFile: Note, CodaTisSetup!\n"); R_FlushConsole();
      }
      CodaTisSetup = 0;
      if (Rf_isNull(_sCodaiTFile) || _sCodaiTFile == NULL ||
        !Rf_isString(_sCodaiTFile)) {
        DDelete(RsCodaiTFile, "RsCodaiTFile");  
        DDelete(RsCodadTFile, "RsCodadTFile");
        DDelete(RsCodaiTLocFile, "RsCodaiTLocFile"); 
        sCodaiTFile = R_NilValue;   
        sCodaiTFile = R_NilValue; 
        return(1);
      }
      //char *MyS =  ;
      char Other[10] = "";
      if (strcmp(CHAR(STRING_ELT(_sCodaiTFile, 0)), Other) == 0) {
        DDelete(RsCodaiTFile, "RsCodaiTFile");  
        DDelete(RsCodadTFile, "RsCodadTFile");
        DDelete(RsCodaiTLocFile, "RsCodaiTLocFile"); 
        sCodaiTFile = R_NilValue;   
        sCodaiTFile = R_NilValue; 
        return(1);
      } 
      Rf_error("set_sCodaTFile: Setup sOnTau first Please!\n");   
    }
    if (_sCodaiTFile == NULL || Rf_isNull(_sCodaiTFile)) {
      if (Verbose >= 2) {
        Rprintf("set_CodaTFile, we're setting file names to NULL. \n"); R_FlushConsole();
      }
      DDelete(RsCodaiTFile, "RsCodaiTFile"); 
      RsCodaiTFile = NULL;  //new AObject(R_NilValue);
        //sCodaiTFile = RsCodaiTFile->asSexp();
      sCodaiTFile = R_NilValue; 
      DDelete(RsCodadTFile, "RsCodadTFile");   
      RsCodadTFile = NULL; //RsCodadTFile = new AObject(R_NilValue);
        //sCodadTFile = RsCodadTFile->asSexp(); 
      sCodadTFile  = R_NilValue;
      DDelete(RsCodaiTLocFile, "RsCodaiTLocFile");
      
      if (Verbose >= 2) {
        Rprintf("set_CodaTFile, successfully set file names to NULL. \n"); R_FlushConsole();
      }
      return(0);            
    }
    if (Rf_isNull(_sCodaiTFile) || !Rf_isString(_sCodaiTFile) ||
      Rf_isNull(_sCodadTFile) || !Rf_isString(_sCodadTFile) ) {
      Rprintf("set_sCodaTFile: Set to a string please!\n");
      return(-2);
    }
    //if (!Rf_isNull(DoRecord) && Rf_length(DoRecord) >= 2 && ANINT(DoRecord, 1) == 0) {
    //  Rprintf("set_sCodaTFile: Do not set sCodaTFile, unless recording in CodaTakes place");
    //  return(-1);
    //}
    if (!Rf_isNull(_sCodaiTLocFile) && Rf_isString(_sCodaiTLocFile)) {
      DDelete(RsCodaiTLocFile, "RsCodaiTLoc");
      RsCodaiTLocFile = new RRObject(_sCodaiTLocFile);
      FFree(CodaTLocjj, "CodaTLocjj");
      if (LengthProbCoda > 0 && LengthProbCoda <= (p > MaxSquaresAllocation ? p : MaxSquaresAllocation)) {
        LengthCodaTLoc = LengthProbCoda;
      } else {
        LengthCodaTLoc = 400;
      }
      OnCodaTLocjj = 0;  TotCodaTLocjj = 0;
      CodaTLocjj = (bufferILocType*) Calloc(LengthCodaTLoc, bufferILocType);
      NewWriteTLoc = 1;
    }
    
    LengthCodaTjj = (int) ceil(((double)MaxCODABeta * (double) Rf_length(sOnTau))/ (double)p); 
    if (LengthCodaTjj < Rf_length(sOnTau) + 50) {
       LengthCodaTjj = Rf_length(sOnTau) + 50;
    }


    OnCodaT = 0; TotOnCodaT = 0;
    DDelete(RsCodaiTFile, "RsCodaiTFile"); 
      RsCodaiTFile = new RRObject(_sCodaiTFile);
    sCodaiTFile = RsCodaiTFile->asSexp(); 
    DDelete(RsCodadTFile, "RsCodadTFile");   
    RsCodadTFile = new RRObject(_sCodadTFile);
    sCodadTFile = RsCodadTFile->asSexp();
    if (CodaTau != NULL) {
      if (Verbose >= 1) {
      Rprintf("-- BayesSpikeCpp.h:: setSCodaTFile() Before allocation, CodaTau is not null, deleting now, tt = %d!\n", tt);
      R_FlushConsole();
      }
      FFree(CodaTau, "CodaTau");
      if (Verbose >= 1) {
        Rprintf("-- BayesSpikeCpp.h:: setSCodaTFile() Successfully Eliminated CodaTau, will allocate length %d.\n", LengthCodaTjj);
        R_FlushConsole();
      }
      LengthWrittenTauCodaBuffer = 0; LengthTotalWrittenTauCodaBuffer = 0;
    } 
    if (Verbose >= 1) {
        Rprintf("-- BayesSpikeCpp.h:: setSCodaTFile() About to Allocate CodaTau, will allocate length %d.\n", LengthCodaTjj);
        R_FlushConsole();
    }
    RMemGetD(CodaTau, "CodaTau", LengthCodaTjj);
    if (Verbose >= 1) {
        Rprintf("-- setSCodaTFile() Successfully allocated LengthCodaTjj =%d for CodaTau!\n",
          LengthCodaTjj);
        R_FlushConsole();
    }
    LengthWrittenTauCodaBuffer = 0; LengthTotalWrittenTauCodaBuffer = 0;
    RMemGetI(CodaTjj, "CodaTjj", LengthCodaTjj);

    if (Verbose >= 1) {
        Rprintf("-- setSCodaTFile() Finished allocating LengthCodaTjj =%d for CodaTau!\n",
          LengthCodaTjj);
        R_FlushConsole();
    }               
    NewTWrite = 1; NewWriteTLoc = 1;
    return(1);
  }
  int get_OnCodaT() {return(OnCodaT);}
  int get_TotCodaT()  { return(TotOnCodaT);}
  int get_TotCodaTLocjj() {return((int) TotCodaTLocjj);}
  int get_OnCodaTLocjj(){ return((int) OnCodaTLocjj);}
  int get_NewWriteTLoc() { return(NewWriteTLoc);}
  int get_NewWriteT() {return(NewTWrite);}
  int get_LengthCodaTjj() { return(LengthCodaTjj); }
  int get_LengthCodaTLoc() { return(LengthCodaTLoc); }
  int get_MaxCODABeta() { return(MaxCODABeta); }
  void set_NewWriteT(int iT) {
    if (iT == 0) { NewTWrite = 0;  NewWriteTLoc = 0;}
    if (iT == 1) { NewTWrite =1;  NewWriteTLoc = 1;}    
  }
  SEXP get_RsCodaiTLocFile() {
    if (RsCodaiTLocFile == NULL) {
      return(R_NilValue);
    }
    return(RsCodaiTLocFile->asSexp());
  }
  void set_SaveDir(SEXP _sSaveDir) {
    if (Verbose >= -1) {
      Rprintf("BayesSpikeCpp.h:set_SaveDir(): Starting. %s\n",
        RsSaveDir==NULL? "RsSaveDir is NULL":"RsSaveDir is not NULL");
    }
    if (Rf_isNull(_sSaveDir) ) {
      if (Verbose > 1) {
        Rprintf("sSaveDir: Setting Savedir to NULL!");  R_FlushConsole();
        DDelete(RsSaveDir, "RsSaveDir"); RsSaveDir = NULL;
        sSaveDir = R_NilValue;
        return;
      }
    }
    if (!Rf_isString(_sSaveDir) ) {
      Rprintf("sSaveDir: Please Set to a string please!\n"); R_FlushConsole();
      return;
    }
    DDelete(RsSaveDir, "RsSaveDir"); sSaveDir = R_NilValue;
      RsSaveDir = new AObject(_sSaveDir);
    if (RsSaveDir == NULL) {
      Rprintf("set_SaveDir: uh oh, RsSaveDir new returned NULL!\n"); R_FlushConsole();
    }  else {
      sSaveDir = RsSaveDir->asSexp();
    }
    return;
  }
  SEXP get_SaveDir() { 
    if (RsSaveDir == NULL) { return(R_NilValue); }
    return(RsSaveDir->asSexp()); 
  }
  SEXP get_sCodaiTFile() {
    if (RsCodaiTFile == NULL) {return(R_NilValue); }
    return(RsCodaiTFile->asSexp());
  }
  SEXP get_sCodadTFile() {
    if (RsCodadTFile == NULL) { return(R_NilValue); }
    return(RsCodadTFile->asSexp());
  }
  SEXP get_sCodaIFile() { 
    if (RsCodaIFile == NULL) { return(R_NilValue); }
    return(RsCodaIFile->asSexp()); }
  SEXP get_sCodaJFile() { 
    if (RsCodaJFile == NULL) { return(R_NilValue); }
    return(sCodaJFile); }
  SEXP get_sCodaOldIFile() { if (RsCodaOldIFile == NULL) {
    return(R_NilValue); } else { return(RsCodaOldIFile->asSexp());} }
  SEXP get_sCodaOldJFile() { if (RsCodaOldJFile == NULL) {
    return(R_NilValue); } else { return(RsCodaOldJFile->asSexp()); } }
  void set_sCodaOldIFile( SEXP isCodaOldIFile) {
    DDelete(RsCodaOldIFile, "RsCodaOldIFile");
    if (Rf_isNull(isCodaOldIFile)) { return; }
    if (!Rf_isString(isCodaOldIFile)) {
      Rprintf("set _sCodaOldIFile, supply string");  R_FlushConsole();
    }
    RsCodaOldIFile = new RRObject(isCodaOldIFile);
  }
  void set_sCodaOldJFile( SEXP isCodaOldJFile) {
    DDelete(RsCodaOldJFile, "RsCodaOldJFile");
    if (Rf_isNull(isCodaOldJFile)) { return; }
    if (!Rf_isString(isCodaOldJFile)) {
      Rprintf("set _sCodaOldJFile, supply string\n"); R_FlushConsole();
    }
    RsCodaOldJFile = new RRObject(isCodaOldJFile);
  }
  SEXP get_AllEigenValues() { if (RAllEigenValues == NULL) { return(R_NilValue); }
    return(RAllEigenValues->asSexp()); }
  SEXP get_AllEigenVectors() { 
    if (RAllEigenVectors == NULL) {
    return(R_NilValue); }
    return(RAllEigenVectors->asSexp()); 
  }
  void set_AllEigenValues(SEXP AllEigenValues_) {
    if (Rf_isNull(tauEndList)) {
      Rf_error("Before Setup AllEigenValues: give tauEndList!");
    }
    if (Rf_isNull(AllEigenValues_)) {
      Rprintf("Hey: AllEigenValues input is null too!"); R_FlushConsole();
    }
    if (Rf_length(AllEigenValues_) != Rf_length(tauEndList)) {
      Rprintf("Error: AllEigenValues_ has length %d, but tauEndList length %d\n",
        Rf_length(AllEigenValues_), Rf_length(tauEndList)); R_FlushConsole();
      Rf_error("Setup AllEigenValues: tauEndList of different length!");
    }
    DDelete(RAllEigenValues, "RAllEigenValues");
    RAllEigenValues = new AObject(AllEigenValues_);
    AllEigenValues = RAllEigenValues->asSexp();}
  void set_AllEigenVectors(SEXP AllEigenVectors_) {
    if (Rf_isNull(tauEndList)) {
      Rf_error("Before Setup AllEigenVectors: give tauEndList!");
    }
    if (Rf_isNull(AllEigenVectors_)) {
      Rf_error("set_AllEigenVectors: you have me null EigenVectors."); 
    }
    if (Rf_length(AllEigenVectors_) != Rf_length(tauEndList)) {
      Rf_error("Setup AllEigenVector: tauEndList of different length!");
    }
    DDelete(RAllEigenVectors, "RAllEigenVectors");
    RAllEigenVectors = new AObject(AllEigenVectors_);
    AllEigenVectors = RAllEigenVectors->asSexp();}  
  void set_Verbose(int Verbose_) { if (Verbose_ < 0) {
    Verbose = 0;} else {
    Verbose = Verbose_;}
    int ProtNum = 0;
    TBSR5 = get_TBSR5();
    if (Verbose > 0  && !Rf_isNull(TBSR5)) {
     //Rprintf("Setting Verbose on TBSRoo as well \n"); R_FlushConsole();
     SEXP sVerbose = R_NilValue;
     Rf_protect(sVerbose = Rf_allocVector(INTSXP, 1));  ProtNum++;
     INTEGER(sVerbose)[0] = Verbose + 1; 
     SEXP sTBSR5SetVerbose = 
       Rf_protect(Rf_findVarInFrame( BayesSpikeNameSpace, Rf_install("TBSR5SetVerbose")));
     ProtNum++;
     //if (Rf_isNull(sTBSRooSetVerbose)) {
     //  Rprintf("BayesSpikeCL: SetVerbose, Error: no sTBSRooSetVerbose");
     //}
     //Rprintf("  Evidently we found sTBSRooSetVerbose \n"); R_FlushConsole();
     SEXP call = Rf_protect(Rf_lcons( sTBSR5SetVerbose, 
       Rf_cons(TBSR5,
         sVerbose))); ProtNum++;
     //Rprintf("  We just created Calling function \n"); R_FlushConsole();
     if (Rf_isNull(call)) {
       Rf_error("BayesSpikeCL Destroy: Failed to delete TBSR5");
     } 
     SEXP RTM = R_NilValue;
     if (Rf_isNull(RTM)) {  
     }
     //Rprintf("  Evaluate Call \n"); R_FlushConsole();
     //Rf_protect(RTM = Rf_eval(call, BayesSpikeNameSpace));
     Rf_unprotect(ProtNum);
     //if (Verbose > 0) {
     //  Rprintf(" BayesSpikeCL:  Finish Verbose set "); R_FlushConsole();
     //} 
    }
  }
  int get_Verbose() { return(Verbose);}
 
  void set_OnPiA(SEXP sOnPiA_) {
    if (Rf_isNull(sOnPiA_) || Rf_length(sOnPiA_) < 0 || Rf_length(sOnPiA_) > 2) {
      Rf_error("set_OnPiA:  Please Give Realistic Length");
    }
    if (Rf_isNull(sOnPiA)) {
      DDelete(RsOnPiA, "RsOnPiA");  RsOnPiA = new AObject(sOnPiA_);
      sOnPiA = RsOnPiA->asSexp();
    } else {
      if (Rf_length(sOnPiA) != Rf_length(sOnPiA_)) {
        Rf_error("set_OnPiA, please give PiA of former length %d \n", Rf_length(sOnPiA));
      }
      for (int ii = 0; ii < Rf_length(sOnPiA); ii++) {
        REAL(sOnPiA)[ii] = REAL(sOnPiA_)[ii];
      }
    }
  }
  void set_PiAPrior(SEXP PiAPrior_) {
    if (PiAPrior_ == NULL || Rf_isNull(PiAPrior_) || Rf_length(PiAPrior_) <= 1) {
      DDelete(RPiAPrior, "RPiAPrior");  
      RPiAPrior = NULL;
      PiAPrior = R_NilValue;
      return;      
      //Rf_error("SetPiAPrior: can't set null prior!");
    }
    if (sOnPiA == NULL || Rf_isNull(sOnPiA) || Rf_length(sOnPiA) <= 0) {
      Rf_error("SetPiAPrior: please set OnPiA first!");
    }
    if (Rf_length(sOnPiA) == 2 && Rf_length(PiAPrior_) == 4) {
      DDelete(RPiAPrior, "RPiAPrior");  
      RPiAPrior = new AObject(PiAPrior_);
      PiAPrior = RPiAPrior->asSexp();
      return;
    }
    if (Rf_length(sOnPiA) == 1 && Rf_length(PiAPrior_) > 2) {
      Rf_error("SetPiAPrior: sOnPiA only has length 1 give valid prior!");
    }
    if (Rf_length(PiAPrior_) <= 1 || Rf_length(PiAPrior_) == 3 ||
      Rf_length(PiAPrior_) >4) {
      Rf_error("SetPiAPrior: Realistic Lengths of PiAPrior Only!");
    }
    if (!Rf_isReal(PiAPrior_)) {
      Rf_error("SetPiAPrior: Give Real Input!");
    }
    if (Rf_isNull(PiAPrior)) {
      DDelete(RPiAPrior, "RPiAPrior");  RPiAPrior = new AObject(PiAPrior_);
      PiAPrior = RPiAPrior->asSexp();
    } else {
      DDelete(RPiAPrior, "RPiAPrior");  RPiAPrior = new AObject(PiAPrior_);
      PiAPrior = RPiAPrior->asSexp();
    }
  }
  int CheckTBSRoo(SEXP PropTBSRoo) {
    SEXP rTBSRoo = get_TBSRoo();
    if ( (void*) rTBSRoo != (void*) PropTBSRoo) {
      Rprintf("CheckTBSRoo gives two pointer errors, NInstance = %d, Uh roh! \n", NInstance); R_FlushConsole();
      Rprintf(" Those two pointers = Connected TBSRoo = %d, Prop TBSRoo = %d \n", 
        rTBSRoo, PropTBSRoo);
      R_FlushConsole();
      return(-1);
    }
    return(1);
  }
  int CheckTBSR5(SEXP PropTBSR5) {
    SEXP rTBSR5 = get_TBSR5();
    if ( (void*) rTBSR5 != (void*) PropTBSR5) {
      Rprintf("CheckTBSR5 gives two pointer errors, NInstance = %d, Uh roh! \n", NInstance); R_FlushConsole();
      Rprintf(" Those two pointers = Connected TBSR5 = %d, Prop TBSR5 = %d \n", 
        rTBSR5, PropTBSR5);
      R_FlushConsole();
      return(-1);
    }
    return(1);
  }
  int CheckX(SEXP _X) {  
    if ( (void*) sX != (void*) _X) {
      Rprintf("CheckX gives error, the _X=%d you tested against is not=%d \n", 
         _X,  sX); R_FlushConsole();  return(-1);
    }
    return(1);
  }
  int CheckY(SEXP _Y) {  
    if ( (void*) sY != (void*) _Y) {
      Rprintf("CheckX gives error, the _Y=%d you tested against is not=%d \n", 
        _Y,  sY); R_FlushConsole();  return(-1);
    }
    return(1);
  }
  int CheckBeta(SEXP _Beta) {  
    if ( (void*) sBeta != (void*) _Beta) {
      Rprintf("CheckX gives error, the _X=%f you tested against is not=%f \n", 
        _Beta, sBeta); R_FlushConsole();  return(-1);
    }
    return(1);
  }  
  SEXP get_OnPiA() { if (RsOnPiA == NULL) { return(NULL) ; } 
    return(RsOnPiA->asSexp());}  SEXP get_PiAPrior() {
    if (RPiAPrior == NULL) { return(R_NilValue); }
    return(PiAPrior);}
  SEXP get_CurrentProb() {
    SEXP sOn = R_NilValue;
    Rf_protect(sOn = Rf_allocVector(REALSXP, 1));
    REAL(sOn)[0] = CurrentProb; Rf_unprotect(1); return(sOn);
  }
  void TestSigma() {
    Rprintf("BayesSpikeCpp.h:TestSigma() called. \n"); R_FlushConsole();
    if (RsOnSigma == NULL) { Rprintf("RsOnSigma is NULL!\n"); return; }
    if (Rf_isNull(RsOnSigma->asSexp())) { Rprintf("RsOnSigma->OnSigma is R_NilValue\n"); return; }
    if (Rf_length(RsOnSigma->asSexp()) <= 0) { Rprintf("RsOnSigma->OnSigma is Zero Length. \n"); return;}
    if (Rf_length(RsOnSigma->asSexp()) >= 2) { Rprintf("RsOnSigma->OnSigma has %d large length?\n",
      Rf_length(RsOnSigma->asSexp()));}
    if (!Rf_isReal(RsOnSigma->asSexp())) { Rf_error("RsOnSigma is not REAL??\n");}
    Rprintf("TestSigma:"); R_FlushConsole();
    for (int ii = 0; ii < Rf_length(RsOnSigma->asSexp()); ii++) {
      Rprintf(" (%d, ", ii); R_FlushConsole();
      REAL(RsOnSigma->asSexp())[ii] = REAL(RsOnSigma->asSexp())[ii] + 1.0;
      REAL(RsOnSigma->asSexp())[ii] = REAL(RsOnSigma->asSexp())[ii] - 1.0;
      Rprintf("PASS)"); R_FlushConsole();
    }
    Rprintf("\n");
    Rprintf("TestSigma: We have manipulated all elements of Sigma. \n");   R_FlushConsole();
    Rprintf("TestSigma: Testing the SigmaPrior. \n"); R_FlushConsole();
    if (RSigmaPrior == NULL) { Rprintf("RSigmaPrior is NULL!\n"); return;}
    if (Rf_isNull(RSigmaPrior->asSexp())) { Rprintf("RSigmaPrior has RNilValue object\n"); return;}
    if (!Rf_isReal(RSigmaPrior->asSexp())) { Rprintf("RSigmaPrior has no real memeber.\n"); return;}
    if (Rf_length(RSigmaPrior->asSexp()) <= 1) { Rprintf("RSigmaPrior only has %d Length!\n",
      Rf_length(RSigmaPrior->asSexp()));  return; }
    Rprintf("Test Sigma Prior: "); R_FlushConsole(); 
    for (int ii = 0; ii < Rf_length(RSigmaPrior->asSexp()); ii++) {
      Rprintf(" (%d, ", ii); R_FlushConsole();
      REAL(RSigmaPrior->asSexp())[ii] = REAL(RSigmaPrior->asSexp())[ii] + 1.0;
      REAL(RSigmaPrior->asSexp())[ii] = REAL(RSigmaPrior->asSexp())[ii] - 1.0;
      Rprintf("PASS)"); R_FlushConsole();
    }
    Rprintf("\n");
  }
  void set_OnSigma(SEXP sOnSigma_) {
    if (Verbose > 2) {
      Rprintf("setOnSigma, assigning OnSigma \n"); R_FlushConsole();
      if (Rf_isNull(sOnSigma_)) {
        Rprintf("  Except, sOnSigma_ is NULL \n"); R_FlushConsole();
      }  else if (!Rf_isReal(sOnSigma_)) {
        Rprintf("  Except, sOnSigma_ is not REAL\n"); R_FlushConsole();
      } else {
        Rprintf("sOnSigma = %f, and is length %d \n", Rf_length(sOnSigma_),
          REAL(sOnSigma_));
      }
    }
    if (Rf_isNull(sOnSigma_) || Rf_length(sOnSigma_) < 0 || Rf_length(sOnSigma_) > 2) {
      Rf_error("set_OnSigma:  Please Give Realistic Length");
    }
    if (!Rf_isReal(sOnSigma_)) { 
      Rf_error("set_OnSigma: Give Numeric OnSigma!\n");
    }
    if (RsOnSigma != NULL && Rf_isReal(RsOnSigma->asSexp()) && Rf_length(RsOnSigma->asSexp()) >= 1) {
      REAL(RsOnSigma->asSexp())[0] = REAL(sOnSigma_)[0];
      return;
    }
    if (RsOnSigma != NULL && Verbose >= 3) {
      Rprintf("set_OnSigma, strange, sOnSigma is not null, overwriting.\n") ; R_FlushConsole();
    }
    DDelete(RsOnSigma, "RsOnSigma");  RsOnSigma = new AObject(sOnSigma_);
    sOnSigma = RsOnSigma->asSexp();
  }
  void set_SigmaPrior(SEXP SigmaPrior_) {
    if (Verbose >= -1) {                                                                               
      Rprintf("set_SigmaPrior: Starting, Check OnSigma \n"); R_FlushConsole();
    }
    if (RsOnSigma == NULL || sOnSigma == NULL || Rf_isNull(sOnSigma) || Rf_length(sOnSigma) <= 0) {
      Rf_error("SetSigma: please set OnSigma first!");
    }
    if (Verbose > 2) {
      Rprintf("set_SigmaPrior: Starting, Check for null prior \n"); R_FlushConsole();
    }
    if (SigmaPrior_ == NULL || Rf_isNull(SigmaPrior_) ) {
      //Rf_error("SetSigma: Give non null prior!!");
      DDelete(RSigmaPrior, "SigmaPrior");  
      RSigmaPrior= NULL;
      SigmaPrior = R_NilValue;
      return;
    }
    
    if (Verbose > 2) {
      Rprintf("set_SigmaPrior: Starting, Check Lengths \n"); R_FlushConsole();
    }
    if (Rf_length(SigmaPrior_) <= 1 || Rf_length(SigmaPrior_) > 2) {
      Rf_error("SetSigmaPrior: SigmaPrior_ only has length %d give valid prior!",
        Rf_length(SigmaPrior_));
    }
    if (Rf_length(SigmaPrior_) <= 1 || Rf_length(SigmaPrior_) == 3 ||
      Rf_length(SigmaPrior_) >4) {
      Rf_error("SetSigmaPrior: Realistic Lengths of SigmaPrior Only!");
    }
    if (!Rf_isReal(SigmaPrior_)) {
      Rf_error("SetSigmaPrior: Give Real numbers for prior!");
    }
    if (Verbose >= -1) {
      Rprintf("set_SigmaPrior: it seems that SigmaPrior_ is valid to set.\n"); R_FlushConsole();
    }
    if (!Rf_isNull(SigmaPrior) && RSigmaPrior == NULL) {
      Rprintf("set_SigmaPrior: it seems that SigmaPrior_ is invalid. \n"); R_FlushConsole();
    }
    if (RSigmaPrior != NULL) {
      DDelete(RSigmaPrior, "SigmaPrior"); SigmaPrior= R_NilValue;
      RSigmaPrior = new AObject(SigmaPrior_); SigmaPrior = RSigmaPrior->asSexp();
      if (Verbose >= -1) {
        Rprintf("Resetted SigmaPrior by deletion\n"); R_FlushConsole();
      }
      return;
    } else {
      RSigmaPrior = new AObject(SigmaPrior_); SigmaPrior = RSigmaPrior->asSexp();
      if (Verbose >= -1) {
        Rprintf("Resetted SigmaPrior by deletion\n"); R_FlushConsole();
      }
      return;    
    }
    if (Rf_isNull(SigmaPrior)) {
      if (RSigmaPrior != NULL) {
        Rprintf("SetSigmaPrior:  Huh, RSigmaPrior is not NULL \n"); R_FlushConsole();
      }
      DDelete(RSigmaPrior, "SigmaPrior");  
      RSigmaPrior= new AObject(SigmaPrior_); SigmaPrior = RSigmaPrior->asSexp();
    } else {
      if (Rf_length(SigmaPrior) != Rf_length(SigmaPrior_)) {
        Rf_error("set_SigmaPrior: Please supply sigmaprior of length %d", 
          Rf_length(SigmaPrior));
      }
      for (int ii = 0; ii < Rf_length(SigmaPrior); ii++) {
        REAL(SigmaPrior)[ii] = REAL(SigmaPrior_)[ii];
      }
    }
    if (Verbose >= -1) {
      Rprintf("set_SigmaPrior: Seems like we are concluding. \n"); R_FlushConsole();
    }
  }
  SEXP BAllLOfX(SEXP SampleStart);
  void set_dfTNoise(double dfTNoise_) { 
    if (Verbose >= 1) {
      Rprintf("  BayesSpikeCpp.h::set_dfTNoise(); \n"); R_FlushConsole();
    }
    if (dfTNoise_ > 0) { dfTNoise = dfTNoise_;}  else if
      (fabs(dfTNoise_ + 3.0) <= .0001) {
      dfTNoise = -3;
      return;
    } else {
      if (Verbose >= 1) {
        Rprintf("  BayesSpikeCpp.h::set_dfTNoise(), dfTNoise_=%f turning off\n", dfTNoise_); R_FlushConsole();
      }
      FFree(iiWeight, "iiWeight");  iiWeight = NULL;
      FFree(iWeightedXtX, "iWeightedXtX"); iWeightedXtX = NULL;
      dfTNoise = -1;
      return;
    }
    if (iiWeight == NULL) {
      RMemGetD(iiWeight, "iiWeight", n);
    }
    if (iWeightedXtX == NULL) {
      RMemGetS(iWeightedXtX, "iWeightedXtX", OnKappaMem+2); 
      for (int iti = 0; iti < OnKappaMem; iti++) {
        iWeightedXtX[iti] = (short int) 0;
      }
    }
    int ii = 0; for (ii = 0; ii < n; ii++) { iiWeight[ii] = 1.0;} 
    nWeight = n;  
    for (ii = 0; ii < OnKappaMem; ii++) {iWeightedXtX[ii] =(short int) 0;}
    if (Verbose >= 1) {
      Rprintf("  BayesSpikeCpp.h::set_dfTNoise(): Finished \n"); R_FlushConsole();
    }
  }
  SEXP get_iWeightedXtX() {
    if (iWeightedXtX == NULL) {
      return(R_NilValue);
    }
    if (OnKappaS == 0) {
      Rprintf("get_iWeightedXtX: There is no OnKappaS because it is zero!\n");
      return(R_NilValue);
    }
    SEXP sOn = R_NilValue;  Rf_protect(sOn = Rf_allocVector(INTSXP, OnKappaS));
    for (int iti = 0; iti < OnKappaS; iti++) {
      INTEGER(sOn)[iti] = iWeightedXtX[iti];
    }
    Rf_unprotect(1);
    return(sOn);
  }
  double get_dfTNoise() { return(dfTNoise); }
  
  SEXP get_RNGState() {
    SEXP sOn = R_NilValue;
    Rf_protect(sOn = Rf_allocVector(INTSXP, 1));
    INTEGER(sOn)[0] = OnRNGState;
    Rf_unprotect(1); return(sOn);
  }
  void set_RNGState(SEXP AnIn) { 
    if (OnRNGState == 1) {
      if (Rf_isLogical(AnIn) && (int) LOGICAL(AnIn)[0] == 0) {
        PutRNGstate();  OnRNGState = 0;
        return;
      } else if (Rf_isInteger(AnIn) && INTEGER(AnIn)[0] == 0) {
        PutRNGstate();  OnRNGState = 0;
      } else if (Rf_isReal(AnIn) && (int) REAL(AnIn)[0] == 0) {
        PutRNGstate();  OnRNGState = 0;
      }
      if (GetFirstInteger(AnIn) == 0) {
        PutRNGstate();  OnRNGState = 0;
        return;
      }
    }
    if (Rf_isNull(AnIn)) {
    
    } else if (Rf_isLogical(AnIn)) {
       if (Rf_isLogical(AnIn) && (int) LOGICAL(AnIn)[0] > 0) {
         GetRNGstate(); OnRNGState = 1;
       }
    }  else if (GetFirstInteger(AnIn) == 1) {
      GetRNGstate(); OnRNGState = 1; 
    } else {
      GetRNGstate(); OnRNGState = 1;
    }
  }
  void set_CodaList(SEXP CodaList_) {
    if (CodaList_ == NULL) {
      Rf_error("setcodaList, CodaList is true NULL pointer!\n");
    }
    if (Rf_isNull(CodaList_)) {
      Rf_error("setCodaList, don't supply NULL CodaList");
    }
    if (Verbose > 1) {
      Rprintf("set_CodaList: Starting \n");
    }
    if (RDoRecord == NULL ||
      DoRecord == NULL || Rf_isNull(DoRecord) || (!Rf_isReal(DoRecord) && !Rf_isInteger(DoRecord)) ||
      Rf_length(DoRecord) < 7) {
      Rprintf("setCodaList, set DoRecord First");
      if (DoRecord == NULL) {
        Rprintf("setCodaList: DoRecord is a NULL Pointer\n");
      } else if (Rf_isNull(DoRecord)) {
        Rprintf("setCodaList: DoRecord is a NULL SEXP \n");
      } else if (!Rf_isReal(DoRecord) && !Rf_isInteger(DoRecord)) {
        Rprintf("setCodaList: DoRecord is not Real or Integer\n"); R_FlushConsole();
      } else if (Rf_length(DoRecord) < 1) {
        Rprintf("setCodaList: DoRecord is length %d less than 1\n", Rf_length(DoRecord)); R_FlushConsole();
      }
      Rprintf("ERROR --------------------------------------------------------");
      Rprintf("set_CodaList: ERROR WAS INADEQUATE DORECORD\n"); R_FlushConsole();
      Rf_error("set_CodaList Error on DoRecord!\n"); 
    }
    SEXP NewDims;
    DDelete(RCodaList, "RCodaList");  
    CodaList = R_NilValue;
    
    if (Verbose > 3) {
      Rprintf("SetCodaList, let's look into NULL CodaList\n"); R_FlushConsole();
    }  
    if (Rf_isNull(CodaList)) {
      DDelete(RCodaList, "RCodaList");  
      RCodaList = new AObject(CodaList_);
      CodaList = RCodaList->asSexp();
      if (Rf_length(CodaList) <= 0) {
        Rf_error("Set RCodaList, CodaList has less than Zero Length!\n");
      }
      if (!Rf_isVector(CodaList)) {
        Rf_error("Set RCodaList, CodaList is not a Vector!\n");
      }
      if (Rf_length(CodaList) <= 0) {
        Rf_error("Set RCodaList, whatever it is, CodaList is length %d!\n",
          Rf_length(CodaList));
      }
      if (Rf_isReal(CodaList)) {
        Rprintf("SetRCodaList, what, CodaList is REAL?\n"); R_FlushConsole();
      }
      if (Verbose > 2) {
        if (Rf_isVector(CodaList)) {
          Rprintf("Set_CodaList, CodaList is definitely a vector, about to pull it.\n");
        }
      }
      SEXP MyS = VECTOR_ELT(CodaList,0);
      if (Verbose > 1) {
      if (Rf_isReal(MyS)) {
        Rprintf("SetCodaList, first element of CodaList is a REAL!\n");
        R_FlushConsole();
      }
      }
      
      if (Rf_length(CodaList) <= 0 ||  
        Rf_length(Rf_getAttrib(VECTOR_ELT(CodaList,0), R_DimSymbol)) != 2) {
        Rf_error("Set RCodaList, error, supplied CodaList_ does not have worthy dimensions");  
      }
      int OnCodaTable;  SEXP OnDimSymbol;
      for (OnCodaTable = 0; OnCodaTable < Rf_length(CodaList); OnCodaTable++) {
        if (Rf_isNull(VECTOR_ELT(CodaList, OnCodaTable))) {
          Rf_error("set_CodaList, one of the set tables is NULL!");
        }
        DDelete(RCodaTable, "RCodaTable"); 
        RCodaTable = new AObject(VECTOR_ELT(CodaList, OnCodaTable));
        CodaTable = RCodaTable->asSexp(); 
        if (Rf_isNull(CodaTable)) {
          Rf_error("set_CodaList, NULL CodaTable %d/%d !\n",
            OnCodaTable, Rf_length(CodaList));
        }
        if (!Rf_isReal(CodaTable)) {
          Rf_error("set_CodaList, CodaTable %d/%d is not REAL!\n",
            OnCodaTable, Rf_length(CodaList)); 
        }
        OnDimSymbol = Rf_getAttrib(CodaTable, R_DimSymbol);
        if (Rf_isNull(OnDimSymbol) || Rf_length(OnDimSymbol) != 2) {
          Rprintf("set_CodaList: Big Error, no Dim to CodaTable %d/%d!\n",
            OnCodaTable, Rf_length(CodaList)); R_FlushConsole();
          Rf_error("set_CodaList, one of the CodaTables[%d] has no dim!", OnCodaTable+1);
        }
        if (INTEGER(OnDimSymbol)[1] != CountRecord()) {
          Rprintf("************************************************************\n");
          Rprintf("** SetCodaList: Error.  \n");
          Rprintf("** CodaList: Error, we have based upon DoRecord settings that total CountRecord should be %d\n",
            CountRecord());
          Rf_error("** set_CodaList, CodaTables[%d] has not right dim %d != %d ",
            OnCodaTable+1, INTEGER(OnDimSymbol)[1], CountRecord());
        } 
      }
      DDelete(RCodaTable, "RCodaTable");
      RCodaTable = new AObject(VECTOR_ELT(CodaList, 0));
      OnCodaTable = 0;
      CodaTable = RCodaTable->asSexp();
      
    } else {
      if (Rf_length(CodaList) != Rf_length(CodaList_)) {
        Rf_error("set_CodaList:  Error, CodaList and CodaList_ are not same length \n");
      }

      for (int ii = 0; ii < Rf_length(CodaList_); ii++) {
        SEXP NewCodaTable = VECTOR_ELT(CodaList_, ii);
        if (Rf_isNull(NewCodaTable)) {
          Rf_error("set_CodaList, CodaTable element %d/%d is NULL \n",ii, Rf_length(CodaList_)); 
        }
        if (!Rf_isReal(NewCodaTable)) {
          Rf_error("set_CodaList, NewCodaTable %d/%d is not REAL!\n", ii, Rf_length(CodaList_));
        }
        NewDims = Rf_getAttrib(NewCodaTable, R_DimSymbol);
        if (INTEGER(NewDims)[1] != INTEGER(Rf_getAttrib(VECTOR_ELT(CodaList,ii),
          R_DimSymbol))[1]) {
          Rprintf("########################################################\n");
          Rprintf("## set_CodaList Error, previous chain %d has dims problem.\n",
            ii); R_FlushConsole();
          Rprintf("## Error in set_CodaList, previous %d chain had dims %d, %d \n",
            ii, INTEGER(Rf_getAttrib(VECTOR_ELT(CodaList,ii),R_DimSymbol))[0],
            INTEGER(Rf_getAttrib(VECTOR_ELT(CodaList,ii), R_DimSymbol))[1]);
          Rf_error("## ERROR: but you gave table that had dims, %d, %d",
            INTEGER(NewDims)[0], INTEGER(NewDims)[1]);
        }
        int CopyLen = Rf_length(CodaList_);
        if (CopyLen > Rf_length(CodaList)) { CopyLen = Rf_length(this->CodaList); }
        int One = 1;
        F77_CALL(dcopy)(&CopyLen, REAL(NewCodaTable), &One, 
          REAL(VECTOR_ELT(CodaList,ii)), &One);
      }
    }
    OnCodaTable = 0;
    DDelete(RCodaTable, "RCodaTable"); 
    RCodaTable = new AObject(VECTOR_ELT(CodaList, 0));
    CodaTable = RCodaTable->asSexp();
    SEXP dCodaTable =  Rf_getAttrib(CodaTable, R_DimSymbol); 
    if (!Rf_isInteger(dCodaTable)) {
      Rprintf("Hey, why is dCodaTable not an INTEGER? \n"); R_FlushConsole();
    }
    if (TemperatureList != NULL && LengthTempDraws != NULL) {
      if (Tempii >= 0 && Tempii < LengthTemperatureList) {
        LengthTempDraws[Tempii] = INTEGER(dCodaTable)[0];
      }
    }
  }

  SEXP get_CodaList() { return(CodaList);}
 // SEXP get_ChainIter() { return(ChainIter);}
 /* void set_ChainIter(int Set) {
    if (Rf_isNull(CodaList)) {
      Rprintf("SetChainIter: Set CodaList first ! "); R_FlushConsole();
    }
    if (Rf_length(CodaList) > Set || Set <= 0) {
       Rprintf("Set Chain Iter, set one from 1 to %d \n", Rf_length(CodaList));
    }
    if (RCodaTable != NULL) {
      DDelete(RCodaTable);
    }
    ChainIter = Set;
    RCodaTable = new AObject(VECTOR_ELT(CodaList, Set-1));
    CodaTable = RCodaTable->asSexp();
  } */
  void rPutRNGstate() { OnRNGState = 1;  PutRNGstate(); }
  SEXP get_OnSigma() { return(sOnSigma);}  
  SEXP get_SigmaPrior() {
    if (RSigmaPrior == NULL) { return(R_NilValue); }
    return(RSigmaPrior->asSexp());
  }
  int getLSigmaPrior() {
    if (RSigmaPrior == NULL) { return(0); }
    return(Rf_length(RSigmaPrior->asSexp()));
  }
  SEXP get_iiWeight() {
    if (Verbose >= 1) {
      if (iiWeight == NULL) {
        Rprintf("get_iiWeight, iiWeight is NULL!\n"); R_FlushConsole();
      } else {
        Rprintf("get_iiWeight, iiWeight is not NULL\n"); R_FlushConsole();
      }
    }
    if (iiWeight == NULL) { return(R_NilValue); }
    SEXP sOut;
    Rf_protect(sOut = Rf_allocVector(REALSXP, n));
    int One = 1;
    F77_CALL(dcopy)(&n, iiWeight, &One, REAL(sOut), &One);
    Rf_unprotect(1);
    return(sOut); 
  }
  SEXP get_TotalOpenedFiles() {
    SEXP sOn= R_NilValue;  Rf_protect(sOn = Rf_allocVector(INTSXP,1));
    INTEGER(sOn)[0]  = TotalOpenedFiles; Rf_unprotect(1); return(sOn);
  }
  SEXP get_TotalClosedFiles() {
    SEXP sOn= R_NilValue;  Rf_protect(sOn = Rf_allocVector(INTSXP,1));
    INTEGER(sOn)[0]  = TotalClosedFiles; Rf_unprotect(1); return(sOn);
  }
  SEXP get_TotalR5OpenedFiles() {
    SEXP sOn= R_NilValue;  Rf_protect(sOn = Rf_allocVector(INTSXP,1));
    INTEGER(sOn)[0]  = TotalR5OpenedFiles; Rf_unprotect(1); return(sOn);
  }
  SEXP get_TotalR5ClosedFiles() {
    SEXP sOn= R_NilValue;  Rf_protect(sOn = Rf_allocVector(INTSXP,1));
    INTEGER(sOn)[0]  = TotalR5ClosedFiles; Rf_unprotect(1); return(sOn);
  }
  void set_TotalR5OpenedFiles(SEXP iIn) {
    if (Rf_isNull(iIn) && !(Rf_isInteger(iIn) || Rf_isReal(iIn))){
      Rf_error("set_TotalR5OpenedFiles: Hey, give us an actuall iIn!\n");
    }
    TotalR5OpenedFiles = GetFirstInteger(iIn);
  }
  void set_TotalR5ClosedFiles(SEXP iIn) {
    if (Rf_isNull(iIn) && !(Rf_isInteger(iIn) || Rf_isReal(iIn))){
      Rf_error("set_TotalR5ClosedFiles: Hey, give us an actuall iIn!\n");
    }
    TotalR5ClosedFiles = GetFirstInteger(iIn);
  }
  int CheckpXtX();
  SEXP get_XtY() {
   SEXP sOut;
   Rf_protect(sOut = Rf_allocVector(REALSXP, p));
   int One = 1;
   F77_CALL(dcopy)(&p, XtY, &One, REAL(sOut), &One);
   Rf_unprotect(1);
   return(sOut);
  }
  SEXP get_CodaTjj() {
    SEXP sOut;
    if (CodaTjj == NULL) { return(R_NilValue); }
    if (OnCodaT > LengthCodaTjj) {OnCodaT = LengthCodaTjj;}
    Rf_protect(sOut = Rf_allocVector(INTSXP, OnCodaT));
    int One = 1;
    for (One = 0; One < OnCodaT; One++) {
      INTEGER(sOut)[One] = Codajj[One];
    }
    Rf_unprotect(1); return(sOut);
  }
  SEXP get_CodaTLocjj() {
    SEXP sOut;
    if (CodaTLocjj == NULL) { return(R_NilValue); }
    if (OnCodaT > LengthCodaTjj) {OnCodaT = LengthCodaTjj;}
    Rf_protect(sOut = Rf_allocVector(INTSXP, OnCodaTLocjj));
    int One = 1;
    for (One = 0; One < OnCodaTLocjj; One++) {
      INTEGER(sOut)[One] = (int) CodaTLocjj[One];
    }
    Rf_unprotect(1); return(sOut);
  }
  SEXP get_Codajj() {
    SEXP sOut;
    if (Codajj == NULL) { return(R_NilValue); }
    Rf_protect(sOut = Rf_allocVector(INTSXP, OnCoda));
    int One = 1;
    for (One = 0; One < OnCoda; One++) {
      INTEGER(sOut)[One] = Codajj[One];
    }
    Rf_unprotect(1); return(sOut);
  }

  SEXP get_CodaBeta() {
    int All = OnCoda;
    SEXP sOut;
    if (CodaBeta == NULL) { return(R_NilValue); }
    Rf_protect(sOut = Rf_allocVector(REALSXP, OnCoda));
    int One = 1;
    F77_CALL(dcopy)(&All, CodaBeta, &One, REAL(sOut), &One);
    Rf_unprotect(1); return(sOut);
  }
  int get_CodaTisSetup() {
    return(CodaTisSetup);
  }
  SEXP get_CodaTau() {
    int All = OnCodaT;
    SEXP sOut;
    if (CodaTau == NULL) { return(R_NilValue); }
    Rf_protect(sOut = Rf_allocVector(REALSXP, OnCodaT));
    int One = 1;
    F77_CALL(dcopy)(&All, CodaTau, &One, REAL(sOut), &One);
    Rf_unprotect(1); return(sOut);
  }
  SEXP get_SubXtX(SEXP sSet) {
    if (Rf_isNull(sSet) || Rf_length(sSet) <= 0 || (!Rf_isInteger(sSet) && !Rf_isReal(sSet))) {
      Rf_error("get_SubXtX Bad Non Integer subset Entry!\n");
    }
    int ATo = 0;
    for (int ii = 0; ii < Rf_length(sSet); ii++) {
      if (ANINT(sSet, ii) <= 0 || ANINT(sSet, ii) > OnKappaS) {
         Rprintf("get_SubStX, set includes %d=%d, which is not good, Note OnKappaS=%d\n",
           ii, ANINT(sSet,ii), OnKappaS); R_FlushConsole(); ATo++;
      }
    }
    if (ATo >= 1) {
      Rf_error("get_SubStX, sorry, you've got to include less than OnKappaS=%d!", OnKappaS);
    }   
    SEXP sOut;
    Rf_protect(sOut = Rf_allocMatrix(REALSXP, p, Rf_length(sSet)));
    int One = 1;
    for (int ii = 0; ii < Rf_length(sSet); ii++) {
      F77_CALL(dcopy)(&p, pXtX[ANINT(sSet,ii)-1], &One, REAL(sOut) + p * ii, &One);
    }
    Rf_unprotect(1);
    return(sOut);
  }
  SEXP get_XtX() {
   int AKappaS = OnKappaS;
   if ( ((double)p)/((double)8388608)  * ((double)OnKappaS) >= 128.0 ) {
     Rprintf("get_XtX: No, not enough space for p=%d, OnKappaS=%d!\n",
       p, OnKappaS); R_FlushConsole();
   }
   AKappaS = (int) floor(((double)8388608) /  ((double)p) * 128.0);
   if (OnKappaS < AKappaS) {AKappaS = OnKappaS;}
   if (Verbose >= 1) {
     Rprintf("BayesSpikeCpp.h:::get_XtX() we are now going to load with p=%d, AKappaS = %d, OnKappaS=%d/%d/%d\n",
       p, AKappaS, OnKappaS, OnKappaMem, p); R_FlushConsole();
   }
   //int All = p * OnKappaS; 
   SEXP sOut;
   Rf_protect(sOut = Rf_allocMatrix(REALSXP, p, AKappaS));
   int One = 1; 
   for (int ii = 0; ii < AKappaS; ii++) {  
     if (pXtX[ii] == NULL) {
       Rprintf("get_XtX, error, ii = %d/%d is NULL. \n", ii, AKappaS); R_FlushConsole();
     } else {
       F77_CALL(dcopy)(&p, pXtX[ii], &One, REAL(sOut) + p*ii, &One);
     }
   }
   Rf_unprotect(1);
   return(sOut);
  }
  SEXP get_PropBeta() {
   int All = NumActive;  SEXP sOut;
   if (PropBeta == NULL) { return(R_NilValue);}
   Rf_protect(sOut = Rf_allocVector(REALSXP, All));
   int One = 1;   
   F77_CALL(dcopy)(&All, PropBeta, &One, REAL(sOut), &One);
   Rf_unprotect(1);
   return(sOut);
  }
  int get_LengthMergeActive() { return(LengthMergeActive); }
  int get_MaxLengthMergeActive() { return(MaxLengthMergeActive); }
  SEXP get_MergeActive() { 
   if (LengthMergeActive <= 0 || MaxLengthMergeActive == 0) {
     return(R_NilValue);
   }
   int All = LengthMergeActive;  SEXP sOut;
   if (MergeActive == NULL) { return(R_NilValue);}
   Rf_protect(sOut = Rf_allocVector(INTSXP, All));   
   for (int ii = 0; ii < LengthMergeActive; ii++) {
     INTEGER(sOut)[ii] = MergeActive[ii];
   }
   Rf_unprotect(1);
   return(sOut);
  }
  SEXP get_MergePropBeta() { 
   if (LengthMergeActive <= 0 || MaxLengthMergeActive == 0) {
     return(R_NilValue);
   }
   int All = LengthMergeActive;  SEXP sOut;
   if (MergePropBeta == NULL) { return(R_NilValue);}
   Rf_protect(sOut = Rf_allocVector(REALSXP, All));
   int One = 1;   
   F77_CALL(dcopy)(&All, MergePropBeta, &One, REAL(sOut), &One);
   Rf_unprotect(1);
   return(sOut);
  }
  SEXP get_XtResid() {
   int All = p;  SEXP sOut;
   Rf_protect(sOut = Rf_allocVector(REALSXP, p));
   int One = 1;   
   F77_CALL(dcopy)(&All, XtResid, &One, REAL(sOut), &One);
   Rf_unprotect(1);
   return(sOut);
  }
  SEXP get_XLC() {
    SEXP sOut;
    Rf_protect(sOut = Rf_allocVector(INTSXP, p));
    int ii = 0;
    for (ii = 0; ii < p; ii++) {
      INTEGER(sOut)[ii] = XLC[ii];
    }
    Rf_unprotect(1);
    return(sOut);
  }
  SEXP get_kXFinder() {
    SEXP sOut;
    Rf_protect(sOut = Rf_allocVector(INTSXP, OnKappaS));
    int ii = 0;
    for (ii = 0; ii < OnKappaS; ii++) {
      INTEGER(sOut)[ii] = kXFinder[ii];
    }
    Rf_unprotect(1);
    return(sOut);
  }
  SEXP get_BFixed() {
    SEXP sOut;  int ii = 0;
    if (sOnTau == NULL || Rf_isNull(sOnTau) || Rf_length(sOnTau) <= 0 ||
      tauEndList == NULL || Rf_isNull(tauEndList) || Rf_length(tauEndList) <= 0) {
      if (BFixed== NULL) {
        Rf_error("get_BFixed: BFixed never Allocated");
      }
      Rf_protect(sOut = Rf_allocVector(INTSXP, p));
      for (ii = 0; ii < p; ii++) {
        INTEGER(sOut)[ii] = (int) BFixed[ii];
      }
      Rf_unprotect(1);
      return(sOut);
    } else if (iFirstRandom <= 0) {
      if (iFirstRandom < 0) {
        Rf_error("get_BFixed: iFirstRandom is equal to zero!");
      }
    } 
    if (BFixed== NULL) {
      Rf_error("get_BFixed: BFixed never Allocated");
    }
    Rf_protect(sOut = Rf_allocVector(INTSXP, iFirstRandom));
    for (ii = 0; ii < iFirstRandom; ii++) {
      INTEGER(sOut)[ii] = (int) BFixed[ii];
    }
    Rf_unprotect(1);
    return(sOut);
  }
  void set_NoShrinkFixed(SEXP _NoShrinkFixed) {
    if (Rf_isNull(_NoShrinkFixed) || Rf_length(_NoShrinkFixed) <= 0) {
      if (NoShrinkFixed != NULL) {
        FFree(NoShrinkFixed, "NoShrinkFixed"); 
        FFree(NoShrinkFixedPrior, "NoShrinkFixedPrior"); 
        LengthNoShrinkFixed = 0;
      }
      return;
    }
    if (!Rf_isReal(_NoShrinkFixed) && !Rf_isInteger(_NoShrinkFixed)) {
      Rf_error("set_NoShrinkFixed: Hey this isn't real or integer, can't read it!\n");
    }
    if (NoShrinkFixed == NULL) {
      RMemGetI(NoShrinkFixed, "NoShrinkFixed", Rf_length(_NoShrinkFixed));
      RMemGetD(NoShrinkFixedPrior, "NoShrinkFixedPrior",  Rf_length(_NoShrinkFixed));
      for (int ii = 0; ii < Rf_length(_NoShrinkFixed); ii++) {
        NoShrinkFixed[ii] = ANINT(_NoShrinkFixed, ii);
        NoShrinkFixedPrior[ii] = 1.0;
      }
      LengthNoShrinkFixed = Rf_length(_NoShrinkFixed);
    } else {
       FFree(NoShrinkFixed, "NoShrinkFixed");  
       FFree(NoShrinkFixedPrior, "NoShrinkFixedPrior");
       LengthNoShrinkFixed = 0;
       RMemGetI(NoShrinkFixed, "NoShrinkFixed", Rf_length(_NoShrinkFixed));
       RMemGetD(NoShrinkFixedPrior, "NoShrinkFixedPrior",  Rf_length(_NoShrinkFixed));
        for (int ii = 0; ii < Rf_length(_NoShrinkFixed); ii++) {
          NoShrinkFixed[ii] = ANINT(_NoShrinkFixed, ii);
          NoShrinkFixedPrior[ii] = 1.0;
        }
        LengthNoShrinkFixed = Rf_length(_NoShrinkFixed);
    }
    return;
  }
  void set_NoShrinkRandom(SEXP _NoShrinkRandom) {
    if (Rf_isNull(_NoShrinkRandom) || Rf_length(_NoShrinkRandom) <= 0) {
      if (NoShrinkRandom != NULL) {
        FFree(NoShrinkRandom, "NoShrinkRandom"); 
        FFree(NoShrinkRandomPrior, "NoShrinkRandomPrior"); 
        LengthNoShrinkRandom = 0;
      }
      return;
    }
    if (!Rf_isReal(_NoShrinkRandom) && !Rf_isInteger(_NoShrinkRandom)) {
      Rf_error("set_NoShrinkRandom: Hey input isn't real or integer, don't know what it is\n");
    }
    if (NoShrinkRandom == NULL) {
      RMemGetI(NoShrinkRandom, "NoShrinkRandom", Rf_length(_NoShrinkRandom));
      RMemGetD(NoShrinkRandomPrior, "NoShrinkRandomPrior", Rf_length(_NoShrinkRandom));
      for (int ii = 0; ii < Rf_length(_NoShrinkRandom); ii++) {
        NoShrinkRandom[ii] = ANINT(_NoShrinkRandom, ii);
        NoShrinkRandomPrior[ii] = 1.0;
      }
      LengthNoShrinkRandom = Rf_length(_NoShrinkRandom);
    } else {
       FFree(NoShrinkRandom, "NoShrinkRandom");  
       FFree(NoShrinkRandomPrior, "NoShrinkRandomPrior");
       LengthNoShrinkRandom = 0;
       RMemGetI(NoShrinkRandom, "NoShrinkRandom", Rf_length(_NoShrinkRandom));
       RMemGetD(NoShrinkRandomPrior, "NoShrinkRandomPrior", Rf_length(_NoShrinkRandom));
        for (int ii = 0; ii < Rf_length(_NoShrinkRandom); ii++) {
          NoShrinkRandom[ii] = ANINT(_NoShrinkRandom, ii);
          NoShrinkRandomPrior[ii] = 1.0;
        }
        LengthNoShrinkRandom = Rf_length(_NoShrinkRandom);
    }
    return;
  }
  int get_LengthNoShrinkFixed() { return(LengthNoShrinkFixed); }
  SEXP get_NoShrinkFixed() {
    SEXP sOut = R_NilValue;
    if (NoShrinkFixed == NULL) { return(R_NilValue); }
    if (LengthNoShrinkFixed <= 0) { return(R_NilValue); }
    Rf_protect(sOut = Rf_allocVector(INTSXP, LengthNoShrinkFixed));
    for (int ii = 0; ii < LengthNoShrinkFixed; ii++) {
      INTEGER(sOut)[ii] = NoShrinkFixed[ii] + 1;
    }
    Rf_unprotect(1);
    return(sOut);
  }
  SEXP get_CNoShrinkFixed() {
    SEXP sOut = R_NilValue;
    if (NoShrinkFixed == NULL) { return(R_NilValue); }
    if (LengthNoShrinkFixed <= 0) { return(R_NilValue); }
    Rf_protect(sOut = Rf_allocVector(INTSXP, LengthNoShrinkFixed));
    for (int ii = 0; ii < LengthNoShrinkFixed; ii++) {
      INTEGER(sOut)[ii] = NoShrinkFixed[ii];
    }
    Rf_unprotect(1);
    return(sOut);  
  }
  SEXP get_NoShrinkFixedPrior() {
    SEXP sOut = R_NilValue;
    if (NoShrinkFixedPrior == NULL) { return(R_NilValue); }
    if (LengthNoShrinkFixed <= 0) { return(R_NilValue); }
    Rf_protect(sOut = Rf_allocVector(REALSXP, LengthNoShrinkFixed));
    for (int ii = 0; ii < LengthNoShrinkFixed; ii++) {
      REAL(sOut)[ii] = NoShrinkFixedPrior[ii] + 1;
    }
    Rf_unprotect(1);
    return(sOut);
  }
  SEXP get_NoShrinkRandom() {
    SEXP sOut = R_NilValue;
    if (NoShrinkRandom == NULL) { return(R_NilValue); }
    if (LengthNoShrinkRandom <= 0) { return(R_NilValue); }
    Rf_protect(sOut = Rf_allocVector(INTSXP, LengthNoShrinkRandom));
    for (int ii = 0; ii < LengthNoShrinkRandom; ii++) {
      INTEGER(sOut)[ii] = NoShrinkRandom[ii] + 1;
    }
    Rf_unprotect(1);
    return(sOut);
  }
  SEXP get_NoShrinkRandomPrior() {
    SEXP sOut = R_NilValue;
    if (NoShrinkRandomPrior == NULL) { return(R_NilValue); }
    if (LengthNoShrinkRandom <= 0) { return(R_NilValue); }
    Rf_protect(sOut = Rf_allocVector(REALSXP, LengthNoShrinkRandom));
    for (int ii = 0; ii < LengthNoShrinkRandom; ii++) {
      REAL(sOut)[ii] = NoShrinkRandomPrior[ii] + 1;
    }
    Rf_unprotect(1);
    return(sOut);
  }
  void set_NoShrinkFixedPrior(SEXP _NoShrinkFixedPrior) {
    if (Rf_isNull(_NoShrinkFixedPrior) || Rf_length(_NoShrinkFixedPrior) <= 0) {
      if (NoShrinkFixedPrior != NULL) {
        FFree(NoShrinkFixedPrior, "NoShrinkFixedPrior");
      }
      return;
    }
    if (NoShrinkFixed == NULL) {
      Rf_error(" set_NoShrinkFixedPrior: First Set No ShrinkFixed");
      R_FlushConsole();
    }
    if (LengthNoShrinkFixed != Rf_length(_NoShrinkFixedPrior)) {
      Rf_error(" set_NoShrinkFixedPrior:  NoShrinkFixedPrior of different length %d than NSV = %d\n",
        Rf_length(_NoShrinkFixedPrior), LengthNoShrinkFixed);
    }
    if (NoShrinkFixedPrior == NULL) {
      RMemGetD(NoShrinkFixedPrior, "NoShrinkFixedPrior", LengthNoShrinkFixed);
      for (int ii = 0; ii < LengthNoShrinkFixed; ii++) {
        NoShrinkFixedPrior[ii] = REAL(_NoShrinkFixedPrior)[ii];
      }
    } else {
      FFree(NoShrinkFixedPrior, "NoShrinkFixedPrior");
      RMemGetD(NoShrinkFixedPrior, "NoShrinkFixedPrior", LengthNoShrinkFixed);
      for (int ii = 0; ii < LengthNoShrinkFixed; ii++) {
        NoShrinkFixedPrior[ii] = REAL(_NoShrinkFixedPrior)[ii];
      }
    }
    return;
  }
  void set_NoShrinkRandomPrior(SEXP _NoShrinkRandomPrior) {
    if (Rf_isNull(_NoShrinkRandomPrior) || Rf_length(_NoShrinkRandomPrior) <= 0) {
      if (NoShrinkRandomPrior != NULL) {
        FFree(NoShrinkRandomPrior, "NoShrinkRandomPrior");
      }
      return;
    }
    if (NoShrinkRandom == NULL) {
      Rf_error(" set_NoShrinkRandomPrior: First Set No ShrinkRandom");
      R_FlushConsole();
    }
    if (LengthNoShrinkRandom != Rf_length(_NoShrinkRandomPrior)) {
      Rf_error(" set_NoShrinkRandomPrior:  NoShrinkRandomPrior of different length %d than NSV = %d\n",
        Rf_length(_NoShrinkRandomPrior), LengthNoShrinkRandom);
    }
    if (NoShrinkRandomPrior == NULL) {
      RMemGetD(NoShrinkRandomPrior, "NoShrinkRandomPrior", LengthNoShrinkRandom);
      for (int ii = 0; ii < LengthNoShrinkRandom; ii++) {
        NoShrinkRandomPrior[ii] = REAL(_NoShrinkRandomPrior)[ii];
      }
    } else {
      FFree(NoShrinkRandomPrior, "NoShrinkRandomPrior");
      RMemGetD(NoShrinkRandomPrior, "NoShrinkRandomPrior", LengthNoShrinkRandom);
      for (int ii = 0; ii < LengthNoShrinkRandom; ii++) {
        NoShrinkRandomPrior[ii] = REAL(_NoShrinkRandomPrior)[ii];
      }
    }
    return;
  }
  int SetupProbCoda(int iLengthProbCoda) {
    if (iLengthProbCoda < 0) {
      LengthProbCoda = MaxGibbsIters;
    } else {
      LengthProbCoda = iLengthProbCoda;
    }
    if (Verbose > 2) {
      Rprintf("SetupProbCoda: Freeing before Load");
    }
    FFree(ProbCoda, "ProbCoda");  FFree(ICoda, "ICoda");
    if (LengthProbCoda <= 0 || LengthProbCoda >= 1000000) {
      Rf_error("Hey, SetupProbCoda, LengthProbCoda = %d, this is not good length %d\n",
        LengthProbCoda);
    }
    ProbCoda = (double*) Calloc(LengthProbCoda, double);
    ICoda = (long int*) Calloc(LengthProbCoda, long int);
    iiOnILocProb = 0;
    return(1);
  }
  SEXP get_OldOrderedActive() {
    SEXP sOut;
    Rf_protect(sOut = Rf_allocVector(INTSXP, NumActive));
    int ii = 0;
    for (ii = 0; ii < NumActive; ii++) {
      INTEGER(sOut)[ii] = (int) OrderedActive[ii];
    }
    Rf_unprotect(1);
    return(sOut);  
  }
  SEXP get_OrderedActive() {
    SEXP sOut;
    if (AllNewCoords > 0) {
      AddAllNewCoords();
    }
    RefreshOrderedActive(1);
    Rf_protect(sOut = Rf_allocVector(INTSXP, NumActive));
    int ii = 0;
    for (ii = 0; ii < NumActive; ii++) {
      INTEGER(sOut)[ii] = (int) OrderedActive[ii];
    }
    Rf_unprotect(1);
    return(sOut);
  }
  SEXP get_TrueOrder() {
    SEXP sOut;
    Rf_protect(sOut = Rf_allocVector(INTSXP, 1));
    INTEGER(sOut)[0] = TrueOrder;
    Rf_unprotect(1);
    return(sOut);
  }
  int get_Onii() {return(MT->Onii);}
  
  int get_NInstance() {return(NInstance);}
  void set_NInstance(int _NInstance) { 
    ////////////////////////////////////////////////////////////
    //  Re-enter R and Kill TBSRoo
    NInstance = 1;
    if (FALSE) {
    if (_NInstance <= 0) { Rf_error("set_NInstance Error: NInstance = %d <= 0", _NInstance);}
    NInstance = _NInstance;
    if (Rf_isNull(BayesSpikeNameSpace)) {
      Rprintf("set_NInstance: BayesSpikeNameSpace is Null, gonna be an error \n"); R_FlushConsole();
    }
    SEXP ListOfBayesSpikeOb = 
      Rf_protect(Rf_findVarInFrame( BayesSpikeNameSpace, Rf_install(".ListBayesSpikeOb")));
    if (Rf_isNull(ListOfBayesSpikeOb)) {
      Rf_error("BayesSpikeCL: Error: no ListOfBayesSPikeOb");
    }
    TBSRoo = VECTOR_ELT(ListOfBayesSpikeOb,_NInstance-1);
    Rf_unprotect(1);
    }
  }
  
  SEXP get_XtXRecorded() {
    if (XtXRecorded == NULL || LengthXtXRecorded == 0) { return(R_NilValue);}
    SEXP sOut = R_NilValue;
    Rf_protect(sOut = Rf_allocVector(REALSXP, LengthXtXRecorded));
    int ii;
    for (ii = 0; ii < LengthXtXRecorded; ii++) {
      REAL(sOut)[ii] = XtXRecorded[ii];
    }
    Rf_unprotect(1);
    return(sOut);
  }
  SEXP get_Beta() { 
    if (RsBeta == NULL) {
      Rf_error("get_Beta: Really bad, rnsBeta is NULL!\n");
    } else if (sBeta == NULL || Rf_isNull(sBeta)) {
      Rf_error("get_Beta: What's up, RsBeta is not null but sBeta is NULL!");
    } else if (Rf_length(sBeta) <= 0) {
      Rprintf("get_Beta: Issue, length sBeta is 0?  p = %d Why?", p);
      Rf_error("  get_Beta: This has got to be an error!\n");
    } else if (Rf_length(sBeta) != p) {
      Rf_error("get_Beta: Awkward, length sBeta=%d, but p = %d!\n", Rf_length(sBeta), p);
    }
    return(sBeta); 
  }
  SEXP get_X() { return(sX); }
  SEXP get_Y() { return(sY); }
  void set_iiWeight(SEXP iiWeight_) {  
    if (Rf_isNull(iiWeight_)) {
      Rprintf("set_iiWeight(): no, iiWeight is NULL. Won't kill. \n"); R_FlushConsole();
      return;
    }
    if (Rf_length(iiWeight_) != n) {
      Rf_error("Please Supply Weight Matrix of correct length!");
    }
    if (!Rf_isReal(iiWeight_)) {
      Rf_error("set_iiWeight(): it is not real!\n");
    }
    if (iiWeight == NULL) {
      RMemGetD(iiWeight, "iiWeight", n);
    }
    int One = 1;
    F77_CALL(dcopy)(&n, REAL(iiWeight_), &One, iiWeight, &One);
    nWeight = F77_CALL(dasum)(&n, iiWeight, &One);
    if (iWeightedXtX == NULL) {
      if (Verbose >= 5) {
        Rprintf("BayesSpikeGibbs.cpp:::UpdateTNoise() iWeightedXtX Is NULL Allocating length %d=OnKappaMem\n", OnKappaMem); R_FlushConsole();
      }
      RMemGetS(iWeightedXtX, "iWeightedXtX", OnKappaMem);
      if (Verbose >= 5) {
        Rprintf("BayesSpikeGibbs.cpp:::UpdateTNoise() Finished Allocating iWeightedXtX length %d=OnKappaMem\n", OnKappaMem); R_FlushConsole();
      }
      for (int iti = 0; iti < OnKappaMem; iti++) {
        iWeightedXtX[iti] = (short int) 0;
      }
    }
    
    if (Verbose >= 1) {
      Rprintf("set_iiWeight(): Completed Copy of iiWeights. Now UpdateXtY\n"); R_FlushConsole();
    }
    UpdateXtY(); 
    if (Verbose >= 1) {
      Rprintf("set_iiWeight(): Now ReorderXtX\n");  R_FlushConsole();
    }
    ReorderXtX(); 
    if (Verbose >= 1) {
      Rprintf("set_iiWeight(): Now RefillXtX\n");  R_FlushConsole();
    }
    RefillXtX(); 
    if (Verbose >= 1) {
      Rprintf("set_iiWeight(): Now UpdateFreshXtResid\n");  R_FlushConsole();
    }
    UpdateFreshXtResid(); 
    if (tauEndList != NULL && !Rf_isNull(tauEndList) && Rf_length(tauEndList)) {
      if (Verbose >= 1) {
       Rprintf("set_iiWeight(): Now ReweightEigenvalues\n");  R_FlushConsole();
      }
      ReweightEigenvalues();
    }
    FullCrazyTestiWeightedXtX((char*) "End set_iiWeight");
    if (Verbose >= 1) {
      Rprintf("set_iiWeight(): Now set_iiWeight Completed.\n");  R_FlushConsole();
    }
  }
  void set_X(SEXP X_) {
    if (INTEGER(Rf_getAttrib(X_, R_DimSymbol))[1] != p ||
      INTEGER(Rf_getAttrib(X_, R_DimSymbol))[0] != n) {
      Rf_error("setX:  Hey, use the same dim as before, %d, %d", n,p);
    }  
    DDelete(RsX, "RsX");  RsX = new AObject(X_);
    sX = RsX->asSexp();
    UpdateXtY(); ReorderXtX(); RefillXtX(); UpdateFreshXtResid();
  }
  void set_TBSRoo(SEXP TBSRoo_) {
    if (Rf_isNull(TBSRoo_)) {
      Rf_error("TBSRoo Given is Null!");
    }
    DDelete(RTBSRoo, "RTBSRoo");
    RTBSRoo = new AObject(TBSRoo_);
    TBSRoo = RTBSRoo->asSexp();
  }
  void set_TBSR5(SEXP TBSR5_) {
    if (Rf_isNull(TBSR5_)) {
      Rf_error("TBSR5 Given is Null!");
    }
    DDeleteA(RTBSR5, "RTBSR5");
    RTBSR5 = new AObject(TBSR5_);
    TBSR5 = RTBSR5->asSexp();
  }
  void ClearPctBetas() {
    if (PctBetas == NULL) {
      if (Verbose >= 0) {
        Rprintf("ClearPctBetas: PctBetas already cleared.\n");
        R_FlushConsole();
      }
      return;
    }
    FFree(PctBetas, "PctBetas");  PctBetas = NULL; return;
  }
  SEXP get_PctBetas() {
    if (PctBetas == NULL) {
      SEXP sStartIter = R_NilValue;
      Rf_protect(sStartIter = Rf_allocVector(INTSXP,1));
      SEXP sEndIter = R_NilValue;
      Rf_protect(sEndIter = Rf_allocVector(INTSXP,1));
      INTEGER(sStartIter)[0] = 0;
      INTEGER(sEndIter)[0] = MaxGibbsIters;
      if (StartIter >= 0) {
        INTEGER(sStartIter)[0] = StartIter; 
      }
      if (EndIter > StartIter && EndIter < MaxGibbsIters) {
        INTEGER(sEndIter)[0] = EndIter;
      } 
      if (Verbose >= 2) {
        Rprintf("About to call ReadSetPct Betas, Start=%d, End=%d\n",
          INTEGER(sStartIter)[0], INTEGER(sEndIter)[0]); 
        R_FlushConsole();
      }
      ReadSetPctBetas(sStartIter, sEndIter);
      Rf_unprotect(2);
    }
    if (PctBetas == NULL) {
      Rprintf("get_PctBetas:  Error, load in did not generate PctBetas");
      return(R_NilValue);
    }
    int One = 1;
    SEXP sOut = R_NilValue;
    Rf_protect(sOut = Rf_allocVector(REALSXP, p));
    F77_CALL(dcopy)(&p, PctBetas, &One, REAL(sOut), &One);
    Rf_unprotect(1);
    return(sOut);
  }
  SEXP get_PctTaus() {
    if (Rf_isNull(tauEndList) || Rf_isNull(sOnTau)) {
      return(R_NilValue);
    }
    if (PctTaus == NULL) {
      SEXP sStartIter = R_NilValue;
      Rf_protect(sStartIter = Rf_allocVector(INTSXP,1));
      SEXP sEndIter = R_NilValue;
      Rf_protect(sEndIter = Rf_allocVector(INTSXP,1));
      INTEGER(sStartIter)[0] = 0;
      INTEGER(sEndIter)[0] = MaxGibbsIters;
      if (StartIter >= 0) {
        INTEGER(sStartIter)[0] = StartIter; 
      }
      if (EndIter > StartIter && EndIter < MaxGibbsIters) {
        INTEGER(sEndIter)[0] = EndIter;
      } 
      ReadSetPctBetas(sStartIter, sEndIter);
      Rf_unprotect(2);
    }
    if (PctBetas == NULL) {
      Rprintf("get_PctTaus:  Error, load in did not generate PctTaus");
      return(R_NilValue);
    }
    int k = Rf_length(tauEndList);
    int One = 1;
    SEXP sOut = R_NilValue;
    Rf_protect(sOut = Rf_allocVector(REALSXP, k));
    F77_CALL(dcopy)(&k, PctTaus, &One, REAL(sOut), &One);
    Rf_unprotect(1);
    return(sOut);
  }
  SEXP get_TBSRoo() {
    if (FALSE) {
    if (NInstance <= 0) { Rf_error("set_NInstance Error: NInstance = %d <= 0", NInstance);}
    SEXP ListOfBayesSpikeOb = 
      Rf_protect(Rf_findVarInFrame( BayesSpikeNameSpace, Rf_install(".ListBayesSpikeOb")));
    if (Rf_isNull(ListOfBayesSpikeOb)) {
      Rf_error("BayesSpikeCL: Error: no ListOfBayesSPikeOb");
    }
    TBSRoo = VECTOR_ELT(ListOfBayesSpikeOb, NInstance-1);
    if (Rf_isNull(TBSRoo)) {
      Rprintf("get_TBSRoo: Returns NULL TBSRoo for NInstance = %d \n", NInstance); R_FlushConsole();
    }
    Rf_unprotect(1); 
    return(TBSRoo);
    }
    if (this->RTBSR5 == NULL) { return(R_NilValue); }
    return(this->RTBSR5->asSexp());
  }
  SEXP get_TBSR5() {
    if (FALSE) {
    if (NInstance <= 0) { Rf_error("set_NInstance Error: NInstance = %d <= 0", NInstance);}
    SEXP ListOfBayesSpikeOb = 
      Rf_protect(Rf_findVarInFrame( BayesSpikeNameSpace, Rf_install(".ListBayesSpikeOb")));
    if (Rf_isNull(ListOfBayesSpikeOb)) {
      Rf_unprotect(1); Rf_error("BayesSpikeCL: Error: no .ListBayesSpikeOb");
    }
    if (Verbose >= 1) {
      Rprintf("BayesSpikeCL: FoundListOfBayesSpikeOb, length %d looking for %d or N-1=%d\n",
        Rf_length(ListOfBayesSpikeOb), NInstance, NInstance-1);
    }
    if (Rf_length(ListOfBayesSpikeOb) >= NInstance-1) {
      Rf_unprotect(1); Rf_error("BayesSpikeCL: Sorry, Ninstance is %d but .ListBayesSpikeOb Length %d \n",
        NInstance, Rf_length(ListOfBayesSpikeOb));
    }
    SEXP TBSR5 = VECTOR_ELT(ListOfBayesSpikeOb, NInstance-1);
    if (Rf_isNull(TBSR5)) {
      Rprintf("get_TBSR5: Returns NULL TBSR5 for NInstance = %d \n", NInstance); R_FlushConsole();
    }
    Rf_unprotect(1); 
    }
    if (this->RTBSR5 == NULL) { return(R_NilValue); }
    if (Verbose >= 1) {
      Rprintf("-- RTSBR5, Acpp Location is %d \n",RTBSR5->getInVectorLoc());
    }
    return(this->RTBSR5->asSexp());
  }
  void set_Y(SEXP Y_) { 
    if (Rf_length(Y_) != n) {
      Rf_error("setY: Hey, use the same length as before %d",n);
    } 
    if (Rf_isNull(sY)) {
      DDeleteA(RsY, "RsY"); RsY = new AObject(Y_);
      sY = RsY->asSexp();
    } else {
      if (Rf_length(Y_) != Rf_length(sY)) {
        Rf_error("set_Y, sorry, submit vector of length %d\n", Rf_length(sY));
      }
      int One = 1;
      F77_CALL(dcopy)(&n, REAL(Y_), &One, REAL(sY), &One);
    }
    UpdateXtY(); UpdateFreshXtResid();
  }
  //int SetupProbILocOldCoda();
  SEXP get_ProbOldCoda() {
    if (LengthOldCoda <= 0 || ProbOldCoda == NULL) {
      return(R_NilValue);
    }
    SEXP sOn = R_NilValue; int One = 1;
    Rf_protect(sOn = Rf_allocVector(REALSXP, LengthOldCoda));
    One = 1;
    F77_CALL(dcopy)(&LengthOldCoda, ProbOldCoda, &One, REAL(sOn), &One);
    Rf_unprotect(1); return(sOn);
  }
  SEXP get_ProbCoda() {
    if (LengthProbCoda <= 0 || ProbCoda == NULL) {
      return(R_NilValue);
    }
    SEXP sOn = R_NilValue; int One = 1;
    Rf_protect(sOn = Rf_allocVector(REALSXP, LengthProbCoda));
    One = 1;
    F77_CALL(dcopy)(&LengthProbCoda, ProbCoda, &One, REAL(sOn), &One);
    Rf_unprotect(1); return(sOn);
  }
  SEXP get_IOldCoda() {
    if (LengthOldCoda <= 0 || IOldCoda == NULL) {
      return(R_NilValue);
    }
    SEXP sOn = R_NilValue;
    Rf_protect(sOn = Rf_allocVector(INTSXP, LengthOldCoda));
    int ii = 0;
    for (ii = 0; ii < LengthOldCoda; ii++) {
      INTEGER(sOn)[ii] = (int) IOldCoda[ii];
    }
    Rf_unprotect(1); return(sOn);
  }
  SEXP get_RegionWidth() {
    SEXP sOn = R_NilValue;
    Rf_protect(sOn = Rf_allocVector(INTSXP, 1));
    INTEGER(sOn)[0] = RegionWidth;  
    Rf_unprotect(1);
    return(sOn);
  }
  SEXP ShowRegionNeed(SEXP aPoint, SEXP AWidth, SEXP MaxSpace);
  void set_RegionWidth(SEXP iIn) {
    if (Rf_isNull(iIn) || !(Rf_isReal(iIn) || Rf_isInteger(iIn)) ||
      Rf_length(iIn) != 1) {
      Rf_error("Set RegionWidth: no iIn is not valid. \n");
    }
    int IntIn = Rf_isInteger(iIn) ? INTEGER(iIn)[0] : (int) REAL(iIn)[0]; 
    if (IntIn <= 1) {
      Rf_error("Set RegionWidth: No, we can only set RegionWidths to non-zero values, %d\n", IntIn);
    }
    if (IntIn > p) {
      Rf_error("Set RegionWidth, why? p = %d, Why IntIn = %d?\n", p, IntIn);
    }
    RegionWidth = IntIn;
    if (RunProbRegionVector != NULL) { Free(RunProbRegionVector); RunProbRegionVector = NULL; }
    int TTotal = p;
    if (sOnTau != NULL && Rf_length(sOnTau) >= 1) {
      if (iFirstRandom < 0 || iFirstRandom > p) {
        TTotal = Rf_length(sOnTau);
      } else {
        TTotal = iFirstRandom + Rf_length(sOnTau);
      }
    }
    RunProbRegionVector = (double *) Calloc(TTotal, double);
    for (int iti = 0; iti < TTotal; iti++) {
      RunProbRegionVector[iti] = 0.0;
    }
    return;
  }
  SEXP get_RegionMIP(); 
  SEXP get_MIP () {
  if (RMIP == NULL || Rf_isNull(RMIP->asSexp()) ||
      Rf_length(RMIP->asSexp()) <= 0) {
      int PTECTED = 0;
      if (RTBSR5 == NULL) {
        Rprintf("BayesSpikeCL:get_MIP, where is TBSR5?\n");
        R_FlushConsole();
        return(R_NilValue);
      }
      SEXP MIPFunction;
      Rf_protect(MIPFunction = 
         Rf_findVarInFrame( get_BayesSpikeNameSpace(), 
         Rf_install("ATryFillMIP")));
      PTECTED++;
      if (Rf_isNull(MIPFunction)) {
        Rprintf("BayesSpikeCL: Error: no ATryFillMIP.\n");
        Rf_unprotect(PTECTED); PTECTED = 0;
        return(R_NilValue);
      }
      SEXP call;
      Rf_protect(call = Rf_lcons( MIPFunction, 
        Rf_cons(RTBSR5->asSexp(), R_NilValue))); 
      PTECTED++;
      if (Rf_isNull(call)) {
        Rprintf("BayesSpikeCL: MIPFunction, Error, could not genreate call.\n");
        Rf_unprotect(PTECTED); return(R_NilValue);
      }
      SEXP RTM = R_NilValue;
      if (Verbose > 1) {
        Rprintf("get_MIP: Starting to call ATryFillMIP.\n "); R_FlushConsole();
      }
      Rf_protect(RTM = Rf_eval(call, get_BayesSpikeNameSpace())); PTECTED++;
      if (RMIP == NULL) {
        Rprintf("BayesSpikeCL, get_MIP,after all that, still Null List!\n");
        if (!Rf_isNull(RTM)) {
          RMIP = new AObject(RTM);
          Rf_unprotect(PTECTED); PTECTED = 0;
          return(RTM);
        } else {
          Rf_unprotect(PTECTED); PTECTED = 0;
          return(RTM);
        }
      } else {
         Rf_unprotect(PTECTED); PTECTED = 0;
         if (RMIP == NULL) { return(R_NilValue); }
         return(RMIP->asSexp());
      }
      return(R_NilValue);
    }
    if (RMIP == NULL) { return(R_NilValue); }
    return(RMIP->asSexp());  
  }
  void set_MIP(SEXP iMIP) {
    if (Rf_isNull(iMIP)) {
      if (Verbose >= 2) {
      Rprintf("BayesSpikeCL: set MIP to NULL!");
      }
      DDelete(RMIP, "RMIP");
      RMIP = NULL;
      return;
    }
    DDelete(RMIP, "RMIP");
    RMIP = new AObject(iMIP);
  }
  
  SEXP get_RsDeCenteredCodaList() {
    if (RsCenteredColumns == NULL) {
      Rprintf("get_RsDeCenteredCodaList: Sorry, please set RsCenteredColumns first!"); R_FlushConsole();
      return(R_NilValue);
    }
    if ((RCodaList == NULL || Rf_isNull(CodaList)) && RsDeCenteredCodaList == NULL) {
      Rprintf("get_RsDeCenteredCodaList: RsDeCenteredCodaList is a NULL Pointer, you get NULL.  RCodaList is also NULL.\n");
      return(R_NilValue);
    }
    if ((RCodaList == NULL || Rf_isNull(CodaList)) && Rf_isNull(RsDeCenteredCodaList->asSexp())) {
      Rprintf("get_RsDeCenteredCodaList: RsDeCenteredCodaList is a NILValue. RCodaList is also NULL.\n");
      return(R_NilValue);
    }
    if (Verbose >= 3) {
      Rprintf("BayesSpikeCpp.h::get_RsDeCenteredCodaList(): We start\n"); R_FlushConsole();
    }
    if (RsDeCenteredCodaList != NULL && !Rf_isNull(RsDeCenteredCodaList->asSexp()) &&
      Rf_length(RsDeCenteredCodaList->asSexp()) == Rf_length(RCodaList->asSexp())) {
      SEXP sOn = VECTOR_ELT(RsDeCenteredCodaList->asSexp(), 0);
      int nLen = Rf_length(sOn); int One = 1;
      double Fsum = F77_CALL(dasum)(&nLen, REAL(sOn), &One);
      if (Fsum <= .1) {
        if (Verbose >= 3) {
          Rprintf("RsDeCenteredCodaList: It was blank so we are deleting. \n"); R_FlushConsole();
        }
        DDelete(RsDeCenteredCodaList, "RsDeCenteredCodaList");
        RsDeCenteredCodaList = NULL;
      } 
    }
    SEXP MyIn = R_NilValue;
    double MySum = 0.0;
    double MyRCodaSum = 0.0;
    int MyMax = 500;  int One = 1;
    if (RsDeCenteredCodaList != NULL && 
      !Rf_isNull(RsDeCenteredCodaList->asSexp()) && Rf_length(RsDeCenteredCodaList->asSexp()) >= 1) {
      
      if (Verbose >= 3) {
        Rprintf("RsDeCenteredCodaList: Totalling now to count MySum. \n"); R_FlushConsole();
      }
      MyIn = VECTOR_ELT(RsDeCenteredCodaList->asSexp(),
        Rf_length(RsDeCenteredCodaList->asSexp())-1);
      if (MyIn == NULL || Rf_isNull(MyIn) || Rf_length(MyIn) <= 0) {
        Rf_error("RsDeCenteredCodaList: Oh no, RsDeCenteredCodaList, last element is NULL!\n");
      }
      MyMax = Rf_length(MyIn);
      if (MyMax > Rf_length(MyIn)) { MyMax = Rf_length(MyIn); }
      MySum = F77_CALL(dasum)(&MyMax, REAL(MyIn), &One);
      MyMax = 500;
      if (RCodaList != NULL && !Rf_isNull(RCodaList->asSexp())) {
        MyIn = VECTOR_ELT(RCodaList->asSexp(), Rf_length(RCodaList->asSexp())-1);
        if (MyIn == NULL || Rf_isNull(MyIn) || Rf_length(MyIn) <= 0) {
          Rf_error("RsDeCenteredCodaList: Oh no, RsDeCenteredCodaList, last element of RCodaList is NULL!\n");
        }
        if (MyMax > Rf_length(MyIn)) { MyMax = Rf_length(MyIn); }
        MyRCodaSum = F77_CALL(dasum)(&MyMax, REAL(MyIn), &One); 
      }
    }
    if (RCodaList != NULL) {
      MyMax = 500;
      if (RCodaList != NULL && !Rf_isNull(RCodaList->asSexp())) {
        if (Verbose >= 2) {
          Rprintf("RsDeCenteredCodaList: About to Check what the sum of RCodaList is \n"); R_FlushConsole();
        }
        MyIn = VECTOR_ELT(RCodaList->asSexp(), Rf_length(RCodaList->asSexp())-1);
        if (MyIn == NULL || Rf_isNull(MyIn) || Rf_length(MyIn) <= 0) {
          Rf_error("RsDeCenteredCodaList: Oh no, RsDeCenteredCodaList, last element of RCodaList is NULL!\n");
        }
        if (MyMax > Rf_length(MyIn)) { MyMax = Rf_length(MyIn); }
        MyRCodaSum = F77_CALL(dasum)(&MyMax, REAL(MyIn), &One); 
      }
    }
    if (RCodaList != NULL && MyRCodaSum > 0.0 && (RsDeCenteredCodaList == NULL || Rf_isNull(RsDeCenteredCodaList->asSexp()) ||
      Rf_length(RsDeCenteredCodaList->asSexp()) <= 0 || MySum <= 0.0)) {
      int PTECTED = 0;
      if (Verbose >= 3) {
        Rprintf("RsDeCenteredCodaList: RCodaList is not NULL, we will try to retrieve RCodaList");
      }
      if (RTBSR5 == NULL) {
        Rprintf("BayesSpikeCL:get_RsDeCenteredCodaList, where is TBSR5?\n");
        R_FlushConsole();
        return(R_NilValue);
      }
      SEXP RsDeCenteredCodaListFunction;
      Rf_protect(RsDeCenteredCodaListFunction = 
         Rf_findVarInFrame( get_BayesSpikeNameSpace(), 
         Rf_install("ATryFillRsDeCenteredCodaList")));
      PTECTED++;
      if (Rf_isNull(RsDeCenteredCodaListFunction)) {
        Rprintf("BayesSpikeCL: Error: no ATryFillRsDeCenteredCodaList.\n");
        Rf_unprotect(PTECTED); PTECTED = 0;
        return(R_NilValue);
      }
      SEXP call;
      Rf_protect(call = Rf_lcons( RsDeCenteredCodaListFunction, 
        Rf_cons(RTBSR5->asSexp(), R_NilValue))); 
      PTECTED++;
      if (Rf_isNull(call)) {
        Rprintf("BayesSpikeCL: RsDeCenteredCodaList, Error, could not genreate call.\n");
        Rf_unprotect(PTECTED); return(R_NilValue);
      }
      SEXP RTM = R_NilValue;
      if (Verbose >= 3) {
        Rprintf("get_SubCodaList: Starting to call ATryFillRsDeCenteredCodaList.\n "); R_FlushConsole();
      }
      Rf_protect(RTM = Rf_eval(call, get_BayesSpikeNameSpace())); PTECTED++;
      if (RsDeCenteredCodaList == NULL) {
        Rprintf("BayesSpikeCL, get_RsDeCenteredCodaList,after all that, still Null List!\n");
        R_FlushConsole();
        if (!Rf_isNull(RTM)) {
          RsDeCenteredCodaList = new AObject(RTM);
          Rf_unprotect(PTECTED); PTECTED = 0;
          return(RTM);
        } else {
          Rf_unprotect(PTECTED); PTECTED = 0;
          return(RTM);
        }
      } else {
         Rf_unprotect(PTECTED); PTECTED = 0;
         if (RsDeCenteredCodaList == NULL) { return(R_NilValue); }
         return(RsDeCenteredCodaList->asSexp());
      }
      return(R_NilValue);
    }
    if (RsDeCenteredCodaList == NULL || Rf_isNull(RsDeCenteredCodaList->asSexp())) {
      return(R_NilValue);
    }   
    return(RsDeCenteredCodaList->asSexp());
  }
  void set_RsDeCenteredCodaList(SEXP _rsDeCenteredCodaList) {
    if (Rf_isNull(_rsDeCenteredCodaList)) {
      if (Verbose >= 1) {
        Rprintf("BayesSpikeCL: set RsDeCenteredCodaList to NULL!\n");
        R_FlushConsole();
      }
      DDelete(RsDeCenteredCodaList, "RsDeCenteredCodaList");
      RsDeCenteredCodaList = NULL; return;
    }
    DDelete(RsDeCenteredCodaList, "RsDeCenteredCodaList");
    RsDeCenteredCodaList = new AObject(_rsDeCenteredCodaList);
  }
  SEXP get_Copy() {
    int PTECTED=0;
    if (RTBSR5 == NULL) {
        Rprintf("BayesSpikeCL:get_Copy, where is TBSR5?\n");
        R_FlushConsole();
        return(R_NilValue);
    }
    SEXP CopyFunction;
      Rf_protect(CopyFunction = 
         Rf_findVarInFrame( get_BayesSpikeNameSpace(), 
         Rf_install("ATryCopy")));
      PTECTED++;
    SEXP call;
      Rf_protect(call = Rf_lcons( CopyFunction, 
        Rf_cons(RTBSR5->asSexp(), R_NilValue))); 
      PTECTED++;
      if (Rf_isNull(call)) {
        Rprintf("BayesSpikeCL: get_Copy, Error, could not generate call.\n");
        Rf_unprotect(PTECTED); return(R_NilValue);
      }
      SEXP RTM = R_NilValue;
      if (Verbose > 1) {
        Rprintf("get_Copy: Starting to call ATryCopy.\n "); R_FlushConsole();
      }
      Rf_protect(RTM = Rf_eval(call, get_BayesSpikeNameSpace())); PTECTED++;
    if (Rf_isNull(RTM)) {
      Rprintf("Error, ATryCopy returns NULL!\n"); R_FlushConsole();
    }
    Rf_unprotect(PTECTED);
    return(RTM);
  }
  SEXP get_SubCodaList() {
    if (RSubCodaSubSetCoords == NULL) {
      Rprintf("get_SubCodaList: Sorry, please set RSubCodaSubSetCoords first!"); R_FlushConsole();
      return(R_NilValue);
    }
    if (RSubCodaList == NULL || Rf_isNull(RSubCodaList->asSexp()) ||
      Rf_length(RSubCodaList->asSexp()) <= 0) {
      int PTECTED = 0;
      if (RTBSR5 == NULL) {
        Rprintf("BayesSpikeCL:get_SubCodaList, where is TBSR5?\n");
        R_FlushConsole();
        return(R_NilValue);
      }
      SEXP SubCodaListFunction;
      Rf_protect(SubCodaListFunction = 
         Rf_findVarInFrame( get_BayesSpikeNameSpace(), 
         Rf_install("ATryFillSubCodaList")));
      PTECTED++;
      if (Rf_isNull(SubCodaListFunction)) {
        Rprintf("BayesSpikeCL: Error: no ATryFillSubCodaList.\n");
        Rf_unprotect(PTECTED); PTECTED = 0;
        return(R_NilValue);
      }
      SEXP call;
      Rf_protect(call = Rf_lcons( SubCodaListFunction, 
        Rf_cons(RTBSR5->asSexp(), R_NilValue))); 
      PTECTED++;
      if (Rf_isNull(call)) {
        Rprintf("BayesSpikeCL: SubCodaList, Error, could not genreate call.\n");
        Rf_unprotect(PTECTED); return(R_NilValue);
      }
      SEXP RTM = R_NilValue;
      if (Verbose > 1) {
        Rprintf("get_SubCodaList: Starting to call ATryFillSubCodaList.\n "); R_FlushConsole();
      }
      Rf_protect(RTM = Rf_eval(call, get_BayesSpikeNameSpace())); PTECTED++;
      if (RSubCodaList == NULL) {
        Rprintf("BayesSpikeCL, get_SubCodaList,after all that, still Null List!\n");
        R_FlushConsole();
        if (!Rf_isNull(RTM)) {
          RSubCodaList = new AObject(RTM);
          Rf_unprotect(PTECTED); PTECTED = 0;
          return(RTM);
        } else {
          Rf_unprotect(PTECTED); PTECTED = 0;
          return(RTM);
        }
      } else {
         Rf_unprotect(PTECTED); PTECTED = 0;
         return(RSubCodaList->asSexp());
      }
      return(R_NilValue);
    }
    return(RSubCodaList->asSexp());
  }
  void set_SubCodaList(SEXP SSubCodaList) {
    if (Rf_isNull(SSubCodaList)) {
      if (Verbose >= 1) {
        Rprintf("BayesSpikeCL: set SubCodaList to NULL!\n");
        R_FlushConsole();
      }
      DDelete(RSubCodaList, "RSubcodaList");
      RSubCodaList = NULL; return;
    }
    DDelete(RSubCodaList, "RSubCodaList");
    RSubCodaList = new AObject(SSubCodaList);
  }
SEXP get_TauCodaList() {
    if (RSubSetTau == NULL || Rf_length(RSubSetTau->asSexp()) <= 0) {
      Rprintf("get_TauCodaList Sorry, please set RSubSetTau first!"); R_FlushConsole();
      return(R_NilValue);
    }
    if (RTauCodaList == NULL || Rf_isNull(RTauCodaList->asSexp()) ||
      Rf_length(RTauCodaList->asSexp()) <= 0) {
      int PTECTED = 0;
      if (RTBSR5 == NULL) {
        Rprintf("BayesSpikeCL:get_TauCodaList, where is TBSR5?\n");
        R_FlushConsole();
        return(R_NilValue);
      }
      SEXP TauCodaListFunction;
      Rf_protect(TauCodaListFunction = 
         Rf_findVarInFrame( get_BayesSpikeNameSpace(), 
         Rf_install("ATryFillTauCodaList")));
      PTECTED++;
      if (Rf_isNull(TauCodaListFunction)) {
        Rprintf("BayesSpikeCL: Error: no ATryFillTauCodaList.\n");
        Rf_unprotect(PTECTED); PTECTED = 0;
        return(R_NilValue);
      }
      SEXP call;
      Rf_protect(call = Rf_lcons( TauCodaListFunction, 
        Rf_cons(RTBSR5->asSexp(), R_NilValue))); 
      PTECTED++;
      if (Rf_isNull(call)) {
        Rprintf("BayesSpikeCL: TauCodaList, Error, could not generate call.\n");
        Rf_unprotect(PTECTED); return(R_NilValue);
      }
      SEXP RTM = R_NilValue;
      if (Verbose > 1) {
        Rprintf("get_SubCodaList: Starting to call ATryFillTauCodaList.\n "); R_FlushConsole();
      }
      Rf_protect(RTM = Rf_eval(call, get_BayesSpikeNameSpace())); PTECTED++;
      if (RTauCodaList == NULL) {
        Rprintf("BayesSpikeCL, get_TauCodaList, after all that, still Null List!\n");
        R_FlushConsole();
        if (!Rf_isNull(RTM)) {
          RTauCodaList = new AObject(RTM);
          Rf_unprotect(PTECTED); PTECTED = 0;
          return(RTM);
        } else {
          Rf_unprotect(PTECTED); PTECTED = 0;
          return(RTM);
        }
      } else {
         Rf_unprotect(PTECTED); PTECTED = 0;
         return(RTauCodaList->asSexp());
      }
      return(R_NilValue);
    }
    return(RTauCodaList->asSexp());
  }
  void set_TauCodaList(SEXP STauCodaList) {
    if (Rf_isNull(STauCodaList)) {
      if (Verbose >= 1) {
        Rprintf("BayesSpikeCL: set TauCodaList to NULL!\n");
        R_FlushConsole();
      }
      DDelete(RTauCodaList, "RTaucodaList");
      RTauCodaList = NULL; return;
    }
    if (Verbose >= 2) {
      Rprintf("Now Setting TauCodaList!\n"); R_FlushConsole();
    }
    DDelete(RTauCodaList, "RTauCodaList");
    RTauCodaList = new AObject(STauCodaList);
  }
  SEXP get_SubCodaSubSetCoords()  {
    if (RSubCodaSubSetCoords == NULL) { return(R_NilValue); }
    if (Rf_length(RSubCodaSubSetCoords->asSexp()) <= 0) { return(R_NilValue); }
    SEXP sOut = R_NilValue;
    Rf_protect(sOut = Rf_allocVector(INTSXP, 
      Rf_length(RSubCodaSubSetCoords->asSexp())));
    if (Rf_isInteger(RSubCodaSubSetCoords->asSexp())) {
      for (int ii = 0; ii < Rf_length(RSubCodaSubSetCoords->asSexp()); ii++) {
        INTEGER(sOut)[ii] = INTEGER(RSubCodaSubSetCoords->asSexp())[ii]+1;
      }
    } else {
      for (int ii = 0; ii < Rf_length(RSubCodaSubSetCoords->asSexp()); ii++) {
        INTEGER(sOut)[ii] = (int) REAL(RSubCodaSubSetCoords->asSexp())[ii]+1;
      }    
    }
    Rf_unprotect(1);
    return(sOut);
  }
  SEXP get_CSubSetCoords()  {
    if (RSubCodaSubSetCoords == NULL) { return(R_NilValue); }
    if (Rf_length(RSubCodaSubSetCoords->asSexp()) <= 0) { return(R_NilValue); }
    return(RSubCodaSubSetCoords->asSexp());
  }  
  // SubCodaSubSetCoords
  //
  //  If "SubCodaList" is to be loaded which aims to read from files
  //  On disk, but only a subset of the coordinates of Beta, this is the
  //  vector encoding the SubSet Coordinates
  //
  //  Users will usually supply a SEXP with the coordinates+1, CSubSetCoords will not return that.
  //
  void set_CSubSetCoords(SEXP CSubSetCoords) {
    if (Rf_isNull(CSubSetCoords) || Rf_length(CSubSetCoords) <= 0) {
      if (Verbose >= 1) {
        Rprintf("BayesSpikeCL: set SSubCodaVector to NULL!\n");
        R_FlushConsole();
      }
      DDelete(RSubCodaSubSetCoords, "RSubCodaSubSetCoords");
      RSubCodaSubSetCoords = NULL;  
      return;
    }
    SEXP sOn = R_NilValue;
    Rf_protect(sOn=Rf_allocVector(INTSXP, Rf_length(CSubSetCoords)));
    if (Rf_isInteger(CSubSetCoords)) {
      for (int ii = 0; ii < Rf_length(CSubSetCoords); ii++) {
        INTEGER(sOn)[ii] = INTEGER(CSubSetCoords)[ii];
      }
    } else if (Rf_isReal(CSubSetCoords)) {
      for (int ii = 0; ii < Rf_length(CSubSetCoords); ii++) {
        INTEGER(sOn)[ii] = (int) REAL(CSubSetCoords)[ii];
      }
    }
    RSubCodaSubSetCoords = new AObject(sOn); 
    Rf_unprotect(1);
  }
  void set_SubCodaSubSetCoords(SEXP SSubCodaSubSetCoords) {
    if (Rf_isNull(SSubCodaSubSetCoords) || Rf_length(SSubCodaSubSetCoords) <= 0) {
      if (Verbose >= 1) {
        Rprintf("BayesSpikeCL: set SSubCodaVector to NULL!\n");
        R_FlushConsole();
      }
      DDelete(RSubCodaSubSetCoords, "RSubCodaSubSetCoords");
      RSubCodaSubSetCoords = NULL;  
      return;
    }
    DDelete(RSubCodaSubSetCoords, "SubCodaSubSetCoords");
    int PTECTED = 0;
      if (RTBSR5 == NULL) {
        Rprintf("BayesSpikeCL:set_RSubCodaSubSetCoords, where is TBSR5?\n");
        R_FlushConsole();
        return;
      }
      SEXP SubCodaSubSetCoordsFunction = R_NilValue;
      Rf_protect(SubCodaSubSetCoordsFunction = 
         Rf_findVarInFrame( get_BayesSpikeNameSpace(), 
         Rf_install("ATrySetSubCodaSubSetCoords")));
      PTECTED++;
      if (Rf_isNull(SubCodaSubSetCoordsFunction)) {
        Rprintf("BayesSpikeCL: Error: no ATrySetSubCodaSubSetCoords.\n");
        Rf_unprotect(PTECTED); PTECTED = 0;
        return;
      }
      SEXP call;
      Rf_protect(call = Rf_lcons( SubCodaSubSetCoordsFunction, 
        Rf_cons(RTBSR5->asSexp(), 
        Rf_cons(SSubCodaSubSetCoords, R_NilValue))) ); 
      PTECTED++;
      if (Rf_isNull(call)) {
        Rprintf("BayesSpikeCL: Set SubCodaSubSetCoords, Error, could not geneerate call.\n");
        Rf_unprotect(PTECTED); return;
      }
      SEXP RTM = R_NilValue;
      if (Verbose > 1) {
        Rprintf("set_SubCodaVector: Starting to call ATrySetSubCodaSubSetCoords.\n "); R_FlushConsole();
      }
      Rf_protect(RTM = Rf_eval(call, get_BayesSpikeNameSpace())); PTECTED++;
      if (RTM != NULL  && !Rf_isNull(RTM)) {
        if (Verbose > 1) {
          Rprintf("BayesSpikeCL, set_SubCodaSubSetCoords,after all that, ");
          Rprintf(" returned not NULL!\n");
          R_FlushConsole();
        }
        if (!Rf_isNull(RTM)) {
          //RSubCodaSubSetCoords = new AObject(RTM);
          Rf_unprotect(PTECTED); PTECTED = 0;
          return;
        } else {
          RSubCodaSubSetCoords = NULL;
          Rf_unprotect(PTECTED); PTECTED = 0;
          return;
        }
      } 
      if (Verbose >= 1) {
        Rprintf("BayesSpikeCL:  End of set_SubCodaSubSetCoords");
        R_FlushConsole();
      }
      return;
  }
  SEXP get_SubSetTau()  {
    if (RSubSetTau == NULL) { return(R_NilValue); }
    if (Rf_length(RSubSetTau->asSexp()) <= 0) { return(R_NilValue); }
    SEXP sOut = R_NilValue;
    Rf_protect(sOut = Rf_allocVector(INTSXP, 
      Rf_length(RSubSetTau->asSexp())));
    if (Rf_isInteger(RSubSetTau->asSexp())) {
      for (int ii = 0; ii < Rf_length(RSubSetTau->asSexp()); ii++) {
        INTEGER(sOut)[ii] = INTEGER(RSubSetTau->asSexp())[ii]+1;
      }
    } else {
      for (int ii = 0; ii < Rf_length(RSubSetTau->asSexp()); ii++) {
        INTEGER(sOut)[ii] = (int) REAL(RSubSetTau->asSexp())[ii]+1;
      }    
    }
    Rf_unprotect(1);
    return(sOut);
  }
  SEXP get_CSubSetTau()  {
    if (RSubSetTau == NULL) { return(R_NilValue); }
    if (Rf_length(RSubSetTau->asSexp()) <= 0) { return(R_NilValue); }
    SEXP sOut = R_NilValue;
    Rf_protect(sOut = Rf_allocVector(INTSXP, 
      Rf_length(RSubSetTau->asSexp())));
    if (Rf_isInteger(RSubSetTau->asSexp())) {
      for (int ii = 0; ii < Rf_length(RSubSetTau->asSexp()); ii++) {
        INTEGER(sOut)[ii] = INTEGER(RSubSetTau->asSexp())[ii];
      }
    } else {
      for (int ii = 0; ii < Rf_length(RSubSetTau->asSexp()); ii++) {
        INTEGER(sOut)[ii] = (int) REAL(RSubSetTau->asSexp())[ii];
      }    
    }
    Rf_unprotect(1);
    return(sOut);
  }
  void set_CSubSetTau(SEXP SSubSetTau) {
    if (Rf_isNull(SSubSetTau) || Rf_length(SSubSetTau) <= 0) {
      if (Verbose >= 1) {
        Rprintf("BayesSpikeCL: set SubSetTau to NULL!\n");
        R_FlushConsole();
      }
      DDelete(RSubSetTau, "RSubSetTau");
      RSubSetTau = NULL;  
      return;
    }
    SEXP sOn;
    Rf_protect(sOn = Rf_allocVector(INTSXP, Rf_length(SSubSetTau)));  
    if (Rf_isInteger(SSubSetTau)) {
      for (int ii = 0; ii < Rf_length(RSubSetTau->asSexp()); ii++) {
        INTEGER(sOn)[ii] = INTEGER(SSubSetTau)[ii];
      }
    } else {
      for (int ii = 0; ii < Rf_length(RSubSetTau->asSexp()); ii++) {
        INTEGER(sOn)[ii] = (int) REAL(SSubSetTau)[ii];
      }    
    }
    RSubSetTau = new AObject(sOn);
    Rf_unprotect(1);
  }  
  // SubCodaSubSetCoords
  //
  //  If "SubCodaList" is to be loaded which aims to read from files
  //  On disk, but only a subset of the coordinates of Beta, this is the
  //  vector encoding the SubSet Coordinates
  //
  void set_SubSetTau(SEXP SSubSetTau) {
    if (Rf_isNull(SSubSetTau) || Rf_length(SSubSetTau) <= 0) {
      if (Verbose >= 1) {
        Rprintf("BayesSpikeCL: set SubSetTau to NULL!\n");
        R_FlushConsole();
      }
      DDelete(RSubSetTau, "RSubSetTau");
      RSubSetTau = NULL;  
      return;
    }
    DDelete(RSubSetTau, "RSubSetTau");
    int PTECTED = 0;
      if (RTBSR5 == NULL) {
        Rprintf("BayesSpikeCL:set_SubSetTau, where is TBSR5?\n");
        R_FlushConsole();
        return;
      }
      SEXP SetSubSetTauFunction = R_NilValue;
      Rf_protect(SetSubSetTauFunction = 
         Rf_findVarInFrame( get_BayesSpikeNameSpace(), 
         Rf_install("ATrySetSubSetTau")));
      PTECTED++;
      if (Rf_isNull(SetSubSetTauFunction)) {
        Rprintf("BayesSpikeCL: Error: no ATrySetSubSetTau.\n");
        Rf_unprotect(PTECTED); PTECTED = 0;
        return;
      }
      SEXP call;
      Rf_protect(call = Rf_lcons( SetSubSetTauFunction, 
        Rf_cons(RTBSR5->asSexp(), 
        Rf_cons(SSubSetTau, R_NilValue))) ); 
      PTECTED++;
      if (Rf_isNull(call)) {
        Rprintf("BayesSpikeCL: Set SubSetTau, Error, could not geneerate call.\n");
        Rf_unprotect(PTECTED); return;
      }
      SEXP RTM = R_NilValue;
      if (Verbose > 1) {
        Rprintf("set_SubSetTau: Starting to call ATrySetSubSetTau.\n "); R_FlushConsole();
      }
      Rf_protect(RTM = Rf_eval(call, get_BayesSpikeNameSpace())); PTECTED++;
      if (RTM != NULL  && !Rf_isNull(RTM)) {
        if (Verbose > 1) {
          Rprintf("BayesSpikeCL, set_SubSetTau,after all that, ");
          Rprintf(" returned not NULL!\n");
          R_FlushConsole();
        }
        if (!Rf_isNull(RTM)) {
          //RSubSetTau = new AObject(RTM);
          Rf_unprotect(PTECTED); PTECTED = 0;
          return;
        } else {
          RSubSetTau = NULL;
          Rf_unprotect(PTECTED); PTECTED = 0;
          return;
        }
      } 
      if (Verbose >= 1) {
        Rprintf("BayesSpikeCL:  End of set_SubCodaSubSetTau");
        R_FlushConsole();
      }
      return;
  }
  int BlankAllNewCoords() {
    if (Verbose >= -1) {
      Rprintf("BayesSpkeCL:  BlankAllNewCoords(): all we do is set AllNewCoords, NewFixedCoords to zero.\n");
      R_FlushConsole();
    }
    AllNewCoords = 0;  NewFixedCoords = 0;
    return(1);
  }
  int SetAllNewCoords(int AllNewIn, int NewFixedIn) {
    AllNewCoords = AllNewIn;  NewFixedCoords = NewFixedIn;
    return(1);
  }
  int CountAllNewCoords();
  //////////////////////////////////////////////////////////////
  //  PostProbCodaList
  //
  //   PostProbCodaList is a table of posterior probabilities stored on disk.
  //   This function loads back into R to load in those files and attach.
  SEXP get_PostProbCodaList() {
    if (RPostProbCodaList == NULL || Rf_isNull(RPostProbCodaList->asSexp()) ||
      Rf_length(RPostProbCodaList->asSexp()) <= 0) {
      int PTECTED = 0;
      if (RTBSR5 == NULL) {
        Rprintf("BayesSpikeCL:get_PostProbCdoaList, where is TBSR5?\n");
        R_FlushConsole(); return(R_NilValue);
      }
      SEXP PostProbFunction;
      Rf_protect(PostProbFunction = 
         Rf_findVarInFrame( get_BayesSpikeNameSpace(), 
         Rf_install("ATryFillPostProbCoda")));
      PTECTED++;
      if (Rf_isNull(PostProbFunction)) {
        Rprintf("BayesSpikeCL: Error: no ATryFillPostProbCoda.\n");
        Rf_unprotect(PTECTED); PTECTED = 0;
        return(R_NilValue);
      }
      SEXP call;
      Rf_protect(call = Rf_lcons( PostProbFunction, 
        Rf_cons(RTBSR5->asSexp(), R_NilValue))); 
      PTECTED++;
      if (Rf_isNull(call)) {
        Rprintf("BayesSpikeCL: PostProbCodaList, Error, could not genreate call.\n");
        Rf_unprotect(PTECTED); return(R_NilValue);
      }
      SEXP RTM = R_NilValue;
      if (Verbose > 1) {
        Rprintf("get_PostProbCodaList: Starting to call ATryFillPostProbCoda.\n "); R_FlushConsole();
      }
      Rf_protect(RTM = Rf_eval(call, get_BayesSpikeNameSpace())); PTECTED++;
      if (RPostProbCodaList == NULL) {
        Rprintf("BayesSpikeCL, getPostProbCodaList,after all that, still Null List!\n");
        if (!Rf_isNull(RTM)) {
          RPostProbCodaList = new AObject(RTM);
          Rf_unprotect(PTECTED); PTECTED = 0;
          return(RTM);
        } else {
          Rf_unprotect(PTECTED); PTECTED = 0;
          return(RTM);
        }
      } else {
         Rf_unprotect(PTECTED); PTECTED = 0;
         return(RPostProbCodaList->asSexp());
      }
      return(R_NilValue);
    }
    return(RPostProbCodaList->asSexp());
  }
  //////////////////////////////////////////////////////////////
  //  PostProbCodaList
  //
  //   PostProbCodaList is a table of posterior probabilities stored on disk.
  //   This function loads back into R to load in those files and attach.
  SEXP get_SubCodaLongList() {
    if (SubCodaLongList == NULL || Rf_isNull(SubCodaLongList->asSexp()) ||
      Rf_length(SubCodaLongList->asSexp()) <= 0) {
      int PTECTED = 0;
      if (RTBSR5 == NULL) {
        Rprintf("BayesSpikeCL:get_SubCodaLongList, where is TBSR5?\n");
        R_FlushConsole();
        return(R_NilValue);
      }
      SEXP SubCodaLongFunction;
      Rf_protect(SubCodaLongFunction = 
         Rf_findVarInFrame( get_BayesSpikeNameSpace(), 
         Rf_install("ATryFillSubCodaLongList")));
      PTECTED++;
      if (Rf_isNull(SubCodaLongFunction)) {
        Rprintf("BayesSpikeCL: Error: no ATryFillSubCodaLongList.\n");
        Rf_unprotect(PTECTED); PTECTED = 0;
        return(R_NilValue);
      }
      SEXP call;
      Rf_protect(call = Rf_lcons( SubCodaLongFunction, 
        Rf_cons(RTBSR5->asSexp(), R_NilValue))); 
      PTECTED++;
      if (Rf_isNull(call)) {
        Rprintf("BayesSpikeCL: SubCodaLongList, Error, could not genreate call.\n");
        Rf_unprotect(PTECTED); return(R_NilValue);
      }
      SEXP RTM = R_NilValue;
      if (Verbose > 1) {
        Rprintf("get_SubCodaLongList: Starting to call ATryFillSubCodaLongList.\n "); R_FlushConsole();
      }
      Rf_protect(RTM = Rf_eval(call, get_BayesSpikeNameSpace())); PTECTED++;
      if (SubCodaLongList == NULL) {
        Rprintf("BayesSpikeCL, getSubCodaLongList,after all that, still Null List!\n");
        if (!Rf_isNull(RTM)) {
          SubCodaLongList = new RRObject(RTM);
          Rf_unprotect(PTECTED); PTECTED = 0;
          return(RTM);
        } else {
          Rf_unprotect(PTECTED); PTECTED = 0;
          return(RTM);
        }
      } else {
         Rf_unprotect(PTECTED); PTECTED = 0;
         return(SubCodaLongList->asSexp());
      }
      return(R_NilValue);
    }
    return(SubCodaLongList->asSexp());
  }
  void set_SubCodaLongList(SEXP aSet) {
    if (Verbose >= 1) {
      Rprintf("BayesSpikeCpp Class: Setting SubCodaLongList, delete current move on\n"); R_FlushConsole();
    }
    DDelete(SubCodaLongList, "SubCodaLongList");
    if (Verbose >= 1) {
      Rprintf("BayesSpikeCpp Class: Deleted old SubCodaLong List onto Set new\n");
      R_FlushConsole();
    }
    SubCodaLongList = new RRObject(aSet);
  }
  int get_LengthWrittenBetaAllDrawBuffer() {
    return(LengthWrittenBetaAllDrawBuffer);
  }
  int get_LengthTotalWrittenBetaAllDrawBuffer() {
    return(LengthTotalWrittenBetaAllDrawBuffer);
  }
  int get_LengthBetaAllDrawBuffer() {
    return(LengthBetaAllDrawBuffer);
  }
  
  SEXP get_BetaAllDrawBuffer() {
    if (BetaAllDrawBuffer == NULL) { return(R_NilValue); }
    if (LengthWrittenBetaAllDrawBuffer <= 0) {return(R_NilValue); }
    SEXP sOut = R_NilValue;
    int ADraw = LengthWrittenBetaAllDrawBuffer/(p+1);
    int One = 1;
    if (Verbose >= 1) {
      Rprintf("Get_BetaAllDrawBuffer, Written = %d, p+1=%d, ADraw=%d\n",
        LengthWrittenBetaAllDrawBuffer, p+1, ADraw); R_FlushConsole();
    }
    if ( ADraw * (p+1) == LengthWrittenBetaAllDrawBuffer) {
      Rf_protect(sOut  = Rf_allocMatrix(REALSXP, p+1, ADraw));
      F77_CALL(dcopy)(&LengthWrittenBetaAllDrawBuffer, 
        BetaAllDrawBuffer, &One, REAL(sOut), &One);
    } else {
      Rf_protect(sOut  = Rf_allocVector(REALSXP, LengthWrittenBetaAllDrawBuffer));
      F77_CALL(dcopy)(&LengthWrittenBetaAllDrawBuffer, BetaAllDrawBuffer,
        &One, REAL(sOut), &One);
    }
    Rf_unprotect(1); return(sOut);  
  }
  SEXP get_RsCodaBetaAllDrawBuffer() {
    if (RsCodaBetaAllDrawBuffer == NULL) { return(R_NilValue); }
    return(RsCodaBetaAllDrawBuffer->asSexp());
  }
  SEXP get_AllBetaAllDrawBuffer() {
    if (BetaAllDrawBuffer == NULL) { return(R_NilValue); }
    if (LengthWrittenBetaAllDrawBuffer <= 0) {return(R_NilValue); }
    SEXP sOut = R_NilValue; 
    if (LengthBetaAllDrawBuffer <= 0) { return(R_NilValue); }
    Rf_protect(sOut = Rf_allocVector(REALSXP, LengthBetaAllDrawBuffer));
    int One = 1;
    F77_CALL(dcopy)(&LengthBetaAllDrawBuffer, BetaAllDrawBuffer, &One,
      REAL(sOut), &One);
    Rf_unprotect(1); return(sOut); 
  }  
  
  int  set_RsCodaBetaAllDrawBuffer( SEXP _sCodaBetaAllDrawBuffer,
    SEXP sLengthAllDrawBuffer) {
    if (Verbose >= 1) {
      Rprintf("Setting up BetaAllDrawBuffer!\n");
    }

    DDelete(RsCodaBetaAllDrawBuffer, "RsCodaBetaAllDrawBuffer");
    if (Rf_isNull(_sCodaBetaAllDrawBuffer)) {
      LengthBetaAllDrawBuffer = 0;
      return(0);
    }
    int AInt = GetFirstInteger(sLengthAllDrawBuffer);
    if (AInt <= 0 || AInt <= 2*(p+1) ) {
      LengthBetaAllDrawBuffer = 0;
      return(0);
    }
    if (BetaAllDrawBuffer == NULL ||
      LengthBetaAllDrawBuffer != AInt) {
      LengthBetaAllDrawBuffer = AInt;
      FFree(BetaAllDrawBuffer, "BetaAllDrawBuffer");
      if (p >= 100000) {
        LengthBetaAllDrawBuffer = p * 5;
      } else if (LengthBetaAllDrawBuffer > 100000 || LengthBetaAllDrawBuffer <= 0) {
        int LT = (int) floor(100000 / p);
        if (LT <= 4) { LT = 5; }
        LengthBetaAllDrawBuffer = p * LT;
        //Rf_error("Hey, LengthBetaAllDrawBuffer = %d, probably not good!\n", LengthBetaAllDrawBuffer);
      }
      BetaAllDrawBuffer = Calloc(LengthBetaAllDrawBuffer, double);
    }
    RsCodaBetaAllDrawBuffer = new RRObject(_sCodaBetaAllDrawBuffer);
    LengthWrittenBetaAllDrawBuffer= 0;
    LengthTotalWrittenBetaAllDrawBuffer = 0;
    
    return(1);
  } 

  void set_PostProbCodaList(SEXP PostProbCodaList) {
    if (Rf_isNull(PostProbCodaList)) {
      Rprintf("BayesSpikeCL: set PostProbCodaList to NULL!");
      DDelete(RPostProbCodaList, "RPostProbCodaList");
      RPostProbCodaList = NULL;    return;
    }
    DDelete(RPostProbCodaList, "RPostProbCodaList");
    RPostProbCodaList = new AObject(PostProbCodaList);
  }
  // PostProbBuffer is a buffer of Posterior probability values for inclusion of each effect.
  int FillPostProbBuffer(); int WritePostProbBuffer(); 
  SEXP get_sPostProbBufferFile() {
    if (RsPostProbBufferFile == NULL) {
      return(R_NilValue);
    }
    return(RsPostProbBufferFile->asSexp());
  } 
  int set_RsProbBufferFile(SEXP _sPostProbBufferFile, SEXP sLengthPostProbCodaBuffer) {
    if (Rf_isNull(_sPostProbBufferFile)) {
       DDelete(RsPostProbBufferFile, "RsPostProbBufferFile");  
       FFree(PostProbCodaBuffer, "ProbCodaBuffer");
       LengthPostProbCodaBuffer=0;
    }
    if (Rf_isInteger(_sPostProbBufferFile) || Rf_isReal(_sPostProbBufferFile)) {
      Rf_error("set_RsPostProbBufferFile, you gave a numeric for file!\n");
    }
    if (sLengthPostProbCodaBuffer == NULL || Rf_isNull(sLengthPostProbCodaBuffer) || 
      !(Rf_isReal(sLengthPostProbCodaBuffer) || Rf_isInteger(sLengthPostProbCodaBuffer))) {
      Rf_error("set_RsPostProbBufferFile: Invalid sLengthProbCodaBuffer\n");
      
    }
    int _LengthPostProbCodaBuffer = ANINT(sLengthPostProbCodaBuffer, 0);
    if (_LengthPostProbCodaBuffer <= 0) {
      Rf_error("RsPostProbBufferFile: Submit > 0 length!\n"); R_FlushConsole();
      return(0);
    }
    if (Rf_isNull(_sPostProbBufferFile) || !Rf_isString(_sPostProbBufferFile)) {
      Rf_error("RsProbBufferFile: Please submit acceptable ProbBufferFile!");
      return(0);
    }
    //if (LengthProbCoda <= 0 || ProbCoda == NULL) {
    //  Rf_error("RsProbBufferFile: Without ProbCoda set, no setting Up ProbBuffer!\n");
    //}
    DDelete(RsPostProbBufferFile, "RsPostProbBufferFile");  
    FFree(PostProbCodaBuffer, "ProbCodaBuffer");
    LengthPostProbCodaBuffer  = _LengthPostProbCodaBuffer;
    if (iFirstRandom < 0 || RsOnTau == NULL || Rf_isNull(sOnTau) || Rf_length(sOnTau) == 0) {
       LengthPostProb = p;
    } else {
      LengthPostProb = Rf_length(sOnTau) + iFirstRandom;  
    }
    if (LengthPostProb <= 0) {
      Rprintf("ERROR ERROR ERROR ERROR ERROR ERROR ERROR ERROR \n");
      Rprintf("ERROR BayesSpikeCpp.cpp:: set_RsProbBufferFile(): Issue, LengthPostProb = %d!! Big Error \n");      
      Rf_error("Error, PostProbBuffer!");
    }
    LengthWrittenPostProb = 0;  LengthTotalWrittenPostProb = 0;
    //int AGo = 0;
    if (LengthPostProbCodaBuffer * LengthPostProb > 5*TooMuchBuffer) {
      //Rprintf("BayesSpikeCpp.cpp: set_RsProbBufferFile(), Issue: LengthPostProb=%d, ",
      //  LengthPostProb);
      //Rprintf("LengthPostProbCodaBuffer=%d, and!\n",
      //  LengthPostProbCodaBuffer);
      LengthPostProbCodaBuffer = (int) floor(
        ((double) 5)*((double) TooMuchBuffer
        / ((double)LengthPostProb)));
      if (LengthPostProbCodaBuffer <= 5) {
        LengthPostProbCodaBuffer = 5;
      }
      //Rprintf("TooMuchBuffer = %d, setting it back to smaller %d, and continue.\n", TooMuchBuffer); R_FlushConsole();
      //AGo = 1;
    }
    //if (LengthPostProbCodaBuffer * LengthPostProb > 5*TooMuchBuffer)
    //  Rf_error("SetupPostProb, hey, LengthPostProb=%d, LengthPostProbCodaBuffer=%d, probably not good!\n",
    //    LengthPostProb, LengthPostProbCodaBuffer);
    //}
    PostProbCodaBuffer = (double *) Calloc(LengthPostProbCodaBuffer * 
      LengthPostProb, double);
    //if (AGo == 1) {
    //  Rprintf(" -- Still we managed to set PostProbCodaBuffer length \n"); R_FlushConsole();
    //}
    RsPostProbBufferFile = new RRObject(_sPostProbBufferFile);
    NewPostProbWrite = 1;
    return(1);
  }
  //int get_LengthPostProb() { return(LengthPostProb);}
  int get_NewPostProbWrite() { return(NewPostProbWrite); }
  void set_NewPostProbWrite(int _NewPostProbWrite) {
    if (_NewPostProbWrite == 0) {
      NewPostProbWrite = 0;
    } else if (_NewPostProbWrite == 1) {
      NewPostProbWrite = 1;
    } else {
      NewPostProbWrite = 1;
    }
  }
  int get_LengthPostProbCodaBuffer() { return(LengthPostProbCodaBuffer); }
  int get_LengthWrittenPostProb() { return(LengthWrittenPostProb); }
  SEXP get_InitFlags() {
    if (RInitFlags == NULL) {
      return(R_NilValue);
    }
    return(RInitFlags->asSexp());
  }
  // We write to a buffer the posterior probabilities of activation of fixed and random effects.
  SEXP get_CurrentPostProbCodaBuffer() {
    if (LengthPostProbCodaBuffer <= 0 || PostProbCodaBuffer == NULL || LengthPostProb <= 0 || LengthWrittenPostProb <= 0) {
      return(R_NilValue);
    }
    SEXP sOn= R_NilValue;
    Rf_protect(sOn = Rf_allocMatrix(REALSXP, LengthWrittenPostProb, LengthPostProbCodaBuffer));
    int One =1;  
    //int Tot = LengthWrittenPostProb * LengthPostProbCodaBuffer;
    int ii;     
    for (ii = 0; ii < LengthWrittenPostProb; ii++) {
      F77_CALL(dcopy)(&LengthPostProb, PostProbCodaBuffer + ii*LengthPostProb, &One, 
        REAL(sOn)+ii, &LengthWrittenPostProb);
    }
    Rf_unprotect(1); 
    return(sOn);
  }   
  int get_LengthTotalWrittenPostProb() { return(LengthTotalWrittenPostProb); }
  SEXP get_PostProbCodaBuffer() {
    if (LengthPostProbCodaBuffer <= 0 || PostProbCodaBuffer == NULL || LengthPostProb <= 0) {
      return(R_NilValue);
    }
    SEXP sOn= R_NilValue;
    Rf_protect(sOn = Rf_allocMatrix(REALSXP, LengthPostProbCodaBuffer, LengthPostProb));
    int One =1;  
    //int Tot = LengthPostProb * LengthPostProbCodaBuffer;
    int ii;       
    for (ii = 0; ii < LengthPostProbCodaBuffer; ii++) {
      F77_CALL(dcopy)(&LengthPostProb, PostProbCodaBuffer + ii*LengthPostProb, &One, 
        REAL(sOn)+ii, &LengthPostProbCodaBuffer);
    }
    Rf_unprotect(1); 
    return(sOn);
  } 
  int GoFindInOld(double AnEnergy);
  SEXP get_ICoda() {
    if (LengthProbCoda <= 0 || ICoda == NULL) {
      return(R_NilValue);
    }
    SEXP sOn = R_NilValue;  int One = 1;
    Rf_protect(sOn = Rf_allocVector(INTSXP, LengthProbCoda));
    if (Verbose >= 5) {
      Rprintf("get_ICoda, One = %d.\n", One); R_FlushConsole();
    }
    int ii = 0;
    for (ii = 0; ii < LengthProbCoda; ii++) {
      INTEGER(sOn)[ii] = (int) ICoda[ii];
    }
    Rf_unprotect(1); return(sOn);
  }
  int UpdateBFixed();
  void set_NamesBeta(SEXP snBeta) {
     if (Rf_isNull(snBeta)) {
       Rf_error("set_NamesBeta: snBeta is null!\n");
     }
     if (Rf_length(snBeta) != Rf_length(sBeta) || Rf_length(snBeta) != p) {
       Rprintf("set_NamesBeta: No way, snBeta length %d, length Beta = %d!\n",
         Rf_length(snBeta), Rf_length(sBeta)); R_FlushConsole();
        Rf_error("set_NamesBeta, not really. \n");
     }
     
     if (!Rf_isString(snBeta)) {
       Rprintf("set_NamesBeta: snBeta is not a string! \n"); R_FlushConsole();
       Rf_error("set_NamesBeta: Not String. \n");
     }
     if (Verbose >= 1) {
       Rprintf("set_NamesBeta: We're about to run. \n"); R_FlushConsole();
     }
     SEXP MyAllS = R_NilValue;
     Rf_protect(MyAllS = Rf_allocVector(STRSXP, p));
     SEXP sStrings = R_NilValue;
     char *OnName = NULL;
     int ANName = 0;
     for (int iti = 0; iti < p; iti++) {
        OnName = (char*) Calloc(200, char);
        sStrings = STRING_ELT(snBeta, iti);
        ANName = 200;
        //if (Rf_length(sStrings) < 200) {
        //ANName = Rf_length(sStrings);
        //}
        if (ANName >= 1) {
        for (int ss = 0; ss < ANName; ss++) {
          OnName[ss] = (char) CHAR(sStrings)[ss];
          if (CHAR(sStrings)[ss] == '\0') {
            ss = 200; 
            break;
          }
        }
        }
        SET_STRING_ELT(MyAllS, iti, Rf_mkChar(OnName));
        Free(OnName);
     }
     Rf_setAttrib(sBeta, R_NamesSymbol, MyAllS);
     Rf_unprotect(1);
     if (Verbose >= 1) {
       Rprintf("set_NamesBeta: All Set. \n"); R_FlushConsole();
     }
  }
  void set_Beta(SEXP Beta_) {  
    int ii = 0, jj = 0;
    if (Rf_isNull(Beta_) || !Rf_isReal(Beta_)) {
      Rf_error("Error: Beta_ is Null!");
    }
    int Rflen = Rf_length(Beta_);
    if (Rf_length(Beta_) != Rf_length(sBeta)) {
      Rprintf("setBeta Error, Beta_[%d] doesn't have same length as sBeta[%d]",
        Rf_length(Beta_), Rf_length(sBeta)); R_FlushConsole();
      Rflen = Rf_length(Beta_);
    }
    int St = 0;
    if (iFirstRandom < 0 || tauEndList == NULL || Rf_isNull(tauEndList) ||
      Rf_length(tauEndList) <= 0) {
      for (ii = 0; ii < Rflen; ii++) {
        if (REAL(Beta_)[ii] != 0.0 && 
          XLC[ii] < 0) {
          NewFixedCoords++;  AllNewCoords++;
        }
        REAL(sBeta)[ii] = REAL(Beta_)[ii];
      }
    } else {
      if (iFirstRandom > 0) {
        for (ii = 0; ii < iFirstRandom; ii++) {
          if (REAL(Beta_)[ii] != 0.0 && 
            XLC[ii] < 0) {
            NewFixedCoords++;  AllNewCoords++;
          }
          REAL(sBeta)[ii] = REAL(Beta_)[ii];
        }
      }
      St = iFirstRandom;
      for (ii = 0; ii < Rf_length(tauEndList); ii++) {
        for (jj = St; jj <= ANINT(tauEndList, ii); jj++) {
          if (REAL(Beta_)[jj] != 0.0 &&
            XLC[jj] < 0) {
            AllNewCoords +=  ANINT(tauEndList,ii) - St+1;
            if (Verbose >= 2) {
              Rprintf("set_Beta, ii=%d/%d, Beta[jj=%d] has XLC[jj=%d] =%d < 0, add %d coordinates, %d total!\n",
                ii, Rf_length(tauEndList), jj, jj, XLC[jj], ANINT(tauEndList,ii)-St+1,
                AllNewCoords); R_FlushConsole();
            }
            for (; jj <= ANINT(tauEndList,ii); jj++) {
              REAL(sBeta)[jj] = REAL(Beta_)[jj];
              if (REAL(sBeta)[jj] == 0.0) {
                REAL(sBeta)[jj] = 0.0000001;
              }
            }
            if (REAL(sOnTau)[ii] <= 0.0) {
              REAL(sOnTau)[ii] = 1.0;
            }
            jj = ANINT(tauEndList, ii)+1;
            break;  
          } else {
            REAL(sBeta)[jj] = REAL(Beta_)[jj];
          }
        }
        St = ANINT(tauEndList,ii) + 1; 
      }
    }
    
    SEXP MyS = Rf_getAttrib(Beta_, R_NamesSymbol);
    if (!Rf_isNull(MyS) && Rf_length(MyS) == p) {
      SEXP sStrings=R_NilValue;
      char *OnName;
      Rf_protect(sStrings=Rf_allocVector(STRSXP,p));
      for (int iti = 0; iti < p ; iti++) {
        OnName  = (char*) CHAR(STRING_ELT(MyS, iti));
        SET_STRING_ELT(sStrings, iti, Rf_mkChar(OnName));
      }
      Rf_setAttrib(sBeta, R_NamesSymbol, sStrings);
      Rf_unprotect(1);
    }
    
    UpdateBFixed();
    if (AllNewCoords >= 1 && DoAddCoordsOnSetBeta > 0) {
      if (Verbose >= 1) {
        Rprintf("set_Beta: We should add at total of %d New Coords, %d Fixed, OnKappaS=%d/%d\n",
          AllNewCoords, NewFixedCoords, OnKappaS, OnKappaMem); R_FlushConsole();
      }
      AddAllNewCoords();
      if (Verbose >= 1) {
        Rprintf("set_Beta: We have added all of those coords,  OnKappaS=%d/%d",
          OnKappaS, OnKappaMem); R_FlushConsole();
      }

      if (Verbose >= 1) {
        Rprintf("set_Beta: We will run Refresh Orderd Active. \n");
      }
      RefreshOrderedActive(1);
    }  

    SetXtResidNoXtX();  
  }
  void set_FirstRandom(SEXP FirstRandom_) {
    int ti = GetFirstInteger(FirstRandom_);
    if (tauEndList == 0 || Rf_isNull(tauEndList) || Rf_length(tauEndList) <= 0) {
      iFirstRandom = p;
      return;
    }
    if (ti > p || ti < 0) {
      Rf_error("setFirstRandom: Please set first Random in Range!\n");
    }
    if (!Rf_isNull(tauEndList)) {
      if (GetFirstInteger(tauEndList) < ti) {
        Rf_error("Cannot set new FirstRandom[%d]] Past tauEndList[%d]",
          ti, GetFirstInteger(tauEndList));
      }
    }
    iFirstRandom = ti;
  }
  void set_CFirstRandom(int CFirstRandom_) {
    if (sOnTau == 0 || Rf_isNull(sOnTau) || Rf_length(sOnTau) <= 0 ||
      tauEndList == 0 || Rf_isNull(tauEndList) || Rf_length(tauEndList) <= 0) {
      iFirstRandom = p;
      return;
    }
    iFirstRandom = CFirstRandom_;
  }
  int RefreshsBeta(SEXP NewBetaVals) {
    if (Verbose >= 0) {
      Rprintf("--- Refreshing Beta in Cpp \n"); R_FlushConsole();
    }
    if (Rf_isNull(NewBetaVals) || !Rf_isReal(NewBetaVals)) {
      Rf_error("  Error: NewBetaVals is not a valid input");
    }
    int MyLen = Rf_length(NewBetaVals);
    if (p != MyLen) {
      Rf_error("RefreshBeta: Please supply New Beta of same length as old\n");
    }
    if (RsBeta == NULL || sBeta == NULL || Rf_isNull(sBeta)) {
      Rf_error("RefreshBeta: Error, sBeta was never setup it is NULL!\n");
    }
    if (Rf_length(sBeta) == 0) {
      Rf_error("RefreshBeta: error, our length is Zero!\n");
    }
    if (p != Rf_length(sBeta)) {
      Rf_error("RefreshBeta:  Error, sBeta does not have old length of p=%d, but length %d",
        p, Rf_length(sBeta));
    }
    //int One;
    if (Rf_isNull(sBeta)) {
      Rf_error("RefreshBeta: You never supplied a satisfactory first sBeta!");
    }
    if (Verbose > 4) {
      Rprintf("--- Refreshing Beta , performing Copy of p = %d\n", p); R_FlushConsole();
      int ii;
      for (ii = 0; ii < p; ii++) {
        if (ii % 10 == 0) { Rprintf("\n"); }
        Rprintf("%f, ", REAL(sBeta)[ii]); R_FlushConsole();
      }
      Rprintf(" and NewBetaVals: ");
      for (ii = 0; ii < p; ii++) {
        if (ii % 10 == 0) { Rprintf("\n"); }
        Rprintf("%f, ", REAL(NewBetaVals)[ii]); R_FlushConsole();
      }     
      Rprintf("\n  Time for Copy Fun \n"); R_FlushConsole();
      
    }
    int ii;
    if (tauEndList == NULL || Rf_isNull(tauEndList) || Rf_length(tauEndList) <= 0 || iFirstRandom >= p) {
      if (Verbose >= 0) {
        Rprintf("RefreshBeta: About to set to zero. \n"); R_FlushConsole();
      }
    for (ii = 0; ii < p; ii++) {
      //Rprintf(" %d ", ii); 
      //if (ii % 10 == 0) { Rprintf("\n"); R_FlushConsole();  }
      REAL(sBeta)[ii] = REAL(NewBetaVals)[ii];
      if (REAL(NewBetaVals)[ii] == 0.0) {
        BFixed[ii] = 0;
      } else {
        BFixed[ii] = 1;
      }
    }
    } else {
     for (ii = 0; ii < iFirstRandom; ii++) {
      //Rprintf(" %d ", ii); 
      //if (ii % 10 == 0) { Rprintf("\n"); R_FlushConsole();  }
      REAL(sBeta)[ii] = REAL(NewBetaVals)[ii];
      if (REAL(NewBetaVals)[ii] == 0.0) {
        BFixed[ii] = 0;
      } else {
        BFixed[ii] = 1;
      }
     }
     for (ii = iFirstRandom; ii < p; ii++) {
      //Rprintf(" %d ", ii); 
      //if (ii % 10 == 0) { Rprintf("\n"); R_FlushConsole();  }
      REAL(sBeta)[ii] = REAL(NewBetaVals)[ii];
     }   
    }
    double StCalcMe = 0.0;
    int jj,kk;
    if (tauEndList != NULL && !Rf_isNull(tauEndList) && Rf_length(tauEndList) >= 1 &&
      sOnTau != NULL && !Rf_isNull(sOnTau) && Rf_length(sOnTau) >= 1) {
       int St = 0;   //int Go = 0;
       if (iFirstRandom > 0) {
         St = iFirstRandom;
       }
       for (jj = 0; jj < Rf_length(tauEndList); jj++) {
         StCalcMe = 0.0;
         for (kk = St; kk < ANINT(tauEndList, jj)+1; kk++) {
           StCalcMe += REAL(sBeta)[kk] * REAL(sBeta)[kk];
         }
         if (StCalcMe <= 0) {
           REAL(sOnTau)[jj] = 0.0;
         } else {
           REAL(sOnTau)[jj] = StCalcMe / (kk - St);
         }
         St = ANINT(tauEndList,jj)+1;
       }
    }
    if (Verbose >= 2) {
      Rprintf("--- That was one attempt to copy \n"); R_FlushConsole();
    }
    //F77_CALL(dcopy)(&p, REAL(NewBetaVals), &One, REAL(sBeta), &One);
    if (Verbose > 0) {
      Rprintf("--- Now UpdateNonFreshXtResid \n"); R_FlushConsole();
    }
    UpdateNonFreshXtResid();
    return(1);
  }
  void set_tauEndList(SEXP tauEndList_) {
    if (Verbose > 2) {
      Rprintf("set_tauEndList: BayesSpikeCpp: Setting tauEndList"); R_FlushConsole();
    }
    if (tauEndList_ == NULL || Rf_isNull(tauEndList_) ||
      GetFirstInteger(tauEndList_) < 0) {
      if (!Rf_isNull(tauEndList)) {
        Rf_error("Error: set_tauEndList, you try to set an end list null that was not?  weird Stuff!");
      }
      if (Verbose > 2) {
        Rprintf("set_tauEndList: set tauEndList to NULL. \n"); R_FlushConsole();
      }
      iFirstRandom = p;
      DDelete(RtauEndList, "RtauEndList"); 
      tauEndList = R_NilValue; MaxTauList = 0; 
      FFree(SmallXtResid, "SmallXtResid");
      FFree(SmallRVec, "SmallRVec"); return;
    }
    if (Rf_isNull(tauEndList_)) {
      if (Verbose > 0) {
        Rprintf("set_tauEndList warning: turning tauEndList off \n"); R_FlushConsole();
      }
      DDelete(RtauEndList, "RtauEndList"); tauEndList = R_NilValue;
      FFree(SmallXtResid, "SmallXtResid");
      FFree(SmallRVec, "SmallRVec");
      MaxTauList = 0;
      return;
    } else if (Rf_isNull(tauEndList)) {
      DDelete(RtauEndList, "RtauEndList");
      if (Rf_isReal(tauEndList_)) {
        SEXP sNT = R_NilValue;

        Rf_protect(sNT = Rf_allocVector(INTSXP, Rf_length(tauEndList_)));
        int iti; for (iti = 0; iti < Rf_length(tauEndList_); iti++) {
          INTEGER(sNT)[iti] = (int) REAL(tauEndList_)[iti];
        }
        RtauEndList = new AObject(sNT);
        tauEndList = RtauEndList->asSexp();  Rf_unprotect(1);
      } else {
        RtauEndList = new AObject(tauEndList_);
        tauEndList = RtauEndList->asSexp();
      }
    } else {
      DDelete(RtauEndList, "RtauEndList");
      if (Rf_isReal(tauEndList_)) {
        SEXP sNT = R_NilValue;
      
        Rf_protect(sNT = Rf_allocVector(INTSXP, Rf_length(tauEndList_)));
        int iti; for (iti = 0; iti < Rf_length(tauEndList_); iti++) {
          INTEGER(sNT)[iti] = (int) REAL(tauEndList_)[iti];
        }
        RtauEndList = new AObject(sNT);
        tauEndList = RtauEndList->asSexp();  Rf_unprotect(1);
      } else {
        RtauEndList = new AObject(tauEndList_);
        tauEndList = RtauEndList->asSexp();
      }
    }
    if (Verbose > 2) {
      Rprintf("set_tauEndList: At End, here is tauEndList: "); R_FlushConsole();
      for (int ii = 0; ii < Rf_length(tauEndList); ii++) {
        Rprintf(" %d", INTEGER(tauEndList)[ii]);
        if (ii < Rf_length(tauEndList)-1) { Rprintf(", "); } else {
          Rprintf("\n"); 
        }
        R_FlushConsole();
      }
    }
    if (!Rf_isNull(Rf_getAttrib(tauEndList_, R_NamesSymbol))) {
      SEXP sStrings = R_NilValue;
      char OnName[400];
      Rf_protect(sStrings=Rf_allocVector(STRSXP,Rf_length(tauEndList)));
      for (int iti = 0; iti < Rf_length(tauEndList); iti++) {
        for (int jtj = 0; jtj < 399; jtj++) {
          OnName[jtj] = '\0';
        }
        strcpy((char*)CHAR(STRING_ELT(Rf_getAttrib(tauEndList_, R_NamesSymbol), iti)),
          OnName);
        SET_STRING_ELT(sStrings, iti, Rf_mkChar(OnName));
      } 
      Rf_setAttrib(tauEndList, R_NamesSymbol, sStrings);
      Rf_unprotect(1);
    }
    int St; int M = 0;  int ii = 0; int D = 0;
    int StiFirstRandom = 0;
    if (iFirstRandom > 0 && iFirstRandom < p) {
      StiFirstRandom = iFirstRandom;
    }
    if (tauEndList != NULL && !Rf_isNull(tauEndList) 
      && Rf_length(tauEndList) > 0) {
      for (ii = 0; ii < Rf_length(tauEndList); ii++) {
        if (ii == 0) { St = StiFirstRandom; } else {
          St = ANINT(tauEndList, (ii-1));
        }
        D = ANINT(tauEndList, ii) - St+1;
        if (D > M) { M = D;}
      }
    }
    MaxTauList = M;
    FFree(SmallXtResid, "SmallXtResid");
    FFree(SmallRVec, "SmallRVec");
    RMemGetD(SmallXtResid, "SmallXtResid", MaxTauList+1); 
    RMemGetD(SmallRVec, "SmallRVec", MaxTauList+1);
    if (this->Verbose >= 2) {
      Rprintf("-- set TauEndList, setup tauEnd list with MaxTauList = %d\n",
      MaxTauList);  R_FlushConsole();
    }
  }
  void set_OnTau(SEXP sOnTau_) {
    if (Rf_isNull(tauEndList) || Rf_length(tauEndList) <= 0) {
      if ( (sOnTau == NULL || sOnTau == R_NilValue) && Rf_isNull(sOnTau_)) {
         sOnTau = R_NilValue; DDelete(RsOnTau, "RsOnTau");
         return;
      }
      Rf_error("Before Setting OnTau, Please set tauEndList");
    }
    if (Rf_isNull(sOnTau)) {
      DDelete(RsOnTau, "RsOnTau");  RsOnTau = new AObject(sOnTau_);
      sOnTau = RsOnTau->asSexp();
    } else {
      if (Rf_length(sOnTau_) != Rf_length(tauEndList)) {
        Rf_error("set_OnTau, only set to same length of tauEndList = %d", Rf_length(tauEndList));        
      }
      DDelete(RsOnTau, "RsOnTau");  RsOnTau = new AObject(sOnTau_);
      sOnTau = RsOnTau->asSexp();
    }
    if (!Rf_isNull(Rf_getAttrib(sOnTau_, R_NamesSymbol))) {
      SEXP sStrings = R_NilValue;
      char* OnName;
      Rf_protect(sStrings=Rf_allocVector(STRSXP,Rf_length(sOnTau)));
      for (int iti = 0; iti < Rf_length(sOnTau); iti++) {
        OnName = (char*) CHAR(STRING_ELT(Rf_getAttrib(sOnTau_, R_NamesSymbol), iti));
        SET_STRING_ELT(sStrings, iti, Rf_mkChar(OnName));
      } 
      Rf_setAttrib(sOnTau, R_NamesSymbol, sStrings);
      Rf_unprotect(1);
    }
  }
 void SetupPiACodaBuffer(SEXP RsPiACodaFile_, SEXP LengthPiACodaBuffer_) {
    if (RsPiACodaFile_ == NULL || Rf_isNull(RsPiACodaFile_)) {
      if (Verbose >= 0) {
        Rprintf("BayesSpikeCpp.h: NULL Names inputted: SetupPiACodaFile: Deleting PiACodaBuffer\n"); R_FlushConsole();
      }
      DDelete(RsPiACodaFile, "RsPiACodaFile"); FFree(PiACodaBuffer, "PiACodaBuffer");
      return;
    }
    if (Verbose >= -1) {
      Rprintf("BayesSpikeCpp.h: SetupPiACodaBuffer: Testing RsPiACodaFile is a String. \n");
      R_FlushConsole();
    }
    if (!Rf_isString(RsPiACodaFile_)) {
      Rf_error("SetupPiACodaBuffer: Error, RsYFile is not a string!\n");
    }
    //if (dfRobit >= 0) {
    //  Rf_error("SetupPiACodaBuffer: Why copy PiA for binomial data?\n"); 
    //}
    if (Rf_isNull(LengthPiACodaBuffer_) || LengthPiACodaBuffer_ == NULL ||
      (!Rf_isInteger(LengthPiACodaBuffer_) && !Rf_isReal(LengthPiACodaBuffer_)) ||
      GetFirstInteger(LengthPiACodaBuffer_) <= 0) {
      Rf_error("SetupPiACodaBuffer: Supplied LengthYBuffer_ is invalid!\n");  
    }
    if (GetFirstInteger(LengthPiACodaBuffer_) <= 0) {
      Rprintf("SetupPiACodaBuffer: Hey, you gave LengthPiACodaBuffer_ = %d!!!\n",
        GetFirstInteger(LengthPiACodaBuffer_)); R_FlushConsole();
    }
    if (Verbose >= -1) {
      Rprintf("BayesSpikeCpp.h:SetupPiACodaBuffer, we have LengthPiACodaBuffer is valid is = %d. \n",
        GetFirstInteger(LengthPiACodaBuffer_));
    }
    if (Verbose >= 2) {
      Rprintf("SetupPiACodaBuffer, Given LengthPiACodaBuffer = %d\n",
        GetFirstInteger(LengthPiACodaBuffer_)); R_FlushConsole();
    }
    if (GetFirstInteger(LengthPiACodaBuffer_) <= 0) {
      LengthPiACodaBuffer = 100;
    } else {
      LengthPiACodaBuffer = GetFirstInteger(LengthPiACodaBuffer_);
    }
    int LsOnPiA = 1;
    if (Rf_length(sOnPiA) > 1) {
      LsOnPiA = 2;
    }
    if (Verbose >= -1) {
      Rprintf("BayesSpikeCpp.h:SetupPiACodaBuffer: About to delete and reattach RsPiACodaFile"); R_FlushConsole();
    }

    DDelete(RsPiACodaFile, "RsPiACodaFile");
    if (Verbose >= -1) {
      Rprintf("BayesSpikeCpp.h:SetupPiACodaBuffer: Now for a new attachment of RsPiACodaFile"); R_FlushConsole();
    }
    RsPiACodaFile = new RRObject(RsPiACodaFile_);
    if (RsPiACodaFile == NULL) {
      Rprintf("SetupPiACodaBuffer: We did not allocate RsPiACodaFile_ into AObject\n"); R_FlushConsole();
      Rf_error("SetupPiACodaBuffer: Error. \n");
    }
    if (Verbose >= -1) {
      Rprintf("BayesSpikeCpp.hSetupPiACodaBufffer: Checking for reasonable Length"); R_FlushConsole();
    } 
    if (LsOnPiA * LengthPiACodaBuffer >= 2*TooMuchBuffer) {
      Rprintf("Hey Error SetupPiACodaBuffer: LsOnPiA=%d, LengthPiACodaBuffer=%d \n",
        LsOnPiA, LengthPiACodaBuffer); R_FlushConsole();
      LengthPiACodaBuffer = (int) floor( ((double)2.0* TooMuchBuffer) / ((double)LsOnPiA) );
      if (LengthPiACodaBuffer <= 5) {
         LengthPiACodaBuffer = 5;
      }
    }
    //if (LsOnPiA * LengthPiACodaBuffer > 2*TooMuchBuffer) {
    //  Rf_error("Hey Error SetupPiACodaBuffer: LsOnPiA=%d, LengthPiACodaBuffer=%d \n",
    //    LsOnPiA, LengthPiACodaBuffer); R_FlushConsole();
    //}
    if (Verbose >= -1) {
      Rprintf("BayesSpikeCpp.h:SetupPiACodaBuffer: About to FFree current PiACodaBuffer. \n"); R_FlushConsole();
    }
    FFree(PiACodaBuffer, "PiACodaBuffer");
    if (Verbose >= -1) {
      Rprintf("BayesSpikeCpp.h:SetupPiACodaBuffer: Allocating PiACodaBuffer length %d. \n",
        LsOnPiA * LengthPiACodaBuffer+4); R_FlushConsole();
    }
    PiACodaBuffer = (double*)Calloc(LsOnPiA * 
      LengthPiACodaBuffer+4, double);
    if (PiACodaBuffer == NULL) {
      Rprintf("SetupPiACodaBuffer: Error MEMISSUE: PiACodaBuffer is now NULL, it was not allocated. \n");
      R_FlushConsole();
      NewWritePiACodaBuffer = -1;
      return;
    }
    if (Verbose >= -1) {
      Rprintf("BayesSpikeCpp.h:SetupPiACodaBuffer, PiACodaBuffer Calloced, onto setting up parts. \n"); R_FlushConsole();
    }
    NewWritePiACodaBuffer = 1;  LengthWrittenPiACodaBuffer = 0;
    LengthTotalWrittenPiACodaBuffer = 0;
    if (Verbose > 1) {
      Rprintf("SetupPiACodaBuffer: LengthPiACodaBuffer = %d\n", LengthPiACodaBuffer);
      R_FlushConsole();
      Rprintf("SetupPiACodaBuffer, Setup PiACodaBuffer to File %s, with length %d\n",
        CHAR(STRING_ELT(RsPiACodaFile->asSexp(), 0)), LengthPiACodaBuffer);
    }
    if (Verbose >= -1) {
      Rprintf("BayesSpikeCpp.h:SetupPiACodaBuffer: Successful setup of PiACodaBuffer. \n"); R_FlushConsole();
    }
    return;
  }
  void set_NewPiACodaBufferWrite(int ANW) {
    if (ANW == 0) { NewWritePiACodaBuffer = 0;}
    NewWritePiACodaBuffer = 1;
  }
  int get_NewPiACodaBufferWrite() { return(NewWritePiACodaBuffer); }
  SEXP get_RsPiACodaFile() {
    if (RsPiACodaFile == NULL) { return(R_NilValue); }
    return(RsPiACodaFile->asSexp());
  }
  int get_LengthTotalWrittenPiACodaBuffer() { return(LengthTotalWrittenPiACodaBuffer); }
  int get_LengthPiACodaBuffer() { return(LengthPiACodaBuffer); }
  int get_LengthWrittenPiACodaBuffer() { return(LengthWrittenPiACodaBuffer); }
  SEXP get_CurrentPiACodaBuffer() {
    if (PiACodaBuffer == NULL) { Rf_error("Hey: PiACodaBuffer is not setup!\n");}
    if (LengthWrittenPiACodaBuffer <= 0) { Rprintf("LengthWrittenPiACodaBuffer = %d\n",
      LengthWrittenPiACodaBuffer); return(R_NilValue); }
    if (LengthPiACodaBuffer <= 0) { 
      Rf_error("Hey: LengthPiACodaBuffer = %d\n", LengthPiACodaBuffer);}
    SEXP sOut = R_NilValue;
    int Tot = LengthWrittenPiACodaBuffer;
    if (Rf_length(sOnPiA) == 2) {
      Rf_protect(sOut = Rf_allocMatrix(REALSXP, 2, LengthWrittenPiACodaBuffer));
      Tot = LengthWrittenPiACodaBuffer*2;
    } else {
      Rf_protect(sOut = Rf_allocVector(REALSXP, LengthWrittenPiACodaBuffer));    
    }
    int One = 1;
    F77_CALL(dcopy)(&Tot, PiACodaBuffer, &One, REAL(sOut), &One);
    Rf_unprotect(1); return(sOut);
  }
  SEXP get_AllPiACodaBuffer() {
    if (PiACodaBuffer == NULL) { Rf_error("Hey: PiACodaBuffer is not setup!\n");}
    if (LengthWrittenPiACodaBuffer < 0) { Rprintf("LengthWrittenPiACodaBuffer = %d\n",
      LengthWrittenPiACodaBuffer); return(R_NilValue); }
    if (LengthPiACodaBuffer <= 0) { 
      Rf_error("Hey: LengthPiACodaBuffer = %d\n", LengthPiACodaBuffer);}
    SEXP sOut = R_NilValue;
    int Tot = LengthPiACodaBuffer;
    if (Rf_length(sOnPiA) == 2) {
      Rf_protect(sOut = Rf_allocMatrix(REALSXP, 2, LengthPiACodaBuffer));
      Tot = LengthPiACodaBuffer*2;
    } else {
      Rf_protect(sOut = Rf_allocVector(REALSXP, LengthPiACodaBuffer));    
    }
    int One = 1;
    F77_CALL(dcopy)(&Tot, PiACodaBuffer, &One, REAL(sOut), &One);
    Rf_unprotect(1); return(sOut);
  }
  void SetupSigCodaBuffer(SEXP RsSigCodaFile_, SEXP LengthSigCodaBuffer_) {
    if (Verbose >= 0) {
      Rprintf("BayesSpikeCpp.h: SetupSigCodaBuffer() start. \n"); R_FlushConsole();
    }
    if (RsSigCodaFile_ == NULL || Rf_isNull(RsSigCodaFile_)) {
      if (Verbose > 0) {
        Rprintf("SetupSigCodaFile: Deleting SigCodaBuffer\n"); R_FlushConsole();
      }
      DDelete(RsSigCodaFile, "RsSigCodaFile"); FFree(SigCodaBuffer, "SigCodaBuffer");
      return;
    }
    if (!Rf_isString(RsSigCodaFile_)) {
      Rf_error("SetupSigCodaBuffer: Error, RsYFile is not a string!\n");
    }
    if (dfRobit >= 0 && AlterWeightdfRobit < 0) {
      Rf_error("SetupSigCodaBuffer: Why copy sigma for binomial data?\n"); 
    }
    if (Rf_isNull(LengthSigCodaBuffer_) || LengthSigCodaBuffer_ == NULL ||
      (!Rf_isInteger(LengthSigCodaBuffer_) && !Rf_isReal(LengthSigCodaBuffer_)) ||
      GetFirstInteger(LengthSigCodaBuffer_) <= 0) {
      Rf_error("SetupSigCodaBuffer: Supplied LengthYBuffer_ is invalid!\n");  
    }
    if (LengthSigCodaBuffer_==NULL || Rf_isNull(LengthSigCodaBuffer_) ||
      (!Rf_isReal(LengthSigCodaBuffer_) && !Rf_isInteger(LengthSigCodaBuffer_))
      || GetFirstInteger(LengthSigCodaBuffer_) <= 0) {
      Rprintf("BayesSpikeCpp.h: SetupSigCodaBuffer(): Hey, something is wrong with LengthSigCodaBuffer!\n");
      R_FlushConsole();
      if (!Rf_isNull(LengthSigCodaBuffer_) && Rf_length(LengthSigCodaBuffer_) >= 1) {
        Rprintf("BayesSpikeCpp.h: SetupSigCodaBuffer: Hey, you gave LengthSigCodaBuffer_ = %d!!!\n",
          GetFirstInteger(LengthSigCodaBuffer_)); R_FlushConsole();
      } else {
        Rprintf("BayesSpikeCpp.h: SetupSigCodaBuffer: you set it to null. \n");
      }
    }
    if (Verbose >= 2) {
      Rprintf("SetupSigCodaBuffer, Given LengthSigCodaBuffer = %d\n",
        GetFirstInteger(LengthSigCodaBuffer_)); R_FlushConsole();
    }
    if (LengthSigCodaBuffer_==NULL || Rf_isNull(LengthSigCodaBuffer_) ||
      (!Rf_isReal(LengthSigCodaBuffer_) && !Rf_isInteger(LengthSigCodaBuffer_))
      || GetFirstInteger(LengthSigCodaBuffer_) <= 0) {
      LengthSigCodaBuffer = 100;
    } else {
      LengthSigCodaBuffer = GetFirstInteger(LengthSigCodaBuffer_);
    }
    DDelete(RsSigCodaFile, "RsSigCodaFile");
    RsSigCodaFile = new RRObject(RsSigCodaFile_);
   
    if (LengthSigCodaBuffer > TooMuchBuffer) {
      Rprintf("Hey: SigCodaBuffer = %d, probably too long!\n", SigCodaBuffer);
      LengthSigCodaBuffer = TooMuchBuffer;
      if (LengthSigCodaBuffer <= 2) {
        LengthSigCodaBuffer = 4;
      }
    }
    //if (LengthSigCodaBuffer > TooMuchBuffer) {    
    //  Rf_error("Hey: SigCodaBuffer=%d is too long!\n", SigCodaBuffer);
    //} 
    SigCodaBuffer = (double*)Calloc(LengthSigCodaBuffer+2, double);
    if (SigCodaBuffer == NULL) {
      Rprintf("Hey: SigCodaBuffer returned NULL. We have a mem issue error, length attempt was %d\n",
        LengthSigCodaBuffer);
      Rprintf("  -- No mem setup for file %s \n",
        CHAR(STRING_ELT(RsSigCodaFile->asSexp(), 0)));
      DDelete(RsSigCodaFile, "RsSigCodaFile");
      LengthSigCodaBuffer = 0;
    }

    LengthWrittenSigCodaBuffer = 0; NewWriteSigCodaBuffer = 1; 
    LengthTotalWrittenSigCodaBuffer = 0;
    if (Verbose >= -1) {
      Rprintf("SetupSigCodaBuffer, LengthSigCodaBuffer = %d", LengthSigCodaBuffer);
      R_FlushConsole();
      Rprintf("SetupSigCodaBuffer, Setup SigCodaBuffer to File %s, with length %d\n",
        CHAR(STRING_ELT(RsSigCodaFile->asSexp(), 0)), LengthSigCodaBuffer);
    }
  }
  void set_NewSigCodaBufferWrite(int ANW) {
    if (ANW == 0) { NewWriteSigCodaBuffer = 0;}
    NewWriteSigCodaBuffer = 1;
  }
  int get_NewSigCodaBufferWrite() { return(NewWriteSigCodaBuffer); }
  SEXP get_RsSigCodaFile() {
    if (RsSigCodaFile == NULL) { return(R_NilValue); }
    return(RsSigCodaFile->asSexp());
  }
  int get_LengthTotalWrittenSigCodaBuffer() { return(LengthTotalWrittenSigCodaBuffer); }
  int get_LengthSigCodaBuffer() { return(LengthSigCodaBuffer); }
  int get_LengthWrittenSigCodaBuffer() { return(LengthWrittenSigCodaBuffer); }
  SEXP get_CurrentSigCodaBuffer() {
    if (SigCodaBuffer == NULL) { Rf_error("Hey: SigCodaBuffer is not setup!\n");}
    if (LengthWrittenSigCodaBuffer <= 0) { Rprintf("LengthWrittenBuffer = %d\n",
      LengthWrittenSigCodaBuffer); return(R_NilValue); }
    if (LengthSigCodaBuffer <= 0) { 
      Rf_error("Hey: LengthSigCodaBuffer = %d\n", LengthSigCodaBuffer);}
    SEXP sOut = R_NilValue;
    Rf_protect(sOut = Rf_allocVector(REALSXP, LengthWrittenSigCodaBuffer));
    int Tot = LengthWrittenSigCodaBuffer; int One = 1;
    F77_CALL(dcopy)(&Tot, SigCodaBuffer, &One, REAL(sOut), &One);
    Rf_unprotect(1); return(sOut);
  }
  SEXP get_AllSigCodaBuffer() {
    if (SigCodaBuffer == NULL) { Rf_error("Hey: SigCodaBuffer is not setup!\n");}
    if (LengthWrittenSigCodaBuffer < 0) { Rprintf("LengthWrittenSigCodaBuffer = %d\n",
      LengthWrittenSigCodaBuffer); return(R_NilValue); }
    if (LengthSigCodaBuffer <= 0) { 
      Rf_error("Hey: LengthSigCodaBuffer = %d\n", LengthSigCodaBuffer);}
    SEXP sOut = R_NilValue;
    Rf_protect(sOut = Rf_allocVector(REALSXP, LengthSigCodaBuffer));
    int Tot = LengthSigCodaBuffer; int One = 1;
    F77_CALL(dcopy)(&Tot, SigCodaBuffer, &One, REAL(sOut), &One);
    Rf_unprotect(1); return(sOut);
  }
  
  void SetupYBuffer(SEXP RsYFile_, SEXP LengthYBuffer_) {
    if (RsYFile_ == NULL || Rf_isNull(RsYFile_)) {
      if (Verbose > 0) {
        Rprintf("SetupYFile: Deleting YBuffer\n"); R_FlushConsole();
      }
      DDelete(RsYFile, "RsYFile"); FFree(YBuffer, "YBuffer");
      return;
    }
    if (!Rf_isString(RsYFile_)) {
      Rf_error("SetupYFile: Error, RsYFile is not a string!\n");
    }
    if (dfRobit < 0) {
      Rf_error("SetupYFile: Why copy Y for regular regression?\n"); 
    }
    if (Rf_isNull(LengthYBuffer_) || LengthYBuffer_ == NULL ||
      (!Rf_isInteger(LengthYBuffer_) && !Rf_isReal(LengthYBuffer_)) ||
      GetFirstInteger(LengthYBuffer_) <= 0) {
      Rf_error("SetupYFile: Supplied LengthYBuffer_ is invalid!\n");  
    }
    if (GetFirstInteger(LengthYBuffer_) <= 0) {
      Rprintf("SetupYFile: Hey, you gave LengthYBuffer_ = %d!!!\n",
        GetFirstInteger(LengthYBuffer_)); R_FlushConsole();
    }
    if (Verbose >= 2) {
      Rprintf("SetupYFile, Given LengthYBuffer = %d\n",
        GetFirstInteger(LengthYBuffer_)); R_FlushConsole();
    }
    if (GetFirstInteger(LengthYBuffer_) <= 0) {
      LengthYBuffer = 100;
    } else {
      LengthYBuffer = GetFirstInteger(LengthYBuffer_);
    }
    DDelete(RsYFile, "RsYFile");
    RsYFile = new AObject(RsYFile_);
    
    if (LengthYBuffer * n >= 4 * TooMuchBuffer) {
      Rprintf("Hey: LengthYBuffer = %d is too long for current memory!\n",
        LengthYBuffer);  R_FlushConsole();
      LengthYBuffer = (int) floor((double) 4.0 * ((double)TooMuchBuffer) / ((double)n));
      if (LengthYBuffer <= 2) { LengthYBuffer = 4; }
    }
    //if (LengthYBuffer * n > 3* TooMuchBuffer) {
    //  Rprintf("Hey: LengthYBuffer = %d is too long for current memory!\n",
    //    LengthYBuffer);  R_FlushConsole();
    //  Rf_error("Hey: LengthYBuffer=%d, is too long.\n", LengthYBuffer);
    //}
    YBuffer = (double*)Calloc(n * LengthYBuffer, double);

    LengthWrittenYBuffer = 0;  NewYBufferWrite = 1; 
    LengthTotalWrittenYBuffer = 0;
    if (Verbose > 1) {
      Rprintf("SetupYBuffer, Setup YBuffer to File %s, with length %d\n",
        CHAR(STRING_ELT(RsYFile->asSexp(), 0)), LengthYBuffer);
    }
  }
  void set_NewYBufferWrite(int ANW) {
    if (ANW == 0) { NewYBufferWrite = 0;}
  }
  int get_NewYBufferWrite() { return(NewYBufferWrite); }
  SEXP get_RsYFile() {
    if (RsYFile == NULL) { return(R_NilValue); }
    return(RsYFile->asSexp());
  }
  int get_LengthTotalWrittenYBuffer() { return(LengthTotalWrittenYBuffer); }
  int get_LengthYBuffer() { return(LengthYBuffer); }
  int get_LengthWrittenYBuffer() { return(LengthWrittenYBuffer); }
  SEXP get_CurrentYBuffer() {
    if (YBuffer == NULL) { Rf_error("Hey: YBuffer is not setup!\n");}
    if (LengthWrittenYBuffer <= 0) { Rprintf("LengthWrittenYBuffer = %d\n",
      LengthWrittenYBuffer); return(R_NilValue); }
    if (LengthYBuffer <= 0) { 
      Rf_error("Hey: LengthYBuffer = %d\n", LengthYBuffer);}
    SEXP sOut = R_NilValue;
    Rf_protect(sOut = Rf_allocMatrix(REALSXP, n, LengthWrittenYBuffer));
    int Tot = n * LengthWrittenYBuffer; int One = 1;
    F77_CALL(dcopy)(&Tot, YBuffer, &One, REAL(sOut), &One);
    Rf_unprotect(1); return(sOut);
  }
  SEXP get_AllYBuffer() {
    if (YBuffer == NULL) { Rf_error("Hey: YBuffer is not setup!\n");}
    if (LengthWrittenYBuffer < 0) { Rprintf("LengthWrittenYBuffer = %d\n",
      LengthWrittenYBuffer); return(R_NilValue); }
    if (LengthYBuffer <= 0) { 
      Rf_error("Hey: LengthYBuffer = %d\n", LengthYBuffer);}
    SEXP sOut = R_NilValue;
    Rf_protect(sOut = Rf_allocMatrix(REALSXP, n, LengthYBuffer));
    int Tot = n * LengthYBuffer; int One = 1;
    F77_CALL(dcopy)(&Tot, YBuffer, &One, REAL(sOut), &One);
    Rf_unprotect(1); return(sOut);
  }
 void SetupWeightBuffer(SEXP RsWeightFile_, SEXP LengthWeightBuffer_) {
    if (RsWeightFile_ == NULL || Rf_isNull(RsWeightFile_)) {
      if (Verbose > 0) {
        Rprintf("SetupYFile: Deleting WeightBuffer\n"); R_FlushConsole();
      }
      DDelete(RsWeightFile, "RsWeightFile"); FFree(WeightBuffer, "WeightBuffer");
      return;
    }
    if (!Rf_isString(RsWeightFile_)) {
      Rf_error("SetupWeightFile: Error, RsYFile is not a string!\n");
    }
    if (iiWeight == NULL) {
      Rf_error("SetupWeightFile: Why copy Weight when no weights!\n"); 
    }
    if (Rf_isNull(LengthWeightBuffer_) || LengthWeightBuffer_ == NULL ||
      (!Rf_isInteger(LengthWeightBuffer_) && !Rf_isReal(LengthWeightBuffer_)) ||
      GetFirstInteger(LengthWeightBuffer_) <= 0) {
      Rf_error("SetupWeightFile: Supplied LengthWeightBuffer_ is invalid!\n");  
    }
    LengthWeightBuffer = GetFirstInteger(LengthWeightBuffer_);
    DDelete(RsWeightFile, "RsWeightFile");
    RsWeightFile = new AObject(RsWeightFile_);
    if (LengthWeightBuffer * n >= 4 * TooMuchBuffer) {
      Rprintf("Hey: LengthWeightBuffer=%d is probably too long\n", LengthWeightBuffer);
      LengthWeightBuffer = (int) floor(4.0 * ((double)TooMuchBuffer) / ((double)n));
      if (LengthWeightBuffer <= 4) {
        LengthWeightBuffer = 4;
      }
    }
    //if (LengthWeightBuffer * n >= 4 * TooMuchBuffer) {    
    //  Rf_error("Make LengthWeightBuffer=%d Shorter!\n", LengthWeightBuffer);
    //}
    WeightBuffer = (double*)Calloc(n * LengthWeightBuffer, double);
    LengthWrittenWeightBuffer = 0;   
    NewWeightBufferWrite = 1;  LengthTotalWrittenWeightBuffer = 0;
    if (Verbose > 1) {
      Rprintf("SetupWeightBuffer, Setup WeightBuffer to File %s, with length %d\n",
        CHAR(STRING_ELT(RsWeightFile->asSexp(), 0)), LengthWeightBuffer);
    }
  }
  int get_LengthTotalWrittenWeightBuffer() {
    return(LengthTotalWrittenWeightBuffer);
  }
  SEXP get_SPrIChol() {
    if (rIChol == NULL || LastStF < 0 || LastStF >= p || LastGoFor <= 0) {
      Rprintf("get_SPrIChol: Sorry, LastStF=%d, LastGoFor=%d but no setup enough for rIChol.\n", LastStF, LastGoFor);
    }
    SEXP sOn = R_NilValue;
    Rf_protect(sOn=Rf_allocMatrix(REALSXP, LastGoFor, LastGoFor));
    int O3 = LastGoFor * LastGoFor;
    for (int jj = 0; jj < O3; jj++) { 
      REAL(sOn)[jj] = 0.0;
    }
    O3 = 0;
    for (int jj = 0; jj < LastGoFor; jj++) {
      for (int ii = jj; ii < LastGoFor; ii++) {
        REAL(sOn)[jj * LastGoFor + ii] = rIChol[O3]; 
        O3++;
      }
    }
    Rf_unprotect(1);  return(sOn);
  }
  SEXP get_SPrI() {
    if (rI == NULL || LastStF < 0 || LastStF >= p || LastGoFor <= 0) {
      Rprintf("get_SPrIChol: Sorry, LastStF=%d, LastGoFor=%d but no setup enough for rIChol.\n", LastStF, LastGoFor);
    }
    SEXP sOn = R_NilValue;
    Rf_protect(sOn=Rf_allocMatrix(REALSXP, LastGoFor, LastGoFor));
    int O3 = 0;
    for (int ii = 0; ii < LastGoFor; ii++) {
      for (int jj = ii; jj < LastGoFor; jj++) {
        REAL(sOn)[ii * LastGoFor + jj] = rI[O3];
        REAL(sOn)[jj * LastGoFor + ii] = rI[O3]; 
        O3++;
      }
    }
    Rf_unprotect(1); 
    return(sOn);
  }
  SEXP get_RsWeightFile() {
    if (RsWeightFile == NULL) { return(R_NilValue); }
    return(RsWeightFile->asSexp());
  }
  void set_NewWeightBufferWrite(int ANW) {
    if (ANW == 0) { NewWeightBufferWrite = 0;}
  }
  int get_NewWeightBufferWrite() { return(NewWeightBufferWrite); }
  int get_LengthWeightBuffer() { return(LengthWeightBuffer); }
  int get_LengthWrittenWeightBuffer() { return(LengthWrittenWeightBuffer); }
  SEXP get_CurrentWeightBuffer() {
    if (WeightBuffer == NULL) { Rf_error("Hey: YBuffer is not setup!\n");}
    if (LengthWrittenWeightBuffer <= 0) { Rprintf("LengthWrittenWeightBuffer = %d\n",
      LengthWrittenWeightBuffer); return(R_NilValue); }
    if (LengthWeightBuffer <= 0) { 
      Rf_error("Hey: LengthYBuffer = %d\n", LengthYBuffer);}
    SEXP sOut = R_NilValue;
    Rf_protect(sOut = Rf_allocMatrix(REALSXP, n, LengthWrittenWeightBuffer));
    int Tot = n * LengthWrittenWeightBuffer; int One = 1;
    F77_CALL(dcopy)(&Tot, WeightBuffer, &One, REAL(sOut), &One);
    Rf_unprotect(1); return(sOut);
  }
  int get_MaximumAllocation() { return(MaximumAllocation); }
  void set_MaximumAllocation(int iIn) {
    if (iIn > 0 && iIn >= OnKappaMem) {
      MaximumAllocation = iIn;  return;
    } else {
      Rprintf("set_MaximumAllocation, OnKappaMem = %d, we cannot lower MaximumAllocation below to %d with %d\n",
        OnKappaMem, MaximumAllocation, iIn);
      MaximumAllocation = OnKappaMem;  R_FlushConsole();
      return;
    }
  }
  SEXP get_AllWeightBuffer() {
    if (WeightBuffer == NULL) { Rf_error("Hey: WeightBuffer is not setup!\n");}
    if (LengthWrittenWeightBuffer < 0) { Rprintf("LengthWrittenWeightBuffer = %d\n",
      LengthWrittenWeightBuffer); return(R_NilValue); }
    if (LengthWeightBuffer <= 0) { 
      Rf_error("Hey: LengthWeightBuffer = %d\n", LengthWeightBuffer);}
    SEXP sOut = R_NilValue;
    Rf_protect(sOut = Rf_allocMatrix(REALSXP, n, LengthWeightBuffer));
    int Tot = n * LengthWeightBuffer; int One = 1;
    F77_CALL(dcopy)(&Tot, WeightBuffer, &One, REAL(sOut), &One);
    Rf_unprotect(1); return(sOut);
  }
  double get_AlterWeightTemperature() {
    return(AlterWeightTemperature);
  }
  int get_AlterInitRun() { return(AlterInitRun); }
  void set_AlterInitRun(int iIn) {
    if (iIn <= 0) {AlterInitRun = 0; return;}
    if (iIn >= 1) {AlterInitRun = 1; return;}
  }
  int get_MaxLargeSizeSquares() {return(MaxLargeSizeSquares); }
  void set_AlterWeightTemperature(double aWNT) {
    if (aWNT > 0) {
      AlterWeightTemperature = aWNT;
    } else if (aWNT == 0) {
      AlterWeightTemperature = 1.0;
    } else if (aWNT <0 )  {
      Rprintf("AlterWeightTemperature?  Why, Sigma = %d!\n");
      R_FlushConsole();
    }
    Temperature = AlterWeightTemperature;
    return;
  }
 void SetupAlterWeightBuffer(SEXP RsAlterWeightFile_, 
   SEXP LengthAlterWeightBuffer_, SEXP AlterWeightdfTNoise_, SEXP AlterWeightdfRobit_) {
   if (Verbose >= 1) {
     Rprintf("Setting up AlterWeightBuffer first time. \n"); R_FlushConsole();
   }
    if (RsAlterWeightFile_ == NULL || Rf_isNull(RsAlterWeightFile_)) {
      if (Verbose > 0) {
        Rprintf("SetupYFile: Deleting AlterWeightBuffer\n"); R_FlushConsole();
      }
      DDelete(RsAlterWeightFile, "RsAlterWeightFile"); FFree(AlterWeightBuffer, "AlterWeightBuffer");
      return;
    }
    if (!Rf_isString(RsAlterWeightFile_)) {
      Rprintf("Issue: Please set AlterWeightFile to be a String!\n"); R_FlushConsole();
    }
    if (AlterWeightdfTNoise_ != NULL && !Rf_isNull(AlterWeightdfTNoise_) &&
      Rf_isReal(AlterWeightdfTNoise_)) {
      if (Verbose >= 2) {
        Rprintf("Setting AlterWeightdfTNoise = %f.\n",REAL(AlterWeightdfTNoise_)[0]);
        R_FlushConsole();
      }
      AlterWeightdfTNoise = REAL(AlterWeightdfTNoise_)[0];  
    }
    if (AlterWeightdfRobit_ != NULL && !Rf_isNull(AlterWeightdfRobit_) &&
      Rf_isReal(AlterWeightdfRobit_)) {
      if (Verbose >= 2) {
        Rprintf("Setting AlterWeightdfRobit = %f.\n",REAL(AlterWeightdfRobit_)[0]);
        R_FlushConsole();
      }
      AlterWeightdfRobit = REAL(AlterWeightdfRobit_)[0];  
    }
    if (Verbose >= 1) {
      Rprintf("SetupAlterWeightBuffer, Alter dfRobit = %f, Alter dfTNoise is %f\n",
        AlterWeightdfRobit, AlterWeightdfTNoise); R_FlushConsole();
    }
    if (AlterWeightdfRobit < 0.0 && AlterWeightdfTNoise < 0.0) {
      Rprintf("AlterWeight Setup, both dfRobit=%f and dfTNoise = %f, clearing!\n",
        AlterWeightdfRobit, AlterWeightdfTNoise);  R_FlushConsole();
      FFree(AlterWeightBuffer, "AlterWeightBuffer");
      DDelete(RsAlterWeightFile, "RsAlterWeightFile");
    }
    if (!Rf_isString(RsAlterWeightFile_)) {
      Rf_error("SetupAlterWeightFile: Error, RsYFile is not a string!\n");
    }
    if (!(LengthAlterWeightBuffer_ != NULL && !Rf_isNull(LengthAlterWeightBuffer_) &&  
       Rf_length(LengthAlterWeightBuffer_) >= 1 && 
      (Rf_isInteger(LengthAlterWeightBuffer_) || Rf_isReal(LengthAlterWeightBuffer_)) &&
      GetFirstInteger(LengthAlterWeightBuffer_) > 0)) {
      Rprintf("SetupAlterWeightFile: Supplied LengthAlterWeightBuffer_ is weird!, we do 100\n");  
      R_FlushConsole();
      LengthAlterWeightBuffer = 100;
    } else {
      LengthAlterWeightBuffer = GetFirstInteger(LengthAlterWeightBuffer_);
    }
    if (AlterWeightBuffer != NULL) {
      Rprintf("SetupAlterWeightFile: Having to clear AlterWeightBuffer. \n"); R_FlushConsole();
      FFree(AlterWeightBuffer, "AlterWeightBuffer");
    }
    if (Verbose >= 1) {
      Rprintf("SetupAlterWeightFile: Consider delete RsAlterWeightFile. \n"); R_FlushConsole();
    }
    DDelete(RsAlterWeightFile, "RsAlterWeightFile");
    RsAlterWeightFile = new AObject(RsAlterWeightFile_);
    if (LengthAlterWeightBuffer > 1) {
      AlterWeightBuffer = (double*)Calloc(LengthAlterWeightBuffer, double);
    } else {
      Rf_error("SetupAlterWeightFile: Not if LengthAlterWeightBuffer = %d\n",
        LengthAlterWeightBuffer); 
    }
    LengthWrittenAlterWeightBuffer = 0;   
    NewAlterWeightBufferWrite = 1;  LengthTotalWrittenAlterWeightBuffer = 0;
    if (Verbose > 1) {
      Rprintf("SetupAlterWeightBuffer, SetupAlterWeightBuffer, length = %d\n",
        LengthAlterWeightBuffer); R_FlushConsole();
      Rprintf("SetupAlterWeightBuffer, Setup AlterWeightBuffer to File %s, with length %d\n",
        CHAR(STRING_ELT(RsAlterWeightFile->asSexp(), 0)), LengthAlterWeightBuffer);
      R_FlushConsole();
    }
  }
  int get_DoLogitNonePostPreProb() {
    return(DoLogitNonePostPreProb);
  }
  void set_DoLogitNonePostPreProb(int iIn) {
    if (iIn == 1) {
      //if (AlterWeightBuffer == NULL) {
      //  Rf_error("set_DoLogitNonePostPreProb: Sorry, fill AlterWeightBuffer first. \n"); R_FlushConsole();
      //}
      DoLogitNonePostPreProb = 1;
      return;
    }
    if (iIn >= 2) {
      DoLogitNonePostPreProb = 2;
      return;
    }
    if (iIn < 0) {
      DoLogitNonePostPreProb = -1;
      return;
    }
    if (iIn == 0) { DoLogitNonePostPreProb = 0; }
    return;
  }
  int CheckkXFinderError(char *aS);
  int DoCheckkXFinderError() {
    return(CheckkXFinderError((char*) "Default"));
  } 
  int CheckkXToXLC(char *aS) {
    if (OnKappaS <= 0 || OnKappaMem <= 0) {
      return(6);
    }
    for (int iti = 0; iti < OnKappaS; iti++) {
      if (kXFinder[iti] < 0 || kXFinder[iti] >= p) {
        Rf_error("CheckkXToXLC[%s], we have iti=%d/%d but kXFinder[iti=%d] = %d!\n",
          aS, iti, OnKappaS, iti, kXFinder[iti]);
      }
      if (XLC[kXFinder[iti]] != iti) {
        Rf_error("CheckkXToXLC[%s], we have iti=%d/%d but kXFinder[iti=%d] = %d! but XLC[kXFinder[iti=%d]=%d]=%d\n",
          aS, iti, OnKappaS, iti, kXFinder[iti],
          iti, kXFinder[iti], XLC[kXFinder[iti]]);      
      }
    }
    for (int iti = 0; iti < p; iti++) {
      if (XLC[iti] >= 0) {
        if (XLC[iti] >= OnKappaS) {
            Rf_error("CheckkXToXLC[%s], we have iti=%d/%d but XLC[iti=%d] = %d!\n",
             aS, iti, p, iti, XLC[iti]);
        }
        if (kXFinder[XLC[iti]] != iti) {
           Rf_error("CheckkXToXLC[%s], we have iti=%d/%d but kXFinder[iti=%d] = %d! but XLC[kXFinder[iti=%d]=%d]=%d\n",
            aS, iti, p, iti, XLC[iti],
            iti, XLC[iti], kXFinder[XLC[iti]]);      
        }        
      }
    }
    return(1);
  }
  int DoCheckkXToXLC() {
    return(CheckkXToXLC((char*)"Default"));
  }
  int get_LengthTotalWrittenAlterWeightBuffer() {
    return(LengthTotalWrittenAlterWeightBuffer);
  }
  SEXP get_RsAlterWeightFile() {
    if (RsAlterWeightFile == NULL) { return(R_NilValue); }
    return(RsAlterWeightFile->asSexp());
  }
  void set_NewAlterWeightBufferWrite(int ANW) {
    if (ANW == 0) { NewAlterWeightBufferWrite = 0;}
  }
  void set_CenteredColumns(SEXP _rSCenteredColumns) {
    if (Rf_isNull(_rSCenteredColumns)) {
      DDelete(RsCenteredColumns, "RsCenteredColumns"); RsCenteredColumns = NULL;
      DDelete(RsDeCenteredCodaList, "RsDeCenteredCodaList");
      return;
    }
    if (Rf_isNull(tauEndList) || Rf_length(tauEndList) <= 0) {
      Rprintf("Error: setCenteredColumns: Tell me centered columns only if tauEndList is set!\n");    R_FlushConsole();
      return;
    }
    if (Rf_length(_rSCenteredColumns) <= 0) {
      Rf_error("Error: setCenteredColumns: must be a length of members longer than 0.\n");
    }
    if (Rf_isInteger(_rSCenteredColumns)) {
      for (int ii = 0; ii < Rf_length(_rSCenteredColumns); ii++) {
        if (INTEGER(_rSCenteredColumns)[ii] < 0 ||
          INTEGER(_rSCenteredColumns)[ii] >= Rf_length(tauEndList)) {
          Rf_error("Error: setCenteredColumns: Error, tauEndList has length %d but bad rSCenteredColumns given!\n"); 
        }
      }
    } else {
      for (int ii = 0; ii < Rf_length(_rSCenteredColumns); ii++) {
        if (REAL(_rSCenteredColumns)[ii] < 0.0 ||
          (int) REAL(_rSCenteredColumns)[ii] >= Rf_length(tauEndList)) {
          Rf_error("Error: setCenteredColumns: Error, tauEndList has length %d but bad rSCenteredColumns given!\n"); 
        }
      }    
    }
    RsCenteredColumns = new AObject(Rf_allocVector(INTSXP, Rf_length(_rSCenteredColumns)));
    if (Rf_isInteger(_rSCenteredColumns)) {
      for (int ii = 0; ii < Rf_length(_rSCenteredColumns); ii++) {
        if (ii > 0) {
          for (int jj = 0; jj < ii; jj++) {
            if (INTEGER(RsCenteredColumns->asSexp())[jj] ==
              INTEGER(_rSCenteredColumns)[ii]) {
              DDelete(RsCenteredColumns, "RsCenteredColumns");
              Rf_error("Error: setCenteredColumns: No Duplicates!!!\n");  
            }
          }
        }
        INTEGER(RsCenteredColumns->asSexp())[ii] = 
          INTEGER(_rSCenteredColumns)[ii];
      }
    } else {
      for (int ii = 0; ii < Rf_length(_rSCenteredColumns); ii++) {
        if (ii > 0) {
          for (int jj = 0; jj < ii; jj++) {
            if (INTEGER(RsCenteredColumns->asSexp())[jj] ==
              (int) REAL(_rSCenteredColumns)[ii]) {
              DDelete(RsCenteredColumns, "RsCenteredColumns");
              Rf_error("Error: setCenteredColumns: No Duplicates!!!\n");  
            }
          }
        }
        INTEGER(RsCenteredColumns->asSexp())[ii] = 
          (int) REAL(_rSCenteredColumns)[ii];
      }    
    }
    DDelete(RsDeCenteredCodaList, "RsDeCenteredCodaList");
    if (Verbose >= 1) {
      Rprintf("Finished setCenteredColumns!\n");
    }
  }
  SEXP get_CenteredColumns() {
    if (RsCenteredColumns == NULL) { return(R_NilValue); }
    return(RsCenteredColumns->asSexp());
  }
  int get_NewAlterWeightBufferWrite() { return(NewAlterWeightBufferWrite); }
  int get_LengthAlterWeightBuffer() { return(LengthAlterWeightBuffer); }
  int get_LengthWrittenAlterWeightBuffer() { return(LengthWrittenAlterWeightBuffer); }
  SEXP get_CurrentAlterWeightBuffer() {
    if (AlterWeightBuffer == NULL) { Rf_error("Hey: AlterWeightBuffer is not setup!\n");}
    if (LengthWrittenAlterWeightBuffer <= 0) { Rprintf("LengthWrittenAlterWeightBuffer = %d\n",
      LengthWrittenWeightBuffer); return(R_NilValue); }
    if (Verbose >= 1) {
      Rprintf("get_CurrentAlterWeightBuffer, with length %d\n",
        LengthWrittenAlterWeightBuffer); R_FlushConsole();
    }
    if (LengthAlterWeightBuffer <= 0) { 
      Rf_error("Hey: LengthAlterWeightBuffer = %d\n", LengthAlterWeightBuffer);}
    SEXP sOut = R_NilValue;
    Rf_protect(sOut = Rf_allocVector(REALSXP, LengthWrittenAlterWeightBuffer));
    int Tot = LengthWrittenAlterWeightBuffer; int One = 1;
    F77_CALL(dcopy)(&Tot, AlterWeightBuffer, &One, REAL(sOut), &One);
    Rf_unprotect(1); return(sOut);
  }
  SEXP get_AllAlterWeightBuffer() {
    if (AlterWeightBuffer == NULL) { Rf_error("Hey: AlterWeightBuffer is not setup!\n");}
    if (LengthWrittenAlterWeightBuffer < 0) { Rprintf("LengthWrittenAlterWeightBuffer = %d\n",
      LengthWrittenAlterWeightBuffer); return(R_NilValue); }
    if (LengthAlterWeightBuffer <= 0) { 
      Rf_error("Hey: LengthAlterWeightBuffer = %d\n", LengthWeightBuffer);}
    SEXP sOut = R_NilValue;
    Rf_protect(sOut = Rf_allocVector(REALSXP, LengthAlterWeightBuffer));
    int Tot = LengthAlterWeightBuffer; int One = 1;
    F77_CALL(dcopy)(&Tot, AlterWeightBuffer, &One, REAL(sOut), &One);
    Rf_unprotect(1); return(sOut);
  }
  int FillYBuffer();  int FillWeightBuffer();
  int WriteYBuffer();  int WriteWeightBuffer();
  
  int WriteSigCodaBuffer();  int WritePiACodaBuffer();
  int FillSigCodaBuffer(); int FillPiACodaBuffer();
  double CalcNewProb(int LenBeta, double *Beta, 
    int *IDBeta, int LenTau,  double *Tau, int *IDTau, 
    double AlterWeightdfRobit, double AlterWeightdfTNoise, 
    double PiA1, double PiA2, double *EY, double Sigma);
  int DeriveAlterProbabilityStartCheck(int *pLenFixed, int *pNumGroups);
  int DeriveAlterProbability(int AInput);
  double get_AlterWeightdfRobit() {
    return(AlterWeightdfRobit);
  }
  double get_AlterWeightdfTNoise() {
    return(AlterWeightdfTNoise);
  }
  void set_AlterWeightdfRobit(double aWdf) {
    AlterWeightdfRobit = aWdf;
  }
  void set_AlterWeightdfTNoise(double aWtdf) {
    AlterWeightdfTNoise = aWtdf;
  }
  
  SEXP get_tauEndList() { 
    if (Verbose >= 1) {
      Rprintf("get_tauEndList(): Hoping to supply a non null version. \n"); R_FlushConsole();
    }
     if (RtauEndList == NULL) {
       if (Verbose >= 1) {
         Rprintf("get_tauEndList(): RtauEndList is NULL so we are supplying NULL back. \n");
         R_FlushConsole();
       }
      return(R_NilValue); }
     if (tauEndList == NULL || Rf_isNull(tauEndList)) { return(R_NilValue); }
     if (Rf_length(tauEndList) <= 0) { return(R_NilValue); }
     return(RtauEndList->asSexp()); }  
  SEXP get_OnTau() { return(sOnTau); } 
  SEXP get_FirstRandom() {
    SEXP sOut = R_NilValue;
    Rf_protect(sOut = Rf_allocVector(INTSXP, 1));
    if (sOnTau == NULL && Verbose >= 3) {
      Rprintf("No First Random because it is sOnTau is NULL. \n"); R_FlushConsole();
    }
    if (Rf_isNull(sOnTau) && Verbose >= 3) {
      Rprintf("No First Random because it is sOnTau is Nil Object. \n"); R_FlushConsole();
    }
    if (Rf_length(sOnTau) <= 0  && Verbose >= 3) {
      Rprintf("No First Random because length sOnTau is 0. \n"); R_FlushConsole();
    }
    if (tauEndList == NULL && Verbose >= 3) {
      Rprintf("No First Random because it is tauEndList is NULL. \n"); R_FlushConsole();
    }
    if (Rf_isNull(tauEndList) && Verbose >= 3) {
      Rprintf("No First Random because it is tauEndList is Nil Object. \n"); R_FlushConsole();
    }
    if (Rf_length(tauEndList) <= 0 && Verbose >= 3) {
      Rprintf("No First Random because length tauEndList is 0. \n"); R_FlushConsole();
    }
    if (sOnTau == NULL || Rf_isNull(sOnTau) || Rf_length(sOnTau) <= 0 ||
      tauEndList == NULL || Rf_isNull(tauEndList) || Rf_length(tauEndList) <= 0) {
      INTEGER(sOut)[0] = p;  
    } else {
      INTEGER(sOut)[0] = iFirstRandom+1;  
    } 
    Rf_unprotect(1);
    return(sOut);
  }
  int get_RawiFirstRandom() {
    return(iFirstRandom);
  }
  void set_RawiFirstRandom(int iFirstRandom_) {
    if (Verbose >= 2) {
      Rprintf("Setting iFirstRandom Raw!\n");  R_FlushConsole();
    }
    iFirstRandom = iFirstRandom_;  return;
  }
  int get_CFirstRandom() {
    if (sOnTau == NULL  && Verbose >= 3) {
      Rprintf("get_CFirstRandom: No First Random because it is sOnTau is NULL. \n"); R_FlushConsole();
    }
    if (Rf_isNull(sOnTau) && Verbose >= 3) {
      Rprintf("get_CFirstRandom: No First Random because it is sOnTau is Nil Object. \n"); R_FlushConsole();
    }
    if (Rf_length(sOnTau) <= 0 && Verbose >= 3) {
      Rprintf("get_CFirstRandom: No First Random because length sOnTau is 0. \n"); R_FlushConsole();
    }
    if (tauEndList == NULL && Verbose >= 3) {
      Rprintf("get_CFirstRandom: No First Random because it is tauEndList is NULL. \n"); R_FlushConsole();
    }
    if (Rf_isNull(tauEndList) && Verbose >= 3) {
      Rprintf("get_CFirstRandom: No First Random because it is tauEndList is Nil Object. \n"); R_FlushConsole();
    }
    if (Rf_length(tauEndList) <= 0 && Verbose >= 3) {
      Rprintf("get_CFirstRandom: No First Random because length tauEndList is 0. \n"); R_FlushConsole();
    }
    if (iFirstRandom == p  && Verbose >= 3)  {
      Rprintf("get_CFirstRandom: Strange, iFirstRandom == p!\n");
    }
    if (sOnTau == NULL || Rf_isNull(sOnTau) || Rf_length(sOnTau) <= 0 ||
      tauEndList == NULL || Rf_isNull(tauEndList) || Rf_length(tauEndList) <= 0) {
      return(p);  
    }
    return(iFirstRandom);
  }
  int get_OnKappaS() { return(OnKappaS);}
  int get_OnKappaMem() { return(OnKappaMem);}
  
  SEXP get_MaxInvert() {
    SEXP sOut = R_NilValue; Rf_protect(sOut = Rf_allocVector(INTSXP,1));
    INTEGER(sOut)[0] = MaxInvert; Rf_unprotect(1); return(sOut);
  }
  void set_MaxInvert(SEXP sIn) {
    if (Rf_isNull(sIn) || Rf_length(sIn) <= 0 || !(Rf_isInteger(sIn) || Rf_isReal(sIn))) {
      Rf_error("set_MaxInvert(): Invalid input. \n"); return;
    }
    int MM = -1;
    if (Rf_isInteger(sIn)) { MM = INTEGER(sIn)[0]; 
    } else if (Rf_isReal(sIn)) { MM = (int) REAL(sIn)[0]; }
    if (MM <= 0) {
      Rf_error("set_MaxInvert(): will not set to %d \n", MM); return;
    }
    MaxInvert = MM;
  }
  ////////////////////////////////////////////////////
  // set/get taubardf/nu
  //
  //   taubardf is number of degrees of freedom for tau^2 inverse gamma prior
  //       It should be stored as SEXP in TBSRoo, and is put into object "MT"
  //
  //   taubarnu is the expected mean from the inverse gamma prior
  void set_staubarnu(SEXP staubarnu_) { 
    if (Rf_isNull(staubarnu)) {
      DDelete(Rstaubarnu, "Rstaubarnu"); Rstaubarnu = new AObject(staubarnu_);
      staubarnu = Rstaubarnu->asSexp(); return;
    }
    if (Rf_length(staubarnu) < Rf_length(staubarnu_)) {
      Rf_error("set_staubarnu, supply vector of length %d\n", Rf_length(staubarnu));
    }
    DDelete(Rstaubarnu, "Rstaubarnu"); Rstaubarnu = new AObject(staubarnu_);
    staubarnu = Rstaubarnu->asSexp(); 
    return;
  }
  SEXP get_staubarnu() { return(staubarnu); }
  void set_staupriordf(SEXP staupriordf_) { 
    if (Rf_isNull(staupriordf)) {
      DDelete(Rstaupriordf, "Rstaupriordf"); Rstaupriordf = new AObject(staupriordf_);
      staupriordf= Rstaupriordf->asSexp(); return;
    }
    if (Rf_length(staupriordf_) < Rf_length(staupriordf)) {
      Rf_error("set_staupriordf, supply vector of length %d\n", Rf_length(staupriordf));
    }
    DDelete(Rstaupriordf, "Rstaupriordf"); Rstaupriordf = new AObject(staupriordf_);
      staupriordf= Rstaupriordf->asSexp();
    return; }
  
  SEXP get_staupriordf() { return(staupriordf); }
  void set_sPriorXScaling(SEXP sPriorXScaling_) { 
    DDelete(RsPriorXScaling, "RsPriorXScaling"); RsPriorXScaling = new AObject(sPriorXScaling_);
    sPriorXScaling= RsPriorXScaling->asSexp(); return;
    return;
  }
  SEXP get_sPriorXScaling() { return(sPriorXScaling); }
  void set_sPriorOftau(SEXP sPriorOftau_) {
    if (Rf_isNull(sPriorOftau)) {
      DDelete(RsPriorOftau, "RsPriorOftau"); RsPriorOftau = new AObject(sPriorOftau_);
      sPriorOftau = RsPriorOftau->asSexp(); return;
    } else {
      if (Rf_length(sPriorOftau) != Rf_length(sPriorOftau_)); {
        Rprintf("sPriorOftau, old length = %d, but new is %d \n",
         Rf_length(sPriorOftau), Rf_length(sPriorOftau_));  R_FlushConsole();  
      }
    DDelete(RsPriorOftau, "RsPriorOftau"); RsPriorOftau = new AObject(sPriorOftau_);
      sPriorOftau = RsPriorOftau->asSexp(); return;
    }
  } 
  SEXP get_sPriorOftau() { return(sPriorOftau); }
  
  SEXP get_LengthOldCoda() { SEXP sOn = R_NilValue;
    Rf_protect(sOn = Rf_allocVector(INTSXP, 1)); 
    INTEGER(sOn)[0] = LengthOldCoda; Rf_unprotect(1);
    return(sOn); }
  SEXP get_LengthProbCoda() { SEXP sOn = R_NilValue;
    Rf_protect(sOn = Rf_allocVector(INTSXP,1));
    INTEGER(sOn)[0] = LengthProbCoda;
    Rf_unprotect(1); return(sOn); }
  SEXP get_LengthPostProb() { SEXP sOn = R_NilValue;
    Rf_protect(sOn = Rf_allocVector(INTSXP,1));
    INTEGER(sOn)[0] = LengthPostProb;
    Rf_unprotect(1); return(sOn); }
  void set_YCodaList(SEXP YCodaList_) {
    DDelete(RsYCodaList, "YCodaList");
    RsYCodaList = new AObject(YCodaList_);
    return;
  }
  SEXP get_YCodaList() {
    if (RsYCodaList == NULL || Rf_isNull(RsYCodaList->asSexp()) ||
      Rf_length(RsYCodaList->asSexp()) <= 0) {
      int PTECTED = 0;
      if (RTBSR5 == NULL) {
        Rprintf("BayesSpikeCL:get_YCodaList, where is TBSR5?\n");
        R_FlushConsole();
        return(R_NilValue);
      }
      SEXP YCodaFunction;
      Rf_protect(YCodaFunction = 
         Rf_findVarInFrame( get_BayesSpikeNameSpace(), 
         Rf_install("ATryFillYCodaList")));
      PTECTED++;
      if (Rf_isNull(YCodaFunction)) {
        Rprintf("BayesSpikeCL: Error: no ATryFillYCoda.\n");
        Rf_unprotect(PTECTED); PTECTED = 0;
        return(R_NilValue);
      }
      SEXP call;
      Rf_protect(call = Rf_lcons( YCodaFunction, 
        Rf_cons(RTBSR5->asSexp(), R_NilValue))); 
      PTECTED++;
      if (Rf_isNull(call)) {
        Rprintf("BayesSpikeCL: get_RsYCodaList, Error, could not genreate call.\n");
        Rf_unprotect(PTECTED); return(R_NilValue);
      }
      SEXP RTM = R_NilValue;
      if (Verbose > 1) {
        Rprintf("get_RsYCodaList: Starting to call ATryFillYCodaList.\n "); R_FlushConsole();
      }
      Rf_protect(RTM = Rf_eval(call, get_BayesSpikeNameSpace())); PTECTED++;
      if (RsYCodaList == NULL) {
        Rprintf("BayesSpikeCL, getYCodaList,after all that, still Null List!\n");
        if (!Rf_isNull(RTM)) {
          RPostProbCodaList = new AObject(RTM);
          Rf_unprotect(PTECTED); PTECTED = 0;
          return(RTM);
        } else {
          Rf_unprotect(PTECTED); PTECTED = 0;
          return(RTM);
        }
      } else {
         Rf_unprotect(PTECTED); PTECTED = 0;
         return(RsYCodaList->asSexp());
      }
      return(R_NilValue);
    }
    return(RsYCodaList->asSexp());
  }

  void set_WeightCodaList(SEXP WeightCodaList_) {
    DDelete(RsWeightCodaList, "WeightCodaList");
    RsWeightCodaList = new AObject(WeightCodaList_);
    return;
  }
  SEXP get_WeightCodaList() {
    if (RsWeightCodaList == NULL || Rf_isNull(RsWeightCodaList->asSexp()) ||
      Rf_length(RsWeightCodaList->asSexp()) <= 0) {
      int PTECTED = 0;
      if (RTBSR5 == NULL) {
        Rprintf("BayesSpikeCL:get_WeightCodaList, where is TBSR5?\n");
        R_FlushConsole();
        return(R_NilValue);
      }
      SEXP WeightCodaFunction;
      Rf_protect(WeightCodaFunction = 
         Rf_findVarInFrame( get_BayesSpikeNameSpace(), 
         Rf_install("ATryFillWeightCodaList")));
      PTECTED++;
      if (Rf_isNull(WeightCodaFunction)) {
        Rprintf("BayesSpikeCL: Error: no ATryFillWeightCoda.\n");
        Rf_unprotect(PTECTED); PTECTED = 0;
        return(R_NilValue);
      }
      SEXP call;
      Rf_protect(call = Rf_lcons( WeightCodaFunction, 
        Rf_cons(RTBSR5->asSexp(), R_NilValue))); 
      PTECTED++;
      if (Rf_isNull(call)) {
        Rprintf("BayesSpikeCL: get_RsWeightCodaList, Error, could not genreate call.\n");
        Rf_unprotect(PTECTED); return(R_NilValue);
      }
      SEXP RTM = R_NilValue;
      if (Verbose > 1) {
        Rprintf("get_RsWeightCodaList: Starting to call ATryFillWeightCodaList.\n "); R_FlushConsole();
      }
      Rf_protect(RTM = Rf_eval(call, get_BayesSpikeNameSpace())); PTECTED++;
      if (RsWeightCodaList == NULL) {
        Rprintf("BayesSpikeCL, getWeightCodaList,after all that, still Null List!\n");
        if (!Rf_isNull(RTM)) {
          RsWeightCodaList = new AObject(RTM);
          Rf_unprotect(PTECTED); PTECTED = 0;
          return(RTM);
        } else {
          Rf_unprotect(PTECTED); PTECTED = 0;
          return(RTM);
        }
      } else {
         Rf_unprotect(PTECTED); PTECTED = 0;
         return(RsWeightCodaList->asSexp());
      }
      return(R_NilValue);
    }
    return(RsWeightCodaList->asSexp());
  }
  int SSave() {
      int PTECTED = 0;
      if (RTBSR5 == NULL) {
        Rprintf("BayesSpikeCL:Save, where is TBSR5?\n");
        R_FlushConsole();
        return(-4);
      }
    SEXP SaveFunction = R_NilValue;
    Rf_protect(SaveFunction = 
         Rf_findVarInFrame( get_BayesSpikeNameSpace(), 
         Rf_install("ATrySave")));
     PTECTED++;
    if (Rf_isNull(SaveFunction)) {
        Rprintf("BayesSpikeCL: Error: no ATrySave.\n");
        Rf_unprotect(PTECTED); PTECTED = 0;
        return(-4);
    }
    SEXP call;
    Rf_protect(call = Rf_lcons( SaveFunction, 
        Rf_cons(RTBSR5->asSexp(), R_NilValue))); 
    PTECTED++;
    if (Rf_isNull(call)) {
        Rprintf("BayesSpikeCL: SSave, Error, could not genreate call.\n");
        Rf_unprotect(PTECTED); return(-4);
    }
    SEXP RTM = R_NilValue;
    if (Verbose > 1) {
        Rprintf("SSave: Starting to call ATrySave.\n "); R_FlushConsole();
    }
    Rf_protect(RTM = Rf_eval(call, get_BayesSpikeNameSpace())); PTECTED++;
     
    Rf_unprotect(PTECTED);
    if (Verbose > 1) {
      Rprintf("SSave: Saved And now trying to quit.!"); R_FlushConsole();
    }
    return(1);
  }
  int UnSave(SEXP InputWant) {
      int PTECTED = 0;
      if (RTBSR5 == NULL) {
        Rprintf("BayesSpikeCL:UnSave, where is TBSR5?\n");
        R_FlushConsole();
        return(-4);
      }
      if (Rf_isNull(InputWant)) {
        Rprintf("BayesSpikeCL:UnSave: no Input given!");
      }
    SEXP UnSaveFunction = R_NilValue;
    Rf_protect(UnSaveFunction = 
         Rf_findVarInFrame( get_BayesSpikeNameSpace(), 
         Rf_install("ATryUnSave")));
     PTECTED++;
    if (Rf_isNull(UnSaveFunction)) {
        Rprintf("BayesSpikeCL: Error: no ATryUnSave.\n");
        Rf_unprotect(PTECTED); PTECTED = 0;
        return(-4);
    }
    SEXP call;
    Rf_protect(call = Rf_lcons( UnSaveFunction, 
        Rf_cons(RTBSR5->asSexp(), Rf_cons(InputWant, R_NilValue)))); 
    PTECTED++;
    if (Rf_isNull(call)) {
        Rprintf("BayesSpikeCL: UnSave, Error, could not generate call.\n");
        Rf_unprotect(PTECTED); return(-4);
    }
    SEXP RTM = R_NilValue;
    if (Verbose > 1) {
        Rprintf("UnSave: Starting to call ATryUnSave.\n "); R_FlushConsole();
    }
    Rf_protect(RTM = Rf_eval(call, get_BayesSpikeNameSpace())); PTECTED++;
     
    Rf_unprotect(PTECTED);
    if (Verbose > 1) {
      Rprintf("UnSave: Saved And now trying to quit.!"); R_FlushConsole();
    }
    return(1);
  }
  int RunRegression() {
      int PTECTED = 0;
      if (RTBSR5 == NULL) {
        Rprintf("BayesSpikeCL:UnSave, where is TBSR5?\n");
        R_FlushConsole();
        return(-4);
      }
    SEXP RunRegressionFunction = R_NilValue;
    Rf_protect(RunRegressionFunction = 
         Rf_findVarInFrame( get_BayesSpikeNameSpace(), 
         Rf_install("ATryRunRegression")));
     PTECTED++;
    if (Rf_isNull(RunRegressionFunction)) {
        Rprintf("BayesSpikeCL: Error: no ATryRunRegression.\n");
        Rf_unprotect(PTECTED); PTECTED = 0;
        return(-4);
    }
    SEXP call;
    Rf_protect(call = Rf_lcons( RunRegressionFunction, 
        Rf_cons(RTBSR5->asSexp(), R_NilValue))); 
    PTECTED++;
    if (Rf_isNull(call)) {
        Rprintf("BayesSpikeCL: RunRegression, Error, could not generate call.\n");
        Rf_unprotect(PTECTED); return(-4);
    }
    SEXP RTM = R_NilValue;
    if (Verbose > 1) {
        Rprintf("RunRegression: Starting to call ATryRunRegression.\n "); R_FlushConsole();
    }
    Rf_protect(RTM = Rf_eval(call, get_BayesSpikeNameSpace())); PTECTED++;
     
    Rf_unprotect(PTECTED);
    if (Verbose > 1) {
      Rprintf("RunRegression: Saved And now trying to quit.!"); R_FlushConsole();
    }
    return(1);
  }

  void DeliberatelySetAlterWeightCodaList(SEXP AlterWeightCodaList_) {
    DDelete(RsAlterWeightCodaList, "RsAlterWeightCodaList");
    RsAlterWeightCodaList = new AObject(AlterWeightCodaList_);
    return;
  }
 SEXP get_AlterWeightCodaList() {
    if (RsAlterWeightCodaList == NULL || Rf_isNull(RsAlterWeightCodaList->asSexp()) ||
      Rf_length(RsAlterWeightCodaList->asSexp()) <= 0) {
      int PTECTED = 0;
      if (RTBSR5 == NULL) {
        Rprintf("BayesSpikeCL:get_AlterWeightCodaList, where is TBSR5?\n");
        R_FlushConsole();
        return(R_NilValue);
      }
      SEXP AlterWeightCodaFunction;
      Rf_protect(AlterWeightCodaFunction = 
         Rf_findVarInFrame( get_BayesSpikeNameSpace(), 
         Rf_install("ATryFillAlterWeightCodaList")));
      PTECTED++;
      if (Rf_isNull(AlterWeightCodaFunction)) {
        Rprintf("BayesSpikeCL: Error: no ATryFillAlterWeightCoda.\n");
        Rf_unprotect(PTECTED); PTECTED = 0;
        return(R_NilValue);
      }
      SEXP call;
      Rf_protect(call = Rf_lcons( AlterWeightCodaFunction, 
        Rf_cons(RTBSR5->asSexp(), R_NilValue))); 
      PTECTED++;
      if (Rf_isNull(call)) {
        Rprintf("BayesSpikeCL: RsAlterWeightCodaList, Error, could not genreate call.\n");
        Rf_unprotect(PTECTED); return(R_NilValue);
      }
      SEXP RTM = R_NilValue;
      if (Verbose > 1) {
        Rprintf("RsAlterWeightCodaList: Starting to call AlterWeightCodaList.\n "); R_FlushConsole();
      }
      Rf_protect(RTM = Rf_eval(call, get_BayesSpikeNameSpace())); PTECTED++;
      if (RsAlterWeightCodaList == NULL) {
        Rprintf("BayesSpikeCL, RsAlterWeightCodaList,after all that, still Null List!\n");
        if (!Rf_isNull(RTM)) {
          RsAlterWeightCodaList = new AObject(RTM);
          Rf_unprotect(PTECTED); PTECTED = 0;
          return(RTM);
        } else {
          Rf_unprotect(PTECTED); PTECTED = 0;
          return(RTM);
        }
      } else {
         Rf_unprotect(PTECTED); PTECTED = 0;
         return(RsAlterWeightCodaList->asSexp());
      }
      return(R_NilValue);
    }
    return(RsAlterWeightCodaList->asSexp());
  }
  
  void set_SigCodaList(SEXP SigCodaList_) {
    DDelete(RsSigCodaList, "SigCodaList");
    RsSigCodaList = new AObject(SigCodaList_);
    return;
  }
  void set_PiACodaList(SEXP PiACodaList_) {
    DDelete(RsPiACodaList, "PiACodaList");
    RsPiACodaList = new AObject(PiACodaList_);
    return;
  }
  SEXP get_PiACodaList() {
    if (RsPiACodaList == NULL || Rf_isNull(RsPiACodaList->asSexp()) ||
      Rf_length(RsPiACodaList->asSexp()) <= 0) {
      int PTECTED = 0;
      if (RTBSR5 == NULL) {
        Rprintf("BayesSpikeCL:get_PiACodaList, where is TBSR5?\n");
        R_FlushConsole();
        return(R_NilValue);
      }
      SEXP PiACodaFunction;
      Rf_protect(PiACodaFunction = 
         Rf_findVarInFrame( get_BayesSpikeNameSpace(), 
         Rf_install("ATryFillPiACodaList")));
      PTECTED++;
      if (Rf_isNull(PiACodaFunction)) {
        Rprintf("BayesSpikeCL: Error: no ATryFillPiACoda.\n");
        Rf_unprotect(PTECTED); PTECTED = 0;
        return(R_NilValue);
      }
      SEXP call;
      Rf_protect(call = Rf_lcons( PiACodaFunction, 
        Rf_cons(RTBSR5->asSexp(), R_NilValue))); 
      PTECTED++;
      if (Rf_isNull(call)) {
        Rprintf("BayesSpikeCL: RsPiACodaList, Error, could not genreate call.\n");
        Rf_unprotect(PTECTED); return(R_NilValue);
      }
      SEXP RTM = R_NilValue;
      if (Verbose > 1) {
        Rprintf("get_RsPiACodaList: Starting to call ATryFillPiACodaList.\n "); R_FlushConsole();
      }
      Rf_protect(RTM = Rf_eval(call, get_BayesSpikeNameSpace())); PTECTED++;
      if (RsPiACodaList == NULL) {
        Rprintf("BayesSpikeCL, getPiACodaList,after all that, still Null List!\n");
        if (!Rf_isNull(RTM)) {
          RsPiACodaList = new AObject(RTM);
          Rf_unprotect(PTECTED); PTECTED = 0;
          return(RTM);
        } else {
          Rf_unprotect(PTECTED); PTECTED = 0;
          return(RTM);
        }
      } else {
         Rf_unprotect(PTECTED); PTECTED = 0;
         return(RsPiACodaList->asSexp());
      }
      return(R_NilValue);
    }
    return(RsPiACodaList->asSexp());
  }
  SEXP get_SigCodaList() {
    if (RsSigCodaList == NULL || Rf_isNull(RsSigCodaList->asSexp()) ||
      Rf_length(RsSigCodaList->asSexp()) <= 0) {
      int PTECTED = 0;
      if (RTBSR5 == NULL) {
        Rprintf("BayesSpikeCL:get_SigCodaList, where is TBSR5?\n");
        R_FlushConsole();
        return(R_NilValue);
      }
      SEXP SigCodaFunction;
      Rf_protect(SigCodaFunction = 
         Rf_findVarInFrame( get_BayesSpikeNameSpace(), 
         Rf_install("ATryFillSigCodaList")));
      PTECTED++;
      if (Rf_isNull(SigCodaFunction)) {
        Rprintf("BayesSpikeCL: Error: no ATryFillSigCoda.\n");
        Rf_unprotect(PTECTED); PTECTED = 0;
        return(R_NilValue);
      }
      SEXP call;
      Rf_protect(call = Rf_lcons( SigCodaFunction, 
        Rf_cons(RTBSR5->asSexp(), R_NilValue))); 
      PTECTED++;
      if (Rf_isNull(call)) {
        Rprintf("BayesSpikeCL: get_RsSigCodaList, Error, could not genreate call.\n");
        Rf_unprotect(PTECTED); return(R_NilValue);
      }
      SEXP RTM = R_NilValue;
      if (Verbose > 1) {
        Rprintf("get_RsWeightCodaList: Starting to call ATryFillSigCodaList.\n "); R_FlushConsole();
      }
      Rf_protect(RTM = Rf_eval(call, get_BayesSpikeNameSpace())); PTECTED++;
      if (RsSigCodaList == NULL) {
        Rprintf("BayesSpikeCL, getSigCodaList,after all that, still Null List!\n");
        if (!Rf_isNull(RTM)) {
          RsSigCodaList = new AObject(RTM);
          Rf_unprotect(PTECTED); PTECTED = 0;
          return(RTM);
        } else {
          Rf_unprotect(PTECTED); PTECTED = 0;
          return(RTM);
        }
      } else {
         Rf_unprotect(PTECTED); PTECTED = 0;
         return(RsSigCodaList->asSexp());
      }
      return(R_NilValue);
    }
    return(RsSigCodaList->asSexp());
  }
  //int get_LengthPostProbCodaBuffer() { return(LengthPostProbCodaBuffer); }
  int get_iiOnILocProb() { return(iiOnILocProb); }
  void set_iiOnILocProb(int iiiOnILocProb) {
    if (iiiOnILocProb >= 0) {
      iiOnILocProb = iiiOnILocProb;
    } else {
      Rprintf("Give a larger iiOnILocProb");
    }
  }
  int FullCrazyTestiWeightedXtX(char *iIn) {
    if (iWeightedXtX == NULL) {
      return(6);
    }
    Rprintf(" -- Full CraztyTest %s - iWeighted XtX ", iIn);
    R_FlushConsole();
    int ii; 
    for (ii = 0; ii < OnKappaMem; ii++) {
      iWeightedXtX[ii] = iWeightedXtX[ii]+0;
    }
    Rprintf(" -- PASS 1 -- "); R_FlushConsole();
    short int *NewiWeighted = NULL;
    RMemGetS(NewiWeighted, "NewiWeighted", OnKappaMem);
    for (ii = 0; ii < OnKappaMem; ii++) {
      NewiWeighted[ii] = iWeightedXtX[ii];
    }
    Rprintf(" -- Yes now Free --"); R_FlushConsole();
    Free(iWeightedXtX);  iWeightedXtX = NewiWeighted;
    Rprintf(" -- Total Pass \n"); R_FlushConsole();
    return(1);
    
  }

};


/////////////////////////////////////////////////////////////////////////
//  Locate these functions in BayesSpliceSampling.cpp
// Slice Sampling of Taus

inline int InsertSmallResid(int n, int p,
  double *XtResid, double **pXtX, double *sBeta, int iFirstRandom, 
  SEXP sOnTau, SEXP tauEndList, int iOnii, double *SmallXtResid, int Verbose);

///////////////////////////////////////////////////////////
// In BayesSpike Helpers
//  If we believe that some effects should not be on if other are good
int CheckDependencies(SEXP DependenceList, int iFirstRandom,
  short *BFixed, SEXP sOnTau, int iOnii);
    
  
//////////////////////////////////////////////////////////
//  Maximization and slice sampling are necessary for grouped factors.
// 
//////////////////////////////////////////////////////////////////
// Functions for Maximizing Tau/~Beta posterior:  
double MaximizeMe(double StartTau, TauOfFContainer *MT,
  int MaximizeMeCauchyTotalIters, double CauchyEpsilon, int Verbose);
double BinFindX(double (*GiveMeFOfX)(double, void*), double StartX,
 double GoalF, double MinX, double MaxX, int MaxInt, double MoveDist,
 void *HelpStruct,
 double Epsilon);
double GetSomeMax(double (*GiveMeFOfX)(double, void*), double StartX,
 double MinX, double MaxX, int MaxInt, double MoveDist,
 void *HelpStruct, double Epsilon);

/////////////////////////////////////////
//  Finalizer for BayesSpikeCL 
void FinalizerForBayesSpikeCL(BayesSpikeCL *AcL);



 /////////////////////////////////////////
 // Fucntions for Splice Sampling:   
 // (AlsSwitchingSampler within the class is another component)
double SpliceSample( double(*GiveMeFOfX)(double, void *), double StartX,
  double SpliceMove, double MinX, double MaxX, double Unif1, double Unif2,
  int MaxInt, void *HelpStruct, double Epsilon);
double SampleTauSpliced(double MaxTauEst, 
  TauOfFContainer *MT, int MaximizeMeCauchyTotalIters, 
  int MaxIters, double CauchyEpsilon,
  double SpliceMove,  int NumSamsPer, int DoMax, 
  int Verbose, double Oldtau);
int FullIntegrate (double (*lFOfX)(double, void*), TauOfFContainer* MT, 
  int NumIntegrands, int iOnii, double *ProbList, double UsePiA, double TweakCons,
  double *pIntegratedDensity, int Verbose, double *PressMe, double *pProbScore);
double FunInt(double (*GiveMelFOfX)(double, void*), 
  TauOfFContainer *MT, int MaxN, double Tweak, double Stop);  
inline int InsertSmallResid(int n, int p,
  double *XtResid, double **pXtX, int *XLC, SEXP sBeta, int iFirstRandom, 
  SEXP sOnTau, SEXP tauEndList, int iOnii, double *SmallXtResid,
  SEXP MyEigenValues, SEXP MyEigenVectors);
#endif
