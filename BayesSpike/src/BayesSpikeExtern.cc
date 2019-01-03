/* ========================================================================== */
/*                                                                            */
/*   BayesSpikeExtern.c                                                               */
/*   (c) 2011 Alan Lenarcic                                                          */
/*                                                                            */
/*   BayesSpike Extern Functions                                                              */
/*                                                                            */
/*   These functions are largely independent of BayesSpikeCpp Rcpp based functions.
/* Predominantly this code is tasked with file reading both compressed Beta
/* And Compressed Tau draws.
/* ========================================================================== */

#ifndef BAYESSPIKECPPH
  #include "BayesSpikeCpp.h"
  #define BAYESSPIKECPPH 0
#endif

#ifndef UNISTDH
  #include <unistd.h>
  #include <fcntl.h>
  #include <sys/stat.h>
  #include <sys/types.h>
  #define UNISTDH 0
#endif


extern "C" {
SEXP TestEntriesForSampler(SEXP sVerbose,
  SEXP sY, SEXP sX, SEXP sXtY, SEXP sXtResid, SEXP sXtCResid,
  SEXP sYSq, SEXP sYResidSq,
  SEXP sXtX, SEXP NewQ, SEXP CholQI, SEXP InvQ,  
  SEXP sOnBeta, SEXP sPropBeta,
  SEXP sOnPiA, SEXP sPiAPrior,
  SEXP sOnSigma, SEXP sSigmaPrior,
  SEXP BFixed, SEXP tauFixed, SEXP TypeFixed,
  SEXP sOnTau,
  SEXP FirstRandom, SEXP tauEndList,
  SEXP sRecordMatrix, SEXP DoRecord,
  SEXP EigenValuesList, SEXP EigenVectorsList, SEXP sRVals,
  SEXP sPriorOftau, SEXP staubarnu, SEXP snu, SEXP sPriorXScaling,
  SEXP sMaximizeMeCauchyEpsilonIters, SEXP sCauchyEpsilon, SEXP sSpliceMove, 
  SEXP NumSpliceSampleReps,
  SEXP DoMax, SEXP DoRNGState, SEXP FakeOuttau, SEXP MaxTauEst,
  SEXP DependenciesBetaFixed, SEXP DependenciesTau,
  //SEXP sOutMatrix, 
  SEXP sNumberNew,
  SEXP stt, SEXP ArtOnii, SEXP OutYResidSq,
  SEXP ProbtauList, SEXP ProbFixedList,
  SEXP sBSMT, SEXP sHowSample) {
  Rprintf("TestEntriesForSampler  \n"); R_FlushConsole();
  PrintDetailOfRsprec(   sY, (char*) "sY");
  PrintDetailOfRsprec(   sX,(char *) "sX");
  PrintDetailOfRsprec(   sXtY,(char *) "sXtY");
  PrintDetailOfRsprec(   sXtResid,(char *) "sXtResid");
  PrintDetailOfRsprec(   sXtCResid,(char *) "sXtCResid");
  PrintDetailOfRsprec(   sYSq,(char *) "sYSq");
  PrintDetailOfRsprec(   sYResidSq,(char *) "sYResidSq");
  PrintDetailOfRsprec(   sXtX,(char *) "sXtX");
  PrintDetailOfRsprec(   NewQ,(char *) "NewQ");
  PrintDetailOfRsprec(   CholQI,(char *) "CholQI");
  PrintDetailOfRsprec(   InvQ,(char *) "InvQ"); 
  PrintDetailOfRsprec(   sOnBeta,(char *) "sOnBeta");
  PrintDetailOfRsprec(   sPropBeta, (char*) "sPropBeta");
  PrintDetailOfRsprec(   sOnPiA,(char *) "sOnPiA");
  PrintDetailOfRsprec(   sPiAPrior,(char *) "sPiAPrior");
  PrintDetailOfRsprec(   sOnSigma,(char *) "sOnSigma");
  PrintDetailOfRsprec(   sSigmaPrior,(char *) "sSigmaPrior");
  PrintDetailOfRsprec(   BFixed,(char *) "BFixed");
  PrintDetailOfRsprec(   tauFixed,(char *) "tauFixed");
  PrintDetailOfRsprec(   TypeFixed,(char *) "TypeFixed");
  PrintDetailOfRsprec(   sOnTau,(char *) "sOnTau");
  PrintDetailOfRsprec(   FirstRandom,(char *) "FirstRandom");
  PrintDetailOfRsprec(   tauEndList,(char *) "tauEndList");
  PrintDetailOfRsprec(   sRecordMatrix,(char *) "sRecordMatrix");
  PrintDetailOfRsprec(   DoRecord,(char *) "DoRecord");
  PrintDetailOfRsprec(   EigenValuesList,(char *) "EigenValueList");
  PrintDetailOfRsprec(   EigenVectorsList,(char *) "EigenVectorsList");
  PrintDetailOfRsprec(   sRVals, (char*) "sRVals");
  PrintDetailOfRsprec(   sPriorOftau,(char *) "sPriorOftau");
  PrintDetailOfRsprec(   staubarnu,(char *) "staubarnu");
  PrintDetailOfRsprec(   snu,(char *) "snu");
  PrintDetailOfRsprec(   sPriorXScaling, (char *) "sPriorXScaling");
  PrintDetailOfRsprec(   sMaximizeMeCauchyEpsilonIters,(char *) "sMaximizeMeCauchyEpsilonIters");
  PrintDetailOfRsprec(   sCauchyEpsilon,(char *) "sCauchyEpsilon"); 
  PrintDetailOfRsprec(   sSpliceMove,(char *) "sSpliceMove");
  PrintDetailOfRsprec(   NumSpliceSampleReps,(char*) "NumSpliceSampleReps");
  PrintDetailOfRsprec(   DoMax,(char *) "DoMax");
  PrintDetailOfRsprec(   DoRNGState,(char *) "DoRNGState");
  PrintDetailOfRsprec(   FakeOuttau,(char *) "FakeOuttau"); 
  PrintDetailOfRsprec(   DependenciesBetaFixed,(char *) "DependenciesBetaFixed");
  PrintDetailOfRsprec(   DependenciesTau,(char *) "DependenciesTau");
  //PrintDetailOfRsprec(   sOutMatrix,(char *) "sOutMatrix");
  PrintDetailOfRsprec(   sNumberNew,(char *) "sNumberNew");
  PrintDetailOfRsprec(   stt,(char *) "stt");
  PrintDetailOfRsprec(   ArtOnii,(char *) "ArtOnii");
  PrintDetailOfRsprec(   OutYResidSq,(char *) "OutYResidSq");
  PrintDetailOfRsprec(   ProbtauList,(char *) "ProbtauList");
  PrintDetailOfRsprec(   ProbFixedList,(char *) "ProbFixedList");
  
  PrintDetailOfRsprec(   sBSMT,(char *) "sBSMT");
   PrintDetailOfRsprec(   sHowSample,(char *) "sHowSample");
   PrintDetailOfRsprec(   sVerbose,(char *) "sVerbose");
   return(sVerbose);
}

SEXP CallWeightedHPD(SEXP SortedXtX, SEXP cumSum, SEXP alphaTarget) {
  if (Rf_isNull(SortedXtX) || !Rf_isReal(SortedXtX)) {
    Rf_error("WeightedHPD, we can't do that.! SortedXtX\n");
  }
  if (Rf_isNull(cumSum) || !Rf_isReal(cumSum)) {
    Rf_error("WeightedHPD, we can't do that.! cumSum\n");
  }
  if (Rf_isNull(alphaTarget) || !Rf_isReal(alphaTarget)) {
    Rf_error("WeightedHPD, we can't do that.! alphaTarget\n");
  }
  if (Rf_length(cumSum) != Rf_length(SortedXtX)) {
    Rf_error("WeightedHPD, cumSum != length SortedXtX!");
  }
  double alpha = REAL(alphaTarget)[0];
  double MinLength = 0.0;     SEXP sOut = R_NilValue;
  int BigMin = 0; int BigMax = 0;
  int OnMin = 0; int OnMax = 0;
  if (alpha >= 1.0) {
    Rf_error("WeightedHPD, won't do this for alpha=%f\n", alpha);
  }
  OnMin = 0; OnMax = 1;
  int n = Rf_length(SortedXtX);
  if (n <= 1) {
    Rf_error("WeightedHPD, no n=%d. Not enough for intervals. \n", n);
  }
  while(OnMax < n &&                         
    REAL(cumSum)[OnMax] -REAL(cumSum)[OnMin] < alpha) {
    OnMax++;  
  } 
  if (OnMax >= n) {
    Rf_protect(sOut = Rf_allocVector(INTSXP, 2));
    INTEGER(sOut)[0] = 0+1;
    INTEGER(sOut)[1] = n;
    Rf_unprotect(1);  return(sOut);
    Rf_error("WeightedHPD, ERROR, Oh No, We somehow have OnMax=%d, n=%d and we didn't hit Alpha=%f\n", OnMax, n, alpha);
  }
  MinLength = REAL(SortedXtX)[OnMax] - REAL(SortedXtX)[OnMin];
  BigMin = OnMin;  BigMax = OnMax;
  //Rprintf("Test Max we have (OnMin=%d,OnMax=%d) to start coverage is %f, and alpha=%f\n",
  //  OnMin, OnMax, REAL(cumSum)[OnMax]-REAL(cumSum)[OnMin], alpha); R_FlushConsole(); 

  int TCount = 0;
  while(OnMax < n && OnMin < n && REAL(cumSum)[OnMax]-REAL(cumSum)[OnMin] > alpha && TCount < n) {
    if (REAL(cumSum)[OnMax]-REAL(cumSum)[OnMin+1] > alpha) {
      OnMin++;
      if (REAL(SortedXtX)[OnMax] - REAL(SortedXtX)[OnMin] < MinLength) {
         BigMin = OnMin;  BigMax = OnMax;  
         MinLength = REAL(SortedXtX)[OnMax]-REAL(SortedXtX)[OnMin];
      }
    } else if (OnMax < n-1) {
      OnMax++;
    } else {
      OnMin= n;  OnMax = n;
    }
    TCount++;
  }
  if (TCount >= n) {
    Rprintf("WeightedHPD: Error, we have TCount=%d, n=%d we hit max!\n", TCount, n); R_FlushConsole();
  }
  Rf_protect(sOut = Rf_allocVector(INTSXP, 2));
  INTEGER(sOut)[0] = BigMin+1;
  INTEGER(sOut)[1] = BigMax+1;
  Rf_unprotect(1);
  return(sOut); 
}


SEXP GiveCodaSubset(SEXP SubSetCoords, SEXP sCodaIFile, SEXP sCodaDFile,
    SEXP sStartIter, SEXP sEndIter, SEXP sVerbose) {
  int Verbose = GetFirstInteger(sVerbose);
  int TotalOpenedFiles = 0;
  int TotalClosedFiles = 0;
  
  if (!Rf_isString(sCodaIFile)) {
    Rprintf("GiveCodaSubSet: Unfortunately, sCodaIFile is not A String!\n"); R_FlushConsole();
  }
  if (!Rf_isString(sCodaDFile)) {
    Rprintf("GiveCodaSubSet: Unfortunately, sCodaDFile is not A String!\n"); R_FlushConsole();
  }
  if (Verbose > 0) {
    Rprintf("GiveCodaSubSet: About to read from IFile %s, DFile %s \n",
      CHAR(STRING_ELT(sCodaIFile, 0)),CHAR(STRING_ELT(sCodaDFile, 0)));
  }

  int fd1, fd2;
  struct stat stbuf; 
  fd1 = open( CHAR(STRING_ELT(sCodaIFile, 0)), O_RDONLY);
  if (fd1 == -1) {
    Rf_error("Setup GiveCodaSubset: We can't open sCodaIFile"); 
  } 
  TotalOpenedFiles++;

  fd2 = open( CHAR(STRING_ELT(sCodaDFile, 0)), O_RDONLY);
  if (fd2 == -1) {
    close(fd1); TotalClosedFiles++;
    Rf_error("Setup GiveCodaSubset: We can't open sCodaJFile"); 
  } 
  TotalOpenedFiles++;
  FILE *fCodaIFile = fdopen( fd1, "rb");
  FILE *fCodaDFile = fdopen( fd2, "rb");
  if (fCodaIFile != NULL) { TotalOpenedFiles++;
  } else { close(fd1); fd1 = -1; TotalClosedFiles++; }
  if (fCodaDFile != NULL) { TotalOpenedFiles++;
  } else { close(fd2); fd2 = -1; TotalClosedFiles++; }
  if (fd1 == -1) {
    Rprintf("*******************************************************\n");
    Rprintf("GiveCodasubSet Initial Error \n");
    Rprintf("Setup GiveCodaSubSet: we had a open fd1 or fd2 error and are quitting!\n");
    R_FlushConsole();
    Rf_error(" Bad GiveCodaSubSet on fd1, fd2 open ");
  } 
  
  if (fstat(fd1, &stbuf) == -1) {
    Rf_error("GiveCodaSubSet sCodaIFile: could not use fstat");
  }
  long lSizeI = stbuf.st_size / sizeof(int);
  if (fstat(fd2, &stbuf) == -1) {
    Rf_error("GiveCodaSubSet sCodaDFile: could not use fstat");
  }
  fclose(fCodaDFile); fCodaDFile = NULL; TotalClosedFiles++;
  fclose(fCodaIFile); fCodaIFile = NULL; TotalClosedFiles++;
  fd1 = -1; TotalClosedFiles++;
  fd2 = -1; TotalClosedFiles++;
  
  long lSizeD = stbuf.st_size / sizeof(double);
  
  size_t result = 0;
  
  // obtain file size:
  
  if ((lSizeI  ) < 
    (lSizeD )) {
    Rprintf("GiveCodaSubSet:  Error After checking length files %s, and %s\n",
      CHAR(STRING_ELT(sCodaIFile, 0)), CHAR(STRING_ELT(sCodaDFile, 0))  ); R_FlushConsole();
    Rprintf("lSizeI = %d, lSizeD = %d, sizeof(int) = %d, sizeof(char) = %d, sizeof(double) = %d\n",
      lSizeI, lSizeD, sizeof(int), sizeof(char), sizeof(double));
    Rf_error("GiveCodaSubSet: Sorry about the Error");
  }
    
  fCodaIFile = fopen( CHAR(STRING_ELT(sCodaIFile, 0)), "rb");
  fCodaDFile = fopen( CHAR(STRING_ELT(sCodaDFile, 0)), "rb");
  int BLOCKSIZE = 10000;
  if (fCodaIFile==NULL) {
    Rf_error("GiveCodaSubset: sCodaIFile failed to open");
  }  
  TotalOpenedFiles++;  
  if (fCodaDFile==NULL) {
    fclose(fCodaIFile); fCodaIFile = NULL;  TotalClosedFiles++;
    Rf_error("GiveCodaSubset: sCodaDFile failed to open");
  } 
  TotalOpenedFiles++;
  int StartIter = GetFirstInteger(sStartIter);
  int EndIter = GetFirstInteger(sEndIter);
  if (StartIter < 0) {
     fclose(fCodaIFile); fclose(fCodaDFile); TotalClosedFiles+=2;
     if (TotalClosedFiles != TotalOpenedFiles) {
       Rprintf("*****Double Error GiveCodaSubset");
       Rprintf("***** PleaseSupply 0 < StartIter=%d, and TotalClosedFiles=%d, TotalOpenedFiles=%d\n",
         StartIter, TotalClosedFiles, TotalOpenedFiles); R_FlushConsole();
     }
     Rf_error(" GiveCodaSubset: Please supply 0 < StartIter(%d?)",
       StartIter);
  }    
  if (EndIter < StartIter) {
     fclose(fCodaIFile); fclose(fCodaDFile); TotalClosedFiles+=2;
     if (TotalClosedFiles != TotalOpenedFiles) {
       Rprintf("*****Double Error GiveCodaSubset");
       Rprintf("***** PleaseSupply 0 < StartIter=%d < EndIter=%d, and TotalClosedFiles=%d, TotalOpenedFiles=%d\n",
         StartIter, EndIter, TotalClosedFiles, TotalOpenedFiles); R_FlushConsole();
     }
    Rf_error("  GiveCodaSubset: Please set EndIter >= StartIter \n");
  }


  
  int TotalReads = lSizeI * sizeof(char) / sizeof(int);
  
  if (Rf_isNull(SubSetCoords)) {
    Rf_error("Error: SubsEtCoords is Nil Object!");
  }                    
  if (Rf_length(SubSetCoords) <= 0) {
    Rf_error("Please give a nonzero subset of coordinates!\n");
  }
  if (Verbose > 2) {
    Rprintf("Goal is to collect data on coordinates: \n");
    if (Rf_isReal(SubSetCoords)) {
      PrintVector(REAL(SubSetCoords), Rf_length(SubSetCoords));
    } else if (Rf_isInteger(SubSetCoords)) {
      PrintVector(INTEGER(SubSetCoords), Rf_length(SubSetCoords));
    } else {
      fclose(fCodaIFile); fclose(fCodaDFile); TotalClosedFiles+=2;
      if (TotalClosedFiles != TotalOpenedFiles) {
        Rprintf("*****Double Error SubSetCoords was not real or integer\n");
        Rprintf("***** PleaseSupply 0 < StartIter=%d, and TotalClosedFiles=%d, TotalOpenedFiles=%d\n",
          StartIter, TotalClosedFiles, TotalOpenedFiles); R_FlushConsole();
      }
      Rf_error("SubSetCoords is not Real or Integer!");
    }
  }
  if (Verbose > 1) {
    Rprintf("GiveCodaSubSet:  Allocating Matrix of length %d,%d\n", 
      EndIter - StartIter+1, Rf_length(SubSetCoords));
    R_FlushConsole();
  }
  SEXP ReturnMat;
  Rf_protect(ReturnMat = Rf_allocMatrix(REALSXP, 
    EndIter - StartIter+1, Rf_length(SubSetCoords)));
  double ZeroD = 0.0;  int One = 1;
  int RLenR = Rf_length(ReturnMat);
  F77_CALL(dscal)(&RLenR, &ZeroD, REAL(ReturnMat), &One); //Zero Out matrix!
  
  if (TotalReads < BLOCKSIZE) { BLOCKSIZE = TotalReads; }  
  
  if (Verbose > 1) {
    Rprintf("Allocating Buffers. \n"); R_FlushConsole();
  }
  // allocate memory to contain the whole file:
  int *bufferI = NULL;  double *bufferD = NULL;
  bufferI = (int*) Calloc(BLOCKSIZE, int);
  bufferD = (double*) Calloc(BLOCKSIZE, double);  
  if (bufferI == NULL) {
    Rf_error("GiveCodaSubSet: Cannot allocate bufferI");
  }  
  if (bufferD == NULL) {
    Rf_error("GiveCodaSubSet: Cannot allocate bufferD");
  } 
  int ReadBlocks = 0;
  int Finished = 0; 
  int ReadBuffer = 0;
  int LastCoord = -1; int OnIter = 0; int OnSexp = 0;
  while(Finished == 0) {
    if (ReadBlocks + BLOCKSIZE >= lSizeI) {
       BLOCKSIZE = lSizeI - ReadBlocks;  Finished = 1;
    }
    if (Verbose > 2) {
      Rprintf("New Read in Buffer, ReadBlocks = %d, ReadBuffer = %d, LastCoord = %d, OnIter = %d\n",
        ReadBlocks, ReadBuffer, LastCoord, OnIter); R_FlushConsole();
    }
    result = fread(bufferI, sizeof(int), BLOCKSIZE, fCodaIFile);
    if ((int) result != (int) BLOCKSIZE) {
      Rf_error("GiveCodaSubSet: Could not read %d after %d from I, total %d",
        BLOCKSIZE, ReadBlocks, lSizeI);
    }
    result = fread(bufferD, sizeof(double), BLOCKSIZE, fCodaDFile);
    if ((int) result != (int) BLOCKSIZE) {
      Rf_error("GiveCodaSubSet: Could not read %d after %d from D, total %d",
        BLOCKSIZE, ReadBlocks, lSizeD);
    }
    ReadBuffer = 0; ReadBlocks += BLOCKSIZE;
    while (ReadBuffer < BLOCKSIZE) {
      if (bufferI[ReadBuffer] < 0) {
        LastCoord = -1;  OnIter=abs(bufferI[ReadBuffer])-1;
        if (OnIter > EndIter)  {
           FFree(bufferI, "bufferI"); FFree(bufferD, "bufferD");
           fclose(fCodaIFile); fclose(fCodaDFile);  TotalClosedFiles+=2;
           if (TotalOpenedFiles != TotalClosedFiles) {
             Rprintf("*********************************************************\n");
             Rprintf("*** GiveCodaSubSet: We have within While close at an Error\n");
             Rprintf("*** TotalOpenedFiles=%d, TotalClosedFiles = %d \n",
              TotalOpenedFiles, TotalClosedFiles);  R_FlushConsole();
             Rf_error("***  Well go seek opened/closed error \n");
           }
           Rf_unprotect(1);
           return(ReturnMat);      
        }
        OnSexp = 0;
      } else if( bufferI[ReadBuffer] < LastCoord) {
        if (Verbose > 3) {
          Rprintf("End iter %d, new iter, bufferI[%d] = %d, LastCoord = %d, bufferD[%d] = %f\n",
             OnIter, ReadBuffer, bufferI[ReadBuffer], LastCoord,
               ReadBuffer, bufferD[ReadBuffer]); R_FlushConsole();
        }
        LastCoord = bufferI[ReadBuffer];  OnIter++;
        if (OnIter > EndIter) {
          FFree(bufferI, "bufferI"); FFree(bufferD, "bufferD");
          fclose(fCodaIFile); fclose(fCodaDFile);
          TotalClosedFiles+=2;
           if (TotalOpenedFiles != TotalClosedFiles) {
             Rprintf("*********************************************************\n");
             Rprintf("*** GiveCodaSubSet: We have at OnIter>EndIter \n");
             Rprintf("*** TotalOpenedFiles=%d, TotalClosedFiles = %d \n",
              TotalOpenedFiles, TotalClosedFiles);  R_FlushConsole();
             Rf_error("***  Well go seek opened/closed error \n");
           }          
          Rf_unprotect(1);
          return(ReturnMat);
        }
        if (OnIter > EndIter) { Finished = 1;  ReadBuffer = 2 * BLOCKSIZE;  break;}
        if (OnIter >= StartIter && OnIter <= EndIter) {
          OnSexp = 0;
          while( OnSexp < Rf_length(SubSetCoords) &&
            ANINT(SubSetCoords, OnSexp) < LastCoord) { OnSexp++; }
          if ( ANINT(SubSetCoords, OnSexp) == LastCoord) {
            REAL(ReturnMat)[OnIter - StartIter + 
              OnSexp * INTEGER(Rf_getAttrib(ReturnMat, R_DimSymbol))[0]]  =
            bufferD[ReadBuffer];    
          }    
        }
      } else {
        LastCoord = bufferI[ReadBuffer];
        if (OnIter >= StartIter && OnIter <= EndIter) {
          while( OnSexp < Rf_length(SubSetCoords) &&
            ANINT(SubSetCoords, OnSexp) < LastCoord) { OnSexp++; }
          if ( ANINT(SubSetCoords, OnSexp) == LastCoord) {
              REAL(ReturnMat)[OnIter - StartIter + 
                OnSexp * INTEGER(Rf_getAttrib(ReturnMat, R_DimSymbol))[0]]  =
                bufferD[ReadBuffer];    
          } 
        }
      }  
      ReadBuffer++;
    }
  }

  FFree(bufferI, "bufferI"); FFree(bufferD, "bufferD"); 
  int weClosed = 0;
  weClosed = fclose(fCodaIFile);  
  if (weClosed != 0) {
    Rprintf("GiveCodaSubset: Well, we failed to close sCodaIFile !\n");  R_FlushConsole();
  }
  TotalClosedFiles++;
  weClosed = fclose(fCodaDFile);
  if (weClosed != 0) {
    Rprintf("GiveCodaSubset: Well, we failed to close sCodaJFile !\n");  R_FlushConsole();
  }
  TotalClosedFiles++;
  if (TotalOpenedFiles != TotalClosedFiles) {
             Rprintf("*********************************************************\n");
             Rprintf("*** GiveCodaSubSet: We have out of loop close error\n");
             Rprintf("*** TotalOpenedFiles=%d, TotalClosedFiles = %d \n",
              TotalOpenedFiles, TotalClosedFiles);  R_FlushConsole();
             Rf_error("***  Well go seek opened/closed error \n");
  }
  Rf_unprotect(1);
  return(ReturnMat);
}
 //End Of Extern!
 
 
 

SEXP GiveCountFrequency(SEXP sCodaIFile,
    SEXP sStartIter, SEXP sEndIter, SEXP GiveCounts, SEXP sVerbose) {
  int Verbose = GetFirstInteger(sVerbose);    
    
  int LengthGiveCounts = Rf_length(GiveCounts);
  if (LengthGiveCounts <= 0) {
    Rprintf("GiveCountFreuncy: You didn't give any length vector");
  }
  
  if (!Rf_isString(sCodaIFile)) {
    Rprintf("GiveCountFreuncy: Unfortunately, sCodaIFile is not A String!\n"); R_FlushConsole();
  }

  if (Verbose > 0) {
    Rprintf("GiveCountFrequency:  Welcome \n"); R_FlushConsole();
    Rprintf("GiveCountFreuncy: About to read from IFile %s. \n",
      CHAR(STRING_ELT(sCodaIFile, 0))); R_FlushConsole();
  }
  FILE *fCodaIFile = fopen( CHAR(STRING_ELT(sCodaIFile, 0)), "rb");
  int BLOCKSIZE = 10000;
  if (fCodaIFile==NULL) {
    Rf_error("GiveCountFreuncy: sCodaIFile failed to open");
  }    
  
  int StartIter = GetFirstInteger(sStartIter);
  int EndIter = GetFirstInteger(sEndIter);
  if (StartIter < 0) {
     Rf_error(" GiveCountFreuncy: Please supply 0 < StartIter(%d?)",
       StartIter);
  }    
  if (EndIter < StartIter) {
    Rf_error("  GiveCountFreuncy: Please set EndIter >= StartIter \n");
  }
  if (Verbose > 1) {
    Rprintf("  GiveCountFrequency: Start = %d, End = %d\n", 
      StartIter, EndIter); R_FlushConsole();
  }
  long lSizeI;
  size_t result = 0;

  // obtain file size:
  fseek (fCodaIFile , 0 , SEEK_END);
  lSizeI = ftell (fCodaIFile) / sizeof(int);
  rewind (fCodaIFile);
  
  int TotalReads = lSizeI * sizeof(char) / sizeof(int);

  if (Verbose > 1) {
    Rprintf("GiveCountFreuncy:  lSizeI = %d, TotalReads = %d \n", lSizeI, TotalReads);
    R_FlushConsole();
  }

  //if (Rf_isReal(GiveCounts)) {
  //  double ZeroD = 0.0;  int One = 1;
  //  int RLenR = Rf_length(GiveCounts);
  //  F77_CALL(dscal)(&RLenR, &ZeroD, REAL(GiveCounts), &One); //Zero Out matrix!
  //} else {
  //  for (int iii = 0; iii < Rf_length(GiveCounts); iii+++) {
  //    INTEGER(GiveCounts)[iii] = 0;
  //  }
  //}
  
  if (TotalReads < BLOCKSIZE) { BLOCKSIZE = TotalReads; }  
  
  if (Verbose > 1) {
    Rprintf("Allocating Buffers. \n"); R_FlushConsole();
  }
  // allocate memory to contain the whole file:
  int *bufferI = NULL;
  bufferI = (int*) Calloc(BLOCKSIZE, int);
  if (bufferI == NULL) {
    Rf_error("GiveCountFreuncy: Cannot allocate bufferI");
  }  
 
  int ReadBlocks = 0;
  int Finished = 0; 
  int ReadBuffer = 0;
  int LastCoord = -1; int OnIter = -1;
  while(Finished == 0) {
    if (ReadBlocks + BLOCKSIZE >= lSizeI) {
       BLOCKSIZE = lSizeI - ReadBlocks;  Finished = 1;
    }
    if (Verbose > 2) {
      Rprintf("New Read in Buffer, BLOCKSIZE = %d, ReadBlocks = %d, ReadBuffer = %d, LastCoord = %d, OnIter = %d\n",
        BLOCKSIZE, ReadBlocks, ReadBuffer, LastCoord, OnIter); R_FlushConsole();
    }
    result = fread(bufferI, sizeof(int), BLOCKSIZE, fCodaIFile);
    if ((int) result != (int) BLOCKSIZE) {
      Rf_error("GiveCodaSubSet: Could not read %d after %d from I, total %d",
        BLOCKSIZE, ReadBlocks, lSizeI);
    }
    ReadBuffer = 0; ReadBlocks += BLOCKSIZE;
    while (ReadBuffer < BLOCKSIZE) {
      if (bufferI[ReadBuffer] < 0) {
        if (OnIter > 0 && abs(bufferI[ReadBuffer])-1 <= OnIter) {
          Rprintf("Doh, Chain just backed up on itself, OnIter = %d, but buffer = %d, don't do that!\n",
            OnIter, abs(bufferI[ReadBuffer]));
           FFree(bufferI, "bufferI"); 
           fclose(fCodaIFile); fCodaIFile = NULL;
           return(GiveCounts);      
        }
        OnIter=abs(bufferI[ReadBuffer])-1; LastCoord = -100;
        if (OnIter > EndIter)  {
           FFree(bufferI, "bufferI"); 
           fclose(fCodaIFile); fCodaIFile = NULL;
           return(GiveCounts);      
        }
      } else if( bufferI[ReadBuffer] < LastCoord) {
        if (Verbose > 3) {
          Rprintf("End iter %d, new iter, bufferI[%d] = %d, LastCoord = %d\n",
             OnIter, ReadBuffer, bufferI[ReadBuffer], LastCoord,
               ReadBuffer); R_FlushConsole();
        }
        OnIter++;
        if (OnIter > EndIter) {
          FFree(bufferI, "bufferI");
          fclose(fCodaIFile); 
          return(GiveCounts);
        }
        if (OnIter >= StartIter && OnIter <= EndIter  && Rf_length(GiveCounts) > bufferI[ReadBuffer]) {
          if (Rf_isReal(GiveCounts)) { REAL(GiveCounts)[bufferI[ReadBuffer]]++;
          } else if (Rf_isInteger(GiveCounts)) { INTEGER(GiveCounts)[bufferI[ReadBuffer]]++; }
        }
      } else {
        if (bufferI[ReadBuffer] == LastCoord) {
          Rprintf("Error, we just read from the buffer the same coordinate we just quit %d \n",
            LastCoord); R_FlushConsole();
          int StartOn = ReadBuffer -5; int Len = 10;
          if (ReadBuffer < 5) { StartOn = 0; }
          if (StartOn + Len >= BLOCKSIZE) { Len = BLOCKSIZE - StartOn - 1;}
          Rprintf("ReadBuffer = %d, and the, StartOn = %d, Len = %d, printing \n");
          PrintVector( bufferI + StartOn, Len);
          Rprintf("\n"); R_FlushConsole();
          
        }
        if (OnIter >= StartIter && OnIter <= EndIter && Rf_length(GiveCounts) > bufferI[ReadBuffer]) {
          if (Rf_isReal(GiveCounts)) { REAL(GiveCounts)[bufferI[ReadBuffer]] +=  1.0;
          } else if (Rf_isInteger(GiveCounts)) { INTEGER(GiveCounts)[bufferI[ReadBuffer]]++; }
        }
        LastCoord = bufferI[ReadBuffer];
      }  
      ReadBuffer++;
    }
  }

  FFree(bufferI, "bufferI");
  int weClosed;
  weClosed = fclose(fCodaIFile);
  if (weClosed != 0) {
    Rprintf("GiveCounts: Well, we failed to close sCodaIFile !\n");  R_FlushConsole();
  }
  return(GiveCounts);
}
 //End Of Extern!


}
