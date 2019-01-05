## ABayesSpikeStart.r
##
##  Alan Lenarcic, 2014
##
##  These R files declare a few important functions to be attached to
##   BAYESSPIKENAMESPACE on load of the BayesSpike Package


## LICENSE INFO: R CODE
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  https://www.R-project.org/Licenses/

BAYESSPIKENAMESPACE <- environment()

require(R.oo);  require(coda);

SetupFixedPrior <- function(nFixed = 0, TypePrior = 0, MeanTau = 1, 
  DfTau = 1, MeanLambda =1,
  DfLambda =1, DfX =-1, MeanX =-1, Min = -1000, Max = 1000,
  sigmaMu = 1, FixedRVec = NULL, FixTau = -1) {
  OntauFixed = c();
  if (!is.null(FixedRVec) && length(FixedRVec) == nFixed) {
    OntauFixed = FixedRVec;
  } else {
    OntauFixed = rep(1, nFixed);
  }
  if (nFixed <= 0) { return(NULL); }
  if (TypePrior == 1) {
    if (DfX >1 && MeanX > 1) {
      OnTauFixed = c(OnTauFixed, 1.0, DfX, MeanX, DfTau, MeanTau);   
    } else if (is.null(FixedRVec) || length(FixedRVec) != nFixed) {
      if (FixTau > 0) {
        tauFixed = FixTau;
      } else {
        tauFixed = c(1.0, DfTau, MeanTau);
      }   
    }
  } else if (TypePrior == 2) {
  
  
  }
  return(OntauFixed)
  
}



SetGEqualsText <- function(AText, EqualsVar, envir = "globalenv()", S=0) {
  if (S==0) {
  return(paste("TryO = NULL; try( TryO <- bindingIsActive(\"", AText, "\", ", envir, " ), silent=TRUE);
    if (!is.null(TryO) && (TryO == TRUE || TryO == FALSE)) {
      unlockBinding( \"", AText, "\", ", envir, " ); 
    }
    try(unlockBinding( \"", AText, "\", ", envir, " ), silent=FALSE); 
    assign( \"", AText, "\", ", EqualsVar, ", ", envir, " );
    ", sep=""));
  } else {
   return(paste("TryO = NULL; try( TryO <- bindingIsActive(\"", AText, "\", ", envir, " ), silent=TRUE);
    if (!is.null(TryO) && (TryO == TRUE || TryO == FALSE)) {
      unlockBinding( \"", AText, "\", ", envir, " ); 
    }
    try(unlockBinding( \"", AText, "\", ", envir, " ), silent=TRUE); 
    assign( \"", AText, "\", ", EqualsVar, ", ", envir, " );
    ", sep="")); 
  }
}



SetGText <- function(AText, envir = "globalenv()", S=0) {
  if (S==0) {
  return(paste("TryO = NULL; try( TryO <- bindingIsActive(\"", AText, "\", ", envir, " ), silent=TRUE);
    if (!is.null(TryO) && (TryO == TRUE || TryO == FALSE)) {
      unlockBinding( \"", AText, "\", ", envir, " ); 
    }
    try(unlockBinding( \"", AText, "\", ", envir, " ), silent=FALSE); 
    assign( \"", AText, "\", ", AText, ", ", envir, " );
    ", sep=""));
  } else {
   return(paste("TryO = NULL; try( TryO <- bindingIsActive(\"", AText, "\", ", envir, " ), silent=TRUE);
    if (!is.null(TryO) && (TryO == TRUE || TryO == FALSE)) {
      unlockBinding( \"", AText, "\", ", envir, " ); 
    }
    try(unlockBinding( \"", AText, "\", ", envir, " ), silent=TRUE); 
    assign( \"", AText, "\", ", AText, ", ", envir, " );
    ", sep="")); 
  }
}


GetG0Text <- function(AText, envir = "globalenv()", S=0) {
  if (S== 0) {
    return(paste("TryO = NULL; try( TryO <- bindingIsActive(\"", AText, "\",  ", envir, " ), silent=TRUE);
    if (!is.null(TryO) && (TryO == TRUE || TryO == FALSE)) {
      try(unlockBinding( \"", AText, "\", ", envir, " ), silent=FALSE);
    }
    try(unlockBinding( \"", AText, "\", ", envir, " ), silent=FALSE); 
    if (exists(\"", AText, "\", ", envir, " )) {
    ", AText, " <- get(\"", AText, "\", ", envir, " );
    } else {
      ", AText, " = 0;
    }
    ", sep=""));
  } else {
    return(paste("TryO = NULL; try( TryO <- bindingIsActive(\"", AText, "\",  ", envir, " ), silent=TRUE);
    if (!is.null(TryO) && (TryO == TRUE || TryO == FALSE)) {
      try(unlockBinding( \"", AText, "\", ", envir, " ), silent=TRUE);
    }
    try(unlockBinding( \"", AText, "\", ", envir, " ), silent=TRUE); 
    if (exists(\"", AText, "\", ", envir, " )) {
    ", AText, " <- get(\"", AText, "\", ", envir, " );
    } else {
      ", AText, " = 0;
    }
    ", sep=""));   
  }
}


LockGText <- function(AText, envir="globalenv()", S=0)  {
  paste("TryO = NULL; try( TryO <- bindingIsActive(\"", AText, "\", ", envir, "), silent=TRUE)
    if (!is.null(TryO) && (TryO == TRUE || TryO == FALSE)) {
      unlockBinding( \"", AText, "\", ", envir, "  ) 
    }
    try(unlockBinding( \"", AText, "\", ", envir, "  ), silent=TRUE); 
    assign( \"", AText, "\", ", AText, ", ", envir, " );
    lockBinding( \"", AText, "\", ", envir, " );
    ", sep="");
}

addname <- function(Set = NULL,FN = NULL, Names = NULL) {
  if (!is.null(Set) && !is.null(names(Set))) {
    newname = names(Set);
  } else if (!is.null(Names)) {
    newname = Names;
  } else if (!is.null(FN) && length(FN) == length(Set)) {
    newname = FN;
  } else {
    newname = rep("", length(Set))
  }
  if (!is.null(FN) && length(FN) == 1) {
    FT =(1:length(newname))[substr(newname, 1,nchar(FN)) != FN];
    newname[FT] = paste(FN, ":", newname[FT], sep="");
  } else if (length(FN) == length(newname)) {
    FT =(1:length(newname))[substr(newname, 1,nchar(FN)) != FN];
    newname[FT] = paste(FN[FT], ":", newname[FT], sep="");  
  }
  FT = (1:length(newname))[substr(newname, nchar(newname),nchar(newname)) == ":"]
  newname[FT] = substr(newname,1,nchar(newname)-1);
  return(newname);
}

########################################################################
##   BayesSpikeRoo is an R.oo Object containing Regression information
##     Ordinarily, this info is saved to .ListBayesSpikeOb() which
##     holds the information in the environment so that R's garbage collector
##     doesn't seize it.  This R.oo object will be saved as element
##     "TBSRoo" in the Rcpp object (Aka, query "MyRcpp"$TBSRoo) to
##     get the BayesSpikeRoo object
##
##   A lot of the code here is about checking inputs (making sure X and Y are right dims, etc.)
##    since potentially the user can enter many times of data which one might
##    want to format correctly
##
##   To learn about the Rcpp class, "BayesSpikeCL", see the src code in
##      "BayesSpikeCpp.h" and "BayesSpikeCpp.cpp" which defines a C++ class
##      which is accessible through Rcpp modules to R users.  This is what
##      is customarily returned by R function from BayesSpike
##      Data in this object however is stored in C, which is potentially
##      toast for R's garbage collector, except that the information
##      is also stored in the TBSRoo element which is attached to the environment
##
setConstructorS3("BayesSpikeRoo", function(X=NULL, Y=NULL, BetaStart = NULL,
  tauEndList = NULL, sOnTau = NULL, FirstRandom = 0,
  CodaTable = NULL, tauFixed = 40, TypeFixedPrior = 0,
  HowSample = 3, OnSigma = 1, OnPiA = .5, PiAPrior = c(1,1), SigmaPrior = c(1,1),
  DoRecord = rep(0,7), dfTNoise = -1, InitFlags = c(0,1,2,3),
  MaxGibbsIters = 1000, NumChains = 3, DependenciesTau = NULL,
  DependenciesFixed = NULL, NumSamples=1000,
  CauchyEpsilon = .00001, MaxIters = 100,
  MaximizeMeCauchyTotalIters = 100,
  tauPriordf = 1, tauPriorMean = 1, ttStart = 0,
  EarlyEndStep = -1, EarlyEndtt = -1, FileSave = "MyS",
  SaveDir = get(".DefaultFileDir", BAYESSPIKENAMESPACE), 
    FileName = "ABSChain",  NewWrite = 1,
  DoMax = 1, PutRNG = 1, MyPS = NULL,  CodaTableNames = -1,
  RStart = 0, Verbose = 0, 
  ListSaveFiles = NULL, TempCoda = NULL, TempCodaIDs = NULL,  AFD=NULL,
  NoShrinkFixed = NULL, NoShrinkFixedPrior = NULL, 
  NoShrinkRandom=NULL, NoShrinkRandomPrior=NULL, DoEEProbSort = 0,
  TempertureDecreasingRate = 1.0,
  ...) {
  t1 = proc.time();
  p=0; n=0;
  sBeta = NULL;
  ##print("BayesSpikeRoo Start"); flush.console();
  if (length(dim(X)) %in% c(0,1)) { X = NULL; Y= NULL; BetaStart = NULL;
    return(NULL);
  } else if (length(dim(X)) ==2) {
    p = dim(X)[2]; n = dim(X)[1];
  }
  if (p == 0) {
    print("p is zero, cant get from X, doh! "); flush.console();
    return(NULL);
  }
  sBeta <- rep(0,p);
  if (!is.null(Y) && length(Y) != n) { X= NULL; Y= NULL; 
    print("Error: Please Input BayesSpike Valid X, Y "); flush.console(); 
    return(NULL);
  }
  ##print(paste("TBSRoo generate or, 1: current length of sBeta is ", length(sBeta), sep="")); 
  ##flush.console();
  if (is.null(DoEEProbSort) || DoEEProbSort == FALSE) {
    DoEEProbSort = 0;
  } else if (DoEEProbSort == TRUE) {
    DoEEProbSort = 1;
  } else {
    DoEEProbSort = 1;
  }
    
  if (p > 0) { 
    if (is.null(BetaStart)) {
    } else if (length(BetaStart) != p) {
      print(paste("Input BayesSpike p = ", p, 
        " but BetaStart wrong length ", length(BetaStart), sep=""));
      flush.console();
    } else {
      sBeta = sBeta + BetaStart;
    }
  }
  ##print(paste("TBSRoo generate or, 2: current length of sBeta is ", length(sBeta), sep="")); 
  ##flush.console();
  NsBeta = NULL;
  try(NsBeta <- names(sBeta));
  if (is.null(sBeta)) {
  } else if (p <= 0) {
  
  } else if (is.null(NsBeta)) {
    NsBetaStart = NULL;
    try(NsBetaStart <- names(BetaStart));
    NsColnamesX  <- NULL;
    try(NsColnamesX <- colnames(X));
    if (!is.null(BetaStart) && !is.null(NsBetaStart) &&
      length(names(BetaStart)) == p) {
     names(sBeta) = addname(Names = NsBetaStart, 
        FN=paste("Beta:", 1:length(sBeta), sep=""))   
    } else if (!is.null(X) && !is.null(NsColnamesX) &&
      length(colnames(X)) == p) {
     names(sBeta) = addname(Names = NsColnamesX, 
        FN=paste("Beta:", 1:length(sBeta), sep="")) 
    }  else {
      names(sBeta) = NULL;
      try(names(sBeta) <- paste("Beta:", (1:length(sBeta)), sep=""));
    }
  } else {
    names(sBeta) = addname(Names = NsBeta,
      FN=paste("Beta:", (1:length(sBeta)),  sep=""))
  }
  if (!is.null(SaveDir) && length(list.files(SaveDir)) >= 1 && NewWrite == 1) {
     ITI = list.files(SaveDir);
     for (ii in 1:length(ITI)) {
        try(unlink(paste(SaveDir, "//", ITI[ii], sep="")), silent=TRUE);
     }
  }
  
  NsOnTau = NULL;
  try(NsOnTau <- names(sOnTau));
  NstauEndList <- NULL;
  try(NstauEndList <- names(tauEndList));
  if (is.null(sOnTau)) {
    if (!is.null(tauEndList)) {
      sOnTau = rep(0, length(tauEndList));
      names(sOnTau) = NstauEndList;
    } else { sOnTau = NULL; }
  } else if (is.null(NsOnTau)) {
    if (!is.null(tauEndList) &&
      length(tauEndList) ==  length(sOnTau) ) {
      names(sOnTau) = addname(Names = NstauEndList, 
        FN=paste("tau:", 1:length(sOnTau), sep="")) 
    }  else {
      names(sOnTau) = paste("tau:", 1:length(sOnTau), sep="");
    }
  } else {
    names(sOnTau) = addname( Names=names(sOnTau),
      FN=paste("tau:", 1:length(sOnTau), sep=""))
  }
  
  
  if (is.null(tauEndList)) {
    FirstRandom = length(X[1,]) +1;
  }
  if (is.null(FirstRandom) && is.null(tauEndList)) {
    FirstRandom = p+1;
  } else if (is.null(FirstRandom)) {
    FirstRandom = 1;
  }
  CFirstRandom = FirstRandom -1;  
  if (!is.null(tauEndList)) {
    CtauEndList = tauEndList - 1;
  } else {
    CtauEndList = NULL;
  }
  ##print("BayesSpikeRoo Start extend"); flush.console();
  if (is.null(X)) { print("We Made X NULL!"); flush.console();}
  ##print(paste("X[1:10] = ", X[1:10], sep="")); flush.console();
  if (length(InitFlags) != 4) {
    InitFlags = rep(0,4);
  }
  if (is.null(sOnTau) && !is.null(tauEndList)) {
    sOnTau = rep(1, length(tauEndList));
    names(sOnTau) = names(tauEndList);
  }
  InitFlags[1] = .NumBayesSpikeOb +1;
  if (tauPriordf[1] < 0) { tauPriordf = 1; }
  if (tauPriorMean[1] < 0) { tauPriorMean = 1;}
  ##print("About to Extend!"); flush.console();
  ##return(NULL);
  NInstance = 0;
  try( NInstance <- get(".NumBayesSpikeOb", BAYESSPIKENAMESPACE) + 1);
  ##print(paste("TBSRoo generate or, 3: current length of sBeta is ", length(sBeta), sep="")); 
  ##flush.console();  
  ##print(paste("TBSRoo, EarlyEndStep = ", EarlyEndStep, 
  ##  "EarlyEndtt = ", EarlyEndtt, sep="")); flush.console();
  if (is.null(dfTNoise)) { dfTNoise = -2; }
  if (!exists("Sigma")) { Sigma = 1; }
  MyR <- extend(Object(), "BayesSpikeRoo",  sBeta=sBeta,
    BetaStart = BetaStart,
    AllEigenValues=NULL, AllEigenVectors=NULL, 
    tauEndList = tauEndList, sOnTau = sOnTau, FirstRandom = FirstRandom,
    CFirstRandom = CFirstRandom, CtauEndList = CtauEndList, CodaTable = CodaTable,
    TypeFixedPrior = TypeFixedPrior, tauFixed=tauFixed, HowSample = HowSample,
    OnPiA = OnPiA, OnSigma= OnSigma, PiAPrior = PiAPrior,
    SigmaPrior = SigmaPrior, DoRecord=DoRecord,
    iiWeights = NULL, dfTNoise = dfTNoise, tauPriordf = tauPriordf,
    tauPriorMean=tauPriorMean, PriorXScaling=NULL, PriorOfTau= NULL,
    DependenciesFixed=NULL, DependenciesTau = NULL,
    ABayesSpikeCL = NULL, NInstance = .NumBayesSpikeOb +1,
    InitFlags = InitFlags, .X = X, .Y = Y,
    MaxGibbsIters = MaxGibbsIters, NumChains = NumChains,
    CodaList = NULL, p = p, n=n,DependenciesTau = DependenciesTau,
    DependenciesFixed = DependenciesFixed,
    NumSamples=NumSamples,
    CauchyEpsilon = CauchyEpsilon, CodaTableNames = CodaTableNames,    
    MaxIters = MaxIters, 
    MaximizeMeCauchyTotalIters = MaximizeMeCauchyTotalIters,
    tauPriordf = tauPriordf, tauPriorMean = tauPriorMean,
    ttStart = ttStart, EarlyEndStep = EarlyEndStep, EarlyEndtt = EarlyEndtt,
    SaveDir = SaveDir, FileSave = FileSave, FileName = FileName,
    OnChainIter = 1, NewWrite = NewWrite,  DoMax = DoMax,
    PutRNG = PutRNG, MyPS = MyPS, RStart = RStart, 
    ListSaveFiles = ListSaveFiles, TempCoda = TempCoda, TempCodaIDs = TempCodaIDs,
    Verbose = Verbose, AFD =AFD, .BeingDestroyed = 0,
    NoShrinkFixed = NoShrinkFixed, NoShrinkFixedPrior = NoShrinkFixedPrior,
    NoShrinkRandom = NoShrinkRandom,NoShrinkRandomPrior=NoShrinkRandomPrior,
    CNoShrinkFixed = NoShrinkFixed -1, 
    CNoShrinkRandom = NoShrinkRandom-1, DoEEProbSort = DoEEProbSort,
    TempertureDecreasingRate=TempertureDecreasingRate
  );
  ## print("Generated MyR"); flush.console();
  if (is.null(MyR$.X)) {
    print("Strange: Here MyR$.X is NULL!"); flush.console();
  }
  
  #################################################################
  ##  Code here locks MyR to the BAYESSPIKENAMESPACE environment;
  ## 
  ##BAYESSPIKENAMESPACE <- environment()
  try(eval(parse(text=GetG0Text(".ListBayesSpikeOb", "globalenv()"))));
  .ListBayesSpikeOb[[MyR$NInstance]] <- MyR;
  try(eval(parse(text=LockGText(".ListBayesSpikeOb", "globalenv()"))));

  try(eval(parse(text=GetG0Text(".NumBayesSpikeOb", "globalenv()"))));
  .NumBayesSpikeOb = .NumBayesSpikeOb + 1;
  try(eval(parse(text=LockGText(".NumBayesSpikeOb", "globalenv()"))));

      
   if (is.null(MyR$X)) {
     print("After Extend: X is null!"); flush.console();
   }
   if (is.null(MyR$.X)) {
     print("After Extend: .X is Null!"); flush.console();
   }
  ##print("BayesSpikeRoo Finish"); flush.console();
  if (MyR$Verbose > 0) {
    t2 = proc.time();
    try(print(paste(" BayesSpikeRoo, took (", paste(round(t2-t1, 4)[1:3], collapse=", "),
      ") to finish. ", sep="")));
  }
  
  eval(parse(text=GetG0Text(".ListBayesSpikeOb", "globalenv()")));
  return( .ListBayesSpikeOb[[MyR$NInstance]] );
});

 setMethodS3("SetupSaveFile", "BayesSpikeRoo", function(this, 
   SaveDir=-1, FileName = "ABSChain", chainiter=1,
   Starttt = 0, NewWrite = 1, DoEEProbSort = -1, ...) {
   if (this$DoSave == FALSE) { return(3); }
   if (is.null(DoEEProbSort) || DoEEProbSort == FALSE) {
     this$DoEEProbSort = 0;
   } else if (DoEEProbSort == TRUE) {
     this$DoEEProbSort = 1;
   } else if (is.numeric(DoEEProbSort) && DoEEProbSort[1] ==0) {
     this$DoEEProbSort = 0;
   } else if (is.numeric(DoEEProbSort) && DoEEProbSort[1] > 0) {
     this$DoEEProbSort = 1;
   }
   if (is.null(this$ListSaveFiles)) {
     this$ListSaveFiles = list();
   }
   if (is.null(this$ListChainIters)) {
     this$ListChainIters = list();
   }
   if (is.null(SaveDir) || !is.character(SaveDir) || SaveDir[1] == -1) {
     if (!is.null(this$SaveDir) && is.character(this$SaveDir)) {
       SaveDir = this$SaveDir
     } else {
        if (!exists(".DefaultFileDir", BAYESSPIKENAMESPACE)) {
          print("SetupSaveFile::  Uh Oh, want to use .DefaultFileDir, but it does not exist!"); flush.console();
        }
        this$SaveDir = get(".DefaultFileDir", BAYESSPIKENAMESPACE);
        SaveDir = .DefaultFileDir;
     }
   } else if (is.character(SaveDir)) {
     this$SaveDir = SaveDir;
   }
   if (this$Verbose >= 1) {
     print(paste("  Starting to create SaveDir = ", as.character(SaveDir), sep="")); flush.console();
   }
   RTE = NULL;
   try(RTE <- dir.create(SaveDir, showWarnings=FALSE));
   if (is.null(RTE)) {
     print(paste(" Looks like we had an error trying to create directory: ", as.character(SaveDir)));
     print(paste("  Note that .DefaultFileDir = ", 
       get(".DefaultFileDir", BAYESSPIKENAMESPACE), sep="")); 
     flush.console();
   }
   if (is.null(SaveDir) || (is.numeric(SaveDir) && SaveDir == -1)) {
     FullFileI = NULL
     FullFileD = NULL
     FullFileiT = NULL
     FullFiledT = NULL
     FullFileP = NULL
     FullFileProb = NULL;  FullFileILoc = NULL;
     
     this$FullFileI = NULL;  this$FullFileD = NULL;
     this$FullFileiT = NULL;  this$FullFiledT = NULL;
     this$FullFileP = NULL;
     ##this$ListSaveFiles[[length(this$ListSaveFiles) + 1]] = NULL;
     this$ListChainIters[[length(this$ListChainIters) + 1]] =  chainiter;
     try(this$ABayesSpikeCL$SetCodaFileNames(NULL, NULL, 0));
     
     try(this$ABayesSpikeCL$SetCodaILocAndProbFile(NULL,NULL,0));
   } else {
     if (!is.null(this$ABayesSpikeCL$sCodaDFile) && this$ABayesSpikeCL$Temperature != 1.0) {
       try(this$ABayesSpikeCL$sCodaOldDFile <- 
         paste(this$ABayesSpikeCL$sCodaDFile, sep=""));
     }
     if (!is.null(this$ABayesSpikeCL$sCodaIFile)  && this$ABayesSpikeCL$Temperature != 1.0) {
       try(this$ABayesSpikeCL$sCodaOldIFile <- 
         paste(this$ABayesSpikeCL$sCodaIFile, sep=""));
     }    
     FullFileI = paste(SaveDir, "//", FileName, chainiter, "I.bin", sep="")
     FullFileD = paste(SaveDir, "//", FileName, chainiter, "D.bin", sep="") 
     FullFileiT = paste(SaveDir, "//", FileName, chainiter, "iT.bin", sep="") 
     FullFiledT = paste(SaveDir, "//", FileName, chainiter, "dT.bin", sep="") 
     FullFileP = paste(SaveDir, "//", FileName, chainiter, "P.bin", sep="")  
     FullFileProb = paste(SaveDir, "//", FileName, chainiter, "Prob.bin", sep="") 
     FullFileiLoc = paste(SaveDir, "//", FileName, chainiter, "iLoc.bin", sep="")
     
     this$FullFileI = FullFileI;  this$FullFileD = FullFileD;
     this$FullFileiT = FullFileiT;  this$FullFiledT = FullFiledT;
     this$FullFileP = FullFileP;
     this$FullFileProb = FullFileProb;
     this$FullFileiLoc = FullFileiLoc
     
     this$ListSaveFiles[[length(this$ListSaveFiles) + 1]] = 
       paste(SaveDir, "//", FileName, chainiter, sep="");
     this$ListChainIters[[length(this$ListChainIters) + 1]] =  chainiter;
     try(this$ABayesSpikeCL$SetCodaFileNames(
       paste(this$FullFileI, sep=""), paste(this$FullFileD, sep=""), chainiter));
     if (this$DoRecord[2] == 1  && !is.null(this$tauEndList)) {
       ATry = 0;
       try(ATry <- this$ABayesSpikeCL$setCodaTFile(
         paste(this$FullFileiT, sep=""),
         paste(this$FullFiledT, sep="")));
       if (ATry!=1) {
         print(paste("SetCodaTFile, we have an Error, FullFileiT = ", this$FullFileiT,
           "  and FullFiledT = ", this$FullFiledT, sep=""));
         return(-1);
       }
       try(this$ABayesSpikeCL$NewWrite  <- NewWrite);
     }
     if ( (this$DoRecord[3] == 1  || this$DoRecord[4] == 1) && 
       (!is.null(this$OnPiA) || !is.null(this$OnSigma)) ) {
        this$ABayesSpikeCL$CodaPFile = paste(this$FullFileP, sep="");
     }
     try(this$ABayesSpikeCL$SetCodaILocAndProbFile(paste(this$FullFileiLoc, sep=""),
       paste(this$FullFileProb, sep=""),
       this$ABayesSpikeCL$MaxGibbsIters));
   }

  this$ABayesSpikeCL$tt = Starttt; 

  eval(parse(text=GetG0Text(".ListBayesSpikeOb", "globalenv()")));
  .ListBayesSpikeOb[[this$NInstance]]$FullFileI <- this$FullFileI;
  .ListBayesSpikeOb[[this$NInstance]]$FullFileD <- this$FullFileD;
  .ListBayesSpikeOb[[this$NInstance]]$FullFileiT <- this$FullFileiT;
  .ListBayesSpikeOb[[this$NInstance]]$FullFiledT <- FullFiledT;
  .ListBayesSpikeOb[[this$NInstance]]$FullFileP <- FullFileP;
  .ListBayesSpikeOb[[this$NInstance]]$FullFileiLoc <- FullFileiLoc;
  .ListBayesSpikeOb[[this$NInstance]]$FullFileProb <- FullFileProb;
  eval(parse(text=LockGText(".ListBayesSpikeOb", "globalenv()")));  
  
 });
ReturnSubProbCodaFromFiles <- function(
    TBSR5 = NULL, TBSCpp = NULL, StartIter = 0, EndIter = 2000,
    Verbose = 0, SubSetCoords = NULL,TBSRoo = NULL,...) {
  if (is.null(TBSR5) && !is.null(TBSRoo)) {
    TBSR5 = TBSRoo;
  }
  CodaList = list();
  SaveDir = NULL;
  if (is.null(SaveDir)) {
    if (!is.null(TBSR5)) {
      if (!is.null(TBSR5$SaveDir)) {
        SaveDir = TBSR5$SaveDir;
      }
    }
  }
  if (is.null(Verbose)) {
    print("ABayesSpikeStart.r:::ReturnSubProbCodaFromFiles: Verbose is Null"); flush.console(); Verbose = 3;
  }
  if (Verbose > 0) {
    print("ABayesSpikeStart.r:::ReturnSubProbCodaFromFiles: Starting"); flush.console();
  }
  if (is.null(TBSRoo) && !is.null(TBSR5)) {
    TBSRoo = TBSR5;
  } else if (is.null(TBSRoo) && is.null(TBSR5) && !is.null(TBSCpp) &&
    !is.null(TBSCpp$TBSR5) ) {
    TBSRoo = TBSCpp$TBSR5;
  } else if (is.null(TBSRoo) && is.null(TBSR5) && !is.null(TBSCpp)) {
    TBSRoo = TBSCpp$TBSRoo;
  }
  if (is.null(TBSRoo)) {
    print("ReturnSubProbCodaFromFiles:  Error Please give non null TBSRoo");
    return(NULL);
  }

  if (is.null(SubSetCoords)) {
    SubSetCoords = 1:TBSRoo$p;
  }
  if (Verbose > 0) {
    print(paste("ReturnSubProbCodaFromFiles: SubSetCoords before fun = (",
      paste(SubSetCoords, collapse = ", "), ")", sep="")); flush.console(); Verbose = 3;
  }
  SubSetCoords = SubSetCoords[SubSetCoords > 0 & SubSetCoords <= TBSRoo$p]
  CSubSetCoords = SubSetCoords -1;
  
  if (is.null(SaveDir)) {
    if (Verbose > 0) {
      print("SaveDir is Null, so can't use it"); flush.console();
    }
    return(-1);
  }
  if (is.null(TBSRoo)) {
    if (Verbose > 0) {
      print("TBSRoo is Null so we can't use it"); flush.console();
    }
  }  else if (is.null(TBSRoo$ListSaveFiles)) {
    if (Verbose > 0) {
      print("TBSRoo is not null, but it does have null ListSaveFiles"); flush.console();
    }
  } else if (is.null(SaveDir)) {
    print("Savedir is NULL Not continuing!");
    return(-1);
  } 
  ##Dir = SaveDir;
  MyNFiles = unlist(list.files(SaveDir));
    iTFiles = MyNFiles[substr(MyNFiles, nchar(MyNFiles)-5, nchar(MyNFiles)) ==
      "iT.bin"];               
    dTFiles = MyNFiles[substr(MyNFiles, nchar(MyNFiles)-5, nchar(MyNFiles)) ==
      "dT.bin"]; 
    PFiles = MyNFiles[substr(MyNFiles, nchar(MyNFiles)-4, nchar(MyNFiles)) ==
      "P.bin"];
      
    ILocFiles = MyNFiles[substr(MyNFiles, nchar(MyNFiles)-nchar("ILoc.bin")+1,
      nchar(MyNFiles)) == "ILoc.bin"]; 
    ProbFiles = MyNFiles[substr(MyNFiles, nchar(MyNFiles)-nchar("Prob.bin")+1,
      nchar(MyNFiles)) == "Prob.bin"];
      
} 

GetG0Text <- function(AText, envir = "globalenv()", S=0) {
  if (S== 0) {
    return(paste("TryO = NULL; try( TryO <- bindingIsActive(\"", AText, "\",  ", envir, " ), silent=TRUE);
    if (!is.null(TryO) && (TryO == TRUE || TryO == FALSE)) {
      try(unlockBinding( \"", AText, "\", ", envir, " ), silent=FALSE);
    }
    try(unlockBinding( \"", AText, "\", ", envir, " ), silent=FALSE); 
    if (exists(\"", AText, "\", ", envir, " )) {
    ", AText, " <- get(\"", AText, "\", ", envir, " );
    } else {
      ", AText, " = 0;
    }
    ", sep=""));
  } else {
    return(paste("TryO = NULL; try( TryO <- bindingIsActive(\"", AText, "\",  ", envir, " ), silent=TRUE);
    if (!is.null(TryO) && (TryO == TRUE || TryO == FALSE)) {
      try(unlockBinding( \"", AText, "\", ", envir, " ), silent=TRUE);
    }
    try(unlockBinding( \"", AText, "\", ", envir, " ), silent=TRUE); 
    if (exists(\"", AText, "\", ", envir, " )) {
    ", AText, " <- get(\"", AText, "\", ", envir, " );
    } else {
      ", AText, " = 0;
    }
    ", sep=""));   
  }
}


setMethodS3("getX", "BayesSpikeRoo", function(this,...) {
   if (!is.null(this$ABayesSpikeCL) &&
     !is.null(this$ABayesSpikeCL$X)) {
     return(this$ABayesSpikeCL$X);
   }
   return(this$.X);
});

setMethodS3("GenerateCodaTable", "BayesSpikeRoo", 
  function(this, DoRecord = 0, NumChains = 0,MBS=NULL, Verbose = 0, ...)  {
  if (is.null(this$MaxGibbsIters)) {
    print("GenerateCodaTable: Please set MaxGibbsIters first"); flush.console();
  }
  if (Verbose > 0) {
    print(paste("GenerateCodaTable: BayesSpikeRoo, setting up, Num Chains = ", NumChains,
      ", Number Iters = ", this$MaxGibbsIters, sep="")); flush.console();
  }
  if (is.null(this$DoRecord) && length(DoRecord) != 7) {
    print("Error: Do Record must be length 7 to make table!");
  } else if (length(DoRecord) == 7) {
    this$DoReord = DoRecord;
  } else if (length(this$DoRecord) != 7) {
    print("Error: Do Record must be length 7 to make table!");
  } else {
    DoRecord = this$DoRecord;
  }
  if (NumChains > 0) {
    this$NumChains = NumChains;
  } else if (this$NumChains > 0) { NumChains = this$NumChains; } else {
    this$NumChains = 3;  NumChains = 3;
  }
  if (sum(abs(DoRecord)) == 0) {
    this$CodaTable = NULL;  this$CodaList = NULL;
    this$ABayesSpikeCL$CodaTable = NULL;
    this$ABayesSpikeCL$CodaList = NULL;  
    if (Verbose > 0) {
      print(paste("GenerateCodaTable: Returning NULL")); flush.console();
    }
    return(NULL);
  }
  Count = 0;  CLN = c();
  if (DoRecord[1] > 0) {
    if (Verbose > 0) {
      print(paste("GenerateCodaTable: DoRecord 1 greater than 0."));
        flush.console();
    }
    Count = Count + this$p
    if (!is.null(names(this$sBeta)) && length(names(this$sBeta)) == this$p) {
      CLN = c(CLN, 
      addname(Names=names(this$sBeta), 
        FN=paste("Beta:", 1:length(this$sBeta),  sep="")));
    } else {
      CLN = c(CLN, paste("Beta:", 1:length(this$sBeta), sep=""));
    }
  }
  if (DoRecord[2] > 0) {
    if (Verbose > 0) {
      print(paste("GenerateCodaTable: DoRecord 2 greater than 0."));
        flush.console();
    }
    Count = Count + length(this$tauEndList);
    if (!is.null(names(this$sOnTau))) {
      CLN = c(CLN, 
      addname(Names=names(this$sOnTau), 
        FN=paste("tau:", 1:length(this$sOnTau), sep="")));
    } else {
      CLN = c(CLN, paste("tau:", 1:length(this$sOnTau), sep=""));
    }
  }
  if (DoRecord[3] > 0) {
    Count = Count + length(this$OnSigma);
    CLN = c(CLN, "Sigma:1")
  }
  if (DoRecord[4] > 0) {
    Count = Count + length(this$OnPiA);
    if (!is.null(names(this$OnPiA))) {
      CLN = c(CLN, 
      addname(Names=names(this$OnPiA), 
        FN= paste("PiA:", 1:length(this$OnPiA),  sep="")));
    } else {
      CLN = c(CLN, paste("PiA:", 1:length(this$OnPiA), sep=""));
    }
  }

  if (DoRecord[5] > 0) {
    if (length(this$tauFixed) > 0) {
      Count = Count + length(this$tauFixed);
      if (!is.null(names(this$tauFixed))) {
        CLN = c(CLN, 
        addname(Names=names(this$tauFixed), 
          FN=paste("TF:", 1:length(this$tauFixed), sep="")));
      } else {
        CLN = c(CLN,
          paste("TF:", 1:length(this$tauFixed), sep=""));
      }
    }
  }
  if (DoRecord[6] > 0) {
    if (this$CFirstRandom > 0 && this$CFirstRandom <= this$p) {
      Count = Count + this$CFirstRandom;
      if (!is.null(names(this$sBeta)) && 
        length(names(this$sBeta)) > this$CFirstRandom) {
        CLN = c(CLN,
          addname(Set=this$sBeta[1:this$CFirstRandom],
            FN=paste("ProbFixed:", 1:this$CFirstRandom,  sep="")));  
      } else {
        CLN = c(CLN, 
          paste("ProbFixed", 1:this$CFirstRandom, sep=""));
      }
    }
  }
  if (DoRecord[7] > 0) {
    if (!is.null(this$tauEndList) && length(this$tauEndList) > 0) {
      Count = Count + length(this$tauEndList);
      if (!is.null(names(this$sOnTau)) && 
        length(names(this$sOnTau)) > length(this$sOnTau) ) {
        CLN = c(CLN,
          addname(Set=this$sOnTau,
            FN=paste("ProbTau:", 1:length(this$sOnTau), sep="")));  
      } else {
        CLN = c(CLN, 
          paste("ProbTau", 1:length(this$sOnTau), sep=""));
      }
    }
  }
  if (Verbose >= 1) {
    print("Other Count starting to setup. "); flush.console();
  }
  OtherCount = this$ABayesSpikeCL$NeedRecord;
  if (Count != OtherCount) {
    print(paste("I expect an error, Count = ", Count, " but OtherCount: ", OtherCount, sep=""))
    print(paste("DoRecord: (", paste(DoRecord, collapse =", "), ")", sep=""));
  }
  if (Verbose >= 2) {
    print("CodaList being setup. "); flush.console();
  }
  this$CodaList = list();
  for (ii in 1:NumChains) {
    if (Verbose >= 2) {
      print(paste("  Setting up Chain number ", ii, sep="")); flush.console();
    }
    this$CodaTable = matrix(0, this$MaxGibbsIters, Count);
    if(!is.null(this$CodaTableNames)) {
      colnames(this$CodaTable) = this$CodaTableNames;
    } else {
      colnames(this$CodaTable) = CLN;
    }
    this$CodaTable = as.mcmc(this$CodaTable);
    this$CodaList[[ii]] = this$CodaTable;
  }  
  this$CodaList = as.mcmc(this$CodaList);

  if (Verbose >= 2) {
    print(paste("  ABayesSpikeCL CodaList being attached to this, Verbose = ", 
      Verbose, sep="")); flush.console();
  }
  this$ABayesSpikeCL$CodaList = this$CodaList;
  if (!is.null(this$CodaList) && length(this$CodaList)  > 0) {
  try(this$ABayesSpikeCL$OnCodaTable <- 1, silent=TRUE);
  }   
  ##this$ABayesSpikeCL$CodaTable = this$CodaTable;
  if (!is.null(MBS)) { MBS$CodaList = this$CodaList; }
  
  eval(parse(text=GetG0Text(".ListBayesSpikeOb", "globalenv()")));
  .ListBayesSpikeOb[[this$NInstance]]$CodaTable <- this$CodaTable;
  .ListBayesSpikeOb[[this$NInstance]]$CodaList <- this$CodaList;
  eval(parse(text=LockGText(".ListBayesSpikeOb", "globalenv()")));  

});
setMethodS3("GenerateEigenList", "BayesSpikeRoo", 
  function(this, AnotherBayesSpikeCL = NULL, ...) {
  if (is.null(this$tauEndList) || length(this$tauEndList) <= 0) {
    return(1);
  }
  if (is.null(this$CtauEndList) || length(this$tauEndList) <= 0) {
    print("GenerateEigenList, why try to do this when it does not exist!");
    flush.console();
    return(1);
  }
  if (!is.null(AnotherBayesSpikeCL)) {
    this$ABayesSpikeCL = AnotherBayesSpikeCL;
  }
  if (is.null(this$ABayesSpikeCL$tauEndList)) {
    print("GenerateEigenList: Error, Roo object End list is not null, BayesSpikeCL end list is!");
    print(paste("this[TBSRoo]$tauEndList = (",
      paste( this$tauEndList, collapse=", "), ")", sep=""));
    flush.console(); return(-1);
  }
  ABayesSpikeCL = this$ABayesSpikeCL;
  if (is.null(ABayesSpikeCL)) {
    print("GenerateEigenList:  Error, no BayesSpikeCL in Roo object ");
    flush.console(); return(NULL);
  } else {
    print("GenerateEigenList: Going with AnotherBayesSpikeCL"); flush.console();
    ABayesSpikeCL = AnotherBayesSpikeCL;
  }
  MSize  = 0;
  if (is.null(ABayesSpikeCL$FirstRandom)) {
     print("error: GenerateEigenList has no FirstRandom!");
     return(NULL);
  }
  Start = ABayesSpikeCL$CFirstRandom + 1;
  X = ABayesSpikeCL$X;  tauEndList = ABayesSpikeCL$tauEndList + 1;
  AllEigenVectorList = list();
  AllEigenValueList = list();
  for (jj in 1:length(tauEndList)) {
    End = tauEndList[jj];
    if (End -Start < 0) {
      tryCatch("Error: GenerateEigenList Bad, disordered tauEndList");
    }
    XXSub = t(X[,Start:End]) %*% X[,Start:End];
    EG = eigen(XXSub);
    AllEigenVectorList[[jj]] = EG$vectors;
    AllEigenValueList[[jj]] = EG$values     
    Start = End +1;
  }
  if (this$Verbose > 1) {
    print("TBSRoo: Assigning the AllEigenValues"); flush.console();
  }
  try(this$AllEigenVectors <- AllEigenVectorList);
  try(this$AllEigenValues <- AllEigenValueList);
  

  if (this$Verbose > 1) {
    print("TBSRoo: EigenValues generated"); flush.console();
  }
  #################################################################
  ##  Code here locks MyR to the BAYESSPIKENAMESPACE environment;
  ## 
  ##BAYESSPIKENAMESPACE <- environment()
  TryO = NULL;
  try( TryO <- bindingIsActive(".ListBayesSpikeOb", BAYESSPIKENAMESPACE), silent = TRUE);
  if (!is.null(TryO) && (TryO == TRUE || TryO == FALSE)) {
   unlockBinding( ".ListBayesSpikeOb" , BAYESSPIKENAMESPACE )
  }
  
  eval(parse(text=GetG0Text(".ListBayesSpikeOb", "globalenv()")));
  .ListBayesSpikeOb <- get(".ListBayesSpikeOb", BAYESSPIKENAMESPACE);
  .ListBayesSpikeOb[[this$NInstance]]$AllEigenVectors <- AllEigenVectorList;
  .ListBayesSpikeOb[[this$NInstance]]$AllEigenValues <- AllEigenValueList;
  eval(parse(text=LockGText(".ListBayesSpikeOb", "globalenv()")));
  
  ##ABayesSpikeCL$TBSRoo = this;
  if (length(.ListBayesSpikeOb[[this$NInstance]]$AllEigenValues)
     != length(this$tauEndList)) {
    print("We're going to have an error in eigenvalues!");
    print(paste(" Length All EigenValues  has length ", length(AllEigenValueList), sep=""));
    for (ii in 1:length(AllEigenValueList)) {
      print(paste("AllEigenValueList[", ii, "] = c(", 
        paste(AllEigenValueList[[ii]], collapse=", "),
        ")", sep="")); flush.console();
    }
    print(paste("tauEndList = c(", paste(tauEndList, collapse=", "),
      ")", sep="")); flush.console();      
    return(NULL);         
  }
  if (this$Verbose > 0) {
    print("Assigning EigenValues now to ABayesSpikeCL "); flush.console();
  }
  try(ABayesSpikeCL$AllEigenVectors <- .ListBayesSpikeOb[[this$NInstance]]$AllEigenVectors);
  try(ABayesSpikeCL$AllEigenValues <- .ListBayesSpikeOb[[this$NInstance]]$AllEigenValues);

  if (this$Verbose > 2) {
    print("GenerateEigenValues:  Finishing"); flush.console();
  }
  return(1);
});


setConstructorS3("PriorStructure", function(Chdf = -1, taubarnu = -1,
  PriorFunc = NULL, MaxX = 500,  StartX = 0, NumPos = 100000,
  DoLog = FALSE, DoSqrt = FALSE) {

    PriorXScaling = c( (MaxX - StartX) / NumPos, 1, StartX,0);
    if (DoLog == TRUE) {
      X = exp(StartX + (1:(NumPos-1))/NumPos * MaxX);
      PriorXScaling[2] = -1;
    } else if (DoSqrt == TRUE) {
      X = (StartX + (1:(NumPos-1))/NumPos * MaxX)^2;
      PriorXScaling[2] = 2;
    } else {
      X = StartX + (1:(NumPos-1))/NumPos * MaxX;
    }
    if (is.null(PriorFunc)) {
      if (Chdf == -1) {
        Chdf = 1;  taubarnu = 1;
      }
      PriorOftau = NULL;  X = NULL; 
    }  else {
      PriorOftau = PriorFunc(X);
    }
  

  extend(Object(), "PriorStructure",
    PriorOftau = PriorOftau, PriorXScaling = PriorXScaling,
    Chdf = Chdf, taubarnu = taubarnu,
    PriorFunc = PriorFunc, MaxX = MaxX, 
    StartX = StartX, NumPos = NumPos,
    DoLog = DoLog, DoSqrt = DoSqrt);
});

FuncMaximizeMe <- function(D,R, MyPriorStruct, MaxIters = 100,
  CauchyEpsilon = .0001, Verbose = 0) {
  sJ = length(D);
  if (!is.numeric(Verbose) && Verbose == TRUE) {Verbose = 1;}
  Outtau = rep(0, MaxIters+2);
  BSMT = .Call("GetBSMT");
  .Call("MaximizeMe", Outtau = Outtau, MaxTauEst = -1,
    PriorOftau = MyPriorStruct$PriorOftau, 
    PriorXScaling = MyPriorStruct$PriorXScaling,
    nu = MyPriorStruct$Chdf,
    taubarnu = MyPriorStruct$taubarnu,
    D=D,R=R, MaxIters = MaxIters, CauchyEpsilon=CauchyEpsilon,
    BSMT = BSMT, Verbose = Verbose);
  return(Outtau);
}

setMethodS3("finalize", "BayesSpikeRoo", function(this, ...) {
  if (!is.null(this$Verbose) && this$Verbose > 0) {
    print("BayesSpikeRoo Is about to be removed from memory"); flush.console();
  }
  if (!is.null(this$ABayesSpikeCL)  && this$ABayesSpikeCL$BeingDestroyed == 0) {
    print("TBSRoo:finalize, we're going to try to call MBS's destructor"); flush.console();
    ##rm(this$ABayesSpikeCL);
    this$ABayesSpikeCL$DeleteMe();
    this$ABayesSpikeCL = NULL;
    gc();
  }
  if (!is.null(this$AFD)  && this$AFD$.BeingDestroyed != 1) {
    finalize(this$AFD);  this$AFD= NULL;
  }

});

################################################################################
##  BayesSpikeSnow
##
##    SnowCall for BayesSpike
##
##   Calls Multiple BayesSpikeRegression algorithms in parrallel.
##   
##
##  A separate BayesSpikeCL is allocated on each node, and operations are
##   saved to parralel buffer files.  After completion of regression,
##  Buffers are returned as SaveFile List.  
##  Use Function "ReadBayesSpikeSnow" to retrieve CodaChains on desired subset
##   of variables.  
##    All of the function calls attached to MBS are in .cpp code (in src)
##    Function calls attached to TBSRoo are located in this document "ABayesSpikeStart.R"
BayesSpikeSnow <- function(NumClusters = 2, X=NULL,Y=NULL, BetaStart = NULL,
  IndexFirstRandomEffect = 1, tauEndList=NULL,
  tauStartValues = NULL, PiAStart = .5, SigmaSq = 1,
  dfTauStart = -1, MStart = -1,
  dfTNoise = 0, Verbose = 0,
  DependenciesTau = NULL, DependenciesFixed=NULL,
  DoMax = 1, 
  PiAPrior = c(-1,-1), HowSample = 3,
  NumSpliceSampleReps = 5, SpliceMove= 1, CauchyEpsilon = .00001,
  MyPS = BayesSpike:::PriorStructureR5(.5,1),  tauFixed = 40,
  TypeFixedPrior = 1, DoRNGState = 1,
  SigmaPrior = c(2,1), InitKKs = 0, DoRecord = c(0,0,0,0,0,0,0),
  NumSamples = 1000, MaxIters = 1000, MaxGibbsIters = 1000,
  MaximizeMeCauchyTotalIters = 100,
  ttStart = 0,  NumChains = 3,
  EarlyEndStep = -1, EarlyEndtt = -1, FileSave = "MyS",
  SaveDir = .DefaultFileDir,  FileName = "ABSChain",
  NewWrite = 1, tauPriordf = 1, tauPriorMean = 1, 
  PutRNG = TRUE, Run = TRUE,
  CodaTableNames = NULL, OnChainIter = 1, RStart = 3,
  NoNoiseBetaStart = FALSE, AFD = NULL, NoShrinkFixed = NULL,
  NoShrinkFixedPrior=NULL, NoShrinkRandom=NULL, NoShrinkRandomPrior=NULL,
  StartRunProbVector = -100,
  PriorProbFixed=NULL, PriorProbTau = NULL, dfRobit = -1, TemperatureList = NULL,
  burnin = 1, EEMergeEvery = NULL, DoEEProbSort = 0, 
  EEProbSortWidth = -1, PrintIter = 200, 
  TemperatureDecreasingRate = 1.0, FirstRandom = -1,
  RandomInfoList = NULL, DoLogitPostProb = FALSE, 
  WriteYBuffer = FALSE, NewYBufferWrite = TRUE,
  WriteWeightBuffer=FALSE, NewWeightBufferWrite=TRUE,
  LengthYBuffer=100, LengthWeightBuffer = 100, EEDoSort = TRUE,
  StartOldEELoad = 50,  AlterWeightFlag = FALSE, DoLongCI = FALSE, 
  CpBuffLongCI = -1, LengthBetaLongBuffer = -1,
  ...) {
  
  ## ExecuteClusters
  ClusterSystem <- makeCluster(NumClusters, type = "MPI") 
  
  RText = paste("
  AnIter = ", 1:NumChains, ";
  AMBS = BayesSpikeRegression(X=X,Y=Y, BetaStart = BetaStart,
  IndexFirstRandomEffect = IndexFirstRandomEffect, tauEndList=tauEndList,
  tauStartValues = tauStartValues, PiAStart = PiAStart, dfTauStart = dfTauStart, MStart = MStart,
  dfTNoise = dfTNoise, Verbose = Verbose,
  DependenciesTau = DependenciesTau, DependenciesFixed=DependenciesFixed,
  DoMax = DoMax, 
  PiAPrior = PiAPrior, HowSample = HowSample,  SigmaSq = SigmaSq, 
  NumSpliceSampleReps = NumSpliceSampleReps, SpliceMove= SpliceMove, 
  CauchyEpsilon = CauchyEpsilon,
  MyPS = MyPS,  tauFixed = tauFixed,
  TypeFixedPrior = TypeFixedPrior, DoRNGState = DoRNGState,
  SigmaPrior = SigmaPrior, InitKKs = InitKKs, DoRecord = c(0,0,0,0,0,0,0),
  NumSamples = NumSamples, MaxIters = MaxIters, MaxGibbsIters = MaxGibbsIters,
  MaximizeMeCauchyTotalIters = MaximizeMeCauchyTotalIters,
  ttStart = ttStart,  NumChains = 1,
  EarlyEndStep = EarlyEndStep, EarlyEndtt = EarlyEndtt, FileSave = FileSave,
  SaveDir = SaveDir,  FileName = \"", FileName, 1:NumChains, "\",
  NewWrite = NewWrite, tauPriordf = tauPriordf, tauPriorMean = tauPriorMean, 
  PutRNG = PutRNG, Run = Run,
  CodaTableNames = CodaTableNames, OnChainIter = OnChainIter, RStart = RStart,
  NoNoiseBetaStart = NoNoiseBetaStart, AFD = AFD, NoShrinkFixed = NoShrinkFixed,
  NoShrinkFixedPrior=NoShrinkFixedPrior, NoShrinkRandom=NoShrinkRandom,
  NoShrinkRandomPrior=NoShrinkRandomPrior,
  StartRunProbVector = StartRunProbVector,
  PriorProbFixed=PriorProbFixed, PriorProbTau = PriorProbTau, dfRobit = dfRobit, 
  TemperatureList = TemperatureList,
  burnin = burnin, EEMergeEvery = EEMergeEvery, DoEEProbSort = DoEEProbSort, 
  EEProbSortWidth = EEProbSortWidth, PrintIter = PrintIter, 
  TemperatureDecreasingRate = TemperatureDecreasingRate, FirstRandom = FirstRandom,
  RandomInfoList = RandomInfoList, DoLogitPostProb = DoLogitPostProb, 
  WriteYBuffer = WriteYBuffer, NewYBufferWrite = NewYBufferWrite,
  WriteWeightBuffer=WriteWeightBuffer, NewWeightBufferWrite=NewWeightBufferWrite,
  LengthYBuffer=LengthYBuffer, LengthWeightBuffer = LengthWeightBuffer, EEDoSort = EEDoSort,
  StartOldEELoad = StartOldEELoad, AlterWeightFlag = AlterWeightFlag,
  DoLongCI = DoLongCI, CpBuffLongCI = CpBuffLongCI, LengthBetaLongBuffer = LengthBetaLongBuffer);
  if (AnIter == 1) {
    return(1);
  }  
  ", sep="");
  
  lAMBS <- clusterEvalQ(ClusterSystem,parse(text=RText));
  AMBS  = lAMBS[[1]];
  
  stopCluster(ClusterSystem)
  
  return(AMBS);
}
