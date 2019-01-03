## A Function that attempts to reload a saved BayesSpikeCpp/TBSR5 Object.
ReLoad <- function(MyList = NULL) {
  if (is.null(MyList)) {
    eval(parse(text=GetG0Text("MyList", S=1)));
  }
  if (is.null(MyList) || length(MyList) == 1) {
    print(
      "ReLoad; BayesSpikeRegression, sorry, no List object to reload from!");
  }
  MyMBS <- NULL;
  try(MyMBS <- BayesSpikeRegression(ReInputList = MyList));
  if (is.null(MyMBS)) {
    print("ReLoad: Error in BayesSpikeRegression Call!");flush.console();
    return(NULL);
  }
  return(MyMBS);
}

BayesSpikeFromMyListText <- function() { 
return(
  "
  try(X <- ReInputList$.X, silent=TRUE);  
  try(Y <-  ReInputList$.Y, silent=TRUE);  
  try(BetaStart <- ReInputList$.BetaStart, silent=TRUE);
  try(IndexFirstRandomEffect <- ReInputList$FirstRandom, silent=TRUE);
  try(tauEndList <- ReInputList$tauEndList, silent=TRUE);
  try(tauStartValues <- ReInputList$sOnTau, silent=TRUE);
  try(PiAStart <- ReInputList$OnPiA, silent=TRUE);
 
  try(PiAPrior <- ReInputList$PiAPrior, silent=TRUE); 
  try(SigmaPrior <- ReInputList$SigmaPrior, silent=TRUE);
  try(dfRobit <- ReInputList$dfRobit, silent=TRUE); 
  try(DoRecord <- ReInputList$DoRecord, silent=TRUE);
  try(dfTNoise <- ReInputList$dfTNoise, silent=TRUE);  
  try(InitFlags <- ReInputList$InitFlags, silent=TRUE);
  try(MaxGibbsIters <- ReInputList$MaxGibbsIters, silent=TRUE); 
  try(NumChains <- ReInputList$NumChains, silent=TRUE);
  try(DependenciesTau <- ReInputList$DependenciesTau, silent=TRUE);
  try(DependenciesFixed <- ReInputList$DependenciesFixed, silent=TRUE);
 
  try(tauPriordf <- ReInputList$tauPriordf, silent=TRUE);  
  try(dfTauStart <- ReInputList$tauPriordf, silent=TRUE);
  try(tauPriorMean <- ReInputList$tauPriorMean, silent=TRUE);  
  try(mStart <- ReInputList$tauPriorMean, silent=TRUE);
  
  try(dfTNoise <- ReInputList$dfTNoise, silent=TRUE);  
  try(ttStart <- ReInputList$ttStart, silent=TRUE);
  try(Verbose <- ReInputList$Verbose, silent=TRUE);  
  try(DoMax <- ReInputList$DoMax, silent=TRUE);
  try(HowSample <- ReInputList$HowSample, silent=TRUE);
  try(NumSpliceSampleReps <- ReInputList$NumSampleReps, silent=TRUE); 
  try(SpliceMove <- ReInputList$SpliceMove, silent=TRUE); 
  try(CauchyEpsilon <- ReInputList$CauchyEpsilon, silent=TRUE);
  try(MyPS <- ReInputList$MyPS, silent=TRUE);  
  try(tauFixed <- ReInputlist$MyPS, silent=TRUE);
  try(TypeFixedPrior <- ReInputList$TypeFixedPrior, silent=TRUE); 
  DoRNGState = ReInputList$DoRNGState;
  InitKKs = ReInputList$InitFlags[3]; DoRecord = ReInputList$DoRecord;
  NumSamples = ReInputList$NumSamples;
  MaxIters = ReInputList$MaxIters;
  MaxGibbsIters = ReInputList$MaxGibbsIters;
  MaximizeMeCauchyTotalIters = ReInputList$MaximizeMeCauchyTotalIters;
  ttStart = ReInputList$ttStart;  NumChains = ReInputList$NumChains;
  EarlyEndStep = -1; EarlyEndtt = -1; FileSave = ReInputList$FileSave;
  SaveDir = ReInputList$SaveDir;  FileName = ReInputList$FileName;
  NewWrite = ReInputList$NewWrite; 
  PutRNG = ReInputList$PutRNG; Run = FALSE;
  CodaTableNames = ReInputList$CodaTableNames;
  OnChainIter = 1;  RStart = ReInputList$RStart;
  NoNoiseBetaStart = TRUE;
  AFD = NULL; 
  try(NoShrinkFixed <- ReInputList$NoShrinkFixed);
  try(NoShrinkFixedPrior <- ReInputList$NoShrinkFixedPrior);
  try(NoShrinkRandom <- ReInputList$NoShrinkRandom);
  try(NoShrinkRandomPrior <- ReInputList$NoShrinkRandomPrior);
  try(StartRunProbVector <- ReInputList$StartRunProbVector);
  try(PriorProbFixed <- ReInputList$PriorProbFixed);
  try(PriorProbTau <- ReInputList$PriorProbTau);
  try(dfRobit <- ReInputList$dfRobit);
  TemperatureList = ReInputList$TemperatureList;
  burnin = ReInputList$burnin; 
  EEMergeEvery = ReInputList$EEMergeEvery;
  DoEEProbSort = ReInputList$DoEEProbSort; 
  EEProbSortWidth = ReInputList$EEProbSortWidth;
  PrintIter = ReInputList$PrintIter; 
  TemperatureDecreasingRate = ReInputList$TemperatureDecreasingRate;
  try(FirstRandom <- ReInputList$FirstRandom);
  RandomInfoList = ReInputList$RandomInfoList;
  DoLogitPostProb = ReInputList$DoLogitPostProb; 
  WriteYBuffer = ReInputList$WriteYBuffer;
  NewYBufferWrite = ReInputList$NewYBufferWrite;
  WriteWeightBuffer=ReInputList$WriteWeightBuffer;
  NewWeightBufferWrite=ReInputList$NewWeightBufferWrite;
  LengthYBuffer=ReInputList$LengthYBuffer;
  LengthWeightBuffer = ReInputList$LengthWeightBuffer;
  EEDoSort = ReInputList$EEDoSort;
  DoLongCI = ReInputList$DoLongCI;
  CpBuffLongCI = ReInputList$CpBuffLongCI;
  StartOldEELoad = ReInputList$StartOldEELoad;
  DoChainIters = ReInputList$DoChainIters;
  AlterWeightFlag = ReInputList$AlterWeightFlag;
  DoLogitNonePostPreProb <- 0;
  try(DoLogitNonePostPreProb <- ReInputList$DoLogitNonePostPreProb);
  if (DoLogitNonePostPreProb == 1) {
    DoLogitPostProb = TRUE;  DoLogitPreProb = FALSE;
  } else if (DoLogitNonePostPreProb == 2) {
    DoLogitPostProb = FALSE; DoLogitPreProb = TRUE;
  }
  

   
  if (is.null(MaxIters) || !is.numeric(MaxIters) ||
    length(MaxIters) <= 0 ) { MaxIters = 100; }
  if (is.null(CauchyEpsilon) || !is.numeric(CauchyEpsilon)) {
    CauchyEpislon = .00001;
  } 
  ");
}
  
  
FinalAttachTryCode <- function() {
  return("
    if (!exists(\"Verbose\")) { Verbose = 1; }
    if (Verbose >= 0) {
      print(\"FinalAttachCode: Try to attach old CodaNames.\");
      flush.console();
    }
    try(MBS$OldCodaNames <- ReInputList$OldCodaNames, silent=TRUE);
    try(MBS$TBSR5$OldCodaNames <- ReInputList$OldCodaNames, silent=TRUE);  
    try(MBS$OtherNameCodaList <- ReInputList$OtherNameCodaList, silent=TRUE); 
    try(MBS$TBSR5$FirstCLN <- ReInputList$FirstCLN, silent=TRUE);
    try(MBS$AllEigenValues <- ReInputList$AllEigenValues, silent=TRUE);
    try(MBS$AllEigenVectors <- ReInputList$AllEigenVectors, silent=TRUE);
    try(MBS$DoLogitNonePostPreProb <- ReInputList$DoLogitNonePostPreProb, silent=TRUE);
    if (exists(\"TimePartsList\") && !is.null(TimePartsList)) {
      try(MBS$TimePartsList <- TimePartsList);
    }
    if (Verbose >= 0) {
      print(\"FinalAttachCode: Try to attach old CodaList.\");
      flush.console();
    }
    if (!is.null(ReInputList$CodaList) && length(ReInputList$CodaList) >= 1 &&
      NCOL(ReInputList$CodaList[[1]]) > 0 && is.null(colnames(ReInputList$CodaList)) ) {
      if (!is.null(MBS$OldCodaNames) && 
        length(MBS$OldCodaNames) == NCOL(ReInputList$CodaList[[1]])) {
        if (Verbose >= 0) {
          print(\"We're going to replace in OldCodaNames for CodaList.\"); flush.console();
        }
        ABC <- list();
        library(coda);
        for (tt in 1:length(ReInputList$CodaList)) {
          ABN <- matrix(as.numeric(ReInputList$CodaList[[tt]]), 
            NROW(ReInputList$CodaList[[tt]]), NCOL(ReInputList$CodaList[[tt]]));
          try(colnames(ABN) <- MBS$OldCodaNames);
          try(ABN <- as.mcmc(ABN));
          try(ABC[[tt]] <- ABN);
        } 
        try(MBS$CodaList <- as.mcmc.list(ABC)); 
        try(ReInputList$CodaList <- as.mcmc.list(ABC));
        try(MBS$TBSR5$CodaList <- MBS$CodaList);
      }  else if (!is.null(MBS$OtherNameCodaList) && 
        length(MBS$OtherNameCodaList) == NCOL(ReInputList$CodaList[[1]])) {
        if (Verbose >= 0) {
          print(\"We're going to replace in OtherNameCodaList for CodaList.\"); flush.console();
        }
        ABC <- list();
        library(coda);
        for (tt in 1:length(ReInputList$CodaList)) {
          ABN <- matrix(as.numeric(ReInputList$CodaList[[tt]]), 
            NROW(ReInputList$CodaList[[tt]]), NCOL(ReInputList$CodaList[[tt]]));
          try(colnames(ABN) <- MBS$OldCodaNames);
          try(ABN <- as.mcmc(ABN));
          try(ABC[[tt]] <- ABN);
        } 
        try(ReInputList$CodaList <- as.mcmc.list(ABC)); 
        try(MBS$CodaList <- as.mcmc.list(ABC));
        try(MBS$TBSR5$CodaList <- MBS$CodaList);
      } else {
        try(MBS$CodaList <- ReInputList$CodaList, silent=TRUE);
        try(MBS$TBSR5$CodaList <- MBS$CodaList);      
      }
    } else {
      try(MBS$CodaList <- ReInputList$CodaList, silent=TRUE);
      try(MBS$TBSR5$CodaList <- ReInputList$CodaList, silent=TRUE);
    }
    if (!is.null(ReInputList$gPrior) && length(ReInputList$gPrior) == 1 &&
      is.numeric(ReInputList$gPrior) && ReInputList$gPrior > 0.0) {
      try(MBS$gPrior <- ReInputList$gPrior);    
    }
    try(MBS$NInstance <- ReInputList$NInstance, silent=TRUE); 
    try(MBS$TBSR5$NInstance <- ReInputList$NInstance, silent=TRUE);  
    if (Verbose >= 0) {
      print(\"FinalAttachCode: Try to attach old ProbCodaList.\");
      flush.console();
    }
    if (!is.null(ReInputList$RunProbVector) && !is.null(ReInputList$RunProbVector)) {
      try(MBS$SetupRunProbVector(MBS$burnin));
      try(MBS$FillRunProbVector(ReInputList$RunProbVector, as.integer(ReInputList$TotEveryProbVector)));
      if (!is.null(ReInputList$RegionWidth) && length(ReInputList$RegionWidth) == 1 &&
        as.numeric(ReInputList$RegionWidth) > 1 && !is.null(ReInputList$RunProbRegionVector) &&
        length(ReInputList$RunProbRegionVector) >= 1) {
        MBS$SetRegionWidthAndRunProbRegionVector(ReInputList$RegionWidth, 
          ReInputList$RunProbRegionVector);              
      }
      if (!is.null(ReInputList$RunSumBeta) && length(ReInputList$RunSumBeta) == MBS$p &&
        !is.null(ReInputList$RunMoreZero) && !is.null(ReInputList$RunLessZero)) {
          try(MBS$SetupRunLessZeroMoreZeroSumBeta(ReInputList$RunLessZero,
            ReInputList$RunMoreZero, ReInputList$RunSumBeta));
        }
    }
    try(MBS$TBSR5$.ProbCodaList <- ReInputList$.ProbCodaList, silent=TRUE);
    try(MBS$ProbCodaList <- ReInputList$.ProbCodaList, silent=TRUE);
    try(MBS$TBSR5$.SigmaCodaList <- ReInputList$.SigmaCodaList, silent=TRUE);
    try(MBS$SigmaCodaList <- ReInputList$.SigmaCodaList, silent=TRUE);
    try(MBS$TBSR5$.TauCodaList <- ReInputList$.TauCodaList, silent=TRUE);
    try(MBS$TauCodaList <- ReInputList$.TauCodaList, silent=TRUE);    
    try(MBS$TBSR5$.PiACodaList <- ReInputList$.PiACodaList, silent=TRUE);
    try(MBS$PiACodaList <- ReInputList$.PiACodaList, silent=TRUE);    
    try(MBS$TBSR5$.PiACodaList <- ReInputList$.PiACodaList, silent=TRUE);
    try(MBS$PiACodaList <- ReInputList$.PiACodaList, silent=TRUE); 
    try(MBS$InitiateTime <- ReInputList$InitiateTime, silent=TRUE); 
    try(MBS$SecondInitiateTime <- ReInputList$SecondInitiateTime, silent=TRUE); 
    try(MBS$CompleteTime <- ReInputList$CompleteTime, silent=TRUE); 
    try(MBS$OnSigma <- ReInputList$OnSigma, silent=TRUE);
    try(MBS$RefreshBeta(ReInputList$Beta));
    if (Verbose >= 0) {
      print(\"FinalAttachCode: Try to setup PostProbBuffer.\");
      flush.console();
    }
    if (!is.null(ReInputList$OtherNameCodaList) && length(ReInputList$OtherNameCodaList) >= 1) {
      try(MBS$OtherNameCodaList <- ReInputList$OtherNameCodaList);
    }
    if (ReInputList$LengthPostProbBuffer > 0) {     
      try(MBS$TBSR5$LengthPostProbBuffer 
        <- ReInputList$LengthPostProbBuffer, silent=TRUE);
      try(MBS$LengthPostProbBuffer <- 
        ReInputList$LengthPostProbBuffer, silent=TRUE);   
    } 
    if (!is.null(ReInputList$SubSetCoords) &&
      length(ReInputList$SubSetCoords) > 0 &&
       max(ReInputList$SubSetCoords) >= 0) {     
    try(MBS$TBSR5$SubSetCoords <- ReInputList$SubSetCoords, silent=TRUE);
    try(MBS$SubSetCoords <- ReInputList$SubSetCoords, silent=TRUE); 
    }  
    if (!is.null(ReInputList$SubSetTau) &&
      length(ReInputList$SubSetTau) > 0 &&
       max(ReInputList$SubSetTau) >= 0) {   
    try(MBS$TBSR5$SubSetTau <- ReInputList$SubSetTau, silent=TRUE);
    try(MBS$SubSetTau <- ReInputList$SubSetTau, silent=TRUE);  
    }
    if (Verbose >= 0) {
      print(\"FinalAttachCode: Try to add YCodaList\"); flush.console();
    }
    if (!is.null(ReInputList$.YCodaList) && 
      length(ReInputList$.YCodaList) > 0) {
      try(MBS$TBSR5$.YCodaList <- ReInputList$.YCodaList, silent=TRUE);
      try(MBS$YCodaList <- ReInputList$.YCodaList, silent=TRUE); 
    }
    if (!is.null(ReInputList$SubCodaList) &&
      length(ReInputList$SubCodaList) > 0  &&
      max(ReInputList$SubCodaList) >= 0) {
    try(MBS$TBSR5$SubCodaList <- ReInputList$SubCodaList, silent=TRUE);
    try(MBS$SubCodaList <- ReInputList$SubCodaList, silent=TRUE); 
    }
    if (!is.null(ReInputList$SubYList) &&
      length(ReInputList$SubYList) > 0  &&
      max(ReInputList$SubYList) >= 0) {
    try(MBS$TBSR5$SubYList <- ReInputList$SubYList, silent=TRUE);
    try(MBS$SubYList <- ReInputList$SubYList, silent=TRUE);
    }
    
    if (!is.null(ReInputList$KeepPosteriorQuantiles) &&
      is.numeric(ReInputList$KeepPosteriorQuantiles)) {
      try(MBS$TBSR5$.KeepPosteriorQuantiles <- ReInputList$KeepPosteriorQuantiles);
    }
    if (!is.null(ReInputList$BetaSymmetricQuantiles) &&
      is.matrix(ReInputList$BetaSymmetricQuantiles) &&
      length(ReInputList$BetaSymmetricQuantiles) > 0) {
      try(MBS$TBSR5$.BetaSymmetricQuantiles <- ReInputList$BetaSymmetricQuantiles);
    }
    if (!is.null(ReInputList$BetaSymmetricUnshrinkQuantiles) &&
      is.matrix(ReInputList$BetaSymmetricUnshrinkQuantiles) &&
      length(ReInputList$BetaSymmetricUnshrinkQuantiles) > 0) {
      try(MBS$TBSR5$.BetaSymmetricUnshrinkQuantiles <- ReInputList$BetaSymmetricUnshrinkQuantiles);
    }
    if (!is.null(ReInputList$TauSymmetricQuantiles) &&
      is.matrix(ReInputList$TauSymmetricQuantiles) &&
      length(ReInputList$TauSymmetricQuantiles) > 0) {
      try(MBS$TBSR5$.TauSymmetricQuantiles <- ReInputList$TauSymmetricQuantiles);
    }
    if (!is.null(ReInputList$SigmaSymmetricQuantiles) &&
      is.vector(ReInputList$SigmaSymmetricQuantiles) &&
      length(ReInputList$SigmaSymmetricQuantiles) > 0) {
      try(MBS$TBSR5$.SigmaSymmetricQuantiles <- ReInputList$SigmaSymmetricQuantiles);
    }
    if (!is.null(ReInputList$PiASymmetricQuantiles) &&
      is.matrix(ReInputList$PiASymmetricQuantiles) &&
      length(ReInputList$PiASymmetricQuantiles) > 0) {
      try(MBS$TBSR5$.PiASymmetricQuantiles <- ReInputList$PiASymmetricQuantiles);
    }
    if (!is.null(ReInputList$BetaHPDQuantiles) &&
      is.matrix(ReInputList$BetaHPDQuantiles) &&
      length(ReInputList$BetaHPDQuantiles) > 0) {
      try(MBS$TBSR5$.BetaHPDQuantiles <- ReInputList$BetaHPDQuantiles);
    }
    if (!is.null(ReInputList$BetaHPDUnshrinkQuantiles) &&
      is.matrix(ReInputList$BetaHPDUnshrinkQuantiles) &&
      length(ReInputList$BetaHPDUnshrinkQuantiles) > 0) {
      try(MBS$TBSR5$.BetaHPDUnshrinkQuantiles <- ReInputList$BetaHPDUnshrinkQuantiles);
    }
    if (!is.null(ReInputList$TauHPDQuantiles) &&
      is.matrix(ReInputList$TauHPDQuantiles) &&
      length(ReInputList$TauHPDQuantiles) > 0) {
      try(MBS$TBSR5$.TauHPDQuantiles <- ReInputList$TauHPDQuantiles);
    }
    if (!is.null(ReInputList$SigmaHPDQuantiles) &&
      is.vector(ReInputList$SigmaHPDQuantiles) &&
      length(ReInputList$SigmaHPDQuantiles) > 0) {
      try(MBS$TBSR5$.SigmaHPDQuantiles <- ReInputList$SigmaHPDQuantiles);
    }
    if (!is.null(ReInputList$PiAHPDQuantiles) &&
      is.matrix(ReInputList$PiAHPDQuantiles) &&
      length(ReInputList$PiAHPDQuantiles) > 0) {
      try(MBS$TBSR5$.PiAHPDQuantiles <- ReInputList$PiAHPDQuantiles);
    }
    if (!is.null(ReInputList$CenteredColumns) && length(ReInputList$CenteredColumns) >= 1) {
      try(MBS$CenteredColumns <- ReInputList$CenteredColumns);
      if (!is.null(ReInputList$DeCenteredCodaList) && length(ReInputList$DeCenteredCodaList) >= 1) {
        ##try(MBS$DeCenteredCodaList <- ReInputList$DeCenteredCodaList);
        try(MBS$DeCenteredCodaList <- NULL);
        try(DCC <- MBS$DeCenteredCodaList);
      }
    }
    if (!is.null(ReInputList$SubCodaLongList)) {
      try(MBS$TBSR5$.SubCodaLongList <- ReInputList$SubCodaLongList);
    }
    if (!is.null(ReInputList$BetaSymmetricLongQuantiles) &&
      is.matrix(ReInputList$BetaSymmetricLongQuantiles) &&
      length(ReInputList$BetaSymmetricLongQuantiles) > 0) {
      try(MBS$TBSR5$.BetaSymmetricLongQuantiles <- ReInputList$BetaSymmetricLongQuantiles);
    }
    if (!is.null(ReInputList$BetaHPDLongQuantiles) &&
      is.matrix(ReInputList$BetaHPDLongQuantiles) &&
      length(ReInputList$BetaHPDLongQuantiles) > 0) {
      try(MBS$TBSR5$.BetaHPDLongQuantiles <- ReInputList$BetaHPDLongQuantiles);
    }   
    if (Verbose >= 0) {
      print(\"FinalAttachCode: Finished.\"); flush.console();
    }
    ##try(print(\"FinalAttachCode: Really Finished.\")); flush.console();
   ");
}