BayesSpikeOnChainIterInitiate <- function() {
    return("
   OnChainIter = STT[iti]
      
  ##if (!is.null(MBS$TemperatureList) && OnChainIter > 1) {
  ##  MBS$Tempii = OnChainIter -1;
  ##}
  if (is.null(MBS)) {
    print(\"ERROR: Way Weird BayesSpikeRegression Error, in middle MBS is null!\");
    return(MBS);
  }
  if (is.numeric(MBS$MaxGibbsIters) && MBS$MaxGibbsIters <= 0) {
    print(paste(\"Issue: BayesSpikeRegression: MBS$MaxGibbsIters = \",
      MBS$MaxGibbsIters, sep=\"\"));
    print(\"Seek Why.\")
    return(MBS);
  }
  if (MBS$Verbose > 0) {
    print(paste(\"BayesSpikeRegression,  Temp = \", TTii, \"/\", length(ItTemps),
      \", and iti = \", iti, \"/\", length(STT), \".\", sep=\"\")); flush.console();
  }
  if (MBS$DoLogitNonePostPreProb == 2) {
    try(MBS$dfRobit <- 9.0);  try(MBS$OnSigma <- 7 / 9 * pi^2 /3);
    try(MBS$AlterdfRobit <- -9.0);
  } else if (MBS$DoLogitNonePostPreProb == 1) {
    try(MBS$dfRobit <- 0.0); try(MBS$AlterdfRobit <- 9.0);
    try(MBS$OnSigma <- 7/9 * pi^2 /3);
  }

  ##if (min(c(
  ##  MBS$CheckX(MBS$TBSR5$.X), 
  ##  MBS$CheckY(MBS$TBSR5$.Y), MBS$CheckBeta(MBS$TBSR5$sBeta))) < 0) {
  ##  print(\"Running Algorithm Error,check TBSroo: \"); flush.console();
  ##  return(MBS);  
  ##}
  if (MBS$TBSR5$PutRNG == 1) {
    MBS$RNGState = 1;
  } else if (isNull(TBSR5$PutRNG) || TBSR5$PutRNG == 0) {

  } else {
    MBS$RNGState  = 1;
  }
  try(MBS$TBSR5$OnChainIter <- OnChainIter);
  if (MBS$TBSR5$DoSave == FALSE || 
    (MBS$TBSR5$SaveDir %in% c(\"NoSave\", \"NOSAVE\") ) ) {
    try(MBS$NoSave <- -1); 
    try(MBS$DoSave <- 0); 
  } else {
    try(MBS$NoSave <- 0);
    try(MBS$DoSave <- 1);
  }

  if (sum(abs(MBS$DoRecord)) >= 1) {
    if (is.null(MBS$CodaList)) {
      print(\"DOH, Looks like from here CodaList was never setup though it should be!\");
      flush.console();
    }
    try(MBS$OnCodaTable <- as.integer(OnChainIter));
  }
  if (TBSR5$Verbose > 0) {
    print(paste(\"Starting at beginning, OnChainIter = \", MBS$TBSR5$OnChainIter, sep=\"\"));
    flush.console();
  }
  RfFlag = 0;
  ");
}

BayesSpike0021SetupToGo <- function() {
 return("
     ##if (RfFlag == 0) {
  ##   try(MBS$RefreshsBeta( TBSR5$BetaStart + rnorm(length(TBSR5$BetaStart),
  ##     sd(TBSR5$BetaStart) * 3)));  
  ##}
  ##if (RfFlag == 0) {
  ##  print(\"MBS$RefreshBeta Error: Returning MBS\"); flush.console();
  ##  return(MBS);
  ##}
  if (TBSR5$Verbose > 0) {
    print(paste(\"BayesSpikeRegression: On Chain \", TBSR5$OnChainIter, sep=\"\"));
    flush.console();
  }
  if (MBS$DoLogitNonePostPreProb == 2) {
    MBS$RobitReplace();
    MBS$AlterInitRun <- 0;
  } else if (MBS$DoLogitNonePostPreProb == 1) {
    MBS$dfRobit <- 0.0;
    MBS$AlterWeightdfRobit <- 9.0;
    MBS$RobitReplace();
    MBS$AlterInitRun <- 0;
  } else if (MBS$AlterdfRobit > 0 && is.null(MBS$iiWeight)) {
    MBS$AlterInitRun <- 0.0;
    print(paste(\" BayesSpike0021SetupToGo: No, AlteredRobit is \", 
      MBS$AltereddfRobit, \" but iiWeight is NULL\", sep=\"\"));
    flush.console();
    print(paste(\"   Try Fill AverageIIWeight anyway. \", sep=\"\")); flush.console();
    if (!is.null(MBS$TBSR5$AverageIIWeight)) {
      try(MBS$iiWeight <- MBS$TBSR5$AverageIIWeight)
    }
    try(print(paste(\"   And now average of MBS$iiWeight is \", 
      mean(MBS$iiWeight), sep=\"\")));
      flush.console();
  }

  if (sum(abs(MBS$DoRecord) > 0)) {
    if (TBSR5$Verbose > 1) {
      print(\"Refresh MBS$CodaTable\"); flush.console();
    }
    ##MBS$CodaTable = MBS$TBSR5$CodaList[[OnChainIter]];
    ##MBS$CodaList = MBS$TBSR5$CodaList;
    if (!is.null(MBS$CodaList) && length(MBS$CodaList) >= OnChainIter) {
      try(MBS$OnCodaTable <- OnChainIter);
    }
    MBS$MaxGibbsIters = length(MBS$CodaList[[1]][,1]);
  }
  MBS$tt = MBS$TBSR5$ttStart;
  if (MBS$TBSR5$Verbose > 0) {
    print(paste(\"0014InsideBayesSpikeRegressionLoop:BayesSpikeRegression: Setup Save file Chain \", TBSR5$OnChainIter, sep=\"\"));
    flush.console();
  }
  if (is.null(MBS$TBSR5$ABayesSpikeCL) ||
    !is.numeric(MBS$TBSR5$ABayesSpikeCL$BeingDestroyed) ||
    MBS$TBSR5$ABayesSpikeCL$BeingDestroyed[1] == 1) {
    print(paste(\"BayesSpikeRegression: Before Second setup file, we have \",
      \"MBS$TBSR5$ABayesspikeCL has problems, find a fix, return MBS.\", sep=\"\"));
    flush.console();
    return(MBS);  
  }
  AD = -1;
  try(AD <- MBS$TBSR5$SetupSaveFile(SaveDir=MBS$TBSR5$SaveDir, FileName = MBS$TBSR5$FileName,
    chainiter = OnChainIter));
  if (AD < 0) {
    print(\"BayesSpikeRegression: SetupSaveFile ended in Error, we return MBS!\");
    return(MBS);
  }
  if (MBS$TBSR5$Verbose >= 0) {
    if (length(STT) < OnChainIter) {
      print(paste(\"BayesSpikeRegression: when OnChainIter = \", OnChainIter,
        \" but STT is \", paste(STT, collapse=\", \"), sep=\"\"));
      flush.console();
    }
    print(paste(\"BayesSpikeRegression: Running Algorithm, TTii=\", TTii, \"/\",
      length(ItTemps), \", and Chain: \", 
      OnChainIter, \"/\", length(STT), \".\", sep=\"\"));
    flush.console();
  }
  ##if (OnChainIter > 1) {
  ##  MBS$Verbose = 6;
  ##}
  MyTry = 0;
  
  if (!is.null(MBS$CodaList)) {
    if (!is.null(colnames(MBS$CodaList[[1]]))) {
      BackupColnames <-  colnames(MBS$CodaList[[1]]);
    } else {
      BackupColnames <- NULL;
    }
  } else {
    BackupColnames <- NULL;
  }
  ##try(MBS$UpdateNonFreshXtResid());
 ");   
}

BayesSpike0021AfterRunAlgorithmTest <- function() {
  return("
   if (MyTry == 0) {
    print(\"**************************************************************************** Error Issue! \"); flush.console();
    print(\"*** BayesSpikeRegression: Looks like we caught an error during RunAlgorithm, here is MBS\"); flush.console();
    print(paste(\"*** MyTry = \", MyTry, \" and MBS$tt = \", MBS$tt, sep=\"\")); flush.console();
    return(MBS);
  }

  if (!is.null(MBS$TemperatureList) && !is.null(MBS$sCodaDFile) && MBS$Temperature != 1.0) {
     MBS$sCodaOldDFile = paste(MBS$sCodaDFile, sep=\"\");
  }
  if (!is.null(MBS$TemperatureList) && !is.null(MBS$sCodaIFile) && MBS$Temperature != 1.0) {
     MBS$sCodaOldIFile = paste(MBS$sCodaIFile, sep=\"\");
  }
  if (!is.null(MBS$AlterWeightFile) && is.character(MBS$AlterWeightFile) && MBS$AlterWeightFile  != \"\") {
    ART = 0;
    try(ART <- MBS$DeriveAlterProbability(MBS$TBSR5$LengthAlterWeightBuffer));
    if (!is.numeric(ART) || ART != 1) {
      print(\"BayesSpikeRegression: Fail on DeriveAlterProbability in BayesSpike0021AfterRunAlgorithmTest  Return MBS.\"); flush.console();
      eval(parse(text=SetGText(\"MBS\", \"globalenv()\", S=1)));
      MyText = \"1 = 2\"; eval(parse(text=MyText));
      return(MBS);
    }
  }
  if (Verbose > -1) {
    if (!is.null(MBS$TemperatureList)) {
      if (is.null(MBS$CodaILocFile)) {
        print(\"Doh: MBS$CodaILocFile is null even though no null TemepratureList\"); flush.console();  
        return(-1);
      }
      if (!is.null(MBS$TemperatureList) &&
        length(MBS$TemperatureList) > 1 && is.null(MBS$CodaProbFile)) {
        print(\"Doh: MBS$CodaILocFile:  MBS$CodaProbFile is null despite TemperatureList\"); flush.console();
        return(-1);
        
      }

    }
    
  }
  ")
}

BayesSpike0020ZeroOutNow <- function() {
  return("
  if (!exists(\"ZeroOutBeforeSample\")) { ZeroOutBeforeSample <- 0; }
  if (ZeroOutBeforeSample == 5) {
    print(\"MBS: You want to return before ZeroOutBeforeSample, returning\");
    return(MBS);
  }
  try(MBS$tt <- 0);
  print(\"$$ InZeroOutBeforeSample about to exectue. \"); flush.console();
  try(print(paste(\"$$ We are OnChainIter \",iti,\"/\", length(STT), \" and Temperature \", TTii, \"/\", length(ItTemps), sep=\"\")));
  if (MBS$p >= 20000 || MBS$ZeroOutBeforeSample == 1) {
     print(\"MBS: We are running a ZeroOutBeforeSample because of a RefreshBeta Condition. \"); flush.console();
     MyTry <- 0;
     
     try(MyTry <- ZeroOutFunction(MBS));
     if (MyTry == 0 || MyTry == -1) {
        print(paste(\"MBS: Some failure in ZeroOutFunction: (\", MyTry, 
          \") we return MBS as a fail\", sep=\"\")); flush.console(); 
        eval(parse(text=SetGText(\"MBS\", \"globalenv()\", S=1)));
        1=0;
        return(MBS);
     }
  } else if (!is.null(MBS$TBSR5$BetaStart)) {
    if (TBSR5$Verbose > 1) {
      print(\"Refresh MBS$Beta using non null TBSR5$BetaStart\"); flush.console();
    }
    if (is.logical(NoNoiseBetaStart) && NoNoiseBetaStart == TRUE)  {
      try(MBS$DoAddCoordsOnSetBeta <- 1);
      try(MBS$Beta <- MBS$TBSR5$BetaStart);      
    }  else {
      try(RfFlag <- MBS$RefreshsBeta( MBS$TBSR5$BetaStart + rnorm(length(MBS$TBSR5$BetaStart),
        sd(MBS$TBSR5$BetaStart) * MBS$TBSR5$RStart)));
    }
  } else if (is.null(MBS$TBSR5$BetaStart)) {
    if (MBS$TBSR5$Verbose > 1) {
      print(\"Refresh MBS$Beta using random RStart\"); flush.console();
    }
    try(RfFlag <- MBS$RefreshBeta(MBS$TBSR5$RStart * rnorm(MBS$TBSR5$pCoef) ));
  }
  if (sum(abs(MBS$DoRecord)) > 0 && (is.null(MBS$CodaList) || length(MBS$CodaList) < iti)) {
    print(\"0014InsideBayesSpikeRegressionLoop.r:BayesSpike0020ZeroOutNow() Error \");
    print(paste(\"  Sum DoRecord = \", sum(abs(MBS$DoRecord)), \" but length(CodaList) is \",
      length(MBS$CodaList), sep=\"\")); flush.console();
    ERRORMBS <- MBS;
    eval(parse(text=SetGText(\"ERRORMBS\", \"globalenv()\", S=1)));
    tryCatch(\"0014InsideBayesSpikeRegressionLoop.r!  CodaList not configured.\");
  } 
  ")
}

BayesSpike00XXFirstZeroOut <- function() {
  return("
  
print(\"$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\")
print(\"$$ ZeroOutFunction\"); flush.console();
print(\"$$ We are on 00XXFirstZeroOut. \"); flush.console();
try(print(paste(\"$$ We are OnChainIter \",iti,\"/\", length(STT), \" and Temperature \", TTii, \"/\", length(ItTemps), sep=\"\")));
print(\"$$ Hey: MBS$ is large before a sample we are going to blank up and sample small\");
flush.console();
MBS$DoAddCoordsOnSetBeta <- 0;
RfFlag <- MBS$RefreshBeta(rep(0,MBS$p));
MyTryProbs <- NULL;
MyTryProbTau <- NULL;  MyTryProbFixed <- NULL;
if (MBS$TBSR5$OnChainIter >= 2) {
  try(MyTryProbs <- MBS$ProbVector);
  if (!is.null(MyTryProbs)) {
    if (is.null(MBS$tauEndList)) {
      MyTryProbFixed <- log(MyTryProbs/(1.0-MyTryProbs));  
      MyTryProbFixed[(MyTryProbFixed > 0 & !is.finite(MyTryProbFixed)) | MyTryProbs == 1.0] <- 4;
      MyTryProbFixed[(MyTryProbFixed < 0 & !is.finite(MyTryProbFixed)) | MyTryProbs == 0.0] <- -4;
    } else if (MBS$FirstRandom > 1) {
      MyTryProbFixed <- log(MyTryProbs[1:MBS$iFirstRandom]/(1.0-MyTryProbs[1:MBS$iFirstRandom]))
      MyTryProbFixed[(MyTryProbFixed > 0 & !is.finite(MyTryProbFixed)) | MyTryProbs[1:MBS$iFirstRandom] == 1.0] <- 4;
      MyTryProbFixed[(MyTryProbFixed < 0 & !is.finite(MyTryProbFixed)) | MyTryProbs[1:MBS$iFirstRandom] == 0.0] <- -4;
    }
    if (!is.null(MBS$tauEndList)) {
      MyTryProbTau <- log(MyTryProbs[MBS$FirstRandom:length(MyTryProbs)] /
        (1.0 - MyTryProbs[MBS$FirstRandom:length(MyTryProbs)]));
      MyTryProbTau[(MyTryProbTau > 0 & !is.finite(MyTryProbTau)) | 
        MyTryProbs[MBS$FirstRandom:length(MyTryProbs)] >= 1.0] <- 4;
      MyTryProbTau[(MyTryProbTau < 0 & !is.finite(MyTryProbTau)) |
        MyTryProbs[MBS$FirstRandom:length(MyTryProbs)]  <= 0.0] <- -4;
    }
  }  
}
try(MBS$BlankAllNewCoords());
  if (length(MBS$tauEndList) >= 1) {
    if (!is.null(MyTryProbTau)) {
      MyProb <- MyTryProbTau;
    } else {
      MBS$SampleTausOnly <- 1;
      MBS$SampleNewTaus();
      MyProb <- MBS$ProbTau;
    }
    GoTau  <- MBS$OnTau;
    print(\"$$ ZeroOut before Sample Taus. \"); flush.console();
    try(print(paste(\"$$ We are OnChainIter \",iti,\"/\", length(STT), \" and Temperature \", TTii, \"/\", length(ItTemps), sep=\"\")));
    print(paste(\"$$ ZeroOutFunction: We sample Taus and the number of nonzero are \", length(GoTau[GoTau>0]),
     \"/\", length(GoTau), sep=\"\")); flush.console();
    Fac = 1;
    if (length(MBS$tauEndList) >= 2)   {
      TauLens <- MBS$tauEndList[2:length(MBS$tauEndList)] -
        MBS$tauEndList[1:(length(MBS$tauEndList)-1)];
      TauLens <- c(TauLens, MBS$tauEndList- MBS$FirstRandom+1);
      Fac <- floor(mean(TauLens));
    } else if (length(MBS$tauEndList) == 1) {
      Fac <- MBS$tauEndList[1] - MBS$FirstRandom+1;
    }
    OnPiAT = MBS$OnPiA[1];
    if (length(MBS$OnPiA) == 2) { OnPiAT <- MBS$OnPiA[2]; }
    KeepCount <- round(5*length(GoTau)  * OnPiAT);
    if (KeepCount >= .5 * length(GoTau)) { KeepCount = round(.5* length(GoTau)); }
    if (KeepCount >= length(GoTau[GoTau>0])) { KeepCount <- length(GoTau[GoTau>0]); }
    if (KeepCount <= 0) { KeepCount = 1; }     
    if (KeepCount < length(GoTau[GoTau>0])) {
      Expt <- GiveMeNewBalance(lX=MyProb/2.0, WantTot = KeepCount, StartMove=1.0, MaxTTs = 300);
      if (length(Expt[Expt > 0]) < KeepCount) {
        KeepCount <- length(Expt[Expt > 0]);
      }
      print(paste(\"$$ after calculating Expt we have its sum is \", round(sum(Expt),4),
        \"/\", length(Expt), sep=\"\")); flush.console();
      Keepers <- sample(1:length(GoTau), prob=Expt, replace=FALSE, size=KeepCount);
      print(paste(\"$$ Now we reduce active set to length(Keepers) = \", length(Keepers), sep=\"\"));
      flush.console();
      DoTau <- rep(0, length(GoTau));
      DoTau[Keepers] <- GoTau[Keepers];
      try(MBS$DoAddCoordsOnSetBeta <- 0);
      try(MBS$BlankAllNewCoords());
      print(paste(\"AssignTau: about to Assign \")); flush.console();
      try(MBS$AssignTau(DoTau));
      print(paste(\"$$ Reduced DoTau length(Keepers) = \", length(Keepers), \" is added\", sep=\"\"));
      print(paste(\"$$ After AssignTau, AllNewCoords = \", MBS$AllNewCoords, sep=\"\")); flush.console();
      flush.console();
    } else {
      print(paste(\"$$ KeepCount = \", KeepCount, \" but length(GoTau) = \",
        length(GoTau[GoTau>0]), \"/\", length(GoTau), \" so no change.\", sep=\"\"));
      flush.console(); 
      try(MBS$BlankAllNewCoords());
      print(paste(\"$$ Now Assigning Tau to GoTau\")); flush.console();
      MBS$AssignTau(GoTau);
      print(paste(\"$$ After AssigningTau to GoTau length \", length(GoTau), \" we \",
        \" have AllNewCoords = \", MBS$AllNewCoords, sep=\"\"));
      flush.console();
    }
  }
  print(\"$$ Well Looking into the length of MBS$tauEndList\"); flush.console();
  ")
}