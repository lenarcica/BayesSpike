
################################################################################
##  RunRegression(MBS)
##
##  If MBS is a Preconstructed BayesDiallel Object, this attemps to continue 
##  running the BayesSpikeRegression on the object (including multiple chains)
##
##  MBS$RunAlgorithm (BayesSpikeGibbs.cpp) will only run a chain to its end 
##
##  This is essentially all of the functions from BayesSpikeRegression
##
##

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

RunRegression <- function(MBS=NULL, MBSAFD=NULL, OnChainIter = -1, TBSR5 = NULL,
  SaveDir= NULL, SaveAFD = FALSE, TBSRoo = NULL, ...) {
  if ((!exists("MBS") ||is.null(MBS)) && exists("MBSAFD") && !is.null(MBSAFD)) {
    MBS = AFD$.MBS;
  }
  if (!exists("MBS") || is.null(MBS)) {
    print(paste("RunRegression function; Out Of Luck.  ",
      "no BayesSpikeCLL object supplied.  ", sep="")); flush.console();
    return(-1);
  }
  if (!exists("Verbose")) {
    Verbose = NULL;
    try(Verbose <- MBS$Verbose);
    if (!is.null(Verbose) || (!is.numeric(Verbose) && !is.integer(Verbose))) {
      print(paste("RunRegression: Hey, could not extract Verbose",
        " from MBS!", sep="")); flush.console();
    }
    Verbose = 3;
  }
  TBSR5 = MBS$TBSR5;
  if (OnChainIter < 1 ) {
    OnChainIter = 1;
  } else if (OnChainIter < 1) { OnChainIter = 1; }  
  StChainIter = OnChainIter;

  STT = OnChainIter:(length(MBS$CodaList));

  MBS$CheckX(MBS$TBSR5$.X); MBS$CheckY(MBS$TBSR5$.Y); MBS$CheckBeta(MBS$TBSR5$sBeta);
  if (is.null(MBS$TBSR5)) {
    print("RunRegression: Cannot run with TBSR5 is null!"); flush.console();
  }
  if (!is.null(MBS$DoRecord) && sum(abs(MBS$DoRecord)) > 0 &&
    (is.null(MBS$CodaList) || length(MBS$CodaList) <= 0)) {
    MBS$TBSR5$GenerateCodaTable();
  }
  if (is.null(MBS$TBSR5$ABayesSpikeCL)) {
    print("RunRegression: ERROR: after GenerateCodaTable: M$T$M is NULL!");
    flush.console();
    MBS$TBSR5$ABayesSpikeCL = MBS;
  }
  
if (!is.null(MBS$TemperatureList)) {
  if (Verbose > -1) {
    print("BayesSpikeRegression: Installing a nonNull TemepratureList"); flush.console();
  }
  if (MBS$Tempii <= 0)   {
    ItTemps = 1:length(MBS$TemperatureList);
  } else {
    ItTemps = (MBS$Tempii+1):length(MBS$TemperatureList);
  }
}  else {
  ItTemps = 1;
}
try(MBS$SecondInitiateTime <- proc.time());
for (TTii in ItTemps)  {
  if (TTii > 1) {
    try(MBS$TBSR5$AllTempCodaLists[[TTii]] <- MBS$CodaList);
    try(MBS$AllTempCodaLists <- MBS$TBSR5$AllTempCodaLists);
  }
  if (!is.null(MBS$TemperatureList) && length(MBS$TemperatureList) > 1) {
    try(MBS$Tempii <- TTii - 1);
  }
  if (TTii == length(ItTemps) && sum(MBS$DoRecord[6:7]) < 2) {
    MBS$SetupRunProbVector(MBS$burnin);
    if (MBS$TBSR5$RegionWidth >= 1) {
      MBS$RegionWidth <- MBS$TBSR5$RegionWidth;
    }
  }
  if (MBS$Verbose > 0) {
    print(paste("RunRegression:  On Temp = ", TTii, "/", length(ItTemps), sep=""));
  }
  if (TTii > ItTemps[1]) {
    OnChainIter = 1;
  }
  if (sum(abs(MBS$TBSR5$DoRecord) > 0)) {
    if (Verbose > -1) {
      print("BayesSpikeRegression: DoRecord > 0 so Generating Coda Table"); flush.console();
    }
    BZ = -1;
    ATryText = "
      MBS$TBSR5$GenerateCodaTable(Verbose = Verbose);
      BZ = 1;
      ";
    try(eval(parse(text=ATryText)));
    if (is.null(BZ) || !is.numeric(BZ) || BZ != 1) {
      print("BayesSpikeRegression: TTii=", TTii, "/", length(ItTemps), ": something went wrong trying to generate coda table!.");
      flush.console();
    }
    if (Verbose > -1) {
      print(paste("BayesSpikeRegression: Finished Generating CodaTable length = ", length(MBS$CodaList),
        ":(", paste(dim(MBS$CodaList[[1]]), collapse=", "), "), TTii=", TTii,
        "/", length(ItTemps), sep="")); flush.console();
    }  
  }
if (MBS$TBSR5$DoSave == FALSE) {
  try(MBS$NoSave <- -1);
  try(MBS$DoSave <- 0);
}
try(MBS$SecondInitiateTime <- proc.time())
for (OnChainIter in STT) {
  ATrySaveText <- "
    if (!is.null(MBS$SaveDir) && is.character(MBS$SaveDir) &&
      MBS$SaveDir[1] != \"\" && MBS$TBSR5$DoSaveTBSR5==TRUE) {
      try(MBS$Save(), silent=TRUE);  
    }
  ";
  if (OnChainIter <= 0 || OnChainIter > MBS$TBSR5$NumChains) {
    print("Oh No, OnChainIter, not very good!"); flush.console();
    ERRORSTT = STT;
    eval(parse(text="ERRORSTT", "globalenv()", S=1));
    return(MBS);
  }
  try(eval(parse(text=ATrySaveText)), silent=TRUE);
  
  if (TBSR5$PutRNG == 1) {
    MBS$RNGState = 1;
  } else if (isNull(TBSR5$PutRNG) || TBSR5$PutRNG == 0) {

  } else {
    MBS$RNGState = 1;
  }
  if (!is.null(MBS$DoRecord) && sum(abs(MBS$DoRecord) > 0) &&
    (is.null(MBS$CodaList) || length(MBS$CodaList) == 0)) {
    print(paste("RunRegression: Error, Temp = ", TTii, ", OnChainIter  ", OnChainIter,
      "/", length(STT), " but although DoRecord is On, CodaList is Gone!",
      sep="")); flush.console();
    return(MBS);  
  }
  if (!is.null(MBS$CodaList) && length(MBS$CodaList) >= OnChainIter) {
    try(MBS$OnCodaTable <- OnChainIter);
  }
  ##if (!is.null(MBS$TemperatureList) && OnChainIter > 1) {
  ##  MBS$Tempii = OnChainIter -1;
  ##}
  if (TBSR5$Verbose > 0) {
    print(paste("Starting at beginning, TTii = ", TTii, "/", length(ItTemps), 
      ", OnChainIter = ", MBS$OnCodaTable, sep=""));
    flush.console();
  }
  try(MBS$TBSR5$OnPiA[1:length(MBS$TBSR5$OnPiA)] <-
     MBS$TBSR5$PiAStart[1:length(MBS$TBSR5$PiAStart)]);
  try(MBS$OnPiA <- MBS$TBSR5$OnPiA, silent=TRUE);
  if (length(MBS$OnPiA) == 1 && (MBS$OnPiA <= 0.0 || MBS$OnPiA > 1.0)) {
    try(MBS$OnPiA <- .5, silent=TRUE);
  } 
  RfFlag = 0;
  if (ZeroOutBeforeSample == 5) {
    print(paste("BayesSpikeRegression: You Want to ",
      "ZeroOutBeforeSample, but here is MBS", sep="")); flush.console();
    return(MBS);
  }
  if (MBS$p >= 5000 || MBS$ZeroOutBeforeSample == 1) {
    MyTry <- 0;
    try(MyTry <- ZeroOutFunction(MBS));
    if (MyTry == 0 || MyTry == -1) {
      print(paste("BayesSpikeRegression:: Fail ZeroOutFunction gives Problem. ", 
        MyTry, sep=""));  flush.console();
      return(MBS);
    }
  } else if (!is.null(MBS$TBSR5$BetaStart)) {
    if (TBSR5$Verbose > 1) {
      print("Refresh MBS$Beta using non null TBSR5$BetaStart"); flush.console();
    }
    try(RfFlag <- MBS$RefreshBeta( MBS$TBSR5$BetaStart + rnorm(length(MBS$TBSR5$BetaStart),
      sd(MBS$TBSR5$BetaStart) * MBS$TBSR5$RStart)));
  } else if (!is.null(MBS$TBSR5$BetaStart)) {
    if (MBS$TBSR5$Verbose > 1) {
      print("Refresh MBS$Beta using random RStart"); flush.console();
    }
    try(RfFlag <- MBS$RefreshBeta(MBS$TBSR5$RStart * rnorm(MBS$TBSR5$pCoef) ));
  }
  if (!is.null(MBS$dfTNoise) && MBS$dfTNoise > 0) {
    MBS$iiWeight = rep(1, MBS$n);
  }
  
  if (MBS$p <= 10000) {
    MBS$UpdateNonFreshXtResid();
  } else {

  }
    MBS$UpdateSigma();
  ##if (RfFlag == 0) {
  ##   try(MBS$RefreshsBeta( TBSR5$BetaStart + rnorm(length(TBSR5$BetaStart),
  ##     sd(TBSR5$BetaStart) * 3)));  
  ##}
  ##if (RfFlag == 0) {
  ##  print("MBS$RefreshBeta Error: Returning MBS"); flush.console();
  ##  return(MBS);
  ##}
  if (MBS$TBSR5$Verbose > 0) {
    print(paste("BayesSpikeRegression: TTii = ", TTii, "/", length(ItTemps),
      ", and On Chain ", TBSR5$OnChainIter, "/", length(STT), ".", 
    " STT = (", paste(STT, collapse=", "), ")", sep=""));
    flush.console();
  }
  if (sum(abs(MBS$DoRecord) > 0)) {
    if (MBS$Verbose > 1) {
      print("Refresh MBS$CodaTable"); flush.console();
    }
    if (!is.null(MBS$CodaList)  && length(MBS$CodaList) >= OnChainIter) {
      try(MBS$OnCodaTable <- OnChainIter);
    }
    ##MBS$CodaTable = MBS$TBSR5$CodaList[[OnChainIter]];
    ##MBS$CodaList = MBS$TBSR5$CodaList;
  }
  MBS$tt = MBS$TBSR5$ttStart;
  if (TBSR5$Verbose > 0) {
    print(paste("BayesSpikeRegression: Setup Save file Chain ", TBSR5$OnChainIter, sep=""));
    flush.console();
  }
  if (is.null(MBS)) {
    print(paste("BayesSpikeRegression: error on chain iter ", TBSR5$OnChainIter,
      ",  MBS went NULL!", sep="")); flush.console();
    return(TBSR5);
  }
  if (is.null(MBS$BeingDestroyed) || !is.numeric(MBS$BeingDestroyed)) {
    print(paste("BayesSpikeRegression: error on chain iter ", TBSR5$OnChainIter,
      ",  MBS has null BeingDestroyed!", sep="")); flush.console();
    return(TBSR5);
  }
  if (MBS$BeingDestroyed == 1) {
    print(paste("BayesSpikeRegression: error on chain iter ", TBSR5$OnChainIter,
      ",  MBS has 1 BeingDestroyed!", sep="")); flush.console();
    return(TBSR5);  
  }
  if (is.null(MBS$TBSR5)) {
    print(paste("BayesSpikeRegression: error on chain iter ", TBSR5$OnChainIter,
      "MBS has NULL TBSR5!", sep="")); flush.console();
    return(TBSR5);  
  }
  if (is.null(MBS$TBSR5$ABayesSpikeCL)) {
     print(paste("BayesSpikeRegression: error on chain iter ", TBSR5$OnChainIter,
      ",  MBS has TBSR5 has NULL ABayesSpikeCL!", sep="")); flush.console();
     try(MBS$TBSR5$ABayesSpikeCL <- MBS);
  }
  if (nchar(MBS$TBSR5$FileName) <= 0) {
    print(paste("BayesSpikeRegression: error Filename = ", 
      MBS$TBSR5$FileName, sep="")); flush.console();
  }
  try(MBS$TBSR5$SetupSaveFile(SaveDir=SaveDir, FileName = MBS$TBSR5$FileName,
    chainiter = OnChainIter));

  if (MBS$Verbose >= 0) {
    print(paste("BayesSpikeRegression: Running Algorithm, TTii=", TTii, "/",
      length(ItTemps), ", and Chain: ", 
      OnChainIter, "/", length(STT), ", after SetupSaveFile.", sep=""));
    flush.console();
  }
  ##if (OnChainIter > 1) {
  ##  MBS$Verbose = 6;
  ##}
  MyTry = 0;
  ##print(paste("Beginning of Algorithm, MBS$Verbose = ", MBS$Verbose, " and tt = ", MBS$tt, sep=""));
  ##print(paste("CodaTable[1,] = (", paste(round(MBS$CodaTable[1,],2), collapse=", "), sep=""));
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
  try(MyTry <- MBS$RunAlgorithm());
  if (MyTry == 0) {
    print("Looks Like MBS$RunAlgorithm had an error, take a look"); flush.console();
    return(MBS);
  }
  if (!is.null(MBS$AlterWeightFile) && is.character(MBS$AlterWeightFile) && MBS$AlterWeightFile  != "") {
    MyAT <- -1;
    try(MyAT <- MBS$DeriveAlterProbability(MBS$TBSR5$LengthAlterWeightBuffer));
    if (!is.null(MyAT) && MyAT == -1) {
      print("RunRegression, DeriveAlterProbability returns poorly, return MyAT. "); flush.console();
    }
    return(MBS);
  }
  ##print(paste("End of Algorithm, MBS$Verbose = ", MBS$Verbose, " and tt = ", MBS$tt, sep=""));
  ##print(paste("CodaTable[1,] = (", paste(round(MBS$CodaTable[1,],2), collapse=", "), sep=""));
  if (MBS$Verbose > 1) {
    print("Trying to Set Old Temperature List files "); flush.console();
  }
  if (!is.null(MBS$TemperatureList) && !is.null(MBS$sCodaDFile) && MBS$Temperature != 1.0) {
     MBS$sCodaOldDFile = paste(MBS$sCodaDFile, sep="");
  }
  if (!is.null(MBS$TemperatureList) && !is.null(MBS$sCodaIFile)&& MBS$Temperature != 1.0) {
     MBS$sCodaOldIFile = paste(MBS$sCodaIFile, sep="");
  }
  if (!is.null(MBS$TemperatureList) &&  !is.null(MBS$CodaILocFile)&& MBS$Temperature != 1.0) {
    try(MBS$SetCodaOldLocAndProbFile(MBS$CodaILocFile, MBS$CodaProbFile,
      MBS$burnin, MBS$MaxGibbsIters, OnChainIter, MBS$TBSR5$DoEEProbSort));
  }
  ATrySaveText <- "
    if (!is.null(MBS$SaveDir) && is.character(MBS$SaveDir) &&
      MBS$SaveDir[1] != \"\" && MBS$TBSR5$DoSaveTBSR5==TRUE) {
      try(MBS$Save(), silent=TRUE);  
    }
  ";
  try(eval(parse(text=ATrySaveText)), silent=TRUE);
  if (MyTry == 0) {
    print("MyTry returns early exit.");
    return(MBS);
  }
  if (MBS$TBSR5$Verbose > 0) {
    print(paste("BayesSpikeRegression: Finished Algorithm, Chain: ", 
      TBSR5$OnChainIter, sep="")); flush.console();
  }
  if (MBS$TBSR5$PutRNG == 1) {
    MBS$PutRNGState();
  } else if (!is.null(MBS$TBSR5$PutRNG) || MBS$TBSR5$PutRNG == 0) {
  
  } else {
    MBS$PutRNGState();
  }
}  
  
  if (MBS$TBSR5$PutRNG == 1) {
    MBS$PutRNGState();
  } else if (!is.null(MBS$TBSR5$PutRNG) || MBS$TBSR5$PutRNG == 0) {
  
  } else {
    MBS$PutRNGState();
  }
  ##MBS$TBSR5 = TBSR5; TBSR5$MBS = MBS;


  try(AD <- MBS$CheckTBSR5(TBSR5));
  if (AD < 0) {
    print("RunRegression: We checked TBSR5 and we got some problems!");flush.console();
    A2 = 0;
    try(min(c(MBS$CheckX(MBS$TBSR5$.X), MBS$CheckY(MBS$TBSR5$.Y),
      MBS$CheckBeta(MBS$TBSR5$sBeta))));
    if (A2 < 0) {
      print("RunRegression, we also have a failure in A2 "); flush.console();
    } else {
      print("But no Failure in individual components so we're safe!");
    }
  }
  ## Make sure lock in the background
  ##TryO = NULL;
  ##try( TryO <- bindingIsActive(".ListBayesSpikeOb", BAYESSPIKENAMESPACE), silent = TRUE);
  ##if (!is.null(TryO) && (TryO == TRUE || TryO == FALSE)) {
  ## unlockBinding( ".ListBayesSpikeOb" , BAYESSPIKENAMESPACE )
  ##}
  ##NLL <- get(".ListBayesSpikeOb", BAYESSPIKENAMESPACE);
  ##NLL[[TBSR5$NInstance]] <- TBSR5;
  ##assign( ".ListBayesSpikeOb", NLL, BAYESSPIKENAMESPACE )
  ##lockBinding( ".ListBayesSpikeOb", BAYESSPIKENAMESPACE )
 
   if (FALSE) {
   if (!is.null(MBS$AFD) && SaveAFD == TRUE) { 
   if (exists("WatchCoordinates") && is.null(WatchCoordinates) || length(WatchCoordinates) == 0) {
     for (ii in 1:length(MBS$CodaList)) {
       if (!is.null(MMS$CodaList[[ii]])) {
         CD = length(MMS$CodaList[[ii]][1,]);
         CM = length(MBS$AFD$AllDiallelObs[[1]]$CodaChains[[ii]][1,]);
         CDL = length(MMS$CodaList[[ii]][,1]);
         CML = length(MBS$AFD$AllDiallelObs[[1]]$CodaChains[[ii]][,1]);
         MyCM = min(c(CD, CM));
         MyCML = min(c(CDL, CML));
         AFD$AllDiallelObs[[1]]$CodaChains[[ii]][1:MyCML, 1:MyCM] = 
           MyBMS$CodaList[[ii]][1:MyCML, 1:MyCM];    
       } 
   
     }
   }  else {
     for (ii in 1:MBS$TBSR5$NumChains) {
       if (is.null(MBS$TBSR5)) {  try(MBS$TBSR5  <- TBSR5); }
       if (is.null(MBS$TBSR5$ABayesSpikeCL)) {
         print(paste("BayesSpikeRegression: ERROR -------------- In Loop ii=", ii,
         "/", MBS$TBSR5$NumChains, "  NULL M$T$M!",sep=""));  flush.console();
         try(MBS$TBSR5$ABayesSpikeCL <- MBS);
        }
        MBS$TBSR5$SetupSaveFile(MBS$TBSR5$SaveDir, 
           FileName = MBS$TBSR5$FileName, chainiter = ii);
        RTM = .Call("GiveCodaSubset", WatchCoordinates, 
          paste(MBS$sCodaIFile, sep=""), paste(MBS$sCodaDFile, sep=""), 0, MaxGibbsIters, 0)
        if (!is.null(RTM) && length(dim(RTM)) == 2 && 
         exists("WatchCoordinates") && dim(RTM)[2] == length(WatchCoordinates)) {
         CD = length(RTM[1,]);
         CM = length(MBS$AFD$AllDiallelObs[[1]]$CodaChains[[ii]][1,]);
         CDL = length(RTM[,1]);
         CML = length(MBS$AFD$AllDiallelObs[[1]]$CodaChains[[ii]][,1]);
         MyCM = min(c(CD, CM));
         MyCML = min(c(CDL, CML));
         MBS$AFD$AllDiallelObs[[1]]$CodaChains[[ii]][1:MyCML, 1:MyCM] = 
           RTM[1:MyCML, 1:MyCM];    
       }  else {
         print("Somehow, not all of the WatchCoordinates were recorded!");
         flush.console();  return(MBS);
       }
     }
   }}
   }
  }
  try(MBS$CompleteTime <- proc.time());
  ATrySaveText <- "
    if (!is.null(MBS$SaveDir) && is.character(MBS$SaveDir) &&
      MBS$SaveDir[1] != \"\" && MBS$TBSR5$DoSaveTBSR5==TRUE) {
      try(MBS$Save(), silent=TRUE);  
    }
  ";
  try(eval(parse(text=ATrySaveText)), silent=TRUE);
}


TBSR5SetVerbose <- function(TBSR5, sVerbose) {
  TBSR5$Verbose[1] = sVerbose[1];
}
