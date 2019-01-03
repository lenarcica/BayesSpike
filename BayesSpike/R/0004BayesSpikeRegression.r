## BayesSpikeRegression.r
##
##  Alan Lenarcic, 2011-2019
##
##   The "BayesSpikeRegression()" main R function takes a large list of
##  default and optional parameters setting many possible T/Logit/Robit++ prior
##  details.  To manage defaults for all of these, and determine which parameters
##  should be active when others aren't, there elements of BayesSpikeRegression
##  from setting defaults, to extracting data, to setting up temperatures,
##  to generating an Rcpp Object and R5 object with data,
##  to the actual Gibbs runs, come in functional parts.

## EvalReDeclares()
##  Running eval(parse(text=EvalReDeclarers)) within a function scope, this
##  function should make sure all necessary variables are declared.
EvalReDeclarers <- function() {
  return("##  Redeclarers:
##
##  Admittedly I always start big functions with default declarers so I can
##  work with the interior of the function to develop it
##  It would be nice to develop a smart workaround someday.
  if (!exists(\"MaximizeMeCauchyTotalIters\")) { 
    MaximizeMeCauchyTotalIters = 100; }
  if (!exists(\"tauStartValues\")) { tauStartValues = NULL; }
  if (!exists(\"MStart\")) { MStart = -1; }
  if (!exists(\"DoTimePartsList\")) { DoTimePartsList = TRUE; }
  if (!exists(\"NoShrinkFixed\")) { NoShrinkFixed = NULL; }
  if (!exists(\"NoShrinkRandom\")) { NoShrinkRandom = NULL; }
  if (!exists(\"NoShrinkFixexPrior\")) { NoShrinkFixedPrior = NULL; }
  if (!exists(\"NoShrinkRandomPrior\")) { NoShrinkRandomPrior = NULL; }
  if (!exists(\"InitKKs\")) { InitKKs = 10; }
  if (!exists(\"RegionWidth\")) { RegionWidth = -1; }
  if (!exists(\"ZeroOutBeforeSample\")) { ZeroOutBeforeSample = 0; }
  if (!exists(\"PiAStart\")) { PiAStart = .5; }
  if (!exists(\"tauPriorMean\")) { tauPriorMean = 1; }
  if (!exists(\"SigmaPrior\")) { SigmaPrior= c(2,1); }
  if (!exists(\"tauPriordf\")) { tauPriordf = 1; }
  if (!exists(\"TypeFixedPrior\")) { TypeFixedPrior = 1;}
  if (!exists(\"NewWrite\")) { NewWrite = 1; }
  if (!exists(\"RStart\")) { RStart = 3; }
  if (!exists(\"dfTauStart\")) { dfTauStart = -1; }
  if (!exists(\"DoMax\")) { DoMax = 1; }
  if (!exists(\"HowSample\")) { HowSample = 3; }
  if (!exists(\"PiAPrior\")) { PiAPrior  = c(-1,-1);  }
  if (!exists(\"SaveDir\")) { SaveDir = .DefaultFileDir; }
  if (!exists(\"BetaStart\")) { BetaStart = NULL; }
  if (!exists(\"NoNoiseBetaStart\")) { NoNoiseBetaStart = FALSE; }
  if (!exists(\"ReInputList\")) { ReInputList = -1; }
  if (!exists(\"AlterWeightFlag\")) { AlterWeightFlag = FALSE; }
  if (!exists(\"DoChainIters\")) { DoChainIters = NULL; }
  if (!exists(\"StartOldEELoad\")) { StartOldEELoad = 50;}
  if (!exists(\"EEDoSort\")) { EEDoSort = TRUE; }
  if (!exists(\"LengthWeightBuffer\")) { LengthWeightBuffer = 100; }
  if (!exists(\"LengthYBuffer\")) { LengthYBuffer = 100; }
  if (!exists(\"NewWeightBufferWrite\")) { NewWeightBufferWrite = TRUE; }
  if (!exists(\"WriteWeightBuffer\")) { WriteWeightBuffer = FALSE; }
  if (!exists(\"NewYBufferWrite\")) { NewYBufferWrite = TRUE; }
  if (!exists(\"WriteYBuffer\")) { WriteYBuffer = FALSE; }
  if (!exists(\"DoLogitPostProb\")) { DoLogitPostProb = FALSE; }
  if (!exists(\"DoLogitPreProb\")) { DoLogitPreProb = FALSE; }
  if (!exists(\"RandomInfoList\")) { RandomInfoList = NULL;}
  if (!exists(\"FirstRandom\")) { FirstRandom = -1; }
  if (!exists(\"TemperatureDecreasingRate\")) { TemperatureDecreasingRate = 1.0; }
  if (!exists(\"NumSamples\")) { NumSamples = 1000; }
  if (!exists(\"MaxIters\")) { MaxIters = 1000;}
  if (!exists(\"PreRunMaxGibbsIters\")) { PreRunMaxGibbsIters = 20; }
  if (!exists(\"MaxGibbsIters\")) { MaxGibbsIters = 1000; }
  if (!exists(\"MaximizeMeCauchyTotalIters\")) { MaximizeMeCauchyTotalIters = 100; }  
  if (!exists(\"NumSpliceSampleReps\")) { NumSpliceSampleReps = 5; }
  if (!exists(\"SpliceMove\")) { SpliceMove = 1; }
  if (!exists(\"CauchyEpsilon\")) { CauchyEpsilon = .00001;}
  if (!exists(\"DoLongCI\")) { DoLongCI = TRUE; }
  if (!exists(\"LengthBetaLongBuffer\")) { LengthBetaLongBuffer = -1; }
  if (!exists(\"CpBuffLongCI\")) { CpBuffLongCI = LengthBetaLongBuffer; }
  if (!exists(\"DoSampleBFixedii\")) { DoSampleBFixedii = 1; }
  if (!exists(\"PrintIter\")) { PrintIter = 200; }
  if (!exists(\"StartRecordingiiWeight\")) { StartRecordingiiWeight = 3; }
if (!exists(\"dfTNoise\")) { dfTNoise = -1; }
if (!exists(\"tauEndList\")) { tauEndList = NULL; }
if (!exists(\"DoSave\")) { DoSave = TRUE; }
if (!exists(\"FileSave\")) { FileSave = \"DefaultSaveFile\"; }
if (!exists(\"NumChains\")) { NumChains = 3; }
if (!exists(\"AFD\")) { AFD = NULL; } 
if (!exists(\"PutRNG\")) { PutRNG = TRUE; }
if (!exists(\"DoRecord\")) { DoRecord = c(1,0,0,0,0,0,0); }
if (!exists(\"DependenciesTau\")) { DependenciesTau = NULL; }
if (!exists(\"DependenciesFixed\")) { DependenciesFixed = NULL; }
if (!exists(\"EarlyEndStep\")) { EarlyEndStep = -1; }
if (!exists(\"EarlyEndtt\")) { EarlyEndtt = -1; }
if (!exists(\"ttStart\")) { ttStart = 1; }
if (!exists(\"CodaTableNames\")) { CodaTableNames = NULL; }
if (!exists(\"WriteYBuffer\")){WriteYBuffer=FALSE;}
if (!exists(\"StartRunProbVector\")) { StartRunProbVector <- 20; }
if (!exists(\"NewYBufferWrite\")){NewYBufferWrite=TRUE;}
if (!exists(\"WriteWeightBuffer\")){WriteWeightBuffer=FALSE;}
if (!exists(\"NewWeightBufferWrite\")){NewWeightBufferWrite=FALSE;}
if (!exists(\"LengthWeightBuffer\")){LengthWeightBuffer=FALSE;}
if (!exists(\"LengthWeightBuffer\")){LengthWeightBuffer=FALSE;}
if (!exists(\"LengthWeightBuffer\")){LengthWeightBuffer=FALSE;}
if (!exists(\"DoSampleBFixedii\")){ DoSampleBFixedii = 1; }
if (!exists(\"RevertTemperatureEvery\")) { RevertTemperatureEvery = -1; }
if (!is.integer(LengthYBuffer)) { LengthYBuffer <- 100;}
if (!is.integer(LengthWeightBuffer)) { LengthWeightBuffer <- 100;}
if (!exists(\"EEDoSort\") || !is.logical(EEDoSort)) { EEDoSort = TRUE;}
if (!exists(\"DoSaveTBSR5\")) { DoSaveTBSR5 = TRUE;}
if (!exists(\"StartOldEELoad\") || !is.numeric(StartOldEELoad) ||
  StartOldEELoad < 0) {  StartOldEELoad = 50; }
if (!exists(\"gPrior\")) { gPrior = 0.0; }
if (!exists(\"SaveDir\") || is.null(SaveDir) || is.numeric(SaveDir) ||
  !is.character(SaveDir) || SaveDir == \"NoSave\" ||
  SaveDir == \"NOSAVE\" ||  DoSave == FALSE) {
  SaveDir = \"NOSAVE\";   DoSave = FALSE;
}
  if (is.character(SaveDir) && (SaveDir %in% c(\"NoSave\", \"NOSAVE\")) ) {
    SaveDir = \"NOSAVE\";  DoSave = FALSE;
  } else if (!is.logical(DoSave) && is.numeric(DoSave) && DoSave > 0) {
    DoSave = TRUE;
  } else if (is.numeric(DoSave) && DoSave <= 0) {
    DoSave = FALSE;
  } else if (is.null(DoSave)) { DoSave = TRUE; }
  
  if (!exists(\"MyPS\") || is.null(MyPS) ||
    !is.numeric(MyPS$Chdf) || !is.numeric(MyPS$MaxX) ) { 
    MyPS = BayesSpike:::PriorStructureR5$new(.5,1); }
  if (!exists(\"tauFixed\")) { 
    if (Verbose >= -1) {
      print(\"BayesSpikeRegression: By Default setting tau Fixed to 40\"); flush.console();
    }
    tauFixed = 40; 
  }    

  if (is.list(ReInputList) && length(ReInputList) > 0) {
    DoneByList = TRUE;
    try(eval(parse(text=BayesSpikeFromMyListText())));
  }  else {
    DoneByList = FALSE;
  }
  ");  
}

## EvalXandY: check that X/Y are valid inputs for a regression of desired type.
##  Also checks that Priors are determined.
EvalXandY <- function() {
  return("
  if (is.null(X)) {
   print(\"BayesSpikeRegression: Please Non Null X\");
}
if (is.null(Y)) {
  print(\"BayesSpikeRegression: Please Non Null Y\")
}
mStart = MStart;
if (dfTauStart > 0) { tauPriordf = dfTauStart; }
if (mStart > 0) { tauPriorMean = mStart}
if (!exists(\"Verbose\")) { try(Verbose  <- 1); }
if (Verbose > 0) {
  print(\"BayesSpikeRegression: Start by allocating R5 object\"); flush.console();
}
p = length(X[1,]);
if (PiAPrior[1] > 0 && length(PiAStart) == 1) {
  PiAGuess = as.numeric(PiAPrior[1] / (PiAPrior[1] + PiAPrior[2]) );
} else if (length(PiAPrior) == 4 && length(PiAStart) == 1 && PiAStart == .5) {
  PiAGuess = as.numeric(c( PiAPrior[1] / (PiAPrior[1] + PiAPrior[2]),
     PiAPrior[3] / (PiAPrior[3] + PiAPrior[4]) ));
} else {
  PiAGuess = as.numeric(PiAStart);
}
if (PiAGuess <= 0 || PiAGuess >= 1.0) {
  PiAGuess = .5;
}
if (!is.null(NoShrinkFixed) && length(NoShrinkFixed) >= 2) {
  NoShrinkFixed <- sort(unique(NoShrinkFixed));
}
if (!is.null(NoShrinkRandom) && length(NoShrinkRandom) >= 2) {
  NoShrinkRandom <- sort(unique(NoShrinkRandom));
}
if (Verbose > 1) {
  print(\"BayesSpikeRegression: -- On NoShrinkFixedPrior\"); flush.console();
}
if (is.null(NoShrinkFixedPrior) && !is.null(NoShrinkFixed)) {
  NoShrinkFixedPrior = rep(100, length(NoShrinkFixed));
} else if (!is.null(NoShrinkFixedPrior) && is.null(NoShrinkFixed)) {
   print(\"Error: Supply No ShrinkFixed, if you supply a prior spread for variables\");
   return(-1);
} else if (length(NoShrinkFixedPrior) < length(NoShrinkFixed)) {
  NoShrinkFixedPrior = rep(NoShrinkFixedPrior[1], length(NoShrinkFixed));
} else if (length(NoShrinkFixedPrior) > length(NoShrinkFixed)) {
  NoShrinkFixedPrior = NoShrinkFixedPrior[1:length(NoShrinkFixed)];
}
if (is.null(NoShrinkRandomPrior) && !is.null(NoShrinkRandom)) {
  NoShrinkRandomPrior = rep(100, length(NoShrinkRandom));
} else if (!is.null(NoShrinkRandomPrior) && is.null(NoShrinkRandom)) {
   print(\"Error: Supply No ShrinkRandom, if you supply a prior spread for variables\");
   return(-1);
} else if (length(NoShrinkFixedPrior) < length(NoShrinkRandom)) {
  NoShrinkRandomPrior = rep(NoShrinkRandomPrior[1], length(NoShrinkRandom));
} else if (length(NoShrinkRandomPrior) > length(NoShrinkRandom)) {
  NoShrinkRandomPrior = NoShrinkRandomPrior[1:length(NoShrinkRandom)];
}
if (!exists(\"IndexFirstRandomEffect\") && exists(\"FirstRandom\")) { 
  IndexFirstRandomEffect = FirstRandom; 
}
if (FirstRandom >= 0 && IndexFirstRandomEffect < 0) { 
IndexFirstRandomEffect = FirstRandom; }
if (InitKKs <= 0) {
  if (p < 300) {
    InitKKs = p;
  } else {
    p = max(6,min( PiAGuess * p, 300 ));
  }
}
if (!exists(\"Sigma\")) { Sigma = 1.0; }
if (any(is.nan(Y))) {
  print(\"Error: Lots of Nan Y\");
  print(Y);
  return(-1);
}
if (any(is.nan(X))) {
  print(\"Error: Lots of Nan X\");
  print(X);
  return(-1);
}
if (!exists(\"dfRobit\")) { dfRobit = -1; }

## DoLogitPostProb: Specifically, this is to study Logistic regresssion
##  which is approximated, as per wikipedia with a Robit distribution
##  of df = 9.0 and s^2 = 7/9 * pi^2 / 3;
if (exists(\"DoLogitPreProb\") && !is.null(DoLogitPreProb) &&
  (DoLogitPreProb[1] == TRUE || DoLogitPreProb[1] >= .5)) {
  DoLogitPreProb = TRUE;  DoLogitPostProb <- FALSE;
  if (all(Y %in% c(-1.0,1.0)))  {
    Y[Y == -1.0] = 0.0;
  }  else if (all(Y %in% c(0.0,1.0))) {
  } else {
    print(\"BayesSpikeRegression():::DoLogitPreProb=TRUE, Hey, Not All Y are 0,1, this will not be good. \"); flush.console();
    DoLogitPreProb = FALSE;
  }
  if (is.null(dfRobit) || dfRobit < 0 || dfRobit == 9.0) {
    dfRobit <- 9;
    SigmaRobit <-  sqrt(7 / 9 * pi^2 / 3);
  } else {
    print(\"BayesSpikeRegression():::Hey, dfRobit is already set though DoLogitPreProb is TRUE, this won't work out that well. \"); flush.console();
    SigmaRobit <- sqrt( 7/ 9 * pi^2 / 3);
  }    
} else if (exists(\"DoLogitPostProb\") &&
  !is.null(DoLogitPostProb) && (DoLogitPostProb[1] == TRUE ||
    DoLogitPostProb[1] >= .5)) {
  DoLogitPreProb <- FALSE;
  if (!is.logical(DoLogitPostProb)) {
    DoLogitPostProb <- TRUE;
  }
  if (all(Y %in% c(-1.0,1.0)))  {
    Y[Y == -1.0] = 0.0;
  }  else if (all(Y %in% c(0.0,1.0))) {
  } else {
    print(\"BayesSpikeRegression():::DoLogitPostProb=TRUE, Hey, Not All Y are 0,1, this will not be good. \"); flush.console();
    DoLogitPostProb = FALSE;
  }
  if (is.null(dfRobit) || dfRobit < 0 || dfRobit == 9.0) {
    dfRobit <- 9; AlterWeightdfRobit <- 9;
    SigmaRobit <-  sqrt(7 / 9 * pi^2 / 3);
  } else {
    print(\"BayesSpikeRegression():::Hey, dfRobit is already set though DoLogitPostProb is TRUE, this won't work out that well. \"); flush.console();
    SigmaRobit <- sqrt( 7/ 9 * pi^2 / 3);
  }
}
if ( (!is.null(dfRobit) && dfRobit >= 0) && (!exists(\"Z\") || is.null(Z) ||
  length(Z) == 0) &&
  (all(Y %in% c(0,1.0)))  ) {
  if (Verbose >= 2) {
    print(paste(\"  0004BayesSpikeRegression.r: dfRobit is \", dfRobit, \" and Y is length \", length(Y), sep=\"\"));
    flush.console();
  }
  try(Z <- NULL);
  try(Z <- rep(0, length(Y)) + Y);
  if (is.null(Z)) {
    print(paste(\" 0004BayesSpikeRegression.r: Error, tried to allocate Z but Y is bad length \", 
      length(Y), \" it is bad.\", \"\")); flush.console();
    eval(parse(text=SetGText(\"Y\", \"globalenv()\", S=1)));
    print(\" ERROR return -1 \"); flush.console();
    return(-1);
  }  
  try(X <- cbind(rep(1, NROW(X)), X));
  if (length(tauEndList) >= 1 && tauEndList[1] > 0) {
    IndexFirstRandomEffect <- IndexFirstRandomEffect +1;
    tauEndList <- tauEndList + 1;
  }
  if (is.null(tauFixed) || tauFixed[1] != -666) {
    tauFixed <- c(-1000, tauFixed);
  } else {
    tauFixed <- c(-1000);
  }
  if(Verbose >= 0) {
    print(paste(\"BayesSpikeRegression after another test tauFixed[1] = \", 
      tauFixed[1], \" and length is \", length(tauFixed), sep=\"\"));
    flush.console();
  }
  NoShrinkFixed <- c(1, NoShrinkFixed);
} else {
  Z = NULL;
}
  ")
}

##EvalTemperatureEE
##
## running eval(parse(text=EvalTemperatureEE())) within a scope, this
##  checks that all Temperature relevant data exist and is usable.
EvalTemperatureEE <- function() {
  return("
if (!exists(\"TemperatureList\")) { TemperatureList = NULL; }
  if (!is.null(TemperatureList) && length(TemperatureList) > 0) {
    if (TemperatureList[length(TemperatureList)] != 1.0) {
        TemperatureList = c(TemperatureList,1.0);
    }   
    ##NumChains = length(TemperatureList);
  }
if (!exists(\"DoEEProbSort\")) { DoEEProbSort = TRUE; }
if (is.null(DoEEProbSort) ||
 (is.logical(DoEEProbSort) && length(DoEEProbSort) == 1 &&
    DoEEProbSort == FALSE)) {
  DoEEProbSort = 0;
} else if (DoEEProbSort == TRUE) {
  DoEEProbSort = 1;
} else if (DoEEProbSort == 0) {
  DoEEProbSort = 0;
} else {
  DoEEProbSort = 1;
}
if (is.null(dfRobit)) { dfRobit = -1; }
##print(paste(\"Before TBSR5, EarlyEndtt = \", EarlyEndtt,
##   \" and EarlyEndStep = \", EarlyEndStep, Sep=\"\"));
##return(1);
NewInstance = 1;  NInstance = 1;
if (FALSE) {
eval(parse(text=GetG0Text(\".ListBayesSpikeOb\", \"globalenv()\", S=1)));
NewInstance = length(.ListBayesSpikeOb)+1;
} else {
  NewInstance = 1;
}
if (is.null(PiAPrior) || length(PiAPrior) <= 0) {
  PiAPrior <- c(2,2);
} else if (length(PiAPrior) == 1) {
  PiAPrior <- c(PiAPrior[1], PiAPrior[1]);
}
  ")
}

## Another function to call eval(parse(text=CheckXandY())) for X/Y examination of right dimensions.
CheckXandY <- function() {
return("
if (is.null(dim(X)) || length(dim(X)) != 2) {
  print(\"BayesSpikeRegression: Error X input must be a matrix!\")
  print(\"Submitted X has bad Dim!\");
  return(-1);
}
p = dim(X)[1];
if (length(Y) != p) {
  print(paste(\"BayesSpikeRegression: Error, length(Y) = \",
    length(Y), \" but p = \", p, sep=\"\")); flush.console();
  return(-1);
}
if (!is.null(tauEndList) && !is.null(tauStartValues) && 
  !is.null(names(tauStartValues)) && is.null(names(tauEndList))) {
  names(tauEndList) <- names(tauStartValues);
}
if (!is.null(tauEndList) && !is.null(tauStartValues) && 
  !is.null(names(tauEndList)) && is.null(names(tauStartValues))) {
  names(tauStartValues) <- names(tauEndList);
}
")
}

## CheckSaveDir()
##
##  Checks to see that desired SaveDir location is available.
CheckSaveDir <- function() {
  return("
  if (DoSave == TRUE &&
    exists(\"SaveDir\") && !is.null(SaveDir) && is.character(SaveDir) &&
    SaveDir != \"\" && SaveDir != \"NOSAVE\" && SaveDir != \"NoSave\") {
    dir.create(SaveDir, showWarnings=FALSE, recursive = TRUE)
  }
  ")
}

## CreateTBSR5()
##
##  Running eval(parse(text=CreateTBSR5())) in right scope generates
##    an R5 native object with the right statistical variables.
CreateTBSR5 <- function() {
  return("  
##  TBSR5 is R5 storage class
##
##    It was probably silly to save R code in an environment along with Rcpp
##     But we're stuck with some of this which is helpful in managing all of 
##     The user Interface code.  \"BayesSpikeR5.r\" has the class definition
##     and a host of functions for I/O in there.
if (Verbose >= -1) {
  print(paste(\"BayesSpikeRegression: before setting tauFixed its zero value is: \", 
  tauFixed[1], \" and length \", length(tauFixed), sep=\"\"));
  flush.console();
}
TBSR5 = BayesSpikeR5$new(X=X,Y=Y, Z=Z,
  BetaStart = BetaStart, FirstRandom = IndexFirstRandomEffect,
  tauEndList = tauEndList, sOnTau = tauStartValues, CodaTable = NULL, 
  DoRecord = DoRecord, dfTNoise = dfTNoise, Verbose = Verbose,
  DependenciesTau = DependenciesTau, DependenciesFixed = DependenciesFixed,
  TypeFixedPrior = TypeFixedPrior , MaximizeMeCauchyTotalIters = MaximizeMeCauchyTotalIters,
  PiAPrior = PiAPrior, HowSample = HowSample, NumberSpliceSampleReps = NumberSpliceSampleReps,
  SpliceMove=SpliceMove, CauchyEpsilon = CauchyEpsilon, tauFixed = tauFixed,
  MyPS = MyPS, SigmaPrior = SigmaPrior, InitKKs = InitKKs,
  NumSamples=NumSamples, MaxIters = MaxIters, MaxGibbsIters = MaxGibbsIters,
  tauPriordf = tauPriordf, tauPriorMean = tauPriorMean,
  ttStart = ttStart, NumChains = NumChains,
  EarlyEndStep = EarlyEndStep, EarlyEndtt = EarlyEndtt, DoLogitPostProb = DoLogitPostProb,
  DoSave = DoSave,  DoSaveTBSR5=DoSaveTBSR5,
  FileSave = FileSave, SaveDir = SaveDir, FileName = FileName,
  OnChainIter = 1, NewWrite = NewWrite, DoMax = DoMax,
  PutRNG = PutRNG, Run = Run, CodaTableNames = CodaTableNames,
  RStart = RStart, OnPiA = PiAGuess, PiAStart = PiAGuess, AFD=AFD, NInstance = NewInstance,
  NoShrinkFixed = NoShrinkFixed, NoShrinkFixedPrior=NoShrinkFixedPrior,
  NoShrinkRandom = NoShrinkRandom, NoShrinkRandomPrior=NoShrinkRandomPrior,
  StartRunProbVector=StartRunProbVector, TemperatureList = TemperatureList,
  dfRobit = dfRobit, DoEEProbSort = DoEEProbSort,TemperatureDecreasingRate = TemperatureDecreasingRate,
  RandomInfoList= RandomInfoList,  PreRunMaxGibbsIters = PreRunMaxGibbsIters,
  WriteYBuffer = WriteYBuffer, NewYBufferWrite = NewYBufferWrite,
  WriteWeightBuffer=FALSE, NewWeightBufferWrite= NewWeightBufferWrite,
  LengthYBuffer=LengthYBuffer, LengthWeightBuffer = LengthWeightBuffer,
  EEDoSort = EEDoSort, StartOldEELoad = StartOldEELoad,
  DoChainIters = DoChainIters, AlterWeightFlag = AlterWeightFlag,
  DoLongCI = DoLongCI, CpBuffLongCI = CpBuffLongCI, RegionWidth=RegionWidth);
  ")
}

## FalseLockObjects()
##
## This checks that ".ListBayesSpikeOb" exists in an environment if necessary to
## check that data is being preserved.
##
## Note that the code is locked in an if (FALSE) statement and shouldn't normally
##  get evaluated.
FalseLockObjects <- function() {
  return("
  
if (FALSE) { 
if (!exists(\".ListBayesSpikeOb\", envir=globalenv() ) ) {
  print(paste(\"Error: .ListBayesSpikeOb does not exist in BAYESPSIKENAMESPACE: \",
    as.character(is.environment(BAYESSPIKENAMESPACE)) )); flush.console();
}
if (!exists(\".ListBayesSpikeOb\", envir=globalenv() )) {
  print(paste(\"Error: .ListBayesSpikeOb does not exist in globalenv: \",
    as.character(is.environment(globalenv)) )); flush.console();
}
if (is.null(get(\".ListBayesSpikeOb\", envir=globalenv() ) ) ) {
  print(paste(\"Error: .ListBayesSpikeOb is still NULL in BAYESPSIKENAMESPACE: \",
    as.character(is.environment(BAYESSPIKENAMESPACE)) )); flush.console();
}
if (length(get(\".ListBayesSpikeOb\", envir=globalenv()   )) < TBSR5$NInstance ) {
  print(paste(\"Error: .ListBayesSpikeOb Does not have length compatible with TBSR5=: \",
    TBSR5$NInstance) ); flush.console();
}
}
  ")
}

## BayesSpikeStartEarlyEndSteps()
##
##  Sometimes when the GibbsSampler is started, one will want to burn in for
##  a certain number of steps before writing to disk.  Most critically that
##  is to test that new temperatures are stable.
##  
##  This function text allows a determination of whether to end before completing
##   many iterations, and whether to save the state of the system at a given EarlyEndStep.
##
BayesSpikeStartEarlyEndSteps <- function() {
  return("
  
if (MBS$dfTNoise > 0) {
   print(\"dfTNoise: MBS is greater than zero, hoo boy!\"); flush.console();
   print(paste(\"MBS$dfTNoise = \", MBS$dfTNoise, \" and TBSR5$dfTNoise = \", 
     TBSR5$dfTNoise, sep=\"\")); flush.console();
   ##return(MBS);
}

if (MBS$TBSR5$EarlyEndStep >= 0) {
  if (Verbose > -1) {
    print(\"BayesSpikeRegression: Setting up Early End Step/ Early End tt\"); flush.console();
  }
MBS$EarlyEndStep = MBS$TBSR5$EarlyEndStep;  MBS$EarlyEndtt = MBS$TBSR5$EarlyEndtt;
}

if (DoneByList == FALSE && Run == TRUE) {

if (MBS$Verbose > 1) {
print(paste(\"Starting with STT = \", paste(STT, collapse=\", \"))); flush.console();
}
  try(eval(parse(text=SetGText(\"STT\", S=1))));
##try(MBS$SecondInitiateTime <- proc.time());
  ")
}

## SetupToCreateMBS()
##
## MBS is an Rcpp-RModules RObject which is an encapsulated C++ Class that
## one can play around with in the terminal.
SetupToCreateMBS <- function() {
  return("
InitFlags = c(TBSR5$NInstance, TBSR5$Verbose[1]-1, 
  InitKKs, TBSR5$CFirstRandom); 
if (Verbose > 1) {
  print(paste(\"Our initial InitFlags = \", paste(InitFlags, collapse=\", \"), sep=\"\")); flush.console();
}
if (is.null(TBSR5)) {
  print(\"BayesSpikeRegression: Failed to generate TBSR5, return error!\");
  return(-1);
}
if (is.null(TBSR5$sBeta) || length(TBSR5$sBeta) <= 0) {
  print(\"BayesSpikeRegression: Error, sBeta is NULL, return error!\");
  return(TBSR5);
}
TBSR5$InitFlags = InitFlags;

if (FALSE) {
eval(parse(text=GetG0Text(\".ListBayesSpikeOb\", \"globalenv()\" )));

if (length(.ListBayesSpikeOb) == 1 && is.numeric(.ListBayesSpikeOb) &&
  .ListBayesSpikeOb == 0) {
  print(\"Error: We were not able to find TBSR5 attached to global .ListBayesSpikeOb\" ); flush.console();  
}
if (Verbose > 1) {
  print(\"BayesSpikeRegression: Recovering TBSR5 from the list\"); flush.console();
}

TBSR5 = .ListBayesSpikeOb[[TBSR5$NInstance]];
}
  ")
}

## CreateAndLockInMBS
## Run this to ensure that MBS is locked somewhere as declared by a call to
## BayesSpikeCL$new();
##
## MBS as a pointer is shared by the TBSR5 R5 object as well as the C++ class.
CreateAndLockInMBS <- function() {
  return("
  if (Verbose > 1) {
  print(\"BayesSpikeRegression: Getting BayesSpikeCL\");  flush.console();
}
if (is.null(modBayesSpikeCL)) {
  print(\"****************************************************************\");
  print(\"** BayesSpikeRegression Error\");
  print(\"** At start modBayesSpikeCL is NULL! \"); flush.console();
}
try(eval(parse(text=GetG0Text(\"BayesSpikeCL\", S=1))))
BayesSpikeCL <- NULL;
try(BayesSpikeCL <- modBayesSpikeCL$BayesSpikeCL);
if (is.null(BayesSpikeCL)) {
  print(\"*****************************************************************\");
  print(\"** BayesSpikeRegression Error \");
  print(\"** Tried to get BayesSpikeCL from modBayesSpikeCL  \");
  print(\"** However this produced a NULL.  Can't move forward to NEW\");
  flush.console();
  print(\"** Return TBSR5 structure with error.  \"); flush.console();
  return(TBSR5);
}


if (is.null(TBSR5$.X)) {
  print(\"Error: TBSR5 did not write correctly, TBSR5$.X is NULL!\");
  return(TBSR5);
}

if (!is.null(dfRobit) && dfRobit > 0) {
  intY = as.integer(Y);
  vY = Y * 2 - 1;
  Y = vY;
}  else if (!is.null(dfRobit) && dfRobit == 0.0) {
  intY = as.integer(Y);
  vY = Y * 2 - 1;
  Y = vY;
} else {
  vY = TBSR5$.Y;
}
if (Verbose > 0) {
  print(\"BayesSpikeRegression: Creating Rcpp C++ class\"); flush.console();
}
MBS <- NULL;
try(MBS <- new( BayesSpikeCL, TBSR5$.X, as.double(vY), as.double(TBSR5$Beta), 
  TBSR5$RandomInfoList, as.integer(InitFlags), globalenv() ));
if (is.null(MBS)) {
  print(\"*****************************************************************\");
  print(\"** BayesSpikeRegression: Error on call new BayesSpikeCL \");
  print(\"** We have that MBS is NULL on allocation. \"); flush.console();
  print(\"** Return TBSR5 instead to check types. \"); flush.console();
  return(TBSR5);
}

  
try(MBS$TBSR5 <- TBSR5);
if (Verbose > 1) {
  print(\"BayesSpikeRegression: We Received MBS, setting the NInstance again\");
}

try(TBSR5$ABayesSpikeCL <- MBS);
try(MBS$TBSR5 <- TBSR5);
if (Verbose > 1) {
  print(\"BayesSpikeRegression: We're about to try to program in MBS to TBSR5!\");
  flush.console();
}
try(MBS$TBSR5$ABayesSpikeCL <- MBS);
if (Verbose > 1) {
  print(paste(\"MBS's BayesSpkeNameSpace = \", as.character(is.environment(MBS$BayesSpikeNameSpace)),
    \"  and BAYESSPIKENAMESPACE = \", as.character(is.environment(BAYESSPIKENAMESPACE)), sep=\"\")); flush.console();
  print(paste(\" MBS's NInstance = \", MBS$NInstance, sep=\"\")); flush.console();
  print(paste(\"  and TBSR5$NInstance = \", TBSR5$NInstance,
    \" and MBS$TBSR5$NInstance = \", MBS$TBSR5$NInstance, sep=\"\")); flush.console();
}
##MBS$BayesSpikeNameSpace = BAYESSPIKENAMESPACE;  
##MBS$NInstance = TBSR5$NInstance;
  try(MBS$DoSampleBFixedii <- DoSampleBFixedii);
  ");
}

SetupDoTimePartsList <- function() {
 return(
 "
  if (DoTimePartsList == TRUE) {
  if (Verbose > 0) {
    print(\"BayesSpikeRegression: Setup DoTimePartsList\"); flush.console();
  }
  ATFlag <- 0;
  MyT <- \" ATFlag <- 0;
  ATFlag <- MBS$SetupTimePartsList();
  \"; 
  try(eval(parse(text=MyT)));
  if (ATFlag == 0) {
    print(\"ERROR: BayesSpikeRegression: Error, want to SetupTimePartsList, but do not!\");
    flush.console();
    return(MBS);
  }
  }
  ");
}

### BetaRegressionSetup
### Setup OnPiA and decide to ZeroOut
SetupOnPiAAndHowSample <- function() {
  return("
try(MBS$OnPiA <- MBS$TBSR5$OnPiA, silent=TRUE);
try(MBS$OnSigma <- SigmaSq, silent=TRUE);



#####################################################################
## Sometimes rather than random starts, if Beta is too large
##  we need to start at zero, see what factors are valuable for activation
##  and start from a subset of those activated.   Code for this is further
##  down, search \"ZeroOutBeforeSample\"  -- 10/15/2013
if (!exists(\"ZeroOutBeforeSample\")) { ZeroOutBeforeSample <- 0;  }
if (MBS$p >= 1000 || ZeroOutBeforeSample >= 1) {
  MBS$ZeroOutBeforeSample <- 1;
  if (!exists(\"StartRunProbVector\") || is.null(StartRunProbVector) ||
    StartRunProbVector[1] < 0) {
    StartRunProbVector <- 5;      
  }
}

if (DoneByList == FALSE) {
if (Verbose >= 1) {
print(\"Checking MBS X\"); flush.console();
}

if (is.null(MBS$X)) {
  print(\"Uh Oh Reorge!,  MBS$X is NULL!\"); 
  flush.console();
}
if (Verbose > 1){
  print(paste(\"  Assigning MBS to TBSR5 \", sep=\"\")); flush.console();
}
}
if (Verbose >= 1) {
print(\"BayesSpike--BayesSpikeRegression: Setting Up Features of MBS\");
}
##MBS$TBSR5$ABayesSpikeCL = MBS;
if (Verbose > 1)  {
  print(paste(\"  Assigning MBS$HowSample from TBSR5\", sep=\"\")); flush.console();
}
if (is.null(HowSample) || HowSample < 0) { HowSample = 3; }
MBS$HowSample = HowSample;
if (Verbose > 2) {
  print(paste(\"MBS$HowSample  \", MBS$HowSample, sep=\"\")); flush.console();  
}

if (Verbose > 0) {
  print(\"BayesSpikeRegression: Set up tauEndList/FirstRandom\"); flush.console();
  if (is.null(TBSR5$CtauEndList) || length(TBSR5$CtauEndList) <= 0) {
    print(\"BayesSpikeRegresion: cTauEndList of no length, won't be filled.\");
    flush.console();
  } 
}
  ")
}
SetupTauEndList <- function() {
return("
    if (Verbose > 0) {
  print(\"BayesSpikeRegression: Setting up CtauEndList = \"); flush.console();
  if(length(TBSR5$CtauEndList) <= 10) {
    print(paste(\"(\", paste(TBSR5$CtauEndList, collapse=\", \"), \")\", sep=\"\"))
  } else if (length(TBSR5$CtauEndList) >= 500) {
    print(paste(\"(\", paste(TBSR5$CtauEndList[1:5], collapse=\", \"), \",..., \", 
      TBSR5$CtauEndList[length(TBSR5$CtauEndList)], \")\", sep=\"\")); flush.console();
  } else { 
    LM <- matrix(0, 15, ceiling(length(TBSR5$CtauEndList)/15));
    LM[1:length(TBSR5$CtauEndList)] <- TBSR5$CtauEndList;
    MaxCT <- max(nchar(as.character(TBSR5$CtauEndList) ));
    APt <- function(x, MaxCT) {
      BS <- as.character(x);
      if (any(nchar(BS) < MaxCT)) {
        for (jj in 1:length(BS)) {
          if (nchar(BS[jj]) < MaxCT) {
             BS[jj] <- paste(paste(rep(\" \", MaxCT-nchar(BS[jj])), collapse=\"\"), BS[jj], sep=\"\");
          }
        }
      }
      return(BS);
    }
    for (ii in 1:NCOL(LM)) {
      if (ii == 1) {
        print(paste(\"( \", paste(APt(LM[,ii], MaxCT), collapse=\", \"), \", \", sep=\"\")); flush.console();
      } else if (ii == NCOL(LM)) {
        print(paste(\"  \", paste(APt(LM[LM[,ii] != 0,ii], MaxCT), collapse=\", \"), \") \", sep=\"\")); flush.console();
      } else {
        print(paste(\"  \", paste(APt(LM[,ii], MaxCT), collapse=\", \"), \", \", sep=\"\")); flush.console();
      }
    }
  }
  flush.console();
}
if (!is.null(TBSR5$CtauEndList)  && length(TBSR5$CtauEndList) > 0) {
  try(MBS$tauEndList <- TBSR5$CtauEndList);
}
if (!is.null(MBS$tauEndList)) {
  try(MBS$OnTau <- TBSR5$sOnTau);
}
if ((is.null(MBS$tauEndList) || length(MBS$tauEndList) <= 0) &&
  (is.null(MBS$tauEndList) || length(MBS$sOnTau) <= 0)) {
  MBS$RawiFirstRandom = -1;  
} else {
  try(MBS$RawiFirstRandom <- TBSR5$CFirstRandom);
}
");
}
SetupAlternateNoise <- function() {
return("    
if (TBSR5$AlterWeightFlag == FALSE && !is.null(dfRobit) && dfRobit >= 0) {
  ## If Robit model, Y = {0,1} and we use robit model
  MyTry <- -1;
  try(MyTry <- MBS$SetupIntY(intY, dfRobit));
  
  if (DoLogitPreProb == TRUE) {
    try(MBS$OnSigma <- 7/9 * pi^2 / 3)
    try(MBS$SigmaPrior <- c(-1,-1));
    try(MBS$DoLogitNonePostPreProb <- 2);
    try(dfRobit <- 9.0);
    try(MBS$AlterdfRobit <- -9.0);    
    try(MBS$dfRobit <- 9.0);
  } else if (DoLogitPostProb == TRUE) {
    try(MBS$OnSigma <- 7/9 * pi^2 / 3)
    try(MBS$SigmaPrior <- c(-1,-1));
    try(MBS$DoLogitNonePostPreProb <- 1);
    try(dfRobit <- 9.0);
    try(MBS$AlterdfRobit <- dfRobit);
  }
  if (MyTry < 0) {
      print(\"SetupAlternateNoise: We have an Error on SetupIntY, return to NULL\");
      flush.console();
      eval(parse(text=SetGText(\"MBS\", \"globalenv()\", S=1)));
      return(MBS);
  }  
} else if (TBSR5$AlterWeightFlag == TRUE && !is.null(dfRobit) &&
  dfRobit >= 0) {
  ## In AlterWeight model, we have Y = {0,1} and simulate a \"Z\" such that
  ##  we can do fixed regressions on that Z and map the Beta through
  ##  importance sampling back to the Robit model.
  MyTry <- -1;
  try(MyTry <- MBS$SetupIntY(intY, dfRobit));
  if (DoLogitPreProb == TRUE) {
    try(MBS$OnSigma <- 7/9 * pi^2 / 3);
    try(MBS$SigmaPrior <- c(-1,-1));
    try(MBS$DoLogitNonePostPreProb <- 2);
    try(MBS$AlterdfRobit <- -9.0);
    try(MBS$dfRobit <- 9.0);
  } else if (DoLogitPostProb == TRUE) {
    try(MBS$OnSigma <- 7/9 * pi^2 / 3);
    try(MBS$SigmaPrior <- c(-1,-1));
    try(MBS$DoLogitNonePostPreProb <- 1);
    try(MBS$AlterdfRobit <- 9.0);
    try(MBS$dfRobit <- 0.0);
  }
  if (MyTry < 0) {
      print(\"SetupAlternateNoise: We have an Error on SetupIntY, return to NULL\");
      flush.console();
      eval(parse(text=SetGText(\"MBS\", \"globalenv()\", S=1)));
      return(MBS);
  }
  if (FALSE) {
  AOCM <- 100;
  try(AOCM <- MBS$CurrentMaxGibbsIters);
  try(MBS$CurrentMaxGibbsIters <- 10);
  ATT = 0;
  if (!is.null(MBS$CodaList)) {
    if (!is.null(colnames(MBS$CodaList[[1]]))) {
      BackupColnames <-  colnames(MBS$CodaList[[1]]);
    } else {
      BackupColnames <- NULL;
    }
  } else {
    BackupColnames <- NULL;
  }
  if (!exists(\"Run\")) { Run = TRUE; }   
  if (is.null(Run)) { Run = 0; }
  if ((is.logical(Run) && Run == TRUE) || Run > 0) {
    ##try(MBS$UpdateNonFreshXtResid());
    try(ATT <- MBS$RunAlgorithm());
  }
  if (!is.numeric(ATT) || ATT <= 0 || length(ATT) != 1) {
    if (!is.numeric(ATT) || length(ATT) != 1) {
      ATT <- -1;
    }
    try(paste(print(\"RunBayesRegression: Setup Weight chain has problem ATT = \",
      ATT, sep=\"\"))); flush.console();
    return(MBS);
  }  
  try(MBS$CurrentMaxGibbsIters <- AOCM);
  try(MBS$tt <- 0);
  try(MBS$AlterdfRobit <- MBS$dfRobit);
  try(MBS$dfRobit <- -1);
  try(MBS$dfTNoise <- -1);
  }  
}
 ");
}
SetupInsertPointers <- function() {
return("
if (Verbose > 0) {
  print(\"BayesSpikeRegression: Inserting TBSR5 pointer\"); flush.console();
}

try(MBS$InitiateTime <- proc.time());
##### Inserts info into MBS to be able to destroy TBSR5
try(MBS$BayesSpikeNameSpace <- globalenv());
##  Important to insert this namespace so we can find R-environment from inside C;

if (Verbose > 0) {
  print(\"BayesSpikeRegression: Inserting DoMax\"); flush.console();
}
MBS$DoMax = DoMax;


if (Verbose > 0) {
  print(\"BayesSpikeRegression: GetRNGState?\"); flush.console();
}
if (TBSR5$PutRNG == 1) {
  MBS$RNGState = 1;
} else if (isNull(MBS$TBSR5$PutRNG) || MBS$TBSR5$PutRNG == 0) {

} else {
  MBS$RNGState = 1;
}

if (is.numeric(PrintIter)) {
  MBS$PrintIter = as.integer(round(PrintIter));
}
");
}
BSSetupEigenList <- function() {
return("

XE = NULL;
if (!is.null(TBSR5$sOnTau) && length(TBSR5$sOnTau) > 0 &&
  !is.null(TBSR5$tauEndList) && length(TBSR5$tauEndList) > 0 ) {
  if (Verbose > 0) {
    print(\"BayesSpikeRegression: Generate EigenList\"); flush.console();
  }
  if (!exists(\"ReInputList\") || is.null(ReInputList) ||
   !is.list(ReInputList) || (
     (is.double(ReInputList) || is.integer(ReInputList) || is.double(ReInputList)) &&
      ReInputList[1] == -1) || 
   is.null(ReInputList$AllEigenValues) || is.null(ReInputList$AllEigenVectors)) {
  XE = NULL;
  try(XE <- TBSR5$GenerateEigenList());   
  if (is.null(MBS$TBSR5)) {
    print(\"BayesSpikeRegression: Error after GenerateEigenList, we've got NULL TBSR5 in MBS!\");
  }
  if (is.null(XE) || (is.numeric(XE) && XE[1] != 1)) {
    print(paste(\"BayesSpikeRegression: GenerateEigenList Error, XE = \", XE, sep=\"\")); flush.console();
    return(MBS);
  }
  } else {
     try(MBS$AllEigenValues <- ReInputList$AllEigenValues);
     try(MBS$AllEigenVectors <- ReInputList$AllEigenValues);
  }
  if (Verbose > 0) {
    print(\"BayesSpikeRegression:  Success on Generate EigenList\"); flush.console();
  }
}


if (is.null(MBS)) {
  print(\"BayesspikeRegression: Error after GenerateEigenList: MBS is NULL\");
  flush.console(); return(TBSR5);
}
if (MBS$BeingDestroyed == 1) {
  print(\"BayesSpikeRegression: Error after GenerateEigenList: MBS Being Destroyed was launched.\");
  flush.console(); return(TBSR5);
}
if (Verbose > 0) {
  print(paste(\"BayesSpikeRegression: Finished EigenList with X == \", XE, sep=\"\"));
  flush.console();
}

")
}

BayesSpikeSetupEarlyEndttNow <- function() {
    return("
    
if (MBS$dfTNoise > 0) {
   print(\"dfTNoise: MBS is greater than zero, hoo boy!\"); flush.console();
   print(paste(\"MBS$dfTNoise = \", MBS$dfTNoise, \" and TBSR5$dfTNoise = \", 
     TBSR5$dfTNoise, sep=\"\")); flush.console();
   ##return(MBS);
}

if (MBS$TBSR5$EarlyEndStep >= 0) {
  if (Verbose > -1) {
    print(\"BayesSpikeRegression: Setting up Early End Step/ Early End tt\"); flush.console();
  }
MBS$EarlyEndStep = MBS$TBSR5$EarlyEndStep;  MBS$EarlyEndtt = MBS$TBSR5$EarlyEndtt;
}


    ")
}


BayesSpikeSetupMBSOnSigmaAndStuff <- function() {
    return("
    
try(eval(parse(text=SetGText(\"MBS\", \"globalenv()\", S=1))));
if (Verbose > 0) {
  print(paste(\"BayesSpikeRegression: Assign TypeFixed Prior to default 1 but TBSR5 has \", 
    MBS$TBSR5$TypeFixedPrior, sep=\"\"));
  flush.console();
}
MBS$TypeFixedPrior <- 1;
if (is.numeric(gPrior) && length(gPrior) == 1 && gPrior  > 0.0) {
   MBS$gPrior <- gPrior;   
}
if (is.null(MBS$TBSR5$TypeFixedPrior)) {
  if (is.null(TypeFixedPrior)) {
    try(MBS$TBSR5$TypeFixedPrior <- as.integer(1));
    try(MBS$TypeFixedPrior <- as.integer(1));
  } else {
    try(MBS$TBSR5$TypeFixedPrior <- as.integer(TypeFixedPrior));
    try(MBS$TypeFixedPrior <- as.integer(TypeFixedPrior[1]));
  }
} else {
  try(MBS$TypeFixedPrior <- MBS$TBSR5$TypeFixedPrior[1]);
}
try(MBS$tt <- ttStart);
try(MBS$NumSamples <- NumSamples);
if (Verbose >= -11) {
  print(paste(\"Assigning MBS$TBSR5$tauFixed which has 0 = \", 
    MBS$TBSR5$tauFixed[1], 
    \" which has length \", length(MBS$TBSR5$tauFixed), sep=\"\")); 
    flush.console();
}
MBS$tauFixed = MBS$TBSR5$tauFixed;
if (is.null(MBS$TBSR5$OnPiA) || length(MBS$TBSR5$OnPiA) <= 0) {
  print(\"Hey, MBS$TBSR5$OnPiA is NULL, what is this!\"); flush.console();
}
try(MBS$OnPiA <- MBS$TBSR5$OnPiA, silent=TRUE);
if (length(MBS$OnPiA) == 1 && (MBS$OnPiA < 1.0 && MBS$OnPiA > 0.0)) {
} else if (length(MBS$OnPiA) == 1 && (MBS$OnPiA <= 0.0 || MBS$OnPiA > 1.0)) {
  try(MBS$OnPiA <- .5, silent=TRUE);
} else if (length(MBS$OnPiA) == 2) {
  try(OnPiA <- MBS$OnPiA);
  try(OnPiA[OnPiA < 0] <- .5);
  try(OnPiA[OnPiA > 1.0] <- .5);
  try(MBS$OnPiA <- OnPiA);
  try(MBS$TBSR5$OnPiA <- as.numeric(OnPiA));
} else {
  print(\"Very weird MBS$OnPiA did not set correctly!\"); flush.console();
  return(MBS);
}
if (Verbose > 1) {
  print(\"Assigning Now OnSigma \"); flush.console();
  print(paste(\"Current OnSigma = \", MBS$OnSigma, sep=\"\"));
}
MBS$OnSigma = MBS$TBSR5$OnSigma;
if (MBS$OnSigma <= 0.0) {
  try(MBS$OnSigma <- 1.0, silent=TRUE);
}
if (Verbose > 1) {
  print(\"BayesSpikeRegression: assigned OnSigma\"); flush.console();
}
if (is.null(MBS$OnSigma)) {
  print(\"Uh Oh: after setup, MBS$OnSigma is Null!\"); flush.console();
  if (is.null(TBSR5$OnSigma)) {
    print(\"TBSR5$OnSigma is null too.\"); flush.console();
    return(TBSR5);
  }
  return(MBS)
}
    ")
}

## BayesSpikeMergeAndPriors
##
## This code is run to set up merge steps where a draw from equal likelihood
## at  higher temperature hains is proposed and taken.
BayesSpikeMergeAndPriors <- function() {
 return("   
if (exists(\"EEMergeEvery\") && !is.null(EEMergeEvery) && EEMergeEvery[1] > 0) {
  MBS$EEMergeEvery = ceiling(EEMergeEvery[1]);  
}
if (exists(\"EEProbSortWidth\") && !is.null(EEProbSortWidth) &&
  EEProbSortWidth[1] > 0) {
  try(MBS$EEProbSortWidth <- EEProbSortWidth);    
}
if (exists(\"burnin\") && !is.null(burnin) && burnin[1] > 0) {
  MBS$burnin = ceiling(burnin[1]);
} else if (!exists(\"burnin\")) {
  burnin=1;  MBS$burnin = 1;
} else if (is.null(burnin)) {
  burnin = 1; MBS$burnin = 1;
}
if (Verbose > 1) {
  print(paste(\"BayesSpikeRegression: Assign PiA/Sigma Priors\", sep=\"\"));
  flush.console();
}
MyC = NULL;
if (MBS$TBSR5$PiAPrior[1] > 0) {
MyText = \"MBS$PiAPrior <- MBS$TBSR5$PiAPrior;  MyC = 1;\";
try(eval(parse(text=MyText)), silent=TRUE);
} else {
  MyText = \"MBS$PiAPrior <- NULL;  MyC = 1;\";
  try(eval(parse(text=MyText)), silent=TRUE);
}
if (is.null(MyC)) {
  print(\"BayesSpikeRegression: Error, PiAPrior Set, MyText fails.\");
  flush.console();
  return(MBS);
} 
MyC = NULL;
if (length(MBS$TBSR5$SigmaPrior) >= 2 && MBS$TBSR5$SigmaPrior[2] > 0) {
  MyText = \"MBS$OnSigma <- MBS$TBSR5$OnSigma; MBS$SigmaPrior <- MBS$TBSR5$SigmaPrior;  MyC = 1;\";
} else {
  MyText = \"MBS$SigmaPrior <- NULL;   MyC = 1;\"
}
try(eval(parse(text=MyText)));
if (is.null(MyC)) {
  print(\"BayesSpikeRegression: Error, SigmaPrior Set, MyText fails.\");
  flush.console();
  return(MBS);
}  

if (Verbose > 1) {
  print(paste(\"BayesSpikeRegression: MaxIters, CauchyEpsilon\", sep=\"\"));
  flush.console();
}
if (is.null(MaxIters) || length(MaxIters) <= 0 || 
  !is.numeric(MaxIters) || MaxIters[1] <= 0) {
  MaxIters = 100;  
}
try(MBS$MaxIters <- MaxIters);  MBS$CauchyEpsilon = CauchyEpsilon; 
MBS$MaximizeMeCauchyTotalIters  = MBS$TBSR5$MaximizeMeCauchyTotalIters;

if (DoneByList == FALSE && StartRunProbVector >= -1) {
  MBS$SetupRunProbVector(MBS$burnin)
  if (RegionWidth >= 2) {
    MBS$RegionWidth <- RegionWidth;
  }
} else {
  MBS$SetupRunProbVector(MBS$burnin)
  if (RegionWidth >= 2) {
    MBS$RegionWidth <- RegionWidth;
  }
}
");
}

BayesSpikeSetupTaus <- function() {
  return("
 if (Verbose > 1) {
  print(paste(\"BayesSpikeRegression: tauPriordf = \", TBSR5$tauPriordf,
    \" and tauPriorMean = \", tauPriorMean, sep=\"\"));
  flush.console();
}
if (!is.null(TBSR5$tauEndList) &&  (length(TBSR5$tauEndList) > 0) &&
  (is.null(TBSR5$tauPriordf) || is.null(TBSR5$tauPriorMean)) ) {
  print(\"Please set Priordf and mean for Mixed Selection Regression\");
  flush.console(); return(TBSR5);
}
if (Verbose > 1) {
  print(paste(\"Asking about Null ness of tauEndList\", sep=\"\"));
  flush.console();
}
if (!is.null(MBS$TBSR5$tauEndList) && length(MBS$TBSR5$tauEndList) > 0) {
  if (is.null(MBS) || is.null(MBS$BeingDestroyed) || MBS$BeingDestroyed == 1) {
    print(\"BayesSpikeRegression: Error before set tauPriordf.\");
    return(TBSR5);
  }
  try(MBS$tauPriordf <- MBS$TBSR5$tauPriordf);
  try(MBS$tauPriorMean <- MBS$TBSR5$tauPriorMean);
  if (Verbose > 0) {
    print(\"BayesSpikeRegression: SetupMT()\"); flush.console();
  }
  if (is.null(MBS) || is.null(MBS$BeingDestroyed) || MBS$BeingDestroyed == 1) {
    print(\"BayesSpikeRegression: Error before SetupMT: MBS Being Destroyed was launched.\");
    return(TBSR5);
  }
  try(MBS$SetupMT());
  if (is.null(MBS) || is.null(MBS$BeingDestroyed) || MBS$BeingDestroyed == 1) {
    print(\"BayesSpikeRegression: Error after SetupMT: MBS Being Destroyed was launched.\");
    return(TBSR5);
  }
}
if (is.null(MBS$TBSR5$ABayesSpikeCL)) {
  print(\"BayesSpikeRegression: ERRROR -------------- after SetupMT, NULL M$T$M!\");
  try(MBS$TBSR5$ABayesSpikeCL <- MBS);
}
  ");
}

## BayesSpikeChecksForMemory
##
##  Sets Default values that govern allocation in memory for XtX matrix
BayesSpikeChecksForMemory <- function() {
    return("
    if (!exists(\"KillPct\")) { KillPct = -1; }
if (!exists(\"MinimumProbKeeper\")) { MinimumProbKeeper = -1; }
if (!exists(\"MaximumAllocation\")) { MaximumAllocation = -1; }
if (!exists(\"StartRemoval\")) { StartRemoval = -1; }
if (!exists(\"CheckRemovable\")) { CheckRemovable = -1; }
if (KillPct > 0 && KillPct < 1) { try(MBS$KillPct <- KillPct); }
if (MinimumProbKeeper > 0 && MinimumProbKeeper < 1) { 
  try(MBS$MinimumProbKeeper <- MinimumProbKeeper); }
if (MaximumAllocation > 50 && MaximumAllocation < MBS$p) { 
  try(MBS$MaximumAllocation <- MaximumAllocation); }
if (StartRemoval > 0 && StartRemoval < MBS$MaximumAllocation ) { 
  try(MBS$StartRemoval <- StartRemoval); } 
if (CheckRemovable > 0) { 
  try(MBS$CheckRemovable <- as.integer(CheckRemovable)); }  
  ");
}

BayesSpikeSetupProbs <- function() {
  return("
  
if (Verbose > 1) {
  print(paste(\"BayesSpikeRegression: Assign DoRecord\", sep=\"\"));
  flush.console();
}
try(MBS$DoRecord <- DoRecord);
if (!is.null(MBS$TBSR5$MaxGibbsIters) && MBS$TBSR5$MaxGibbsIters > 0) {
try(MBS$MaxGibbsIters <- MBS$TBSR5$MaxGibbsIters);
} else {
  try(MBS$MaxGibbsIters <- MaxGibbsIters);
} 


if (MBS$CFirstRandom != 0) {
  if (Verbose > 1) {
   print(\"BayesSpikeRegression: Setting up Prob Fixed. \"); flush.console();
  }
  ATT = -1;
  try(ATT <- MBS$SetupProbFixed());
  if (is.null(ATT) || !is.numeric(ATT) || ATT == -1 ) {
    print(paste(\"BayesSpikeRegression: Tried to Run Setup \", 
      \"ProbFixed but got an error, one more time!\", sep=\"\"));
    flush.console();
    MBS$Verbose = 5;
    try(ATT <- MBS$SetupProbFixed());
    if (is.null(ATT) || !is.numeric(ATT) || ATT == -1) {
      print(paste(\"BayesSpikeRegression: Ran ProbFixed twice, \",
        \" failed!\", ATT, sep=\"\")); flush.console();
    }
    return(MBS);
  }
  if (is.null(MBS) || is.null(MBS$BeingDestroyed) || MBS$BeingDestroyed == 1) {
    print(\"BayesSpikeRegression: Error after SetupProbFixed: MBS Being Destroyed was launched.\");
    return(TBSR5);
  }
}
if (!is.null(MBS$TBSR5$tauEndList) && length(MBS$TBSR5$tauEndList) > 0) {
  if (Verbose > 1) {
    print(\"BayesSpikeRegression: SetupProbTau.\"); flush.console();
  }
  if (is.null(MBS) || is.null(MBS$BeingDestroyed) || MBS$BeingDestroyed == 1) {
    print(\"BayesSpikeRegression: Error before SetupProbTotal: MBS Being Destroyed was launched.\");
    return(TBSR5);
  }
  try(MBS$SetupProbTau()); 
}
if (is.null(MBS$TBSR5$ABayesSpikeCL)) {
  print(\"BayesSpikeRegression: ERRROR -------------- after SetupProbTau, NULL M$T$M!\");
  try(MBS$TBSR5$ABayesSpikeCL <- MBS);
}
try(MBS$MaxGibbsIters <- MBS$TBSR5$MaxGibbsIters);
if (is.null(MBS) || is.null(MBS$BeingDestroyed) ||
  (is.numeric(MBS$BeingDestroyed) && MBS$BeingDestroyed == 1)) {
  print(\"BayesSpikeRegression: Error before set MBS$dfTNoise, MBS being destroyed.\");
  flust.console(); return(TBSR5);  
}
if (is.null(MBS$TBSR5)) {
  print(\"BayesSpikeRegression: Error before set MBS$dfTNoise, TBSR5 is inaccessible\");
  flush.console();
  return(TBSR5);
}
if (!is.numeric(MBS$TBSR5$dfTNoise)) {
  print(\"BayeSpikeRegression: Error before set MBS$dfTNoise, TBSR5$dfTNoise is not numeric!\");
  flush.console();
  return(TBSR5);
}
  ")
}
BayesSpikeUpToEarlyEndStep <- function() { 
return("
BZ = 0;
ATryText = \"
if (Verbose > 1) {
  print(paste(\\\"BayesSpikeRegression: MBS$TBSR5$dfTNoise[1] = \\\",
    MBS$TBSR5$dfTNoise[1], sep=\\\"\\\")); flush.console();
}
if (!is.null(MBS) && !is.null(MBS$TBSR5) && 
  is.numeric(MBS$TBSR5$dfTNoise) && length(MBS$TBSR5$dfTNoise) >= 1 &&
  MBS$TBSR5$dfTNoise[1] > 0) {
  if (Verbose > 1) {
    print(\\\"BayesSpikeRegression: Setup dfTNoise\\\");flush.console();
  }
  try(MBS$dfTNoise <- MBS$TBSR5$dfTNoise);
}
BZ = 1;
\";
try(eval(parse(text=ATryText)));
if (BZ != 1) {
  print(\"BayesSpikeRegression: Issue after attempt to setup dfTNoise for MBS\");
  flush.console();
  return(TBSR5);
}

if (Verbose > 0) {
  print(\"BayesSpikeRegression: Setup Dependencies\"); flush.console();
}
if (is.null(MBS) || is.null(MBS$BeingDestroyed) || (is.numeric(MBS$BeingDestroyed) &&
  MBS$BeingDestroyed == 1)) {
  print(\"BayesSpikeRegression: Error before check Early End Steps!\");
  flush.console(); return(TBSR5);  
}
if (!is.numeric(MBS$TBSR5$EarlyEndStep)) {
  print(\"BayesSpikeRegression: Check, EarlyEndStep is not numeric in TBSR5!\");
  flush.console();
}
if (!is.null(MBS$TBSR5$EarlyEndStep) && MBS$TBSR5$EarlyEndStep >= 0)  {
  if (Verbose > 1) {
    print(paste(\"BayesSpikeRegression: Setup EarlyEndStep = \", 
      MBS$TBSR5$EarlyEndStep, \" and EarlyEndtt=\",
      MBS$TBSR5$EarlyEndtt, sep=\"\")); flush.console();
  }
  try(MBS$EarlyEndStep <- MBS$TBSR5$EarlyEndStep);
  try(MBS$EarlyEndtt <- MBS$TBSR5$EarlyEndtt);
}
try(eval(parse(text=SetGText(\"MBS\", S=1))));
if (is.null(MBS$TBSR5$ABayesSpikeCL)) {
  print(paste(\"BayesSpikeRegression: ERRROR -------------- \",
    \"before Dependencies, NULL M$T$M!\", sep=\"\"));
  try(MBS$TBSR5$ABayesSpikeCL <- MBS);
}
");
}

BayesSpikeSetupDependencies <- function() {
  return("
  
if (Verbose > 1) {
  print(\"BayesSpikeRegression: About to insert Dependencies.\");
  flush.console();
}
BZ = 0;
ATryText = \"
if (!is.null(MBS$TBSR5$DependenciesTau)  && length(MBS$TBSR5$DependenciesTau) > 0) {
  print(\\\"BayesSpikeRegression: Setting up MBS$DependenciesTau.\\\");
  flush.console();
  if (length(MBS$TBSR5$DependenciesTau)  != length(MBS$tauEndList)) {
    print(paste(\\\"BayesSpikeRegression: issue MBS$DependenciesTau length = \\\",
      length(MBS$TBSR5$DependenciesTau), \\\" but length tauEndList = \\\",
      length(MBS$tauEndList), sep=\\\"\\\")); flush.console();
    return(TBSR5);
  }
  try(MBS$DependenciesTau <- MBS$TBSR5$DependenciesTau);
  try(MBS$DependenciesFixed <- MBS$TBSR5$DependenciesFixed);
}

BZ = 1;
\"
try(eval(parse(text=ATryText)));
if (BZ != 1) {
  print(\"BayesSpikeRegression: Something was awry in setup DependenciesTau.\");
  flush.console();
}
  ")
}
BayesSpikeSetupItersAndSave <- function() {
  return("
  
if (!is.null(AFD)) { 
  if (Verbose > 1) {
    print(\"BayesSpikeRegression: Inserting AFD to MBS.\");
    flush.console();
  }
  try(MBS$AFD <- AFD); 
}

if (!exists(\"OnChainIter\")) { OnChainIter = 1; }
if (OnChainIter <= 0 || OnChainIter > MBS$TBSR5$NumChains) {
  OnChainIter = 1;
}
try(MBS$TBSR5$OnChainIter <- as.integer(OnChainIter));
if (sum(abs(DoRecord) > 0)) { 
 if (!is.null(MBS$CodaList) && length(MBS$CodaList) >= OnChainIter) {
  try(MBS$OnCodaTable <- OnChainIter);
 }
}
STT <- 1:MBS$TBSR5$NumChains;
if (exists(\"DoChainIters\") && !is.null(DoChainIters)  && length(DoChainIters) >= 1) {
  if (!is.null(STT[STT %in% DoChainIters])) {
    STT <- STT[STT %in% DoChainIters];
  }
}
eval(parse(text=SetGText(\"STT\", \"globalenv()\", S=1)));
if (!is.null(SaveDir)  && !(SaveDir %in% c(\"NoSave\", \"NOSAVE\")) ) {
try(dir.create(SaveDir, showWarnings=FALSE));
} else {
  try(MBS$DoSave <- 0);  try(MBS$NoSave <- -1);
}

if (!is.null(MBS$TBSR5$NoShrinkFixed)) {
  try(MBS$NoShrinkFixed <- unlist(MBS$TBSR5$CNoShrinkFixed));
  try(MBS$NoShrinkFixedPrior <- unlist(MBS$TBSR5$NoShrinkFixedPrior));
}
if (!is.null(MBS$TBSR5$NoShrinkRandom)) {
  try(MBS$NoShrinkRandom <- unlist(MBS$TBSR5$CNoShrinkRandom));
  try(MBS$NoShrinkRandomPrior <- unlist(MBS$TBSR5$NoShrinkRandomPrior));
}

  ")
}
BayesSpikeSetupSaveFileToBeDestroyed <- function() {
  return("
  
if (Verbose > 0) {
  print(paste(\"BayesSpikeRegression: SetupSaveFile: \", TBSR5$FileName, sep=\"\")); flush.console();
}
if (is.null(MBS)) {
  print(\"BayesSpikeRegression, weird, MBS is Null before SetupFile!\");
  return(TBSR5);
}
if (is.null(MBS$BeingDestroyed) || !is.numeric(MBS$BeingDestroyed) ||
  MBS$BeingDestroyed == 1) {
  print(\"BayesSpikeRegression, weird, MBS is being destroyed before Setup SaveFile.\");
  flush.console();  
  return(TBSR5);
}
if (is.null(MBS$TBSR5)) {
  print(\"BayesSpikeRegression: Weird, MBS starts with TBSR5 being NULL!\");
  return(MBS);
}
if (is.null(MBS$TBSR5$ABayesSpikeCL)) {
  print(\"BayesSpikeRegression: Do, before Save File, MBS$TBSR5$ABSCL is null!\");
  try(MBS$TBSR5$ABayesSpikeCL <- MBS);
}
if (is.null(MBS$TBSR5$ABayesSpikeCL$BeingDestroyed) ||
  !is.numeric(MBS$TBSR5$ABayesSpikeCL$BeingDestroyed)) {
  print(\"BayesSpikeRegression: Don't know if TBSR5$MBS is being destroyed\");
  flush.console(); return(MBS);  
}

if (MBS$TBSR5$ABayesSpikeCL$BeingDestroyed[1] == 1) {
  print(paste(\"BayesSpikeRegression: Before setup file, we have \",
    \"MBS$TBSR5$ABayesspikeCL has problems, find a fix, return MBS.\", sep=\"\"));
  flush.console();
  return(MBS);  
}
  ")
}
BayesSpikeCodaListSetter <- function() {
return("
  TempRunFlag = 1;
  if (DoneByList == FALSE &&  DoLogitPreProb == TRUE || 
   (AlterWeightFlag == TRUE && is.numeric(dfRobit) && length(dfRobit) == 1 &&
   dfRobit >= 0.0) || DoLogitPostProb == TRUE) {
  if (Verbose >= 3) {
    print(\"RunRegression: Trying to setup realistic weights for Y .\");
    flush.console();
  }
  MyText = \"
  AOCM <- 100;
  try(AOCM <- MBS$CurrentMaxGibbsIters);
  try(MBS$CurrentMaxGibbsIters <- MBS$TBSR5$PreRunMaxGibbsIters);
  ATT = 0;
  if (MBS$Verbose >= -2) {
    print(paste(\\\"BayesSpikeRegression: We're going to do a \\\",
      \\\"temporary Run Campaign!\n\\\", sep=\\\"\\\")); flush.console();
  }
  if (MBS$TBSR5$DoSave == TRUE) {
    try(MBS$NoSave <- 1);
  } else {
    try(MBS$NoSave <- 1);
  }
  if (DoLogitPreProb == FALSE && DoLogitPostProb == TRUE) {
    if (!exists(\\\"StartRecordingiiWeight\\\") || StartRecordingiiWeight < 0) {
      StartRecordingiiWeight <- 3;   
    }
  }
  if (DoLogitPreProb == TRUE) {
  } else if (exists(\\\"StartRecordingiiWeight\\\") && StartRecordingiiWeight >= 0) {
    if (MBS$Verbose >= -2) {
      print(\\\"BayesSpikeRegression, we will setup SetupSumiiWeight\\\"); flush.console();
    }
    MBS$SetupSumiiWeight(StartRecordingiiWeight);
  } else if (MBS$AlterWeightdfRobit > 0.0) {
    print(paste(\\\"BayesSpikeRegression: Hey, why is MBS$AlterWeightdfRobit = \\\",
      MBS$AlterWeightdfRobit, \\\"  and the StartRecordingiiWeight not set up? \\\", sep=\\\"\\\"));
    flush.console();
    if (StartRecordingiiWeight > 0) {
      MBS$SetupSumiiWeight(StartRecordingiiWeight);
    }
  }
  if (DoLogitPreProb == TRUE) {
    try(MBS$dfRobit <- 9.0);
    try(MBS$AlterWeightdfRobit <- -9.0);
    try(MBS$AlterWeightdfTNoise <- -9.0);
    MBS$OnSigma <- 7/9 * pi^2 /3;
    MBS$SigmaPrior <- c(-1.0,-1.0); DoLogitPostProb = FALSE;
    try(MBS$DoLogitNonePostPreProb <- as.integer(2)); DoLogitPostProb = FALSE;
  } else if (DoLogitPostProb == TRUE) {
    ## Note, logistic regression is close to doing robit regression
    ##  with dfRobit = 9.0 and A spread parameter of 7/9 *pi^2 /3
    if (MBS$Verbose >= -1) {
      print(\\\"BayesSpikeRegression(): DoLogitPostProb is true and so we set AlterWeightdfRobit = 9.0\\\");
      flush.console();        
    }
    MBS$AlterWeightdfRobit <- 9.0; DoLogitPreProb = FALSE;
    MBS$dfRobit <- 9.0;
    try(MBS$OnSigma <- 7/9 * pi^2 / 3);
    try(MBS$SigmaPrior <- c(-1.0,-1.0));
    try(MBS$DoLogitNonePostPreProb <- as.integer(1));
  } else if (DoLogitPreProb == TRUE) {
    MBS$DoLogitNonePostPreProb <- 2.0; DoLogitPostProb = FALSE;
  } else if (DoLogitPostProb == TRUE) {
    MBS$DoLogitNonePostPreProb <- 1.0; DoLogitPreProb = FALSE;
  } else if (MBS$AlterWeightdfTNoise > 0.0 && MBS$AlterWeightdfRobit < 0.0) {
    MBS$dfRobit <- -3; 
    MBS$dfTNoise <- MBS$AlterWeightdfTNoise;
  } else if (MBS$AlterWeightdfRobit > 0.0) {
    MBS$dfTNoise <- -3;
    MBS$dfRobit <- MBS$AlterWeightdfRobit;
  }
  if (MBS$Verbose >= -2) {
    print(paste(\\\"  --- Temp Run set MBS$NoSave to: \\\", MBS$NoSave,
      \\\" - dfRobit = \\\", MBS$dfRobit, sep=\\\"\\\")); flush.console();
  }
  if (!is.null(MBS$CodaList)) {
    if (!is.null(colnames(MBS$CodaList[[1]]))) {
      BackupColnames <-  colnames(MBS$CodaList[[1]]);
    } else {
      BackupColnames <- NULL;
    }
  } else {
    BackupColnames <- NULL;
  }
  if (!exists(\\\"Run\\\")) { Run = TRUE; }
  if (DoLogitPreProb == TRUE) {
  } else if (Run == TRUE || DoLogitPostProb == TRUE) {
    ##try(MBS$UpdateNonFreshXtResid());
    TempRunFlag <- NULL;
    if (Verbose >= -2) {
      print(\\\"*******************************************************************\\\"); flush.console();
      print(paste(\\\"*** MBS$RunAlgorithm: We will execute a TempRun, MBS$tt = 0, run \\\",
        MBS$CurrentMaxGibbsIters, sep=\\\"\\\")); flush.console();
      print(\\\"***  --- Go \\\"); flush.console();
    }
    try(MBS$AlterInitRun <- 1);
    try(MBS$tt <- 0);
    try(TempRunFlag <- MBS$RunAlgorithm());
    if (is.null(TempRunFlag) || (is.numeric(TempRunFlag) && TempRunFlag < 0)) {
      print(\\\" --- ERROR After Temporary Run, we have an MBS$ Fail, return MBS \\\"); flush.console();
      return(MBS);
    }
    if (Verbose >= -2) {
      print(\\\"******************************************************************\\\"); flush.console();
    }
  } else {
    TempRunFlag = 1;
  }
  if (MBS$Verbose >= -2) {
    print(paste(\\\"  --- After Temp Run set MBS$NoSave is: \\\", MBS$NoSave,
      \\\" -- dfRobit = \\\", MBS$dfRobit, \\\" and tt = \\\", MBS$tt, sep=\\\"\\\")); flush.console();
  }
  ##print(\\\"Temporary Quit now\\\"); flush.console();
  ##return(MBS);
  if (MBS$Verbose >= -2 && (DoLogitPreProb != TRUE && (Run==TRUE || DoLogitPostProb == TRUE))) {
     MyWeights <- MBS$iiWeight;
     if (is.null(MyWeights) || length(MyWeights) <= 0) {
        print(\\\"   ---   Hey, iiWeight is NULL, what was the point of NoSave Temp Run? \\\", sep=\\\"\\\");
        flush.console();
     } else {
        print(paste(\\\"  --- iiWeight after run is mean=\\\", mean(MyWeights),
         \\\" and sd \\\", sd(MyWeights), \\\" on tt = \\\", MBS$tt, sep=\\\"\\\")); flush.console();
     }
  }
  if (!is.numeric(ATT) || ATT <= 0 || length(ATT) != 1) {
    if (!is.numeric(ATT) || length(ATT) != 1) {
      ATT <- -1;
    }
    try(print(paste(\\\"RunBayesRegression: Setup Weight chain has problem ATT = \\\",
      ATT, sep=\\\"\\\"))); flush.console();
    return(MBS);
  } 
  if (MBS$TBSR5$DoSave == TRUE) { 
    try(MBS$NoSave <- 0);
  } else {
    try(MBS$NoSave <- -1);
  }
  
  if (DoLogitPreProb == TRUE) {
    if (MBS$TBSR5$DoSave == TRUE) {
      try(MBS$NoSave <- 0);
    }
    MBS$AlterdfRobit <- -9.0;
    MBS$dfRobit <- 9.0;  
  } else if (MBS$dfRobit > 0.0 || MBS$dfTNoise > 0.0 || DoLogitPostProb == TRUE) {
    if (MBS$Verbose >= -2) {
      print(paste(\\\"BayesSpikeRegression: Now calculating AverageIIWeight with CountiiWeight = \\\",
        MBS$CountiiWeight, sep=\\\"\\\")); flush.console();
    }
    try(AverageIIWeight <- MBS$AverageiiWeight);
    if (MBS$Verbose >= -2) {
      print(paste(\\\" --- We got AverageIIWeight = (\\\",
        paste(AverageIIWeight[1:5], collapse=\\\", \\\"), \\\"...)\\\", sep=\\\"\\\")); 
        flush.console();
      print(\\\" ------ \\\"); flush.console();
    }
    if (MBS$dfTNoise > 0.0 && MBS$dfRobit < 0.0) {

      MBS$dfTNoise <- -3;
    } else {
      ##try(MBS$AlterdfRobit <- MBS$dfRobit);
      ##try(MBS$AlterdfTNoise <- MBS$dfTNoise);
      try(MBS$dfRobit <- 0.0);
    }
    
    try(MBS$SetupSumiiWeight(-1));
   
    if (!is.null(AverageIIWeight) && length(AverageIIWeight) == MBS$n) {
      try(MBS$iiWeight <- AverageIIWeight);
      if (MBS$Verbose >= -1) {
        AMM <- min(MBS$n, 5);
        print(paste(\\\"RunBayesRegression: AlteredWeightdfRobit, now we set first av weights: (\\\",
          paste(round(AverageIIWeight,4), collapse=\\\",\\\", sep=\\\" \\\"), \\\")\\\",
          sep=\\\"\\\")); flush.console();
      }
      try(MBS$TBSR5$AverageIIWeight <- AverageIIWeight);
      if (MBS$Verbose >= -1) {
        AMii <- MBS$iiWeight;
        if (!is.null(AMii)) {
          AvWeight = mean(AMii);
          print(paste(\\\"RunBayesRegression: Currently a setting for iiWeight is \\\", AvWeight,
            sep=\\\"\\\")); flush.console();
        } else {
          print(paste(\\\"RunBayesRegression: Sorry we could not set weights correctly AMii does not exist. \\\", sep=\\\"\\\"));
          flush.console();
        }
      }
    } else {
      print(\\\"RunBayesRegression: AlteredWeightdfRobit, Sorry, AverageIIWeight does not exist.\\\"); flush.console();
    }
  }
  if (MBS$Verbose >= -2) {
    print(paste(\\\"  --- After Run set MBS$NoSave to: \\\", MBS$NoSave, \\\" tt was \\\", MBS$tt,
      sep=\\\"\\\")); flush.console();
  }
  try(MBS$CurrentMaxGibbsIters <- MBS$MaxGibbsIters);
  try(MBS$tt <- 0);
  if (DoLogitPreProb == TRUE) {
  } else if (MBS$AlterdfRobit >= 0.0 || MBS$AlterdfTNoise >= 0.0)  {
    MBS$dfTNoise <- -3; MBS$dfRobit = 0.0;
  } else {
    try(MBS$dfTNoise <- -1);  MBS$dfRobit = -1.0;
  }
  if (DoLogitPreProb == TRUE) {
    print(\\\"  -- RunBayesRegression: After SetupCoda, we are using DoLogitPreProb==TRUE\\\"); flush.console();
    try(MBS$DoLogitNonePostPreProb <- 2);
    try(MBS$dfRobit <- 9.0);
    try(MBS$NoSave <- TRUE);
    Run=TRUE;
    try(MBS$AlterdfRobit <- -9.0); 
  }
  \";
  try(eval(parse(text=MyText)));
  if (DoLogitPreProb != TRUE && Verbose >= -1) {
    ATii <- MBS$iiWeight;
    if (!is.null(ATii) && length(ATii) >= 1) {
      AviiWeight <- mean(ATii)
      print(paste(\"BayesSpikeRegression: Done running a few \",
      \"half iterations of BayesSpike, tt = \", MBS$tt, \", true average of iiWeight is \", round(AviiWeight,3), 
       \"  and dfRobit = \", MBS$dfRobit, sep=\"\"));
    } else {
      AviiWeight <- 0;
    }
    print(paste(\"BayesSpikeRegression: Done running a few \",
      \"half iterations of BayesSpike, tt = \", MBS$tt, \", average of iiWeight is \", round(AviiWeight,3), 
      \" and dfRobit = \", MBS$dfRobit, sep=\"\"));
    flush.console();
  }
  try(MBS$CurrentMaxGibbsIters <- MBS$MaxGibbsIters);
  try(MBS$AlterInitRun <- 0.0);  ## Will Run algorithm on Non Alternate dfTNoise if AlterdfRobit/AlterdfTNosie is set
}
");
}
BayesSpikeSetupSaveToWrite <- function() {
  return("
    if (MBS$TBSR5$DoSave == FALSE) {
    try(MBS$NoSave <- -1);
    try(MBS$DoSave <- 0);
  }
TTii = -1; ItTemps = -1;
ABTest = TRUE;
if (DoneByList != FALSE) {ABTest = FALSE; }
if (MBS$TBSR5$DoSave != TRUE)  {ABTest = FALSE; }
if (is.character(MBS$TBSR5$SaveDir) &&
   MBS$TBSR5$SaveDir %in% c(\"NoSave\", \"NOSAVE\"))  {ABTest = FALSE; }
if (!(is.character(MBS$TBSR5$FileName) &&
  nchar(MBS$TBSR5$FileName) > 0))  {ABTest = FALSE; }
if (ABTest == TRUE) {
  AD = -1;
  if (MBS$Verbose >= 2) {
    if (!exists(\"STT\")) { STT = 1:MBS$TBSR5$NumChains; 
      eval(parse(text=SetGText(\"STT\", \"globalenv()\", S=1)));
    }
    print(paste(\"BayesSpikeRegression: Starting SetupSaveFile: TTii=\",
      TTii, \"/\", length(ItTemps), \", Chain = \", OnChainIter, \"/\", length(STT), sep=\"\")); flush.console();
  }
  try(AD <- MBS$TBSR5$SetupSaveFile(SaveDir=MBS$TBSR5$SaveDir, FileName = TBSR5$FileName,
    chainiter = OnChainIter));
  if (is.null(AD) || !is.numeric(AD) || AD < 0) {
    print(\"MBS$TBSR5: Error SetupSaveFile Ended badly!\")
  }
  if (MBS$Verbose >= 2) {
    if (!exists(\"STT\")) { STT = 1:MBS$TBSR5$NumChains; 
      eval(parse(text=SetGText(\"STT\", \"globalenv()\", S=1)));
    }
    print(paste(\"BayesSpikeRegression: Successful Done SetupSaveFile: Tii=\",
      TTii, \"/\", length(ItTemps), \", Chain = \", OnChainIter, \"/\", length(STT), 
      \", AOnFile = \", MBS$TBSR5$OnFile, \", and AOldFile = \", 
      MBS$TBSR5$OldFile, \".\", sep=\"\")); flush.console();
  }
}
if (MBS$Verbose >= 2) {
  print(\"BayesSpikeRegression: Starting NewWrite.\"); flush.console();
}
try(MBS$NewWrite <- MBS$TBSR5$NewWrite);
  ")
}

BayesSpikeSetupPriorProbFixedTau <- function() {
  return("
  
MyEval = NULL;
if (!exists(\"PriorProbTau\")) { PriorProbTau = NULL; }
if (!exists(\"PriorProbFixed\")) { PriorProbFixed = NULL; }
if (MBS$Verbose >= 2) {
  print(\"BayesSpikeRegression: Starting to setup PriorProbTau, PriorProbFixed.\"); flush.console();
}
if (MBS$Verbose >= 1) {
  print(\"BayesSpikeRegression: Wrangling with Priors, PriorProbTau, PriorProbFied\"); flush.console();
}
if (is.null(tauEndList) || length(tauEndList) <= 0) {
  LenFixed <- MBS$p;
} else {
  LenFixed <- MBS$CFirstRandom;
} 
if (!is.null(PriorProbFixed) && length(PriorProbFixed) >= 1) {
  if (LenFixed <= 0) {
  } else {
    try(OldPriorProbFixed <- unlist(PriorProbFixed) + 0);
    if (length(PriorProbFixed) == LenFixed) {  
    } else if (length(PriorProbFixed) == 2 * LenFixed) {
    } else if (length(PriorProbFixed) < 1) {
      PriorProbFixed <- NULL;
    } else if (!is.null(names(PriorProbFixed))) {
      mMNames <- match(names(PriorProbFixed), colnames(MBS$X)[1:LenFixed]);
      fMNames <- match(names(PriorProbFixed), paste(\"Fixed:\", colnames(MBS$X)[1:LenFixed], sep=\"\"));
      nMNames <- match(names(PriorProbFixed), paste(1:LenFixed, sep=\"\"))
      aMNames <- match(names(PriorProbFixed), paste(\"Fixed:\", 1:LenFixed, sep=\"\"))
      aMNames2 <- match(names(PriorProbFixed), paste(\"FixedEffect:\", 1:LenFixed, sep=\"\")) 
      aMNames3 <- match(names(PriorProbFixed), paste(\"FixedEffect:\", colnames(MBS$X)[1:LenFixed], sep=\"\"))            
      mMNames[is.na(mMNames) & !is.na(fMNames)] <- fMNames[is.na(mMNames) & !is.na(fMNames)];
      mMNames[is.na(mMNames) & !is.na(nMNames)] <- nMNames[is.na(mMNames) & !is.na(nMNames)];
      mMNames[is.na(mMNames) & !is.na(aMNames)] <- aMNames[is.na(mMNames) & !is.na(aMNames)];
      mMNames[is.na(mMNames) & !is.na(aMNames2)] <- fMNames[is.na(mMNames) & !is.na(aMNames2)];             
      mMNames[is.na(mMNames) & !is.na(aMNames3)] <- fMNames[is.na(mMNames) & !is.na(aMNames3)]; 
      nPriorProbFixed <- rep(-1, LenFixed);
      try(nPriorProbFixed[mMNames[!is.na(mMNames)]] <- PriorProbFixed[!is.na(mMNames)]);  
      try(PriorProbFixed <- nPriorProbFixed);      
    } else {
      print(\"*****************************************************************\");
      print(\"****  Error Error: BayesSpike Regression Input Error On PriorProbFixed! \"); flush.console();
      if (!exists(\"LenFixed\")) {
        if (length(tauEndList) <= 0 || is.null(tauEndList)) {
          LenFixed <- MBS$p;
        } else {
          LenFixed <- MBS$RawiFirstRandom;
        }
      }
      print(paste(\"**** PriorProbFixed Fitting error: we had LenFixed = \", LenFixed, sep=\"\")); flush.console();
      print(paste(\"**** But there are no names to PriorProbFixed and it is: \", sep=\"\")); flush.console();
      print(PriorProbFixed);
      print(paste(\"**** Old PriorProb Fixed was: \")); flush.console();
      try(print(OldPriorProbFixed)); flush.console();
      print(\"*** We are going to set PriorProbFixed elements to GlobalEnv so you can check\");
      eval(parse(text=SetGText(\"MBS\", \"globalenv()\", S=1)));
      eval(parse(text=SetGText(\"OldPriorProbFixed\", \"globalenv()\", S=1)));
      eval(parse(text=SetGText(\"PriorProbFixed\", \"globalenv()\", S=1)));
      eval(parse(text=SetGText(\"PriorProbTau\", \"globalenv()\", S=1)));
      return(MBS);
    }
    try(PriorProbFixed <- unlist(PriorProbFixed));
    if (!is.null(PriorProbFixed) && length(PriorProbFixed) >= 1) {
      PriorProbFixed[is.na(PriorProbFixed)] <- -1;
      PriorProbFixed[is.nan(PriorProbFixed)] <- -1;
    }
    if (is.null(PriorProbFixed) || all(PriorProbFixed < 0 )) {
      try(PriorProbFixed <- NULL);
      try(MBS$PriorProbFixed <- NULL);
    } else if (!is.null(PriorProbFixed)) { 
      try(MBS$PriorProbFixed <- as.real(PriorProbFixed), silent=TRUE);
    }
  }
} else {
  PriorProbFixed <- NULL;
  if (!exists(\"LenFixed\")) {
    if (length(tauEndList) <= 0 || is.null(tauEndList)) {
      LenFixed <- MBS$p;
    } else {
      LenFixed <- MBS$RawiFirstRandom;
    }
  }
  if (LenFixed >= 1) {  
    try(MBS$PriorProbFixed <- NULL, silent=TRUE);
  }
}
if (!is.null(PriorProbTau)) {
  MyEval = 0;
  OldPriorProbTau <- PriorProbTau;
  MyEval1 = \"
  if (length(tauEndList) >= 1 && (!is.numeric(PriorProbTau) || length(PriorProbTau) < 1)) {
    PriorProbTau <- rep(-1, length(tauEndList));
  } else if (length(tauEndList) >= 1 && is.numeric(PriorProbTau) &&
    length(PriorProbTau) == length(tauEndList)) {
    names(PriorProbTau) <- names(tauEndList);
  } else if (length(tauEndList) >= 1 && length(PriorProbTau) >= 0) {
    nPriorProbTau = rep(-1, length(tauEndList));
    mNPPT <- match(names(PriorProbTau), names(tauEndList));
    fNPPT <- match(names(PriorProbTau), paste(1:length(tauEndList)));
    aNPPT <- match(names(PriorProbTau), paste(\\\"tau:\\\", 1:length(tauEndList), sep=\\\"\\\"));
    aNPPT2 <- match(names(PriorProbTau), paste(\\\"RandomEffect:\\\", 1:length(tauEndList), sep=\\\"\\\"));
    aNPPT3 <- match(names(PriorProbTau), paste(\\\"tau:\\\", names(tauEndList), sep=\\\"\\\"));
    mNPPT[is.na(mNPPT) & !is.na(fNPPT)] <- mNPPT[is.na(mNPPT) & !is.na(fNPPT)];
    mNPPT[is.na(mNPPT) & !is.na(aNPPT)] <- mNPPT[is.na(mNPPT) & !is.na(aNPPT)];
    mNPPT[is.na(mNPPT) & !is.na(aNPPT2)] <- mNPPT[is.na(mNPPT) & !is.na(aNPPT2)];
    mNPPT[is.na(mNPPT) & !is.na(aNPPT3)] <- mNPPT[is.na(mNPPT) & !is.na(aNPPT3)];
    nPriorProbTau[mNPPT[!is.na(mNPPT)]] <- PriorProbTau[!is.na(mNPPT)];
    PriorProbTau <- unlist(nPriorProbTau);
  } else if (length(PriorProbTau) == 2 * length(tauEndList)) {
    
  } else if (length(PriorProbTau) == 1) {
    PriorProbTau <- rep(PriorProbTau[1], length(tauEndList))
  } else if (length(PriorProbTau) == 2) {
    PriorProbTau <- matrix(rep(PriorProbTau, length(tauEndList)), 2, length(tauEndList));
  } else {
    print(\\\"**************************************************************\\\"); flush.console();
    print(\\\"***** Error:  Hey BayesSpikeRegression Trying to Set Prior ProbTau!\\\");
    print(paste(\\\"***** PriorProbTau was given no names and bad length = \\\", 
      length(PriorProbTau), sep=\\\"\\\"));
    print(paste(\\\"*****  Length we really need is \\\", length(tauEndList), sep=\\\"\\\"));
    print(\\\"*****  SaveMBS Error to globalenv and quit! \\\");
    eval(parse(text=SetGText(\\\"tauEndList\\\", \\\"globalenv()\\\", S=1)));
    eval(parse(text=SetGText(\\\"MBS\\\", \\\"globalenv()\\\", S=1)));
    eval(parse(text=SetGText(\\\"PriorProbTau\\\", \\\"globalenv()\\\", S=1))); 
    eval(parse(text=SetGText(\\\"OldPriorProbTau\\\", \\\"globalenv()\\\", S=1))); 
    return(MBS);   
  }
  if (is.list(PriorProbTau)) { PriorProbTau <- unlist(PriorProbTau); }
  if (!is.null(PriorProbTau) && length(PriorProbTau) >= 1) {
    PriorProbTau[is.na(PriorProbTau)] <- -1;
    PriorProbTau[is.nan(PriorProbTau)] <- -1;
  }
  if (is.null(PriorProbTau) || all(PriorProbTau < 0)) {
    PriorProbTau <- NULL;
    try(MBS$PriorProbTau <- NULL);
  } else {
    try(MBS$PriorProbTau <- as.real(PriorProbTau));
  }
  MyEval <- 1;
  \";
  try(eval(parse(text=MyEval1)));
  if (MyEval == 0) {
     print(\"******************************************************************\");
     print(\"**** Error: Hey BayesSpikeRegression Trying to set Prior ProbTau \");
     print(\"**** Failed set, PriorProbTau is now\"); flush.console();
     print(\"**** Return to go find error! \"); flush.console();
    eval(parse(text=SetGText(\"OldPriorProbTau\", \"globalenv()\", S=1)));
    try(eval(parse(text=SetGText(\"tauEndList\", \"globalenv()\", S=1)))); 
    eval(parse(text=SetGText(\"MBS\", \"globalenv()\", S=1))); 
    eval(parse(text=SetGText(\"PriorProbTau\", \"globalenv()\", S=1))); 
    eval(parse(text=SetGText(\"PriorProbFixed\", \"globalenv()\", S=1))); 
    return(MBS);
  }
} else {
  MyEval <- 1;
}
if (is.null(MyEval)) {
  print(\"MBS: Error, attempt to set PriorProbTau, PriorProbFixed ended in fail\"); flush.console(); 
  try(eval(parse(text=SetGText(\"tauEndList\", \"globalenv()\", S=1)))); 
  eval(parse(text=SetGText(\"MBS\", \"globalenv()\", S=1))); 
  eval(parse(text=SetGText(\"PriorProbTau\", \"globalenv()\", S=1))); 
  eval(parse(text=SetGText(\"PriorProbFixed\", \"globalenv()\", S=1))); 
}
if (Verbose >= 2) {
  print(\"BayesSpikeRegression: Finished installing PriorProbTau, PriorProbFixed.\"); flush.console();
}

  
  ")
}
BayesSpikeSetupTemperatureLists <- function() {
  return("
  MBS$HowSample = MBS$TBSR5$HowSample;
##MBS$SetupMTForOnii(2);
##MBS$MaximizeOn();

##XX = (1:400)/10;
##MyF= XX * 0;
##for (ii in 1:length(XX)) { MyF[ii] = MBS$GiveLOfX(XX[ii])}
##MBS$GiveLOfX(1);
##MBS$GiveLOfX(.1);
##MBS$SampleANewTau(0);

if (!is.null(TemperatureList)  && length(TemperatureList) >= 1) {
  if (Verbose > -1) {
    print(\"BayesSpikeRegression: Installing a nonNull TemepratureList\"); flush.console();
  }
  if (is.list(TemperatureList)) {
    print(\"BayesSpikeRegression: Error, TemperatureList should be a vec, converting\");
    if (length(TemperatureList[[1]]) == 1) {
      TemperatureList <- unlist(TemperatureList);
    } else {
       TemperatureList <- TemperatureList[[1]]; 
    }
  }
  try(MBS$TemperatureList <- TemperatureList);
  if (Verbose > -1) {
    print(\"BayesSpikeRegression: TemepratureList Installed\"); flush.console();
  }
  try(MBS$TemperatureDecreasingRate <- TBSR5$TemperatureDecreasingRate); 
  ItTemps = 1:length(TemperatureList);
}  else {
  ItTemps = 1;
}
  MBS$RevertTemperatureEvery <- RevertTemperatureEvery;

if (MBS$Verbose >= 2) {
  print(paste(\"BayesSpikeRegression: About to Start, Total Temps = \", ItTemps,
    sep=\"\")); flush.console();
}
try(MBS$SecondInitiateTime <- proc.time())
  if (MBS$Verbose > 1) {
    print(paste(\"Starting with STT = \", paste(STT, collapse=\", \"))); flush.console();
  }
  try(eval(parse(text=SetGText(\"STT\", S=1))));
  ")
}
BayesSpikeStartSetupOnIterACoda <- function() {
return("
  if (!exists(\"TTii\")) { TTii = 1; }
     if (MBS$Verbose >= 1) {
    print(paste(\"BayesSpikeRegression:  On TTii = \", TTii, \"/\", length(ItTemps),
      sep=\"\")); flush.console();
  }
  
  if (MBS$p >= 20000 ||(TTii == length(ItTemps) && sum(MBS$DoRecord[6:7]) < 2)) {
    MBS$SetupRunProbVector(MBS$burnin);
    if (RegionWidth > 1) {
      MBS$RegionWidth <- RegionWidth;
    }
  }
    try(MBS$TBSR5$AllTempCodaLists[[TTii]] <- MBS$CodaList);
    try(MBS$AllTempCodaLists <- MBS$TBSR5$AllTempCodaLists);

  if (!is.null(MBS$TemperatureList) && length(MBS$TemperatureList) > 1) {
    try(MBS$Tempii <- TTii - 1);
  }
  if (MBS$TBSR5$DoSave == FALSE || 
    (MBS$TBSR5$SaveDir %in% c(\"NoSave\", \"NOSAVE\") ) ) {
    try(MBS$NoSave <- -1);  
  }

");
}
BayesSpikeCodaListSetAndCheck <- function() {
  return("
  ##MBS$Verbose = 10;
##MBS$EarlyEndStep = -1;  MBS$EarlyEndtt = -1;
if (DoneByList == FALSE && sum(abs(MBS$TBSR5$DoRecord)) > 0) {
  if (MBS$Verbose >= --2) {
    print(\"BayesSpikeRegression: DoRecord > 0 so Generating Coda Table\"); flush.console();
  }
  BZ = -1;
  ATryText = \"
  MBS$TBSR5$GenerateCodaTable(Verbose = Verbose);
  MBS$OnCodaTable = 1;
  BZ = 1;
  \";
  try(eval(parse(text=ATryText)));
  if (is.null(BZ) || !is.numeric(BZ) || BZ != 1) {
    print(paste(\"BayesSpikeRegression: TTii=\", 
      TTii, \"/\", length(ItTemps), 
      \": something went wrong trying to generate coda table!.\", sep=\"\"));
    flush.console();
  }
  if (Verbose >= 1) {
    print(paste(\"BayesSpikeRegression: Finished Generating CodaTable length = \", length(MBS$CodaList),
      \":(\", paste(dim(MBS$CodaList[[1]]), collapse=\", \"), \"), TTii=\", TTii,
      \"/\", length(ItTemps), sep=\"\")); flush.console();
  }  
}
if (sum(abs(DoRecord)) > 0 && is.null(MBS$CodaList)) {
  print(\"BayesSpikeRegression: Bad news after Genreate CodaList nothing happened!\");
  flush.console();
}
if (Verbose > 1) {
  print(\"BayesSpikeRgression: checking out CodaList.\"); flush.console();
}
OnChainIter = 1;
if (DoneByList == FALSE) {
if (!is.null(MBS$CodaList)) {
  if (Verbose > 1) {
    print(\"BayesSpikeRegression: Setting STT\"); flush.console();
  }
  eval(parse(text=GetG0Text(\"STT\", S=1)));
  try(STT <- OnChainIter:length(MBS$CodaList));
  if (is.numeric(STT) && min(STT) == 0) {
    try(STT <- 1:MBS$TBSR5$NumChains);
  }
  if (!is.null(DoChainIters)) {
    if (!is.null(STT[STT %in% DoChainIters])) {
      STT <- STT[STT %in% DoChainIters];
    }
  }
  try(eval(parse(text=SetGText(\"STT\", S=1))));
  try(MBS$OnCodaTable <- STT[1]);
  try(MBS$MaxGibbsIters <- length(MBS$CodaList[[1]][,1]));
} else if (MBS$TBSR5$NumChains > 0) {
  eval(parse(text=GetG0Text(\"STT\", S=1)));
  try(STT <- 1:MBS$TBSR5$NumChains);
  eval(parse(text=SetGText(\"STT\", S=1)));
} else {
  eval(parse(text=GetG0Text(\"STT\", S=1)));
  if (Verbose > 1) {
    print(\"Unfortunately as feature of STT, must set STT to NULL\");
  }
  STT = NULL;
  eval(parse(text=SetGText(\"STT\", S=1)));
}
}
  ")
}
################################################################################
##  BayesSpikeRegression
##
##    Main Interface Function for BayesSpike
##
##   Creates TBSRoo (a R.oo object containing Rlevel info), 
##    and then MBS (a RCpp object holding both preserved SEXPs and dynamic memory)
##
##   MBS$RunAlgorithm is called to run algorithm from the RCpp object
##    All of the function calls attached to MBS are in .cpp code (in src)
##    Function calls attached to TBSRoo are located in this document "ABayesSpikeStart.R"
BayesSpikeRegression <- function(X=NULL,Y=NULL, BetaStart = NULL,
  IndexFirstRandomEffect = 1, tauEndList=NULL,
  tauStartValues = NULL, PiAStart = .5, SigmaSq = 1, dfTauStart = -1, MStart = -1,
  dfTNoise = 0, Verbose = 0,
  DependenciesTau = NULL, DependenciesFixed=NULL,
  DoMax = 1, 
  PiAPrior = c(-1,-1), HowSample = 3,
  NumSpliceSampleReps = 5, SpliceMove= 1, CauchyEpsilon = .00001,
  MyPS = NULL,  tauFixed = 40,
  TypeFixedPrior = 1, DoRNGState = 1,
  SigmaPrior = c(2,1), InitKKs = 10, DoRecord = c(0,0,0,0,0,0,0),
  NumSamples = 1000, MaxIters = 1000, MaxGibbsIters = 1000,
  MaximizeMeCauchyTotalIters = 15,
  ttStart = 0,  NumChains = 3,
  EarlyEndStep = -1, EarlyEndtt = -1, DoSave = TRUE, DoSaveTBSR5=TRUE,
  FileSave = "MyS",
  SaveDir = .DefaultFileDir,  FileName = "ABSChain",
  NewWrite = 1, tauPriordf = 1, tauPriorMean = 1, 
  PutRNG = TRUE, Run = TRUE,
  CodaTableNames = NULL, OnChainIter = 1, RStart = 3,
  NoNoiseBetaStart = FALSE, AFD = NULL, NoShrinkFixed = NULL,
  NoShrinkFixedPrior = NULL, NoShrinkRandom=NULL,
  NoShrinkRandomPrior =NULL, StartRunProbVector = -100,
  PriorProbFixed=NULL, PriorProbTau = NULL, dfRobit = -1, 
  TemperatureList = NULL,  RevertTemperatureEvery = -1, 
  burnin = 1, EEMergeEvery = NULL, DoEEProbSort = 0, 
  EEProbSortWidth = -1, PrintIter = 200, 
  TemperatureDecreasingRate = 1.0, FirstRandom = -1,
  RandomInfoList = NULL, DoLogitPostProb = FALSE, DoLogitPreProb = FALSE, 
  WriteYBuffer = FALSE, NewYBufferWrite = TRUE,
  WriteWeightBuffer=FALSE, NewWeightBufferWrite=TRUE,
  LengthYBuffer=100, LengthWeightBuffer = 100, EEDoSort = TRUE,
  StartOldEELoad = 50, DoChainIters = NULL,
  AlterWeightFlag = FALSE,
  ReInputList = -1,  DoTimePartsList = TRUE,
  DoLongCI = FALSE, CpBuffLongCI = -1,
  ZeroOutBeforeSample = 0, KillPct = -1, MinimumProbKeeper = -1,
  MaximumAllocation = -1, StartRemoval = -1,
  CheckRemovable =-1, DoSampleBFixedii = 1, RegionWidth = -1,
  gPrior = 0.0, StartRecordingiiWeight = 3, PreRunMaxGibbsIters = 20,
  ...) {

eval(parse(text=EvalReDeclarers()));
eval(parse(text=EvalXandY()));
eval(parse(text=EvalTemperatureEE()));
eval(parse(text=CheckXandY()));
eval(parse(text=CheckSaveDir()));
eval(parse(text=CreateTBSR5()));
eval(parse(text=SetupToCreateMBS()));

eval(parse(text=CreateAndLockInMBS()));  ## Real Generation of MBS
eval(parse(text=SetupOnPiAAndHowSample()));
eval(parse(text=SetupDoTimePartsList()));
eval(parse(text=SetupTauEndList()));
eval(parse(text=SetupAlternateNoise()));
eval(parse(text=SetupInsertPointers()));
eval(parse(text=BSSetupEigenList()));

eval(parse(text=BayesSpikeSetupMBSOnSigmaAndStuff()));
eval(parse(text=BayesSpikeMergeAndPriors()));
eval(parse(text=BayesSpikeSetupTaus()));
eval(parse(text=BayesSpikeChecksForMemory()));

eval(parse(text=BayesSpikeSetupProbs()));
eval(parse(text=BayesSpikeUpToEarlyEndStep()));
eval(parse(text=BayesSpikeSetupDependencies()));
eval(parse(text=BayesSpikeSetupItersAndSave()));
eval(parse(text=BayesSpikeSetupSaveFileToBeDestroyed()));
TempRunFlag <- 1;
if (Run == FALSE) {
  TTii = 1; if (!exists("ItTemps")) { ItTemps = c(1); }
  eval(parse(text=BayesSpikeCodaListSetAndCheck())); 
  eval(parse(text=BayesSpikeCodaListSetter())); 
  ##return(MBS);
  if (is.null(TempRunFlag) || (is.numeric(TempRunFlag) && TempRunFlag < 0)) {
    print("Run Fail return MBS. "); flush.console(); return(MBS);
  }
} else {
  eval(parse(text=BayesSpikeCodaListSetter()));
  if (is.null(TempRunFlag) || (is.numeric(TempRunFlag) && TempRunFlag < 0)) {
    print("Run Fail return GibbsSS. "); flush.console(); return(MBS);
  }
  ##return(MBS);
}
eval(parse(text=BayesSpikeSetupSaveToWrite()));

eval(parse(text=BayesSpikeSetupPriorProbFixedTau()));
eval(parse(text=BayesSpikeSetupTemperatureLists()));
if (DoneByList == FALSE && Run == TRUE) {

for (TTii in ItTemps)  {
   eval(parse(text=BayesSpikeStartSetupOnIterACoda()));
   eval(parse(text=BayesSpikeCodaListSetAndCheck()));
   eval(parse(text=BayesSpikeSetupEarlyEndttNow()));
for (iti in 1:length(STT)) {
  eval(parse(text=BayesSpikeOnChainIterInitiate()));
  eval(parse(text=BayesSpike0020ZeroOutNow()));
  MBS$UpdateNonFreshXtResid();
  eval(parse(text=BayesSpike0021SetupToGo()));

  try(MyTry <- MBS$RunAlgorithm());
  eval(parse(text=BayesSpike0021AfterRunAlgorithmTest()));
  if (!is.null(MBS$TemperatureList) &&
    length(MBS$TemperatureList) > 1 && is.null(MBS$CodaILocFile)) {
    print("Doh: You want temperatures but MBS$CodaILocFile is null!"); flush.console();
  } 
  if (!is.null(MBS$TemperatureList) &&
    length(MBS$TemperatureList) > 1 && is.null(MBS$CodaProbFile)) {
    print("Even though Doh: MBS$CodaProbFile is NULL"); flush.console();
  }     
  if (!is.null(MBS$TemperatureList) &&  !is.null(MBS$CodaILocFile) && MBS$Temperature != 1.0) {
    try(MBS$SetCodaOldLocAndProbFile(MBS$CodaILocFile, MBS$CodaProbFile,
      MBS$burnin, MBS$MaxGibbsIters, OnChainIter, MBS$TBSR5$DoEEProbSort));
    try(MBS$OldTemperature <- MBS$Temperature);
  }
  ##print(paste("End of Algorithm, MBS$Verbose = ", MBS$Verbose, " and tt = ", MBS$tt, sep=""));
  ##print(paste("CodaTable[1,] = (", paste(round(MBS$CodaTable[1,],2), collapse=", "), sep=""));
  
  if (MyTry == 0) {
    print("MyTry returns early exit.");
    return(MBS);
  }
  if (MBS$Verbose >= 1) {
    print(paste("BayesSpikeRegression: Finished Algorithm, Temperature-Tii=", TTii, "/",
      length(ItTemps), ", Chain: ", 
      MBS$TBSR5$OnChainIter, "/", length(STT), ".", sep="")); flush.console();
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
  } else if (!is.null(TBSR5$PutRNG) || TBSR5$PutRNG == 0) {
  
  } else {
    MBS$PutRNGState();
  }
  ## MBS$TBSR5 = TBSR5; TBSR5$MBS = MBS;
  ATrySaveText <- "
    if (!is.null(MBS$SaveDir) && is.character(MBS$SaveDir) &&
      MBS$SaveDir[1] != \"\" && MBS$TBSR5$DoSaveTBSR5 == TRUE) {
      try(MBS$Save(), silent=TRUE);  
    }
  ";
  try(eval(parse(text=ATrySaveText)), silent=TRUE);
}
}

  try(MBS$CompleteTime <- proc.time())
   if (DoneByList == TRUE && 
     exists("ReInputList") && !is.null(ReInputList) &&
    is.list(ReInputList) && length(ReInputList) > 0) {
    try(eval(parse(text=FinalAttachTryCode())));
    ##print("Finished  RunningFinalAttachTryCode"); flush.console();
    ##return(MBS);
    ATText <- "
    if (Verbose >= 0) {
       print(\"LoadIn: Successfully ReInputList.  FinalAttachTrycode Returns \"); flush.console();
    }";
    try(eval(parse(text=ATText)));
    return(MBS);
  } else {
  if (MBS$Verbose >= 0) {
    print("BayesSpikeRegression: no more lock TBSR5 and check."); flush.console();
  }
  return(MBS);
  ##try(eval(parse(text=GetG0Text(".ListBayesSpikeOb", BAYESSPIKENAMESPACE))));
  if (FALSE) {
  if (length(TBSR5$NInstance) >= length(.ListBayesSpikeOb)) {
    print("BayesSpikeRegression, .ListBayesSpikeOb not as long as saved object!"); flush.console();
  } else {
    ##try(MBS$CheckTBSR5(.ListBayesSpikeOb[[TBSR5$NInstance]]), silent=TRUE); 
  }}
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
  ##try(MBS$TBSR5$MBS <- NULL);
  ATrySaveText <- "
    if (!is.null(MBS$SaveDir) && is.character(MBS$SaveDir) &&
      MBS$SaveDir[1] != \"\" && MBS$TBSR5$DoSaveTBSR5==TRUE) {
      try(MBS$Save(), silent=TRUE);  
    }
  ";
  try(eval(parse(text=ATrySaveText)), silent=TRUE);
  return(MBS); 
}