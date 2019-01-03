.onLoad <- function(libname, pkgname) {
  if (!exists("flush.console")) {
    library(utils);
    try(flush.console <- utils::flush.console);
  }
## load the module and store it in our namespace
  BAYESSPIKENAMESPACE <- environment()
  TryO = NULL;
  try( TryO <- bindingIsActive("modBayesSpikeCL", BAYESSPIKENAMESPACE), silent = TRUE);
  if (!is.null(TryO) && (TryO == TRUE || TryO == FALSE)) {
   unlockBinding( "modBayesSpikeCL" , BAYESSPIKENAMESPACE )
  }
  assign( "modBayesSpikeCL", Module( "modBayesSpikeCL" ), BAYESSPIKENAMESPACE )
  lockBinding( "modBayesSpikeCL", BAYESSPIKENAMESPACE )
  try(eval(parse(text=LockGText("modBayesSpikeCL", "globalenv()"))), silent=TRUE);
  try(eval(parse(text=LockGText("modBayesSpikeCL", "BAYESSPIKENAMESPACE"))), silent=TRUE);
  eval(parse(text=GetG0Text("Verbose", S=1)));
  
  
  .onLoadAcpp(libname, pkgname);
  if (Verbose >= 5) {
    print("ABayesSpikeStart: onLoad: starting"); flush.console();
  }   
  if (is.null(modBayesSpikeCL)) {
    print("*******************************************");
    print("** OnLoad ! BayesSpikeCL error declaring modBayesSpikeCL is NULL at start!");
    flush.console();
  }
  eval(parse(text=GetG0Text("BayesSpikeCL", S=1))); 
  BayesSpikeCL <- modBayesSpikeCL$BayesSpikeCL;
  eval(parse(text=SetGText("BayesSpikeCL", S=1)));
  
  FinalizeCAFD <- function(CAFD, Verbose = 0, ...) {
    print("FinalizeCAFD Start"); flush.console();
    Verbose = 20;
    if (!is.null(Verbose) && Verbose > 0) {
      print("FinalizeCAFD: Calling "); flush.console();
    }
    try(finalize(CAFD));
    return;
  }
  try(eval(parse(text=LockGText("FinalizeCAFD", "globalenv()"))), silent=TRUE);
  try(eval(parse(text=LockGText("FinalizeCAFD", "BAYESSPIKENAMESPACE"))), silent=TRUE);
  
  if (Verbose >= 5) {
    print("Set FinalizeCAFD, on to ATry Functions."); flush.console();
  }

  eval(parse(text=GetG0Text("ATryRunRegression", S=1)));
  ATryRunRegression <- function(TBSR5 = NULL) {
    if (is.null(TBSR5)) {
      print("ATryRunRegression: Faulty TBSR5 Supplied");
      flush.console();
      return(NULL);
    }
    if (is.null(TBSR5$Verbose) || !(is.logical(TBSR5$Verbose) || is.numeric(TBSR5$Verbose)) ) {
      print("ATryRunRegression: No Verbose to TBSR5, it must be false.");
      flush.console();
      return(NULL);
    }
    if (TBSR5$Verbose > 1) {
      print("ATryRunRegression: Try to call RunRegression!"); flush.console();
    }
    if (is.null(TBSR5$MBS)) {
      print("ATryRunRegression: No MBS attached to TBSR5!"); flush.console();
    }
    try(RunRegression(TBSR5$MBS));
    if (TBSR5$Verbose > 1) {
      print("ATryRunRegression: Finished"); flush.console();
    }
    return(1);
  }
  try(eval(parse(text=LockGText("ATryRunRegression", "globalenv()"))), silent=TRUE);
  try(eval(parse(text=LockGText("ATryRunRegression", "BAYESSPIKENAMESPACE"))), silent=TRUE);
  
  eval(parse(text=GetG0Text("ATryFillAlterWeightCodaList", S=1)));
  ATryFillAlterWeightCodaList <- function(TBSR5 = NULL) {
    if (is.null(TBSR5)) {
      print("ATryFillAlterWeightCodaList: Faulty TBSR5 Supplied");
      flush.console();
      return(NULL);
    }
    if (is.null(TBSR5$Verbose) || !(is.logical(TBSR5$Verbose) || is.numeric(TBSR5$Verbose)) ) {
      print("ATryFillAlterWeightCodaList: No Verbose to TBSR5, it must be false.");
      flush.console();
      return(NULL);
    }
    if (TBSR5$Verbose > 1) {
      print("ATryFillAlterWeightCodaList: Try to call FillAlterWeight!"); flush.console();
    }
    try(TBSR5$GenerateAlterWeightCodaList());
    if (TBSR5$Verbose > 1) {
      print("ATryFillAlterWeightCodaList: Finished"); flush.console();
    }
    return(TBSR5$.AlterWeightCodaList);
  }
  try(eval(parse(text=LockGText("ATryFillAlterWeightCodaList", "globalenv()"))), silent=TRUE);
  try(eval(parse(text=LockGText("ATryFillAlterWeightCodaList", "BAYESSPIKENAMESPACE"))), silent=TRUE);
  eval(parse(text=GetG0Text("ATryCopy", S=1)));
  ATryCopy <- function(TBSR5 = NULL) {
    if (is.null(TBSR5)) {
      print("ATryCopy: Faulty TBSR5 Supplied");
      flush.console();
      return(NULL);
    }
    if (is.null(TBSR5$Verbose) || !(is.logical(TBSR5$Verbose) || is.numeric(TBSR5$Verbose)) ) {
      print("ATryCopy: No Verbose to TBSR5, it must be false.");
      flush.console();
      return(NULL);
    }
    if (TBSR5$Verbose > 1) {
      print("ATryCopy: Try to call TBSR%$CopyAll!"); flush.console();
    }
    try(AC <- TBSR5$CopyAll());
    if (TBSR5$Verbose > 1) {
      print("ATryCopy: Finished"); flush.console();
    }
    return(AC);
  }
  try(eval(parse(text=LockGText("ATryCopy", "globalenv()"))), silent=TRUE);
  try(eval(parse(text=LockGText("ATryCopy", "BAYESSPIKENAMESPACE"))), silent=TRUE);
  eval(parse(text=GetG0Text("ATryFillSubCodaList", S=1)));
  ATryFillSubCodaList <- function(TBSR5 = NULL) {
    if (is.null(TBSR5)) {
      print("ATryFillSubCodaList: Faulty TBSR5 Supplied");
      flush.console();
      return(NULL);
    }
    if (is.null(TBSR5$Verbose) || !(is.logical(TBSR5$Verbose) || is.numeric(TBSR5$Verbose)) ) {
      print("ATryFillSubCodaList: No Verbose to TBSR5, must be false.");
      flush.console();
      return(NULL);
    }
    if (TBSR5$Verbose > 1) {
      print("ATryFillSubCodaList: Try to call SubCodaList!"); flush.console();
    }
    EndIter = 0;
    if (TBSR5$MaxGibbsIters >= 0) {
      EndIter = TBSR5$MaxGibbsIters;
    } else {
      EndIter = TBSR5$NumSamples;
    }
    try(TBSR5$LoadSubCodaList(TBSR5$MBS$burnin, EndIter));
    if (TBSR5$Verbose > 1) {
      print("ATryFillSubCodaList: Finished"); flush.console();
    }
    return(TBSR5$.SubCodaList);
  }
  try(eval(parse(text=LockGText("ATryFillSubCodaList", "globalenv()"))), silent=TRUE);
  try(eval(parse(text=LockGText("ATryFillSubCodaList", "BAYESSPIKENAMESPACE"))), silent=TRUE);
  
   eval(parse(text=GetG0Text("ATryFillTauCodaList", S=1)));
  ATryFillTauCodaList <- function(TBSR5 = NULL) {
    if (is.null(TBSR5)) {
      print("ATryFillTauCodaList: Faulty TBSR5 Supplied");
      flush.console();
      return(NULL);
    }
    if (is.null(TBSR5$Verbose) || !(is.logical(TBSR5$Verbose) || is.numeric(TBSR5$Verbose)) ) {
      print("ATryFillTauCodaList: No Verbose to TBSR5, must be false.");
      flush.console();
      return(NULL);
    }
    if (TBSR5$Verbose > 1) {
      print("ATryFillTauCodaList: Try to call SubCodaList!"); flush.console();
    }
    EndIter = 0;
    if (TBSR5$MaxGibbsIters >= 0) {
      EndIter = TBSR5$MaxGibbsIters;
    } else {
      EndIter = TBSR5$NumSamples;
    }
    try(TBSR5$GenerateTauCodaList());
    if (TBSR5$Verbose > 1) {
      print("ATryFillTauCodaList: Finished"); flush.console();
    }
    return(TBSR5$.TauCodaList);
  }
  try(eval(parse(text=LockGText("ATryFillTauCodaList", "globalenv()"))), silent=TRUE);
  try(eval(parse(text=LockGText("ATryFillTauCodaList", "BAYESSPIKENAMESPACE"))), silent=TRUE);
  eval(parse(text=GetG0Text("ATryFillSubCodaLongList", S=1)));
  ATryFillSubCodaLongList <- function(TBSR5 = NULL) {
    if (is.null(TBSR5)) {
      print("ATryFillSubCodaLongList: Faulty TBSR5 Supplied");
      flush.console();
      return(NULL);
    }
    if (is.null(TBSR5$Verbose) || !(is.logical(TBSR5$Verbose) || is.numeric(TBSR5$Verbose)) ) {
      print("ATryFillSubCodaLongList: No Verbose to TBSR5, must be false.");
      flush.console();
      return(NULL);
    }
    if (TBSR5$Verbose > 1) {
      print("ATryFillSubCodaLongList: Try to call SubCodaLongList!"); flush.console();
    }
    EndIter = 0;
    if (TBSR5$MaxGibbsIters >= 0) {
      EndIter = TBSR5$MaxGibbsIters;
    } else {
      EndIter = TBSR5$NumSamples;
    }
    try(TBSR5$LoadSubCodaLongList(TBSR5$MBS$burnin, EndIter));
    if (TBSR5$Verbose > 1) {
      print("ATryFillSubCodaLongList: Finished"); flush.console();
    }
    return(TBSR5$SubCodaLongList);
  }
  try(eval(parse(text=LockGText("ATryFillSubCodaLongList", "globalenv()"))), silent=TRUE);
  try(eval(parse(text=LockGText("ATryFillSubCodaLongList", "BAYESSPIKENAMESPACE"))), silent=TRUE);
  eval(parse(text=GetG0Text("ATrySetSubCodaSubSetCoords", S=1)));
  ATrySetSubCodaSubSetCoords <- function(TBSR5 = NULL, SubSetCoords = NULL) {
    if (is.null(TBSR5)) {
      print("ATrySetSubCodaSubSetCoords: Faulty TBSR5 Supplied");
      flush.console();
      return(NULL);
    }
    if (is.null(SubSetCoords)) {
      print("ATrySetSubCodaSubSetCoords: SubSetCoords Supplied Is NULL!");
      return(NULL);
    }
    if (is.null(TBSR5$Verbose) || !(is.logical(TBSR5$Verbose) || is.numeric(TBSR5$Verbose)) ) {
      print("ATryFillSubCodaList: No Verbose to TBSR5, it mut be false.");
      flush.console();
      return(NULL);
    }
    CSubSetCoords = as.integer(sort(unique(as.integer(SubSetCoords[SubSetCoords >= 1.0 & 
      SubSetCoords <= TBSR5$p]-1))));
    if (is.null(CSubSetCoords) || length(CSubSetCoords) <= 0) {
      print("ATrySetSubCodaSubSetCoords: faulty indices suplied!"); flush.console();
      return(NULL);
    }
    try(TBSR5$.SubSetCoords <- as.integer(CSubSetCoords));
    try(TBSR5$ABayesSpikeCL$CSubSetCoords <- as.integer(CSubSetCoords));
    try(TBSR5$.SubCodaList <- NULL);  
    try(TBSR5$SubCodaLongList <- NULL);
    if (TBSR5$Verbose > 1) {
      print("ATrySetSubCodaSubSetCoords: Tried to set SubSetCoords!"); flush.console();
    }
    if (TBSR5$Verbose > 1) {
      print("ATrySetSubCodaSubSetCoords: Finished"); flush.console();
    }
    return(TBSR5$.SubSetCoords);
  }
  try(eval(parse(text=LockGText("ATrySetSubCodaSubSetCoords", "globalenv()"))), silent=TRUE);
  try(eval(parse(text=LockGText("ATrySetSubCodaSubSetCoords", "BAYESSPIKENAMESPACE"))), silent=TRUE);
  eval(parse(text=GetG0Text("ATrySetSubSetTau", S=1)));
  ATrySetSubSetTau <- function(TBSR5 = NULL, SubSetTau = NULL) {
    if (is.null(TBSR5)) {
      print("ATrySetSubSetTau: Faulty TBSR5 Supplied");
      flush.console();
      return(NULL);
    }
    if (is.null(SubSetTau)) {
      print("ATrySetSubSetTau: SubSetTau Supplied Is NULL!");
      return(NULL);
    }
    if (is.null(TBSR5$Verbose) || !(is.logical(TBSR5$Verbose) || is.numeric(TBSR5$Verbose)) ) {
      print("ATrySetSubSetTau: No Verbose to TBSR5, it mut be false.");
      flush.console();
      return(NULL);
    }
    CSubSetTau = as.integer(sort(unique(as.integer(SubSetTau[SubSetTau >= 1.0 & 
      SubSetTau <= length(TBSR5$OnTau)]-1))));
    if (is.null(CSubSetTau) || length(CSubSetTau) <= 0) {
      print("ATrySetSubSetTau: faulty indices suplied!"); flush.console();
      return(NULL);
    }
    try(TBSR5$.SubSetTau <- as.integer(CSubSetTau));
    try(TBSR5$.TauCodaList <- NULL);  
    try(TBSR5$ABayesSpikeCL$TauCodaList <- NULL);
    try(TBSR5$ABayesSpikeCL$CSubSetTau <- CSubSetTau);
    if (TBSR5$Verbose > 1) {
      print("ATrySetSubSetTau: Tried to set SubSetTau!"); flush.console();
    }
    if (TBSR5$Verbose > 1) {
      print("ATrySetSubSetTau: Finished"); flush.console();
    }
    return(TBSR5$.SubSetTau);
  }
  try(eval(parse(text=LockGText("ATrySetSubSetTau", "globalenv()"))), silent=TRUE);
  try(eval(parse(text=LockGText("ATrySetSubSetTau", "BAYESSPIKENAMESPACE"))), silent=TRUE);
  
  eval(parse(text=GetG0Text("ATrySave", S=1)));
  ATrySave <- function(TBSR5 = NULL) {
    if (is.null(TBSR5)) {
      print("ATrySave: Faulty TBSR5 Supplied");
      flush.console();
      return(NULL);
    }
    if (is.null(TBSR5$Verbose) || 
      !(is.logical(TBSR5$Verbose) || is.numeric(TBSR5$Verbose)) ) {
      print("ATrySave: No Verbose to TBSR5, it must be false.");
      flush.console();
      return(NULL);
    }
    if (TBSR5$Verbose > 1) {
      print("ATrySave: Try to call save!"); flush.console();
    }
    try(TBSR5$save());
    if (TBSR5$Verbose > 1) {
      print("ATrySave: Finished"); flush.console();
    }
    return(1);
  }
  try(eval(parse(text=LockGText("ATrySave", "globalenv()"))), silent=TRUE);
  try(eval(parse(text=LockGText("ATrySave", "BAYESSPIKENAMESPACE"))), silent=TRUE);

  eval(parse(text=GetG0Text("ATryUnSave", S=1)));
  ATryUnSave <- function(TBSR5 = NULL, IReallyWantToDoThis=FALSE) {
    if (is.null(TBSR5)) {
      print("ATryUnSave: Faulty TBSR5 Supplied");
      flush.console();
      return(NULL);
    }
    if (!exists("IReallyWantToDoThis") || is.null(IReallyWantToDoThis)) {
      print("ATryUnSave: Sorry you don't really want to do this!")
    }
    if (is.null(TBSR5$Verbose) || 
      !(is.logical(TBSR5$Verbose) || is.numeric(TBSR5$Verbose)) ) {
      print("ATryUnSave: No Verbose to TBSR5, it must be false.");
      flush.console();
      return(NULL);
    }
    if (TBSR5$Verbose > 1) {
      print("ATryUnSave: Try to call Unsave!"); flush.console();
    }
    try(TBSR5$unsave(IReallyWantToDoThis));
    if (TBSR5$Verbose > 1) {
      print("ATryUnSave: Finished"); flush.console();
    }
    return(1);
  }
  try(eval(parse(text=LockGText("ATryUnSave", "globalenv()"))), silent=TRUE);
  try(eval(parse(text=LockGText("ATryUnSave", "BAYESSPIKENAMESPACE"))), silent=TRUE);  
  
  eval(parse(text=GetG0Text("ATryFillMIP", S=1)));
  ATryFillMIP <- function(TBSR5 = NULL) {
    if (is.null(TBSR5)) {
      print("ATryFillMIP: Faulty TBSR5 Supplied");
      flush.console();
      return(NULL);
    }
    if (is.null(TBSR5$Verbose) || !(is.logical(TBSR5$Verbose) || is.numeric(TBSR5$Verbose)) ) {
      print("ATryFillMIP: No Verbose to TBSR5, it mut be false.");
      flush.console();
      return(NULL);
    }
    if (TBSR5$Verbose > 1) {
      print("ATryFillMIP: Try to call FillMIP!"); flush.console();
    }
    try(TBSR5$GenerateMIP());
    if (TBSR5$Verbose > 1) {
      print("ATryFillMIP: Finished"); flush.console();
    }
    return(TBSR5$.MIP);
  }
  try(eval(parse(text=LockGText("ATryFillMIP", "globalenv()"))), silent=TRUE);
  try(eval(parse(text=LockGText("ATryFillMIP", "BAYESSPIKENAMESPACE"))), silent=TRUE);
    
  
  eval(parse(text=GetG0Text("ATryFillPostProbCoda", S=1)));
  ATryFillPostProbCoda <- function(TBSR5 = NULL) {
    if (is.null(TBSR5)) {
      print("ATryFillPostProbCoda: Faulty TBSR5 Supplied");
      flush.console();
      return(NULL);
    }
    if (is.null(TBSR5$Verbose) || !(is.logical(TBSR5$Verbose) || is.numeric(TBSR5$Verbose)) ) {
      print("ATryFillPostProbCoda: No Verbose to TBSR5, it mut be false.");
      flush.console();
      return(NULL);
    }
    if (TBSR5$Verbose > 1) {
      print("ATryFillPostProbCoda: Try to call FillPostProbCoda!"); flush.console();
    }
    try(TBSR5$GenerateProbCodaList());
    if (TBSR5$Verbose > 1) {
      print("ATryFillpostProbCoda: Finished"); flush.console();
    }
    return(TBSR5$.ProbCodaList);
  }
  try(eval(parse(text=LockGText("ATryFillPostProbCoda", "globalenv()"))), silent=TRUE);
  try(eval(parse(text=LockGText("ATryFillPostProbCoda", "BAYESSPIKENAMESPACE"))), silent=TRUE);
eval(parse(text=GetG0Text("ATryFillAlterWeightProbCoda", S=1)));
  ATryFillAlterWeightProbCoda <- function(TBSR5 = NULL) {
    if (is.null(TBSR5)) {
      print("ATryFillAlterWeightProbCoda: Faulty TBSR5 Supplied");
      flush.console();
      return(NULL);
    }
    if (is.null(TBSR5$Verbose) || !(is.logical(TBSR5$Verbose) || is.numeric(TBSR5$Verbose)) ) {
      print("ATryFillAlterWeightProbCoda: No Verbose to TBSR5, it mut be false.");
      flush.console();
      return(NULL);
    }
    if (TBSR5$Verbose > 1) {
      print("ATryFillAlterWeightProbCoda: Try to call FillAlterWeightCoda!"); flush.console();
    }
    try(TBSR5$GenerateAlterWeightCodaList());
    if (TBSR5$Verbose > 1) {
      print("ATryFillAlterWeightCoda: Finished"); flush.console();
    }
    return(TBSR5$.ProbCodaList);
  }
  try(eval(parse(text=LockGText("ATryFillAlterWeightProbCoda", "globalenv()"))), silent=TRUE);
  try(eval(parse(text=LockGText("ATryFillAlterWeightProbCoda", "BAYESSPIKENAMESPACE"))), silent=TRUE);
  eval(parse(text=GetG0Text("ATryFillYCodaList", S=1)));
  ATryFillYCodaList <- function(TBSR5 = NULL) {
    if (is.null(TBSR5)) {
      print("ATryFillYCodaList: Faulty TBSR5 Supplied");
      flush.console();
      return(NULL);
    }
    if (is.null(TBSR5$Verbose) || !(is.logical(TBSR5$Verbose) || is.numeric(TBSR5$Verbose)) ) {
      print("ATryFillYCodaList: No Verbose to TBSR5, it mut be false.");
      flush.console();
      return(NULL);
    }
    if (TBSR5$Verbose > 1) {
      print("ATryFillYCodaList: Try to call FillPostProbCoda!"); flush.console();
    }
    try(TBSR5$GenerateYCodaList());
    if (TBSR5$Verbose > 1) {
      print("ATryFillYCodaList: Finished"); flush.console();
    }
    return(TBSR5$.YCodaList);
  }
  try(eval(parse(text=LockGText("ATryFillYCodaList", "globalenv()"))), silent=TRUE);
  try(eval(parse(text=LockGText("ATryFillYCodaList", "BAYESSPIKENAMESPACE"))), silent=TRUE);

  try(eval(parse(text=GetG0Text("ATryFillRsDeCenteredCodaList", "globalenv()", S=1))), silent=TRUE);
  ATryFillRsDeCenteredCodaList <- function(TBSR5 = NULL) {
    if (is.null(TBSR5)) {
      print("ATryFillRsDeCenteredCodaList: Faulty TBSR5 Supplied");
      flush.console();
      return(NULL);
    }
    if (is.null(TBSR5$Verbose) || !(is.logical(TBSR5$Verbose) || is.numeric(TBSR5$Verbose)) ) {
      print("ATryFillRsDeCenteredCodaList: No Verbose to TBSR5, it must be false.");
      flush.console();
      return(NULL);
    }
    if (TBSR5$Verbose >= 1) {
      print("ATryFillRsDeCenteredCodaList: About to AFunctionNewUnCenteredCodaList"); flush.console();
    }
    try(ABB <- .AFunctionNewUnCenteredCodaList(TBSR5$MBS));
    if (is.null(ABB)) {
      print("ATryFillRsDeCenteredCodaList: Well, the .AFunctionNew returned a NULL object!"); flush.console();
    }  
  }
  try(eval(parse(text=LockGText("ATryFillRsDeCenteredCodaList", "globalenv()"))), silent=TRUE);
  try(eval(parse(text=LockGText("ATryFillRsDeCenteredCodaList", "BAYESSPIKENAMESPACE"))), silent=TRUE);  
  
  try(eval(parse(text=GetG0Text("ATryFillWeightCodaList", S=1))));
  ATryFillWeightCodaList <- function(TBSR5 = NULL) {
    if (is.null(TBSR5)) {
      print("ATryFillWeightCodaList: Faulty TBSR5 Supplied");
      flush.console();
      return(NULL);
    }
    if (is.null(TBSR5$Verbose) || !(is.logical(TBSR5$Verbose) || is.numeric(TBSR5$Verbose)) ) {
      print("ATryFillWeightCodaList: No Verbose to TBSR5, it mut be false.");
      flush.console();
      return(NULL);
    }
    if (TBSR5$Verbose > 1) {
      print("ATryFillWeightCodaList: Try to call FillPostProbCoda!"); flush.console();
    }
    try(TBSR5$GenerateWeightCodaList());
    if (TBSR5$Verbose > 1) {
      print("ATryFillWeightCodaList: Finished"); flush.console();
    }
    return(TBSR5$.WeightCodaList);
  }
  try(eval(parse(text=LockGText("ATryFillWeightCodaList", "globalenv()"))), silent=TRUE);
  try(eval(parse(text=LockGText("ATryFillWeightCodaList", "BAYESSPIKENAMESPACE"))), silent=TRUE);
  try(eval(parse(text=GetG0Text("ATryFillPiACodaList", S=1))));
  ATryFillPiACodaList <- function(TBSR5 = NULL) {
    if (is.null(TBSR5)) {
      print("ATryFillPiACodaList: Faulty TBSR5 Supplied");
      flush.console();
      return(NULL);
    }
    if (is.null(TBSR5$Verbose) || !(is.logical(TBSR5$Verbose) || is.numeric(TBSR5$Verbose)) ) {
      print("ATryFillPiACodaList: No Verbose to TBSR5, it mut be false.");
      flush.console();
      return(NULL);
    }
    if (TBSR5$Verbose > 1) {
      print("ATryFillPiACodaList: Try to call FillPostProbCoda!"); flush.console();
    }
    try(TBSR5$GeneratePiACodaList());
    if (TBSR5$Verbose > 1) {
      print("ATryFillPiACodaList: Finished"); flush.console();
    }
    return(TBSR5$.PiACodaList);
  }
  try(eval(parse(text=LockGText("ATryFillPiACodaList", "globalenv()"))), silent=TRUE);
  try(eval(parse(text=LockGText("ATryFillPiACodaList", "BAYESSPIKENAMESPACE"))), silent=TRUE);
  try(eval(parse(text=GetG0Text("ATryFillSigCodaList", S=1))));
  ATryFillSigCodaList <- function(TBSR5 = NULL) {
    if (is.null(TBSR5)) {
      print("ATryFillSigCodaList: Faulty TBSR5 Supplied");
      flush.console();
      return(NULL);
    }
    if (is.null(TBSR5$Verbose) || !(is.logical(TBSR5$Verbose) || is.numeric(TBSR5$Verbose)) ) {
      print("ATryFillSigCodaList: No Verbose to TBSR5, it mut be false.");
      flush.console();
      return(NULL);
    }
    if (TBSR5$Verbose > 1) {
      print("ATryFillSigCodaList: Try to call FillSigCoda!"); flush.console();
    }
    try(TBSR5$GenerateSigmaCodaList());
    if (TBSR5$Verbose > 1) {
      print("ATryFillSigCodaList: Finished"); flush.console();
    }
    return(TBSR5$.SigmaCodaList);
  }
  try(eval(parse(text=LockGText("ATryFillSigCodaList", "globalenv()"))), silent=TRUE);
  try(eval(parse(text=LockGText("ATryFillSigCodaList", "BAYESSPIKENAMESPACE"))), silent=TRUE);
  ## FreeTBSRoo deletes TBSRoo from locked background environment when
  ##  the RCpp object for the regressions is deleted.
  ##  TBSRoo will still not be deleted if user keeps pointers to TBSRoo in R
  ##  environment.  
  FreeTBSRoo = function(NInstance,Verbose = 0,...) {
    if (!is.numeric(NInstance)  || NInstance <= 0) {
      print(paste("FreeTBSRoo:: Sorry, supply NInstance"));
    }
    NInstance = NInstance[1];
    try(eval(parse(text=GetG0Text(".ListBayesSpikeOb", "globalenv()"))),silent=TRUE);
    if (is.numeric(.ListBayesSpikeOb) && .ListBayesSpikeOb[1] == 0) {
      if (Verbose > 0) {
        print("FreeTBSRoo, did not get .ListBayesSpikeOb, is already NULL "); flush.console();
      }
      return(NULL);
    }
    if (length(.ListBayesSpikeOb) < NInstance) {
      if (Verbose > 0) {
        print(paste("FreeTBSRoo, could not get .ListBayesSpikeOb of this length, ", NInstance,sep=""));
         flush.console();
      }
      return(NULL);   
    }
    
   TBSRoo = .ListBayesSpikeOb[[NInstance]];
    if (is.null(TBSRoo)) {
      if (Verbose > 0) {
        print(paste("FreeTBSR5, TBSRoo is already NULL", sep="")); flush.console();
      }
      return(NULL);
    }
    print(paste("FreeTBSRoo: We're about to Free ")); flush.console();
    if (is.null(TBSRoo)) {
      print(paste("FreeTBSRoo: Can't Free because we're NULL ")); flush.console();
      return(NULL);}
    if (TBSRoo$Verbose> 0){
      print(paste("Freeing TBSR5 by turning it to null and removing from namespace", sep="")); flush.console();
    }
    if (!is.null(TBSRoo$AFD)) {
      if (TBSRoo$Verbose> 0){
        print(paste("Detatching AFD Connection to TBSRoo ",  sep="")); flush.console();
      }
      try(TBSRoo$AFD <- NULL);
    }
    if(!is.null(TBSRoo$ABayesSpikeCL)) {
      VK <- TBSRoo$ABayesSpikeCL;
      try(VK$TBSRoo$ABayesSpikeCL <- NULL, silent=TRUE);
      try(VK$TBSRoo <- NULL, silent=TRUE)  
      try(TBSRoo$ABayesSpikeCL$TBSRoo <- NULL);
      try(TBSRoo$ABayesSpikeCL<-NULL);
    }
    if (!is.null(TBSRoo$.BeingDestroyed)  && TBSRoo$.BeingDestroyed == 1) {
      return(NULL);
    }
    ##eval(parse(text=GetG0Text(".ListBayesSpikeOb", "globalenv()")));
    ## print(paste("FreeTBSRoo: We just got .ListBayesSpikeOb , NInstance = ", 
    ##   TBSRoo$NInstance, sep="")); flush.console(); 
     if (is.null(.ListBayesSpikeOb) || length(.ListBayesSpikeOb) < TBSR5$NInstance ) {
       print("FreeTBSRoo: We have an error, null ListBayesSpikeOb, or length too short "); flush.console();
     }
    .ListBayesSpikeOb[[TBSRoo$NInstance]] <- NULL;
    try(eval(parse(text=LockGText(".ListBayesSpikeOb", "globalenv()"))),silent=TRUE);  
    print("The .ListBayesSpikeOb has been re attached"); flush.console();
    ##rm(TBSRoo);
    return(NULL);
  }
  try(eval(parse(text=GetG0Text("FreeTBSR5", S=1))));
  FreeTBSR5 = function(NInstance,Verbose = 0,...) {
    if (!is.numeric(NInstance)  || NInstance <= 0) {
      print(paste("FreeTBSR5:: Sorry, supply NInstance"));
    }
    NInstance = NInstance[1];
    try(eval(parse(text=GetG0Text(".ListBayesSpikeOb", "globalenv()"))),silent=TRUE);
    if (is.numeric(.ListBayesSpikeOb) && .ListBayesSpikeOb[1] == 0) {
      if (Verbose > 0) {
        print("FreeTBSR5, did not get .ListBayesSpikeOb, is already NULL "); flush.console();
      }
      return(NULL);
    }
    if (length(.ListBayesSpikeOb) < NInstance) {
      if (Verbose > 0) {
        print(paste("FreeTBSR5, could not get .ListBayesSpikeOb of this length, ", NInstance,sep=""));
         flush.console();
      }
      return(NULL);   
    }  
   
    TBSR5 = .ListBayesSpikeOb[[NInstance]];
    if (is.null(TBSR5)) {
      if (Verbose > 0) {
        print(paste("FreeTBSR5, TBSR5 is already NULL", sep="")); flush.console();
      }
      return(NULL);
    }
    print(paste("FreeTBSR5: We're about to Free ")); flush.console();
    if (is.null(TBSR5)) {
      print(paste("FreeTBSR5: Can't Free because we're NULL ")); flush.console();
      return(NULL);}
    if (TBSR5$Verbose>= 0){
      print(paste("Freeing TBSR5 by turning it to null and removing from namespace", sep="")); flush.console();
    }
    try(eval(parse(text="IKILL = FALSE; TBSR5$KillMyTBSR5(); IKILL = TRUE")));
    if (TBSR5$Verbose >= 0) {
      print(paste("Free TBSR5: Just finnished KillMyTBSR5(), return was IKILL = ", IKILL, sep="")); flush.console();
    }
    if (!is.null(TBSR5$AFD)) {
      if (TBSR5$Verbose> 0){
        print(paste("Detatching AFD Connection to TBSR5 ",  sep="")); flush.console();
      }
      try(TBSR5$AFD <- NULL);
    }
    if(!is.null(TBSR5$ABayesSpikeCL)) { 
      VK <- TBSR5$ABayesSpikeCL;
      if (!is.null(VK$TBSR5)) {
        try(VK$TBSR5$ABayesSpikeCL <- NULL, silent=TRUE);
        try(VK$TBSR5 <- NULL, silent=TRUE);
      }
      try(TBSR5$ABayesSpikeCL$TBSR5 <- NULL, silent=TRUE);
      try(TBSR5$ABayesSpikeCL<-NULL);
    }
    try(TBSR5$ABayesSpikeCL <- NULL);
    if (!is.null(TBSR5$.BeingDestroyed)  && TBSR5$.BeingDestroyed == 1) {
      return(NULL);
    }
    ##eval(parse(text=GetG0Text(".ListBayesSpikeOb", "globalenv()")));
    ## print(paste("FreeTBSRoo: We just got .ListBayesSpikeOb , NInstance = ", 
    ##   TBSRoo$NInstance, sep="")); flush.console(); 
     if (is.null(.ListBayesSpikeOb) || length(.ListBayesSpikeOb) < TBSR5$NInstance ) {
       print("FreeTBSR5: We have an error, null ListBayesSpikeOb, or length too short "); flush.console();
     }
    .ListBayesSpikeOb[[TBSR5$NInstance]] <- NULL;
    try(eval(parse(text=LockGText(".ListBayesSpikeOb", "globalenv()"))),silent=TRUE);  
    print("The .ListBayesSpikeOb has been re attached"); flush.console();
    ##rm(TBSR5);
    return(NULL);
  }
  try(eval(parse(text=LockGText("FreeTBSR5", "globalenv()"))), silent=TRUE);
  try(eval(parse(text=LockGText("FreeTBSRoo", "globalenv()"))), silent=TRUE);
  try(eval(parse(text=LockGText("FreeTBSR5", "BAYESSPIKENAMESPACE"))), silent=TRUE);
  try(eval(parse(text=LockGText("FreeTBSRoo", "BAYESSPIKENAMESPACE"))), silent=TRUE);
    
  try(eval(parse(text=GetG0Text("KillFreeTBSR5", S=1))));
  KillFreeTBSR5 = function(TBSR5,Verbose = 0,...) {
    if (is.null(TBSR5)) {
      if (Verbose > 0) {
        print(paste("KillFreeTBSR5, TBSR5 is already NULL", sep="")); flush.console();
      }
      return(NULL);
    }
    if (Verbose > 0) {
      print(paste("KillFreeTBSR5: We're about to Free ")); flush.console();
    }
    if (TBSR5$Verbose > 0) {
      print(paste("Freeing TBSR5 by turning it to null and removing from namespace", sep="")); flush.console();
    }
    if (!is.null(TBSR5$AFD)) {
      if (TBSR5$Verbose> 0){
        print(paste("Detatching AFD Connection to TBSR5 ",  sep="")); flush.console();
      }
      try(TBSR5$AFD <- NULL);
    }

    if(!is.null(TBSR5$ABayesSpikeCL)) { 
      if (Verbose > 0) {
        print(paste("About to attempt to detach as many ABaysSpikeCL as possible")); flush.console();
      }
      VK <- TBSR5$ABayesSpikeCL;
      if (!is.null(VK$TBSR5)) {
        if (!is.null(VK$TBSR5$ABayesSpikeCL)) {
          try(VK$TBSR5$ABayesSpikeCL <- NULL, silent=TRUE);
        }
        if (!is.null(VK$TBSR5)) {
          try(VK$TBSR5 <- NULL, silent=TRUE);
        }
      }
      if (!is.null(TBSR5$ABayesSpikeCL)) {
        try(TBSR5$ABayesSpikeCL<-NULL, silent=TRUE);
      }
      if (Verbose > 0) {
          print("Did our best to try to detach ABayesSpikeCL from TBSR5");
      }
    }

    if (!is.null(TBSR5$.BeingDestroyed)  && TBSR5$.BeingDestroyed == 1) {
      return(NULL);
    }
    if (TBSR5$Verbose >= 0) {
      print("001BayesSpikeOnLoad.r:KillFreeTBSR5(): About to KillMyTBSR5(). "); flush.console();
    }
    try(eval(parse(text="IKILL = 0; TBSR5$KillMyTBSR5(); IKILL = 1")));
    if (TBSR5$Verbose >= 0) {
      print(paste("001BayesSpikeOnLoad.r:KillFreeTBSR5(): Finished KillMyTBSR5(): IKILL = ", IKILL, sep=""));
    }
    ##eval(parse(text=GetG0Text(".ListBayesSpikeOb", "globalenv()")));
    ## print(paste("FreeTBSRoo: We just got .ListBayesSpikeOb , NInstance = ", 
    ##   TBSRoo$NInstance, sep="")); flush.console(); 
     if (is.null(.ListBayesSpikeOb) || length(.ListBayesSpikeOb) < TBSR5$NInstance ) {
       print("FreeTBSR5: We have an error, null ListBayesSpikeOb, or length too short "); flush.console();
     }
    .ListBayesSpikeOb[[TBSR5$NInstance]] <- NULL;
    try(eval(parse(text=LockGText(".ListBayesSpikeOb", "globalenv()"))),silent=TRUE);  
    print("The .ListBayesSpikeOb has been re attached"); flush.console();
    ##rm(TBSR5);
    return(NULL);
  }
  try(eval(parse(text=LockGText("KillFreeTBSR5", "globalenv()"))), silent=TRUE);
  try(eval(parse(text=LockGText("KillFreeTBSR5", "BAYESSPIKENAMESPACE"))), silent=TRUE);

  TBSRooSetVerbose <- function(TBSRoo, sVerbose) {
    TBSRoo$Verbose[1] = sVerbose[1];
  }
  try(eval(parse(text=LockGText("TBSRooSetVerbose"))));
  eval(parse(text=GetG0Text("TBSR5GetVerbose", S=1)));
  TBSR5SetVerbose <- function(TBSR5, sVerbose) {
    TBSR5$Verbose[1] = sVerbose[1];
  }
  try(eval(parse(text=LockGText("TBSR5SetVerbose"))));
  
  try(eval(parse(text=GetG0Text("ReturnSubCodaFromFiles", "BAYESSPIKENAMESPACE", S=1))));
  if (is.numeric(ReturnSubCodaFromFiles) && length(ReturnSubCodaFromFiles) == 1 && ReturnSubCodaFromFiles == 0) {
     print(".OnLoad, Attempting to do ReturnSubCodaFromFiles, but apparently it sucks!"); flush.console();
  }
  try(eval(parse(text=LockGText("ReturnSubCodaFromFiles", "globalenv()", S=1))));
  ReturnSubCodaFromFiles2 = ReturnSubCodaFromFiles;
  try(eval(parse(text=LockGText("ReturnSubCodaFromFiles2", "globalenv()", S=1))));
  try(eval(parse(text=LockGText("ReturnSubCodaFromFiles2", "BAYESSPIKENAMESPACE", S=1))));  

##############################################################################
##  ReturnNewCountBeta()
##
##  Used in BayesSpikeGibbs.cpp::ReadAndSetPctBetas to read from filess and
##   calculate Posterior inclusion scores based upon compact I Files 
##
##  aBayesSpikeCL$PctBetas activates this function.
ReturnNewCountBeta <- function(
    TBSR5 = NULL,  StartIter = 0, EndIter = 2000,
    Verbose = 0) {
  if (is.null(Verbose)) {
    print("ReturnNewCountBeta:  Issue Verbose is NULL "); Verbose = 4; flush.console();
  }
  if (Verbose > 0) {
    print("ReturnNewCountBeta:  In and we've started \n"); flush.console();
    print(paste("NInstance = ", NInstance, sep="")); flush.console();
  }
  if (is.null(TBSR5)) {
    print("ReturnNewCountBeta: Can't continue, because TBSR5 is NULL!");
    flush.console();
    return(NULL);
  }
  if (is.null(TBSR5$ABayesSpikeCL)) {
    print("ReturnNewCountBeta: Can't continue, TBSR5 has no BayesSpikeCL");
    flush.console();
    return(NULL); 
  }
  MBS = TBSR5$ABayesSpikeCL;
  
  if (is.null(MBS$TBSR5$Verbose)) {
    print("ReturnNewCountBeta: TBSR5$Verbose is NULL not a good sign!");
    flush.console(); Verbose = 3;
  } else if (MBS$TBSR5$Verbose <= 0) {
  } else if (is.null(Verbose) || Verbose <= 0) {
    Verbose = TBSR5$Verbose;
  }
  if (Verbose > 0) {
    print(paste("Trying to Get TBSR5 for NInstance = ", NInstance));
  }

  CodaList = list();
  SaveDir = NULL;
  if (is.null(SaveDir)) {
    if (!is.null(MBS$TBSR5)) {
      if (!is.null(MBS$TBSR5$SaveDir)) {
        SaveDir = MBS$TBSR5$SaveDir;
      }
    }
  }
  if (SaveDir %in% c("NoSave", "NOSAVE")) {
    print("ReturnNewCountBeta: Settting was NOSAVE directory!"); flush.console();
    return(NULL);
  }
  if (is.null(Verbose)) {
    print("ReturnNewCountBeta: Verbose is Null"); flush.console(); Verbose = 3;
  }
  if (Verbose > 0) {
    print("ReturnNewCountBeta: Starting"); flush.console();
  }
  if (is.null(TBSR5)) {
    print("ReturnNewCountBeta:  Error Please give non null TBSRoo");
    return(NULL);
  }
  p = MBS$p;
  NewCountBeta = rep(0, p);
  if (!is.null(MBS$tauEndList) && length(MBS$tauEndList) > 0) {
    NewCountTau = rep(0, length(MBS$TBSR5$tauEndList));
  } else {
    NewCountTau = NULL;
  }

  if (is.null(TBSR5)) {
    if (Verbose > 0) {
      print("TBSR5 is Null so we can't use it"); flush.console();
    }
  }  else if (is.null(TBSR5$ListSaveFiles)) {
    if (Verbose > 0) {
      print("TBSR5 is not null, but it does have null ListSaveFiles"); flush.console();
    }
  }
  IFiles = NULL; DFiles = NULL; iTFiles = NULL; dTFiles = NULL;
  ##if (!is.null(MBS$TBSR5$ListSaveFiles)) {
  ##  MBS$TBSR5$ListSaveFiles <- sort(unique(MBS$TBSR5$ListSaveFiles));
  ##}
  if (!is.null(SaveDir)) {
    Dir = SaveDir;
    if (Verbose > 1) {
      print(paste("We're going to pull files from Save Dir = ", SaveDir, sep=""));
      print(paste("ListSaveFiles = (", paste(unlist(TBSR5$ListSaveFiles), collapse=", "),
         ")", sep="")); flush.console();
      flush.console();
    }
    MyNFiles = unlist(list.files(SaveDir));
    MyNFiles <- MyNFiles[substr(MyNFiles, 1, nchar(TBSR5$FileName)) == 
      TBSR5$FileName];
    if (Verbose > 2) {
      print(paste("Long List of Files in ", SaveDir, " is : ", sep=""));
      print(paste("(", paste(MyNFiles, collapse=", "),
         ")", sep="")); flush.console();
      flush.console();
    }
    IFiles = MyNFiles[substr(MyNFiles, nchar(MyNFiles)-4, nchar(MyNFiles)) ==
      "I.bin"];
    if (length(MBS$TBSR5$ListSaveFiles) >= 1) {
      ATO <- rep(0, length(IFiles));
      AFTT <- sort(unique(unlist(MBS$TBSR5$ListSaveFiles)));
      APT <- strsplit(AFTT, "/");
      NewList <- rep("", length(AFTT));
      for (jj in 1:length(AFTT)) {
        NewList[jj] = APT[[jj]][length(APT[[jj]])];
      }
      APT <- strsplit(NewList, "\\\\");
      NewList <- rep("", length(APT));
      for (jj in 1:length(APT)) {
        NewList[jj] = APT[[jj]][length(APT[[jj]])];
      }
      NewList <- sort(unique(NewList));
      
      for (jj in 1:length(MBS$TBSR5$ListSaveFiles)) {
         ALT <- substr(IFiles, 1, nchar(NewList[jj]));
         AM <-(1:length(ALT))[ ];
         ATO[ATO == 0 & ALT == NewList[jj]] = jj;
      }
      IFiles = IFiles[ATO > 0];
    }
    DFiles = MyNFiles[substr(MyNFiles, nchar(MyNFiles)-4, nchar(MyNFiles)) ==
      "D.bin"];  
    # IFiles = c("ABC1I.bin", "A2CFEI.bin", "ADC2I.bin", "ADEC2I.bin", "AC3I.bin");
    # DFiles = c("ABC1D.bin", "ADCXD.bin", "ADEC2D.bin", "A2CFED.bin");
    MDI = match(substr(DFiles, 1,nchar(DFiles)-5), 
      substr(IFiles, 1,nchar(IFiles)-5));
    MDI[is.na(MDI)] = -10000;
    sMDI = sort(MDI, index=TRUE);
    if (length(sMDI$x[sMDI$x != -10000]) > 0) {
      DFiles = DFiles[sMDI$ix[sMDI$x != -10000]];
      IFiles = IFiles[sMDI$x[sMDI$x != -10000]];      
    }
    
    iTFiles = MyNFiles[substr(MyNFiles, nchar(MyNFiles)-5, nchar(MyNFiles)) ==
      "iT.bin"];               
    dTFiles = MyNFiles[substr(MyNFiles, nchar(MyNFiles)-5, nchar(MyNFiles)) ==
      "dT.bin"]; 
      
    MTIT = match(substr(dTFiles, 1, nchar(dTFiles)-6), 
      substr(iTFiles, 1, nchar(iTFiles)-6));
    MTIT[is.na(MTIT)] = -10000;
    sMTIT = sort(MTIT, index=TRUE);
    if (length(sMTIT$x[sMTIT$x != -10000]) > 0) {
      dTFiles = dTFiles[sMTIT$ix[sMTIT$x != -10000]];
      iTFiles = iTFiles[sMTIT$x[sMTIT$x != -10000]];
    }
      
    MITI = match(substr(iTFiles, 1, nchar(iTFiles)-6), 
      substr(IFiles, 1,nchar(IFiles)-5));
    MITI[is.na(MITI)] = -10000;
    sMITI = sort(MITI, index=TRUE);
    if (length(sMITI$x[sMITI$x != -10000]) > 0) {
      dTFiles = dTFiles[sMITI$ix[sMITI$x != -10000]];
      iTFiles = iTFiles[sMITI$ix[sMITI$x != -10000]];
    }
  }  else if (!is.null(TBSR5$ListSaveFiles)) {
    if (Verbose > 1) {
      print(paste("We're going to work with TBSR5$ListSaveFiles", sep=""));
      print(paste("ListSaveFiles = (", paste(unlist(TBSR5$ListSaveFiles), collapse=", "),
         ")", sep="")); flush.console();
      flush.console();
    }
    Dir = TBSRoo$SaveDir;
    DFiles =  unique(paste(unlist(TBSR5$ListSaveFiles), "D.bin", sep=""));
    IFiles =  unique(paste(unlist(TBSR5$ListSaveFiles), "I.bin", sep=""));
    dTFiles =  unique(paste(unlist(TBSR5$ListSaveFiles), "dT.bin", sep=""));
    iTFiles =  unique(paste(unlist(TBSR5$ListSaveFiles), "iT.bin", sep=""));
  } else {
    print(paste("UhRoh:  TBSR5$ListSaveFiles = ",
      paste(TBSR5$ListSaveFiles, collapse=", "), sep=""));   flush.console();
    print(paste("  And SaveDir = ", SaveDir, sep=""));
    print("That Leaves us pretty S out of luck!"); flush.console();
  }
  if (Verbose > 1) {
    print(paste("The IFiles are: ", paste(IFiles, collapse=", "), sep=""));
    flush.console();
    print(paste("The DFiles are: ", paste(DFiles, collapse=", "), sep=""));
    flush.console();
    print(paste("The dTFiles are: ", paste(dTFiles, collapse=", "), sep=""));
    flush.console();
    print(paste("The iTFiles are: ", paste(iTFiles, collapse=", "), sep=""));
    flush.console();
  }

  if (Verbose > 0) {
    print(" GetNewCountBetas:  Well, we're going into the IFiles"); flush.console();
  }
  if (length(IFiles) > 0) {
  NewCountBeta = NewCountBeta * 0;
  for (ii in 1:length(IFiles)) {  
    NewCoda = NULL;
    if (Verbose > 0) {
      print(paste("  Going to GetNewcountBetas for IFiles[", ii, "] = ",
        IFiles[ii], sep="")); flush.console();
    }
    if (Verbose >= 2) {
      print("Call Statement is:");
      IOnFile = paste(Dir, "/", IFiles[ii], sep="");  
      ATText <- paste("
      IOnFile = \"", IOnFile,"\";
      StartIter = ", StartIter, ";
      EndIter = ", EndIter, ";
      NewCountBeta = rep(0, ", TBSR5$p, ");
      Verbose = ", TBSR5$Verbose -1, ";
      try(AN <- .Call(\"GiveCountFrequency\", IOnFile, StartIter, EndIter, 
        NewCountBeta, Verbose-1), silent=FALSE);
      ", sep="");
      print(cat(ATText));
    }
    IOnFile = paste(Dir, "/", IFiles[ii], sep="");  
    
    try( AN <- .Call("GiveCountFrequency", IOnFile, 
       StartIter, EndIter, NewCountBeta, Verbose-1),  silent=FALSE);
    if (Verbose > 2) {
      print(paste(" First Few NewCountBeta Now : ", 
        paste(NewCountBeta[1:8], collapse=", "), sep=""));
    }
  }
  RTL = (EndIter);
  if (EndIter > TBSR5$ABayesSpikeCL$MaxGibbsIters) { RTL = TBSR5$ABayesSpikeCL$MaxGibbsIters; }
  if (StartIter > 0) { RTL = RTL - StartIter + 1; }
  NewCountBeta = NewCountBeta / ( RTL * length(IFiles));
  }
  if ((length(iTFiles) > 0) && !is.null(NewCountTau)) {
  if (Verbose > 1) {
    print(paste("The iTFiles are: ", paste(iTFiles, collapse=", "), sep=""));
    flush.console();
    print(paste("The dTFiles are: ", paste(dTFiles, collapse=", "), sep=""));
    flush.console();
  }  
  for (ii in 1:length(iTFiles)) {  
    if (Verbose > 0) {
      print(paste("  Going to GiveCodaTauSubSet for iTFiles[", ii, "] = ",
        iTFiles[ii], sep="")); flush.console();
    }
    IOnFile = paste(Dir, "/", iTFiles[ii], sep="");  
    
    try( AN <- .Call("GiveCountFrequency", IOnFile, 
      StartIter, EndIter, NewCountTau, Verbose),  silent=FALSE);
  }
  NewCountTau = NewCountTau / ( RTL * length(iTFiles) );
  }  

  return(list(NewCountBeta=NewCountBeta, NewCountTau=NewCountTau));
}  
  BAYESSPIKENAMESPACE <- environment()
  try(eval(parse(text=LockGText("BAYESSPIKENAMESPACE"))));
  try(eval(parse(text=LockGText("BAYESSPIKENAMESPACE", "BAYESSPIKENAMESPACE"))));  
  try(eval(parse(text=LockGText("ReturnNewCountBeta"))));
  try(eval(parse(text=LockGText("ReturnNewCountBeta", "BAYESSPIKENAMESPACE"))));      
  try(eval(parse(text=LockGText("ReturnSubCodaFromFiles"))));
  try(eval(parse(text=LockGText("ReturnSubCodaFromFiles", "BAYESSPIKENAMESPACE"))));  
  
  DefaultTemperatureList = c(10,5,2,1);
  try(eval(parse(text=LockGText("DefaultTemperatureList", "globalenv()"))));

  MAXDIMXKeep  = 4;
  try(eval(parse(text=LockGText("MAXDIMXKeep", "globalenv()"))));
 TryO = NULL;
 
.ListBayesSpikeOb = list();
.NumBayesSpikeOb = 0;
  try(eval(parse(text=GetG0Text(".DefaultDiallelFileDir","BAYESSPIKENAMESPACE", S=1))));
  .DefaultDiallelFileDir = "";
  try(eval(parse(text=LockGText(".DefaultDiallelFileDir"))));
  try(eval(parse(text=LockGText(".DefaultDiallelFileDir", "BAYESSPIKENAMESPACE"))));
  try(eval(parse(text=LockGText(".ListBayesSpikeOb"))));
  try(eval(parse(text=LockGText(".ListBayesSpikeOb", "BAYESSPIKENAMESPACE"))));  
  try(eval(parse(text=LockGText(".NumBayesSpikeOb"))));
  try(eval(parse(text=LockGText(".NumBayesSpikeOb", "BAYESSPIKENAMESPACE"))));  
  

eval(parse(text=GetG0Text(".DefaultFileDir", "BAYESSPIKENAMESPACE", S=1)));  
for (ii in 1:length(.libPaths())) {
  RL = unlist(list.files(.libPaths()[[ii]]));
  if (any(RL == "BayesSpike")) {
    dir.create(.libPaths()[[ii]], "//BayesSpike", showWarnings=FALSE);
    .DefaultFileDir = paste( .libPaths()[[ii]], "//", "BayesSpike//data//", sep="");
    dir.create(.DefaultFileDir, showWarnings=FALSE);
    break;
  }
}
if(.DefaultFileDir == "") {
  print("Well, No places to save, BayesSpike will require SaveDir!");
} 
  try(eval(parse(text=LockGText(".DefaultFileDir"))), silent=TRUE);
  try(eval(parse(text=LockGText(".DefaultFileDir", "BAYESSPIKENAMESPACE"))),silent=TRUE);


  try(eval(parse(text=GetG0Text(".DefaultFileDir", "BAYESSPIKENAMESPACE"))),silent=TRUE);
  eval(parse(text=GetG0Text(".DefaultDiallelFileDir", "BAYESSPIKENAMESPACE")));
  .DefaultDiallelFileDir = paste(.DefaultFileDir, "//", "Diallel", sep=""); flush.console();
  dir.create(.DefaultDiallelFileDir, showWarnings=FALSE);
  eval(parse(text=LockGText(".DefaultFileDir")));
  eval(parse(text=LockGText(".DefaultFileDir", "BAYESSPIKENAMESPACE")));
  eval(parse(text=LockGText(".DefaultDiallelFileDir")));
  try(eval(parse(text=LockGText(".DefaultDiallelFileDir", "BAYESSPIKENAMESPACE"))),silent=TRUE);  

  eval(parse(text=GetG0Text("Verbose", "BAYESSPIKENAMESPACE", S=1)));
  eval(parse(text=GetG0Text("Verbose", "globalenv()", S=1)));  
  Verbose = 0;
  eval(parse(text=LockGText("Verbose", "BAYESSPIKENAMESPACE", S=1)));
  eval(parse(text=LockGText("Verbose", "globalenv()", S=1)));  
     
  try(library(AlanDirectories), silent=TRUE)
  try(AlanDirectores:::SetSaveHomes(),silent=TRUE);
   
  print("Completed function .OnLoad ");flush.console();
}
