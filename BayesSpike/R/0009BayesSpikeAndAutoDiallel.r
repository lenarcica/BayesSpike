

GetAFDAllCentNames <- function(AFD) {
  nX <- colnames(AFD$X);
  AV <- AFD$AllCenteredRandomVariables;
  CentNKeep <- rep(0, length(AV));
  for (ii in 1:length(AV)) {
    if (any(substr(nX,1,nchar(AV[ii])) == AV[ii])) {
      CentNKeep[ii] <- 1;
    }
  } 
  AllRandom  <- AFD$.AllRandomVariables;
  OtherRandom <- AllRandom[!(AllRandom %in% AV)];
  if (length(OtherRandom) >= 1) {
    OtherNKeep <- rep(0, length(OtherRandom));
    for (ii in 1:length(OtherRandom)) {
      if (any(substr(nX,1,nchar(OtherRandom[ii])) == OtherRandom[ii])) {
        OtherNKeep[ii] <- 1;
      }      
    }
  } else {
    OtherNKeep <- NULL;
  }
  FixedNames <- nX[nX %in% AFD$.AllFixedVariables];
  AFDAllNames <- FixedNames;
  if (!is.null(OtherNKeep) && sum(OtherNKeep) >= 1) {   
    for (ii in 1:length(OtherRandom)) {
      if (OtherNKeep[ii] == 1) {
        AFDAllNames <- c(AFDAllNames,
          nX[substr(nX,1,nchar(OtherRandom[ii])) == OtherRandom[ii]]);
      }
    } 
  }
  if (length(CentNKeep) >= 1) {
    for (ii in 1:length(CentNKeep)) {
      if (CentNKeep[ii] == 1) {
       BP <- nX[substr(nX,1,nchar(AV[ii])) == AV[ii]];
       if (AV[ii] %in% c("SymCrossjk", "Gender:SymCrossjk",
         "ASymCrossjkDkj", "Gender:ASymCrossjkDkj", "AllCrossjk",
         "Gender:AllCrossjk")) {
         jj <- nchar(BP[length(BP)]);
         while(jj > 0 && substr(BP[length(BP)],jj,jj) != ":") {
            jj <- jj-1;
         }   
         AF <- as.numeric(substr(BP[length(BP)], jj+1, nchar(BP[length(BP)])));
         ASS <- paste(substr(BP[length(BP)],1,jj), AF+1, sep="");
         BP <- c(BP, ASS);
       } else {
         AF <- length(BP);
         BP <- c(BP, paste(AV[ii], ":", AF+1, sep=""));
       }
       AFDAllNames <- c(AFDAllNames, BP);
      }
    }
  }
  if (!is.null(OtherNKeep) && sum(OtherNKeep) >= 1) {
    AFDAllNames <- c(AFDAllNames, paste("tau:", OtherRandom[OtherNKeep==1], sep=""));
  }
  if (!is.null(CentNKeep) && sum(CentNKeep) >= 1) {
    AFDAllNames <- c(AFDAllNames, paste("tau:", AV[CentNKeep==1], sep=""));
  }
  AFDAllNames <- c(AFDAllNames, "Sigma:1");
  return(AFDAllNames);
}
################################################################################
##  RenameDeCenteredBSAFDCoda(BSAFD, )
##
## Renaming utility for BayesDiallel companion, essential for when BayesSpike
##  runs on a "Centered" set of covariates, and an uncentering procedure is done.
##
## Decenterd have a few problems comparing to the centered chains.
RenameDeCenteredBSAFDCoda <-   function(BSAFD, TrueParametersList = NULL,  ...) {
  DSList <- BSAFD$DeCenteredCodaList;
  if (is.null(DSList) || length(DSList) <= 0) {
    if (BSAFD$Verbose >= 2) {
       print("0009BayesSpikeAndAutoDiallel.r:: RenameDeCenteredBSAFDCoda can't do because DeCenteredCodaList is not configured.");
       flush.console();
    }
    return(-1);
  }
  if (is.null(BSAFD$CodaList) || length(BSAFD$CodaList) <= 0) {
    if (Verbose >= 1) {
      print("RenameDeCenteredBSAFDCoda: No Rename if there is no CodaList!");
      flush.console();
    }
    return(-1);
  }
  if (!exists("TrueParametersList") || is.null(TrueParametersList) || length(TrueParametersList) <= 1) {
    NamesAFD = GetAFDAllCentNames(BSAFD$AFD)
  } else {                                                             
    NamesAFD = names(TrueParametersList);
  }
  ColBSAFDNot <- colnames(BSAFD$CodaList[[1]]);
  NamesBSAFD = BSAFD$OldDeCenterCodaNames;
  if (sum(
    substr(NamesBSAFD[substr(ColBSAFDNot,1,nchar("Beta")) == "Beta"], 1, nchar("Beta")) ==
    "Beta") <= 3 && !is.null(BSAFD$OldDeCenterCodaNames)) {
    BSAFD$OtherDeCenterCodaNames <- BSAFD$OldDeCenterCodaNames;
    return(1);  
  }
  if (is.null(NamesAFD) || min(nchar(as.character(NamesAFD))) <= 0) {
    print("RenameBSAFDCoda: Sorry, true parameters list is null."); 
      flush.console();
    return(-1);
  }
  

  if (is.null(NamesBSAFD) || length(NamesBSAFD) <= 1 ||
    min(nchar(as.character(NamesBSAFD))) <= 0) {
    print("RenameBSAFDCoda: No Names to colnames(BSAFD$CodaList[[1])");
    flush.console();  
    NamesBSAFD <- paste(BSAFD$OldCodaNames, sep="");
  }
  if (is.null(NamesBSAFD) || length(NamesBSAFD) <= 1 ||
    min(nchar(as.character(NamesBSAFD))) <= 0) {
    print("RenameBSAFDCoda: No Names to OldCodaNames");
    flush.console();  
    MTT <-  -1;
    try(MTT <- BSAFD$TBSR5$GenerateCLN());
    if (is.numeric(MTT) && MTT[1] == -1) {
      print("RenameBSAFDCoda: GenerateCLN failed!\n"); flush.console();
      return(-1);
    }
    NamesBSAFD <- paste(MTT$CLN, sep="");
  }
  if (is.null(NamesBSAFD)) {
    print("RenameDeCenteredBSAFDCoda: We give up, not possible! NamesBSAFD is still NULL!");
    flush.console();
    return(-1);
  }
  if (is.null(NamesAFD)) {
    print("RenameDeCenteredBSAFDCoda: We give up, NamesAFD is NULL!");
    flush.console();
    return(-1);
  }
  MyMatch = NULL;
  try(MyMatch <- MatchBSAFDNames(NamesBSAFD, NamesAFD));
  if (all(is.na(MyMatch))) {
    print("RenameDeCenteredBSAFDCoda: Oh no, we have nothing but NAs!"); flush.console();
    return(-1);
  }
  try(NewColNames <- colnames(BSAFD$DeCenteredCodaList[[1]]));
  if (is.null(NewColNames) || length(NewColNames) <= 1 ||
    min(nchar(as.character(NewColNames))) <= 0) {
    try(NewColNames <- BSAFD$OldCodaNames);  
  }
  try(NewColNames[(1:length(MyMatch[,1]))[!is.na(MyMatch[,1])]] <- 
    NamesAFD[(1:length(MyMatch[,1]))[!is.na(MyMatch[,1])]]);
  ##NewColNames[is.na(MyMatch[,1])] = NamesBSAFD[is.na(MyMatch[,1])];
  try(NMatchProbTau <- MatchProbTau(NamesBSAFD, NamesAFD));
  try(NewColNames[NMatchProbTau$OrderedTau] <-
    NMatchProbTau$NewNames);
  if (any(substr(NamesBSAFD, 1, nchar("ProbFixed")) %in% c("ProbFixed",
    "probfixed", "Probfixed"))) {
     MySii = (1:length(NewColNames))[
      substr(NewColNames, 1, nchar("ProbFixed:")) %in% 
      c("ProbFixed:", "Probfixed", "probfixed") |
      substr(NewColNames, 1, nchar("ProbFixed")) %in%
      c("ProbFixed", "Probfixed", "probfixed")]
     ANum = as.numeric(substr(NewColNames[MySii], nchar("ProbFixed:")+1,
       nchar(NewColNames[MySii])));
     ANum = ANum[!is.na(ANum)]
     MySet = c(NewColNames[NewColNames %in% c("Mu", "mu", "Beta:Mu", "Beta:mu",
       "beta:Mu", "beta:mu", "Intercept", "intercept",
       "Gender:Av", "gender:av", "Beta:Gender:Av", "BetaGenderAv",
       "beta:Gender:Av", "Beta:Gender:Av", "Gender", "gender",
       "female", "Female", "Female-Male",
       "BetaHybrid", "betahybrid", "BetaHybrid:Av", 
       "betahybrid:Av", "Inbred:Av", "Inbreed:Av", "inbreed:av",
       "Inbred", "Inbreed", "inbred", "inbreed",
       "BetaInbred", "BetaInbreed", "BetaInbred:Av", "BetaInbreed:Av",
       "Beta:Inbred", "Beta:Inbred:Av", 
       "Gender:BetaHybrid", "BetaHybridTimesSex",
       "BetaHybrid:Gender:Av", "BetaHybridGenderAv",
       "betahybrid:gender:av", "Hybrid:Gender:Av",
       "Inbred:Gender:Av", "inbred:gender:av","inbredgenderav",
       "inbredgender", "Inbred:Gender")],
       NewColNames[substr(NewColNames, 1, nchar("FixedEffect")) %in%
         c("fixedeffect", "FixedEffect") |
         substr(NewColNames, 1, nchar("Fixed")) %in%
         c("fixed", "Fixed")  |
         substr(NewColNames, 1, nchar("Beta:Fixed")) %in%
         c("Beta:Fixed", "beta:fixed")
         ]);
       AF = (1:length(ANum))[ANum > 0 & ANum <= length(MySet)]  
       try(NewColNames[MySii[AF]] <-  
         paste("Prob:Fixed:",MySet[ANum[AF]], sep=""));
  }
  NewDim = list(NULL,NewColNames);
  if (is.null(NewColNames)) {
    print("Well Nothings good on rename because NewColNames is NULL!");
    flush.console();
  }
  if (BSAFD$Verbose > 0) {
    print("About to call Rename CodaList!"); flush.console();
    print(paste("NewDim = ", paste(NewDim[[2]], collapse=", "), sep=""));
    flush.console();
  }
  try(BSAFD$OtherDeCenterCodaNames <- NewColNames);
  NewDeCenteredCodaList <- list();
  for (ii in 1:length(BSAFD$DeCenteredCodaList)) {
    NewDeCenteredCodaList[[ii]] = matrix(BSAFD$DeCenteredCodaList[[ii]],
      1:NROW(BSAFD$DeCenteredCodaList[[1]]),
      1:NCOL(BSAFD$DeCenteredCodaList[[1]]));
    try(colnames(NewDeCenteredCodaList[[ii]]) <- NewColNames);
    try(NewDeCenteredCodaList[[ii]] <- as.mcmc(NewDeCenteredCodaList[[ii]]));
  }
  try(NewDeCenteredCodaList < as.mcmc.list(NewDeCenteredCodaList));
  try(BSAFD$DeCenteredCodaList <- NewDeCenteredCodaList);
  ##try(BSAFD$RenameCodaList(NewDim));
  if (BSAFD$Verbose > 0) {
    print("After Rename Coda List, names of column 1 is");
    print(paste(colnames(BSAFD$DeCenteredCodaList[[1]]), collapse=", "));
  }
  return(1);
}

################################################################################
##  RenameBSAFDCoda()
##
##   See BayesSpikeFromAutoDiallel(AFD) function.
##
##  Once a BayesDiallel type design has been run by BayesSpike, this
##   is code to try to rename the CodaList object of BSAFD (centered or not)
##   to the conventions typical for the BayesDiallel user.  This includes
##   insertion of strain identifiers into columns.
##
RenameBSAFDCoda <- function(BSAFD, TrueParametersList = NULL, DoDeCentered=TRUE, ...) {
  if (!exists("TrueParametersList")) { TrueParametersList <- NULL; }
  if (!exists("DoDeCentered")) { DoDeCentered <- FALSE; }
  if (is.null(BSAFD)) { return(NULL);}
  if (is.null(BSAFD$AFD)) { print("Sorry, no AFD in BSAFD"); flush.console();}
  if (is.null(BSAFD$CodaList) || length(BSAFD$CodaList) <= 0) {
    if (Verbose >= 1) {
      print("RenameDeCenteredBSAFDCoda: No Rename if there is no CodaList!");
      flush.console();
    }
    return(-1);
  }
  if (DoDeCentered==TRUE && !is.null(BSAFD$CenteredColumns)) {
     RenameDeCenteredBSAFDCoda(BSAFD);
  }
  NamesBSAFD = BSAFD$OldCodaNames;

  if (is.null(NamesBSAFD) || length(NamesBSAFD) <= 1 ||
    min(nchar(as.character(NamesBSAFD))) <= 0) {
    print("RenameBSAFDCoda: No Names to colnames(BSAFD$CodaList[[1])");
    flush.console();  
    NamesBSAFD <- paste(BSAFD$OldCodaNames, sep="");
  }
  if (is.null(NamesBSAFD) || length(NamesBSAFD) <= 1 ||
    min(nchar(as.character(NamesBSAFD))) <= 0) {
    print("RenameBSAFDCoda: No Names to OldCodaNames");
    flush.console();  
    MTT <-  -1;
    try(MTT <- BSAFD$TBSR5$GenerateCLN());
    if (is.numeric(MTT) && MTT[1] == -1) {
      print("RenameBSAFDCoda: GenerateCLN failed!\n"); flush.console();
      return(-1);
    }
    NamesBSAFD <- paste(MTT$CLN, sep="");
  }

  if (is.null(NamesBSAFD) || min(nchar(as.character(NamesBSAFD)))  <= 0) {
    print("RenameBSAFDCoda: At start Names BSAFD is NULL names!");
    flush.console();
    return(-1);
  }
  if (is.null(TrueParametersList) || length(TrueParametersList) <= 1) {
    nX <- colnames(BSAFD$AFD$X);
    NKeep <- rep(0, length(BSAFD$AFD$.AllRandomVariables));
    for (ij in 1:length(NKeep)) {
      if (any(substr(nX,1,nchar(BSAFD$AFD$.AllRandomVariables[ij])) == 
         BSAFD$AFD$.AllRandomVariables[ij])) {
        NKeep[ij] = 1;  
      }
    }
    if (sum(NKeep) >= 1) {
      NamesAFD <- c(colnames(BSAFD$AFD$X), paste("tau:", 
        BSAFD$AFD$.AllRandomVariables[NKeep==1], sep=""), "Sigma:1");
    } else {
      NamesAFD <- c(colnames(BSAFD$AFD$X), "Sigma:1")        
    }
  } else {                                                             
    NamesAFD = names(TrueParametersList);
  }
  if (is.null(NamesAFD) || min(nchar(as.character(NamesAFD))) <= 0) {
    print("RenameBSAFDCoda: Sorry, true parameters list is null."); 
      flush.console();
    return(-1);
  }
  if (any(substr(NamesBSAFD, nchar(NamesBSAFD) - nchar("SymCrossjk:j:1;k:1")+1,
    nchar(NamesBSAFD)) == "SymCrossjk:j:1;k:1") &&
    !any(NamesAFD == "SymCrossjk:j:1;k:1")) {
    if (any(NamesAFD == "SymCrossjk:j:2;k:1")) {
      idi = (1:length(NamesAFD))[NamesAFD == "SymCrossjk:j:2;k:1"]
      STI = NamesAFD[substr(NamesAFD,1,nchar("aj:")) == "aj:"];
      numj = max(as.numeric(substr(STI, nchar("aj:") +1, nchar(STI))))
      ATP = paste("SymCrossjk:j:", (1:numj), ";k:", 1:numj, sep="");
      NamesAFD = c(NamesAFD[1:(idi-1)], ATP, NamesAFD[idi:length(NamesAFD)])
    }
  }
  if (any(substr(NamesBSAFD, 
    nchar(NamesBSAFD) - 
      nchar("Gender:SymCrossjk:j:1;k:1")+1, nchar(NamesBSAFD))
    == "Gender:SymCrossjk:j:1;k:1") &&
     !any(NamesAFD == "Gender:SymCrossjk:j:1;k:1")) {
    if (any(NamesAFD == "SymCrossjk:j:2;k:1")) {
      idi = (1:length(NamesAFD))[NamesAFD == "Gender:SymCrossjk:j:2;k:1"]
      STI = NamesAFD[substr(NamesAFD,1,nchar("aj:")) == "aj:"];
      numj = max(as.numeric(substr(STI, nchar("aj:") +1, nchar(STI))))
      ATP = paste("Gender:SymCrossjk:j:", (1:numj), ";k:", 1:numj, sep="");
      NamesAFD = c(NamesAFD[1:(idi-1)], ATP, NamesAFD[idi:length(NamesAFD)])
    }
  }
  if (any(substr(NamesBSAFD, 
    nchar(NamesBSAFD) - nchar("AllCrossjk:j:1;k:1")+1, nchar(NamesBSAFD))
    == "AllCrossjk:j:1;k:1") &&
    !any(NamesAFD == "AllCrossjk:j:1;k:1")){
    if (any(NamesAFD == "AllCrossjk:j:2;k:1")) {
      idi = (1:length(NamesAFD))[NamesAFD == "AllCrossjk:j:2;k:1"]
      STI = NamesAFD[substr(NamesAFD,1,nchar("aj:")) == "aj:"];
      numj = max(as.numeric(substr(STI, nchar("aj:") +1, nchar(STI))))
      ATP = paste("AllCrossjk:j:", (1:numj), ";k:", 1:numj, sep="");
      NamesAFD = c(NamesAFD[1:(idi-1)], ATP, NamesAFD[idi:length(NamesAFD)])
    }
  }
  if (any(substr(NamesBSAFD, 
    nchar(NamesBSAFD) - nchar("Gender:AllCrossjk:j:1;k:1")+1, 
    nchar(NamesBSAFD))
    == "Gender:AllCrossjk:j:1;k:1") &&
    !any(NamesBSAFD == "Gender:AllCrossjk:j:1;k:1")) {
    if (any(NamesAFD == "Gender:AllCrossjk:j:2;k:1")) {
      idi = (1:length(NamesAFD))[NamesAFD == "Gender:AllCrossjk:j:2;k:1"]
      STI = NamesAFD[substr(NamesAFD,1,nchar("aj:")) == "aj:"];
      numj = max(as.numeric(substr(STI, nchar("aj:") +1, nchar(STI))))
      ATP = paste("Gender:AllCrossjk:j:", (1:numj), ";k:", 1:numj, sep="");
      NamesAFD = c(NamesAFD[1:(idi-1)], ATP, NamesAFD[idi:length(NamesAFD)])
    }
  }

  if (is.null(NamesBSAFD) || length(NamesBSAFD) <= 1 ||
    min(nchar(as.character(NamesBSAFD))) <= 0) {
    print("RenameBSAFDCoda: No Names to colnames(BSAFD$CodaList[[1])");
    flush.console();  
    NamesBSAFD <- paste(BSAFD$OldCodaNames, sep="");
  }
  if (is.null(NamesBSAFD) || length(NamesBSAFD) <= 1 ||
    min(nchar(as.character(NamesBSAFD))) <= 0) {
    print("RenameBSAFDCoda: No Names to OldCodaNames");
    flush.console();  
    MTT <-  -1;
    try(MTT <- BSAFD$TBSR5$GenerateCLN());
    if (is.numeric(MTT) && MTT[1] == -1) {
      print("RenameBSAFDCoda: GenerateCLN failed!\n"); flush.console();
      return(-1);
    }
    NamesBSAFD <- paste(MTT$CLN, sep="");
  }
  MyMatch = NULL;
  try(MyMatch <- MatchBSAFDNames(NamesBSAFD, NamesAFD));
  try(NewColNames <- colnames(BSAFD$CodaList[[1]]));
  if (is.null(NewColNames) || length(NewColNames) <= 1 ||
    min(nchar(as.character(NewColNames))) <= 0) {
    try(NewColNames <- BSAFD$OldCodaNames);  
  }
  try(NewColNames[(1:length(MyMatch[,1]))[!is.na(MyMatch[,1])]] <- 
    NamesAFD[(1:length(MyMatch[,1]))[!is.na(MyMatch[,1])]]);
  ##NewColNames[is.na(MyMatch[,1])] = NamesBSAFD[is.na(MyMatch[,1])];
  try(NMatchProbTau <- MatchProbTau(NamesBSAFD, NamesAFD));
  try(NewColNames[NMatchProbTau$OrderedTau] <-
    NMatchProbTau$NewNames);
  if (any(substr(NamesBSAFD, 1, nchar("ProbFixed")) %in% c("ProbFixed",
    "probfixed", "Probfixed"))) {
     MySii = (1:length(NewColNames))[
      substr(NewColNames, 1, nchar("ProbFixed:")) %in% 
      c("ProbFixed:", "Probfixed", "probfixed") |
      substr(NewColNames, 1, nchar("ProbFixed")) %in%
      c("ProbFixed", "Probfixed", "probfixed")]
     ANum <- suppressWarnings(as.numeric(substr(NewColNames[MySii], nchar("ProbFixed:")+1,
       nchar(NewColNames[MySii]))));
     ANum = ANum[!is.na(ANum)]
     MySet = c(NewColNames[NewColNames %in% c("Mu", "mu", "Beta:Mu", "Beta:mu",
       "beta:Mu", "beta:mu", "Intercept", "intercept",
       "Gender:Av", "gender:av", "Beta:Gender:Av", "BetaGenderAv",
       "beta:Gender:Av", "Beta:Gender:Av", "Gender", "gender",
       "female", "Female", "Female-Male",
       "BetaHybrid", "betahybrid", "BetaHybrid:Av", 
       "betahybrid:Av", "Inbred:Av", "Inbreed:Av", "inbreed:av",
       "Inbred", "Inbreed", "inbred", "inbreed",
       "BetaInbred", "BetaInbreed", "BetaInbred:Av", "BetaInbreed:Av",
       "Beta:Inbred", "Beta:Inbred:Av", 
       "Gender:BetaHybrid", "BetaHybridTimesSex",
       "BetaHybrid:Gender:Av", "BetaHybridGenderAv",
       "betahybrid:gender:av", "Hybrid:Gender:Av",
       "Inbred:Gender:Av", "inbred:gender:av","inbredgenderav",
       "inbredgender", "Inbred:Gender")],
       NewColNames[substr(NewColNames, 1, nchar("FixedEffect")) %in%
         c("fixedeffect", "FixedEffect") |
         substr(NewColNames, 1, nchar("Fixed")) %in%
         c("fixed", "Fixed")  |
         substr(NewColNames, 1, nchar("Beta:Fixed")) %in%
         c("Beta:Fixed", "beta:fixed")
         ]);
       AF = (1:length(ANum))[ANum > 0 & ANum <= length(MySet)]  
       try(NewColNames[MySii[AF]] <-  
         paste("Prob:Fixed:",MySet[ANum[AF]], sep=""));
  }
  NewDim = list(NULL,NewColNames);
  if (is.null(NewColNames)) {
    print("Well Nothings good on rename because NewColNames is NULL!");
    flush.console();
  }
  if (BSAFD$Verbose > 0) {
    print("About to call Rename CodaList!"); flush.console();
    print(paste("NewDim = ", paste(NewDim[[2]], collapse=", "), sep=""));
    flush.console();
  }
  if (any(is.na(NewColNames))) {
    print("Error in RenameCodaList: about to set OtherNameCodaList but we have NAs!");
    print("We're saying to level"); flush.console();
    eval(parse(text=SetGText("NewColNames", "globalenv()", S=1)));
    return(-1);
  }
  try(BSAFD$OtherNameCodaList <- NewColNames);
  try(BSAFD$RenameCodaList(NewDim));
  if (BSAFD$Verbose > 0) {
    print("After Rename Coda List, names of column 1 is");
    print(paste(colnames(BSAFD$CodaList[[1]]), collapse=", "));
  }
  return(1);
}
  
  
#####################################################################
##   ReFillBSAFDAFDCodaList  Helper function or "BayesDiallel" duties
##
##   Code that helps to attach correct names to a coda list for combination
##    of BayesSpike with BayesDiallel
##
##           Alan Lenarcic 10/30/2013 2013
##  
ReFillBSAFDAFDCodaList <- function(BSAFD,...) {
  if (is.null(BSAFD$AFD)) {
    print("Error: ReFillBSAFDAFDCodaList Sorry No AFD");
  }
  for (ii in 1:length(BSAFD$AFD$AllDiallelObs[[1]]$CodaChains)) {
    if (length(BSAFD$CodaList) >= ii) {
      MyM = match( colnames(BSAFD$CodaList[[ii]]),
        colnames(BSAFD$AFD$AllDiallelObs[[1]]$CodaChains[[ii]]));
        
    }
  }
  return(1);
}
MatchProbTau <- function( NewColNames, NamesAFD) {
  TauNames = NamesAFD[substr(NamesAFD,1,4) == "tau:"]
  OrderedTau = (1:length(NewColNames))[
    substr(NewColNames,1,nchar("ProbTau:")) == "ProbTau:"];  
  list( OrderedTau = OrderedTau, NewNames = paste("Prob:", TauNames, sep="")); 

}
#######################################################################
##  Defaullt BayesSpike names do not match FullDiallelNAmes
##  This matches what names are valid.  

#######################################################################
##  Defaullt BayesSpike names do not match FullDiallelNAmes
##  This matches what names are valid.  
MatchBSAFDNames <- function(NamesBS, NamesAFD,...) {
  if (is.null(NamesBS) ||  length(NamesBS) <= 0 ||
    min(nchar(as.character(NamesBS))) <= 0) {
      print("Hey Issue MatchBSAFDNames, NameBS to start is NULL!");
      flush.console();
    }
    NamesBSBetaIter <- (1:length(NamesBS))[substr(NamesBS,1,4) == "Beta"];
    NBs <- strsplit(NamesBS[NamesBSBetaIter], "Beta:[[:digit:]]*:");
    try(NBSBeta <- rep("", length(NamesBS[NamesBSBetaIter])));
    for (ii in 1:length(NBs)) {
      try(NBSBeta[ii] <- NBs[[ii]][length(NBs[[ii]])]);
    }
    NamesAFDBetaIter = (1:length(NamesAFD))[
      substr(NamesAFD,1,3) != "tau" &
      substr(NamesAFD,1,5) != "Sigma"];
    NAFDBeta = NamesAFD[NamesAFDBetaIter];
    try(NABs <- strsplit(NAFDBeta, "Beta:[[:digit:]]*:"));
    for (ii in 1:length(NAFDBeta)) {
      try(NAFDBeta[ii] <- NABs[[ii]][length(NABs[[ii]])]);
    }
    try(NABs <- strsplit(NAFDBeta, "Beta:"));
    for (ii in 1:length(NAFDBeta)) {
      try(NAFDBeta[ii] <- NABs[[ii]][length(NABs[[ii]])]);
    }
    matchBeta = match(NBSBeta, NAFDBeta);
    if (all(is.na(matchBeta)) && length(NBSBeta) == length(NAFDBeta)) {
       matchBeta <- 1:length(NBSBeta);
    }
    matchAllBeta = NamesAFDBetaIter[matchBeta];

    NamesBStauIter = (1:length(NamesBS))[substr(NamesBS,1,3) == "tau"];
    NBs = strsplit(NamesBS[NamesBStauIter], "tau:[[:digit:]]*:");
    try(NBStau <- rep("", length(NamesBS[NamesBStauIter])));
    for (ii in 1:length(NBs)) {
      try(NBStau[ii] <- NBs[[ii]][length(NBs[[ii]])]);
    }
    NBs = strsplit(NBStau, "tau:");
    for (ii in 1:length(NBs)) {
      try(NBStau[ii] <- NBs[[ii]][length(NBs[[ii]])]);
    }
    NamesAFDtauIter = (1:length(NamesAFD))[
      substr(NamesAFD,1,3) == "tau"];
    NAFDtau = NamesAFD[NamesAFDtauIter];
    try(NABs <- strsplit(NAFDtau, "tau:[[:digit:]]*:"));
    for (ii in 1:length(NAFDtau)) {
      try(NAFDtau[ii] <- NABs[[ii]][length(NABs[[ii]])]);
    }
    try(NABs <- strsplit(NAFDtau, "tau:"));
    for (ii in 1:length(NAFDtau)) {
      try(NAFDtau[ii] <- NABs[[ii]][length(NABs[[ii]])]);
    }
    matchtau = match(NBStau, NAFDtau); 
    if (all(is.na(matchtau)) && length(NBStau) == length(NAFDtau)) {
      matchtau <- 1:length(NBStau);
    }   
    matchAlltau = NamesAFDtauIter[matchtau];
    
    NamesBSSigmaIter = (1:length(NamesBS))[substr(NamesBS,1,5) == "Sigma"];
    NBs = strsplit(NamesBS[NamesBSSigmaIter], "Sigma:[[:digit:]]*:");
    try(NBSSigma <- rep("", length(NamesBS[NamesBSSigmaIter])));
    for (ii in 1:length(NBs)) {
      try(NBSSigma[ii] <- NBs[[ii]][length(NBs[[ii]])]);
    }
    NBs = strsplit(NBSSigma, "Sigma:");
    for (ii in 1:length(NBs)) {
      try(NBSSigma[ii] <- NBs[[ii]][length(NBs[[ii]])]);
    }
    NamesAFDSigmaIter = (1:length(NamesAFD))[
      substr(NamesAFD,1,5) == "Sigma"];
    NAFDSigma = NamesAFD[NamesAFDSigmaIter];
    NABs = strsplit(NAFDSigma, "Sigma:[[:digit:]]*:");
    for (ii in 1:length(NAFDSigma)) {
     try(NAFDSigma[ii] <- NABs[[ii]][length(NABs[[ii]])]);
    }
    try(NABs <- strsplit(NAFDSigma, "Sigma:"));
    for (ii in 1:length(NAFDSigma)) {
      try(NAFDSigma[ii] <- NABs[[ii]][length(NABs[[ii]])]);
    }
    matchSigma = match(NBSSigma, NAFDSigma);
    matchAllSigma = NamesAFDSigmaIter[matchSigma];
    BSList = c(NamesBSBetaIter, NamesBStauIter, NamesBSSigmaIter);
    MatchAFDList = c(matchAllBeta, matchAlltau, matchAllSigma);
    return(cbind(BSList, MatchAFDList));  
}


####################################################################
##  BayesSpikeFromAutoDiallel()
##
##  Create a Bayes Spike object from a "DiallelObject", 
## from package "BayesDiallel_x.x.tar.gz".  This extracts information
## So that it can be supplied to "BayesSpikeRegression()" function
##
##
##       
BayesSpikeFromAutoDiallel <- function(AFD, Y=NULL, X=NULL,
  PiAStart = c(.5,.5), 
  PiAPrior = c(20,20,20,20),
  dfTNoise = 0, Verbose = 0,
  DependenciesTau =NULL,DependenciesFixed=NULL, HowSample = 3,
  LengthChains = -1, NumChains = 3,
  BetaStart = NULL,
  IndexFirstRandomEffect = 1, tauEndList=NULL,
  tauStartValues = NULL, dfTauStart = -1, MsStart = -1,
  DoMax = 1,
  NumSpliceSampleReps = 5, SpliceMove= 1, CauchyEpsilon = .00001,
  MyPS = PriorStructureR5$new(.5,1),  tauFixed = 20,
  TypeFixedPrior = 1, DoRNGState = 1,
  SigmaPrior = c(2,1), InitKKs = 5, DoRecord = c(1,1,1,0,0,1,1),
  NumSamples = 1000, MaxIters = 1000, MaxGibbsIters = 1000,
  MaximizeMeCauchyTotalIters = 100,
  ttStart = 0, 
  EarlyEndStep=-1, EarlyEndtt=-1, DoSave = TRUE,
  FileSave = "MyS",
  SaveDir = NULL,  FileName = "AutoDiallelChain",
  NewWrite = 1, tauPriorMean = -1, tauPriordf = -1,
  WatchCoordinates = NULL, Run = FALSE,
  PriorProbTau = NULL, PriorProbFixed = NULL, ZeroOutBeforeSample = FALSE,
  TestEndBayesSpikeAutoDiallel = -1) {
  
  if (!exists("tauFixed")) { tauFixed <- 20; }
  if (!exists("MaximizeMeCauchyTotalIters")) { MaximizeMeCauchyTotalIters <- 100; }
  if (!exists("NumChains")) { NumChains <- 3; }
  if (!exists("CauchyEpsilon")) { CauchyEpsilon <- .00001; }
  if (!exists("DependenciesFixed")) { DependenciesFixed <- NULL; }
  if (!exists("DependenciesTau")) { DependenciesTau <- NULL; }
  if (!exists("HowSample")) { HowSample <- 3; }
  if (!exists("ttStart")) { ttStart <- 0; }
  if (!exists("dfTauStart")) { dfTauStart <- -1; }
  if (!exists("MaxIters")) { MaxIters <- 1000; }
  if (!exists("MaxGibbsIters")) { MaxGibbsIters <- MaxIters; }
  if (!exists("SigmaPrior")) { SigmaPrior <- c(2,1); }
  if (!exists("SpliceMove")) { SpliceMove <- 1; }
  if (!exists("PiAStart")) {   PiAStart = c(.5,.5); }
  if (!exists("PiAPrior")) {   PiAPrior = c(20,20,20,20); }
  if (!exists("PriorProbTau")) { PriorProbTau <- NULL }
  if (!exists("MyPS")) { MyPS <- PriorStructureR5$new(.5,1); }
  if (!exists("dfTauStart")) { dfTauStart <- -1; }
  if (!exists("tauPriorMean")) { tauPriorMean <- -1; }
  if (!exists("tauPriordf")) { tauPriordf <- -1; }
  if (!exists("TypeFixedPrior")) { TypeFixedPrior <- 1; } 
  if (!exists("DoRNGState")) { DoRNGState <- 1; } 
  if (!exists("DoRecord")) { DoRecord <- c(1,1,1,0,0,1,1) }
  if (!exists("InitKKs")) { InitKKs <- 5; } 
  if (!exists("NewWrite")) { NewWrite <- 1; }
  if (!exists("tauStartValues")) { tauStartValues <- NULL; }
  if (!exists("MsStart")) { MsStart <- -1; }
  if (!exists("BetaStart")) { BetaStart <- NULL; }
  if (!exists("MaxGibbsIters")) { MaxGibbsIters <- 1000; }
  if (!exists("FileSave")) { FileSave <- "MyS"; }
  if (!exists("SaveDir")) { SaveDir <- NULL; }
  if (!exists("dfTNoise")) { dfTNoise <- -1;}
  if (!exists("DoSave")) { DoSave = FALSE; }
  if (!exists("DoMax")) { DoMax = 1; }

  library(BayesDiallel);
   if (dfTauStart > 0) { tauPriordf = dfTauStart;}  else if (tauPriordf[1] > 0) {
      tauPriordf = tauPriordf;
   } else if (tauPriordf[1] < 0) {
     tauPriordf = mean(ListRList(AFD$df));
   }
   if (MsStart > 0) { tauPriorMean = MsStart;}  else if (tauPriorMean[1] > 0) {
      tauPriorMean = tauPriorMean;
   } else {
     tauPriorMean = mean(ListRList(AFD$m));
   }
   if (!exists("AFD") || is.null(AFD)) {
    print("Hey BayesSpikeFromAutoDiallel, we were not given AFD"); flush.console();
    return(-1);   
   }

   if (!exists("LengthChains") || is.null(LengthChains) || LengthChains < 0) { LengthChains = 3000; }
   if (exists("LengthChains") && !is.null(LengthChains) && LengthChains > 0) {
      MaxGibbsIters = LengthChains; NumSamples = LengthChains;
   }
   namesX = colnames(AFD$X);
   FixedVariableIds = (1:length(namesX))[namesX %in% AFD$.AllFixedVariables];
   StList = c();  EndList = c();   namesList = c();
   for (ii in 1:length(AFD$.AllRandomVariables)) {
     MyVar = (1:length(namesX))[ substr(namesX,1,nchar(AFD$.AllRandomVariables[ii])) ==
       AFD$.AllRandomVariables[ii]];
     if (length(MyVar) > 0) {
       StList = c(StList, min(MyVar)); EndList = c(EndList, max(MyVar));
       namesList = c(namesList, AFD$.AllRandomVariables[ii]);
     }
   }
   names(EndList) = namesList;
   tauEndList = EndList;
   FirstRandom = min(StList);
   dfTauStart = mean(ListRList(AFD$df));
   if (MsStart < 0) {
     MsStart = mean(ListRList(AFD$m));
   }
   if (dfTNoise > 0) {
   } else if (!is.null(AFD$dfTNoise) && AFD$dfTNoise > 0) {
     dfTNoise = AFD$dfTNoise;
   } else if (!is.null(AFD$AllDiallelObs[[1]]$dfTNoise) && 
     AFD$AllDiallelObs[[1]]$dfTNoise > 0) {
     dfTNoise = AFD$dfTNoise;
   }  else {
     dfTNoise = -1;
   }
   try(eval(parse(text=GetG0Text(".DefaultDiallelFileDir", "globalenv()", S=1))));
   if (is.null(SaveDir) && is.character(.DefaultDiallelFileDir) &&
     .DefaultDiallelFileDir != "") {
     try(SaveDir <- paste(.DefaultDiallelFileDir, "//BayesSpike", sep=""));
     flush.console();
     try(dir.create(SaveDir, showWarnings=FALSE, recursive=TRUE));
     try(eval(parse(text=LockGText(".DefaultDiallelFileDir", "globalenv()", S=1))));
   } else  if (is.null(SaveDir)) {
     DidFail = TRUE;
     .DefaultDiallelFileDir = "";
     for (ii in 1:length(.libPaths())) {
       RL = unlist(list.files(.libPaths()[[ii]]));
       if (any(RL == "BayesDiallelWorkData")) {
         print(paste("BayesDiallel:  Trying to find a directory to save ",
           "work data to! -- Going to BayesDiallelWorkData ", sep="")); 
         flush.console();
         try(dir.create(.libPaths()[[ii]], 
           "//BayesDiallelWorkData", showWarnings=FALSE));
         DidFail = FALSE;
         .DefaultDiallelFileDir = paste( .libPaths()[[ii]], "//", 
           "BayesDiallelWorkData//DiallelChainRecords//", sep="");
         EvalText = "dir.create(.DefaultDiallelFileDir, 
           showWarnings=TRUE);  DidFail = TRUE;";
         try(eval(parse(text=EvalText)));   ## Run Eval Text
         if (DidFail == FALSE) {
           print(paste("Well: we can't create to : ", 
             .DefaultDiallelFileDir, ", please supply BSSaveDir/SaveDir or we won't quite continue.", sep=""));
         }
          break;
       }
     }
     if (.DefaultDiallelFileDir == "") {
     for (ii in 1:length(.libPaths())) {
       RL = unlist(list.files(.libPaths()[[ii]]));
       if (any(RL == "BayesDiallel")) {
         print(paste("BayesDiallel:  Trying to find a directory to ",
           "save work data to! -- Going to BayesDiallel ", sep="")); 
         flush.console();
         dir.create(.libPaths()[[ii]], "//BayesDiallel", showWarnings=TRUE);
         .DefaultDiallelFileDir = paste( .libPaths()[[ii]], "//", 
           "BayesDiallel//DiallelChainRecords//", sep="");
         DidFail = FALSE;
         EvalText = paste( "dir.create(.DefaultDiallelFileDir, ",
           "showWarnings=TRUE); DidFail = TRUE;", sep="");
         try(eval(parse(text=EvalText)));
         if (DidFail == FALSE) {
           print(paste("Well: we can't create to : ", .DefaultDiallelFileDir, 
             ", please supply BSSaveDir/SaveDir or we won't quite continue.", 
             sep=""));
         }
         break;
       }
     }        
     }
     if (.DefaultDiallelFileDir == "") {
     try(eval(parse(text=GetG0Text(".DefaultDiallelFileDir", "globalenv()", S=1)))); 
     for (ii in 1:length(.libPaths())) {
       RL = unlist(list.files(.libPaths()[[ii]]));
       if (any(RL == "BayesDiallel")) {
         dir.create(.libPaths()[[ii]], "//BayesSpike", showWarnings=FALSE);
         try(.DefaultDiallelFileDir <- paste( .libPaths()[[ii]], "//", "BayesSpike//data//", sep=""));
          dir.create(.DefaultDiallelFileDir, showWarnings=FALSE);
          break;
       }
     }
     try(eval(parse(text=LockGText(".DefaultDiallelFileDir", "globalenv()", S=1))));         
     } 
     
      
     SaveDir = .DefaultDiallelFileDir;
     try(eval(parse(text=LockGText(".DefaultDiallelFileDir", "globalenv()", S=1))));  
   }
   
   if (!exists("DependenciesTau") || is.null(DependenciesTau)) {
     print("For the Diallel you should REALLY Set DependenciesTau"); 
     flush.console();
     DependenciesTau <- NULL;
   }
   if (!exists("Run")) { Run <- 1;}
   if (Run > 0  && (is.null(DoRecord) || length(DoRecord != 7) || sum(abs(DoRecord))==0)) {
     DoRecord = c(1,1,1,0,0,1,1)
   }
   
   FixedParam = colnames(AFD$X)[colnames(AFD$X) %in% AFD$AllFixedVariables]
   NoShrinkFixed = (1:length(AFD$X[1,]))[colnames(AFD$X) == "Mu"]
   NoShrinkFixedPrior = 1000;
   NoShrinkRandom <- NULL;  NoShrinkRandomPrior <- NULL;

   if (any(is.nan(AFD$Y))) {
     print("BayesSpike From AutoDiallel:  Error, nan Y!");
     print(AFD$Y);
     flush.console();
   }
   print("Copying AFD");
   CAFD=Copy(AFD, MyFullDiallelName =
     paste("Bayes Spike Copy: ", AFD$MyFullDiallelName, sep=""));
   CAFD$MyFullDiallelName =
     paste("Bayes Spike Copy: ", AFD$MyFullDiallelName, sep="");
  if (tauPriordf[1] < 0 || tauPriorMean[1] < 0) {
    tauPriordf = rep(2, length(tauEndList));
    tauPriorMean = rep(2, length(tauEndList));
    names(tauPriordf) = names(EndList);
    for (ii in 1:length(ListRList(AFD$m))) {
      if (names(ListRList(AFD$m))[ii] %in% names(tauPriorMean)) {
        tauPriorMean[names(tauPriorMean) == names(ListRList(AFD$m))[ii]] = 
          ListRList(AFD$m)[names(ListRList(AFD$m))[ii] == names(ListRList(AFD$m))]
      }
    }
    for (ii in 1:length(ListRList(AFD$df))) {
      if (names(ListRList(AFD$df))[ii] %in% names(tauPriordf)) {
        tauPriordf[names(tauPriordf) == names(ListRList(AFD$df))[ii]] = 
          ListRList(AFD$df)[names(ListRList(AFD$df))[ii] == names(ListRList(AFD$df))]
      }
    }
  }

  ## This is a big fight to automatrically find PriorProbTau which are 
  ##   priors for groups.
  if (is.null(PriorProbTau)) {
    PriorProbTau <- NULL;
  } else if (length(PriorProbTau) != length(tauEndList)) {
     if (is.null(names(PriorProbTau))) {
       if (length(PriorProbTau) < length(tauEndList)) {
         try(PriorProbTau <- c(PriorProbTau, rep(PriorProbTau[length(PriorProbTau)],
           length(EndList) - length(PriorProbTau)) ));
       } else {
         try(PriorProbTau <- PriorProbTau[1:length(tauEndList)]);
       } 
       try(names(PriorProbTau) <- AFD$.AllRandomVariables[1:length(tauEndList)], silent=TRUE);  
     } else {
       MyN = AFD$.AllRandomVariables;
       if (length(MyN) != length(tauEndList)) {
         ##print("-------------Curious-----------------------------------------");
         ##print(paste("Length of tauEndList = ", length(tauEndList),
         ##  "  and length(MyN) = ", length(MyN), sep=""));
         ##print(paste("names MyN = (", 
         ##  paste(MyN, collapse=", "), ")", sep=""));
         ##print(paste("names tauEndList = (", paste(names(tauEndList),
         ##  collapse=", "), ")", sep="")); flush.console();
         if (is.null(names(tauEndList))) {
           print("FF!!, no names to tauEndList!");
           return(AFD);
         } else {
           MyN = names(tauEndList);
         }
       }
       APTM <- "
       NTM = match(MyN, names(PriorProbTau));
       NewPriorProbTau = rep(0, length(tauEndList));
       NewPriorProbTau[(1:length(NTM))[!is.na(NTM)]] = 
         PriorProbTau[NTM[!is.na(NTM)]];
       NewPriorProbTau[NewPriorProbTau == 0.0] = mean(PriorProbTau);
       PriorProbTau = NewPriorProbTau;
       ";
       try(eval(parse(text=APTM)));
     }
   }
   if (!exists("PriorProbFixed") || is.null(PriorProbFixed) ||
     length(PriorProbFixed) <= 0) {
     APTM <- "
     if (FirstRandom > 1) {
       PriorProbFixed = rep(-666, FirstRandom-1);   
     } else if (FirstRandom == 1) {
       PriorProbFixed = NULL;
     } else {
       PriorProbFixed = rep(-666, NCOL(AFD$X));   
     }";
     try(eval(parse(text=APTM)));
   }
   if (!is.null(AFD$DoFirstCenter) && AFD$DoFirstCenter == TRUE && length(AFD$AllCenteredRandomVariables) >= 0 && length(tauEndList) >= 1) {
     print("Trying to identify the AllCenteredRandomVariables Set!");
     AllCenteredRandomVariables = AFD$AllCenteredRandomVariables;
     WhatRandoms <- rep(0, length(tauEndList));
     if (FirstRandom >= 1 && FirstRandom <= NCOL(AFD$X)) {
       ABS <- FirstRandom:tauEndList[1];
     }  else {
       ABS <- 1:tauEndList[1];
     }
     AST <- colnames(AFD$X)[ABS[1]];
     for (gg in 1:length(AllCenteredRandomVariables)) {
       if (substr(AST,1, nchar(AllCenteredRandomVariables[gg])) ==  AllCenteredRandomVariables[gg]) {
         try(WhatRandoms[1] <- 1); 
         break;
       }
     }
     if (length(tauEndList) >= 2) {
     for (iti in 2:length(tauEndList)) {
       ABS <- (tauEndList[iti-1]+1):tauEndList[iti];
       try(AST <- colnames(AFD$X)[ABS[1]]);
       for (gg in 1:length(AllCenteredRandomVariables)) {
         if (substr(AST, 1, nchar(AllCenteredRandomVariables[gg])) == AllCenteredRandomVariables[gg]) {
           WhatRandoms[iti] = iti;
           break;
         }
       }
     }
     }
     CenteredColumns <- WhatRandoms[WhatRandoms >= 1];
   } else {
     CenteredColumns = NULL;
   }
   tauStartValues <- rep(1, length(tauEndList));
   try(names(tauStartValues) <- names(tauEndList));
   
    nX <- colnames(AFD$X);
    NKeep <- rep(0, length(AFD$.AllRandomVariables));
    for (ij in 1:length(NKeep)) {
      if (any(substr(nX,1,nchar(AFD$.AllRandomVariables[ij])) == 
         AFD$.AllRandomVariables[ij])) {
        NKeep[ij] = 1;  
      }
    }
    if (DoRecord[1] == 1) {
      CodaTableNames <- nX
    } else {  CodaTableNames <- NULL; }
    if (DoRecord[2] == 1 && sum(NKeep) >= 1) {
      CodaTableNames <- c(CodaTableNames, paste("tau:",
        AFD$.AllRandomVariables[NKeep==1], sep=""));
    }
    if (DoRecord[3] == 1) {
      CodaTableNames <- c(CodaTableNames, "Sigma:1");
    }
    if (DoRecord[4] == 1) {
      if (length(PiAStart) == 2 && sum(NKeep) >= 1) {
        CodaTableNames <- c(CodaTableNames, "PiA:1", "PiA:2");
      } else {
        CodaTableNames <- c(CodaTableNames, "PiA:1");
      }
    }
    if (DoRecord[5] == 1 && length(tauFixed) >= 1) {
      CodaTableNames <- c(CodaTableNames, 
        paste("tauFixed:", 1:length(tauFixed), sep=""));
    }
    if (DoRecord[6] == 1 && FirstRandom > 1) {
      CodaTableNames <- c(CodaTableNames,
        paste("ProbFixed:", nX[1:(FirstRandom-1)], sep="")
      );  
    }
    if (DoRecord[7] == 1 && sum(NKeep) >= 1) {
      CodaTableNames <- c(CodaTableNames,
        paste("ProbTau:", AFD$.AllRandomVariables[NKeep==1], sep="")
      );  
    }

   ##try(names(tauPriordf) <- names(tauEndList));
   ##try(names(tauPriorMean) <- names(tauEndList));
   MyBMS = BayesSpikeRegression(Y=AFD$Y,X=AFD$X,
     IndexFirstRandomEffect = FirstRandom, tauEndList=tauEndList,
     tauStartValues = tauStartValues, 
     PiAStart = PiAStart, PiAPrior = PiAPrior,
     tauPriordf = tauPriordf, tauPriorMean = tauPriorMean,
     dfTNoise = dfTNoise, Verbose = Verbose,
     DependenciesTau = DependenciesTau, HowSample =HowSample,
     BetaStart = BetaStart, DependenciesFixed=DependenciesFixed,
     DoMax = DoMax,  NumSpliceSampleReps = NumSpliceSampleReps, 
     SpliceMove= SpliceMove, CauchyEpsilon = CauchyEpsilon,
     MyPS = MyPS,  tauFixed = tauFixed,
     TypeFixedPrior = TypeFixedPrior, DoRNGState = DoRNGState,
     SigmaPrior = SigmaPrior, InitKKs = InitKKs, DoRecord = DoRecord,
     NumSamples = NumSamples, MaxIters = MaxIters, 
     MaxGibbsIters = MaxGibbsIters,
     MaximizeMeCauchyTotalIters = MaximizeMeCauchyTotalIters,
     ttStart = ttStart,  NumChains = NumChains,
     EarlyEndStep = -2, EarlyEndtt = -2, DoSave=DoSave,
     FileSave = FileSave,
     SaveDir = SaveDir,  FileName = FileName,
     NewWrite = NewWrite, 
     Run = Run, CodaTableNames = CodaTableNames,
     AFD = CAFD, NoShrinkFixed=NoShrinkFixed, NoShrinkFixedPrior = NoShrinkFixedPrior,
     NoShrinkRandom=NoShrinkRandom, NoShrinkRandomPrior=NoShrinkRandomPrior,
     PriorProbTau = unlist(PriorProbTau), PriorProbFixed=unlist(PriorProbFixed),
     ZeroOutBeforeSample = ZeroOutBeforeSample);

   if (TestEndBayesSpikeAutoDiallel == 1) {
     print("********************************************************************");
     print("*** TestEndBayesSpikeAutoDiallel Is 1 so we quit. "); flush.console();
     print("*** BayesSpikeAutoDiallel Early Exit 1 "); flush.console();
     print("********************************************************************");
     flush.console();
     eval(parse(text=SetGText("MyBMS", "globalenv()", S=1)));
     eval(parse(text=SetGText("CenteredColumns", "globalenv()", S=1)));
     return(MyBMS);
   }
   ##Link in elements of CAFD
   if (is.null(names(MyBMS$Beta))) {
     if (Verbose >= 2) {
       print(paste("BayesSpikeAutoDiallel: We have to set names for Beta", sep=""));
       flush.console();
     }
     try(MyBMS$SetNamesBeta(colnames(MyBMS$X)));
     if (Verbose >= 1) {
        print(paste("BayesSpikeAutoDiallel: Successfully Inserted Beta", sep="")); flush.console();
     }
   }
   CAFD$.MBS=MyBMS;
   CAFD$.TBSR5=MyBMS$TBSR5
   MBS <- MyBMS;
  if (is.null(MBS$tauEndList) || length(MBS$tauEndList) <= 0) {
    LenFixed <- MBS$p;
  } else {
    LenFixed <- MBS$CFirstRandom;
  } 
 
   if (Verbose >= 1) {
      print("BayesSpikeAutoDiallel: Checking about CenteredColumns"); flush.console();
   }  
   if (!is.null(CenteredColumns)  && length(CenteredColumns) >= 0) {
     try(MBS$CenteredColumns <- CenteredColumns -1);
   }
   if (TestEndBayesSpikeAutoDiallel == 2) {
     print("TestEndBayesSpikeAutoDiallel Is 2 so we quit after Center Columns "); flush.console();
     eval(parse(text=SetGText("MyBMS", "globalenv()", S=1)));
     eval(parse(text=SetGText("CenteredColumns", "globalenv()", S=1)));
     return(MyBMS);
   }
   if (!is.null(MyBMS$OtherNameCodaList)) {
     if (Verbose >= 1) {
        print("BayesSpikeAutoDiallel: We are checking OtherNameCodaList"); flush.console();
     }
     ART <- list(NULL, paste(MyBMS$OtherNameCodaList, sep=""));
     try(MyBMS$RenameCodaList(ART));
   } else {
     if (Verbose >= 1) {
        print("BayesSpikeAutoDiallel: We will have to conduct RenameBSAFDCoda"); flush.console();
     }
     try(BayesSpike:::RenameBSAFDCoda(BSAFD=MyBMS, 
       TrueParametersList=NULL ));  
   }
   if (Verbose >= 1) {
     print("BayesSpikeAutoDiallel: About to check OtherNameCodaList"); flush.console();
   }
   if (!is.null(MyBMS$OtherNameCodaList) &&
     is.null(colnames(MyBMS$CodaList[[1]])) &&
     NCOL(MyBMS$CodaList[[1]]) == length(MyBMS$OtherNameCodaList)) {
     NewList <- list();
     if (Verbose >= 1) {
        print("BayesSpikeAutoDiallel: About to Fill OtherNameCodaList");
        flush.console();
     }
     for (ii in 1:length(MyBMS$CodaList)) {
       try(AOb <- as.matrix(MyBMS$CodaList[[ii]])); 
       try(colnames(AOb) <- MyBMS$OtherNameCodaList);
       try(NewList[[ii]] <- as.mcmc(AOb));
     } 
     try(NewList <- as.mcmc.list(NewList));
     try(MyBMS$CodaList <- NewList);
   }
   MyTry = NULL;
   if (Verbose >= 1) {
     print("BayesSpikeAutoDiallel: About to Run PriorProbFixed"); flush.console();
   }
   if (!is.null(PriorProbFixed) && length(PriorProbFixed) >= 1) {
     if (MyBMS$RawiFirstRandom > 0 || (is.null(tauEndList) || length(tauEndList)  <= 0)) {
       MyBMS$RawiFirstRandom;
     }  else {
       LenFixed = MyBMS$p;
     }
     if (length(PriorProbFixed) == LenFixed) {
     } else if (length(PriorProbFixed) == 2 * LenFixed) {
     } else if (length(PriorProbFixed) > 2*LenFixed) {
       PriorProbFixed <- PriorProbFixed[1:(2*LenFixed)];
     } else if (length(PriorProbFixed) > LenFixed) {
        PriorProbFixed <- PriorProbFixed[1:LenFixed];
     } else if (!is.null(names(PriorProbFixed))) {
        NN <- colnames(MyBMS$X);
        if (is.null(NN)) {
          NN <- names(MyBMS$Beta);
        }
        if (!is.null(NN)) {
          mMM <- match(names(PriorProbFixed), colnames(NN));
          fMM <- match(names(PriorProbFixed), paste(1:length(tauEndList), sep=""));
          mMM[is.na(mMM) && !is.na(fMM)] <- fMM[is.na(mMM) && !is.na(fMM)];
          nPriorProbFixed <- rep(-1, LenFixed);
          nPriorProbFixed[mMM[!is.na(mMM)]] <- PriorProbFixed[!is.na(mMM)];
          PriorProbFixed <- unlist(nPriorProbFixed);
        } else {
          PriorProbFixed <- NULL;
        }
     } else if (length(PriorProbFixed) == 1) {
        PriorProbFixed <- rep(PriorProbFixed, LenFixed);
     }
     try(MyBMS$PriorProbFixed <- unlist(PriorProbFixed));
   } else {
     PriorProbFixed <- NULL;
     try(MyBMS$PriorProbFixed <- NULL);
   }
   if (Verbose >= 1) {
     print("BayeSpikeFromAutoDiallel: About to run PriorProbTau"); flush.console();
   }
   if (is.null(PriorProbTau) || length(PriorProbTau) <= 0) {
     PriorProbTau <- NULL;
     try(MyBMS$PriorProbTau <- unlist(PriorProbTau));
   } else if (!is.null(PriorProbTau)) {
     if (length(PriorProbTau) == length(tauEndList)) {
     } else if (length(PriorProbTau) == 2 * length(tauEndList)) {     
     } else if (!is.null(names(PriorProbTau))) {
       NN <- names(tauEndList);
       if (is.null(NN)) {
         NN <- names(tauStartValues);
       }
       if (!is.null(NN)) {
         mMM <- match(names(PriorProbTau), NN);
         fMM <- match(names(PriorProbTau), paste(1:length(tauEndList), ""));
         tMM <- match(names(PriorProbTau), paste("tau:", 1:length(tauEndList), ""));
         TTMM <- match(names(PriorProbTau), paste("tau:", NN, "")); 
         mMM[is.na(mMM) && !is.na(fMM)] <- fMM[is.na(mMM) && !is.na(fMM)];
         mMM[is.na(mMM) && !is.na(tMM)] <- tMM[is.na(mMM) && !is.na(tMM)];
         mMM[is.na(mMM) && !is.na(TTMM)] <- TTMM[is.na(mMM) && !is.na(TTMM)];
         nPriorProbTau <- rep(-1, length(tauEndList));
         nPriorProbTau[mMM[!is.na(mMM)]] <- PriorProbTau[!is.na(mMM)];
         PriorProbTau <- unlist(nPriorProbTau);
       } else {
         print("BayesSpikeRegression: PriorProbTau we forced to set NULL!"); flush.console();
         PriorProbTau <- NULL;
       }
     }
     try(MyBMS$PriorProbTau <- unlist(PriorProbTau));
   }

   if (Verbose >= 1) {
     print("BayesSpikeFromAutoDiallel: Now going into WatchCoordinates"); flush.console();
   }   
   if (Run> 0 && FALSE) { 
   if (exists("WatchCoordinates") && 
     is.null(WatchCoordinates) || length(WatchCoordinates) == 0) {
     for (ii in 1:length(MyBMS$CodaList)) {
       if (!is.null(MyBMS$CodaList[[ii]])) {
         CD = length(MyBMS$CodaList[[ii]][1,]);
         CM = length(AFD$AllDiallelObs[[1]]$CodaChains[[ii]][1,]);
         CDL = length(MyBMS$CodaList[[ii]][,1]);
         CML = length(AFD$AllDiallelObs[[1]]$CodaChains[[ii]][,1]);
         MyCM = min(c(CD, CM));
         MyCML = min(c(CDL, CML));
         AMatchNames = match(colnames(MyBMS$CodaList[[ii]]),
           colnames(AFD$AllDiallelObs[[1]]$CodaChains[[ii]]));
         match123 = (1:length(AMatchNames))[!is.na(AMatchNames)];
         match456 = AMatchNames[!is.na(AMatchNames)];
         AFD$AllDiallelObs[[1]]$CodaChains[[ii]][1:MyCML, match123] = 
           MyBMS$CodaList[[ii]][1:MyCML, match456];    
       } 
   
     }
   }  else {
     for (ii in 1:MyBMS$TBSR5$NumChains) {
        MyBMS$TBSR5$SetupSaveFile(MyBMS$TBSR5$SaveDir, 
           FileName = MyBMS$TBSR5$FileName, chainiter = ii);
        RTM = .Call("GiveCodaSubset", WatchCoordinates, 
          paste(MyBMS$sCodaIFile, sep=""), 
          paste(MyBMS$sCodaDFile, sep=""), 0, MaxGibbsIters, 0)
        if (!is.null(RTM) && length(dim(RTM)) == 2 && 
          exists("WatchCoordinates") && 
          dim(RTM)[2] == length(WatchCoordinates)) {
         CD = length(RTM[1,]);
         CM = length(AFD$AllDiallelObs[[1]]$CodaChains[[ii]][1,]);
         CDL = length(RTM[,1]);
         CML = length(AFD$AllDiallelObs[[1]]$CodaChains[[ii]][,1]);
         MyCM = min(c(CD, CM));
         MyCML = min(c(CDL, CML));
         AFD$AllDiallelObs[[1]]$CodaChains[[ii]][1:MyCML, 1:MyCM] = 
           RTM[1:MyCML, 1:MyCM];    
       }  else {
         print("Somehow, not all of the WatchCoordinates were recorded!");
         flush.console();  return(MyBMS);
       }
     }
   }}
  if (Verbose >= 1) {
    print("BayesSpikeFromAutoDiallel: Rename CodaList again"); flush.console();
  }
  if (!is.null(MyBMS$OtherNameCodaList) &&
     is.null(colnames(MyBMS$CodaList[[1]])) &&
     NCOL(MyBMS$CodaList[[1]]) == length(MyBMS$OtherNameCodaList)) {
     NewList <- list();
     for (ii in 1:length(MyBMS$CodaList)) {
       try(AOb <- as.matrix(MyBMS$CodaList[[ii]])); 
       try(colnames(AOb) <- MyBMS$OtherNameCodaList);
       try(NewList[[ii]] <- as.mcmc(AOb));
     } 
     try(NewList <- as.mcmc.list(NewList));
     try(MyBMS$CodaList <- NewList);
  }
  if (Verbose >= 1) {
    print("BayesSpikeFromAutoDiallel: About to return MyBMS");
    flush.console();
  }
  return(MyBMS); 
}


################################################################################
##  BayesSpikeAutoDiallelPosterior(BSAFD, start...)
##
##    R package "BayesDiallel", a companion package designed for fitting
##  the genetic experiment of the diallel was an inspiration for this package.
##  "AFD" is the full Diallel Object and is attached to BayesSpikeCL BSAFD object.
##  This tried to find the largest model from a list of Models in AFD.
##  Then the design matrix for this largest model is constructed and the
##  BayesSpikeRegression algorithm is run.
BayesSpikeAutoDiallelPosterior <- function(BSAFD, start =1, end = -1, thin =1,
  WithNoise = FALSE) {
   if (end < start) {
      end = length(BSAFD$CodaList[[1]][,1]);
   }
   if (is.null(BSAFD$AFD)) {
    print("BayesSpikeAutoDiallelPosterior: Sorry not with null AFD");
    return(-1);   
   }
   MaxOverAll = FindMaximizerTable(BSAFD$AFD$ModelsList, Verbose = BSAFD$Verbose);
   CrossModelsList = sort(unique(BSAFD$AFD$ModelsList[,3]));
   numj=BSAFD$AFD$.numj;
   FakeListjk = cbind(rep(1:numj,numj), rep(1:numj, each=numj));
   if (MaxOverAll[4] > 0) {
     GenderOn = 0;
     FakeSexVector =  c(rep(0, numj^2), rep(1,,numj^2));
     FakeListjk = rbind(FakeListjk, FakeListjk);
   } else {
     FakeSexVector = NULL;
   }

   FakeX = ConstructXJK(numj, ajModel = MaxOverAll[1],
    MotherModel = MaxOverAll[2],
    CrossModelsList = CrossModelsList,
    CrossModel=MaxOverAll[3], SexModel=MaxOverAll[4], 
    BetaHybridLevel=MaxOverAll[5], 
    BetaHybridTimesSex=MaxOverAll[6], SexVector = FakeSexVector, 
    Listjk=FakeListjk, 
    CO = BSAFD$AFD$ACrossLocations, FixedEffects=BSAFD$AFD$FixedEffects,
    RandomEffectsGroups = BSAFD$AFD$RandomEffectsGroups,
    RandomEffects=BSAFD$AFD$RandomEffects,
    verbose = BSAFD$Verbose, ALLINCLUSIVE = TRUE);
    nX <- colnames(FakeX);
    NKeep <- rep(0, length(BSAFD$AFD$.AllRandomVariables));
    for (ij in 1:length(NKeep)) {
      if (any(substr(nX,1,nchar(BSAFD$AFD$.AllRandomVariables[ij])) == 
         BSAFD$AFD$.AllRandomVariables[ij])) {
        NKeep[ij] = 1;  
      }
    }
    if (sum(NKeep) >= 1) {
      NamesAFD <- c(colnames(FakeX), paste("tau:", 
        BSAFD$AFD$.AllRandomVariables[NKeep==1], sep=""), "Sigma:1");
    } else {
      NamesAFD <- c(colnames(FakeX), "Sigma:1")        
    }    

  if (any(substr(colnames(BSAFD$CodaList[[1]]),1,nchar("Beta:")) == "Beta:")) {
    TrueInVector = rep(0, length(NamesAFD));
    names(TrueInVector) = NamesAFD;
    if (!exists("TrueInVector")) {
      print("Weird right before RenameBSAFDCoda: TrueInVector does not exist!");
      flush.console();
    }
    try(BayesSpike:::RenameBSAFDCoda(BSAFD=BSAFD, TrueParametersList=TrueInVector ));  
  }
  FakeCodaList = list();
  for (ii in 1:length(BSAFD$CodaList)) {
    NewFakeCodaOb = 0;
    if (end < start) { end = length(BSAFD$CodaList[[ii]][,1]); }
    idd = match(colnames(BSAFD$CodaList[[ii]]), colnames(FakeX));
    Cidd = (1:length(colnames(BSAFD$CodaList[[ii]])))[!is.na(idd)];
    NewFakeCodaOb = BSAFD$CodaList[[ii]][start:end,Cidd] %*% t(FakeX[,idd[!is.na(idd)]]);
    if (is.null(SexVector)) {
      colnames(NewFakeCodaOb) = paste("i:", FakeListjk[,1], ";k:", FakeListjk[,2], sep=""); 
    } else {
      colnames(NewFakeCodaOb) = paste("S:", FakeSexVector, ";i:", 
        FakeListjk[,1], ";k:", FakeListjk[,2], sep="");
    }
    if (WithNoise == TRUE) {
      dfTNoise = -1;  try(dfTNoise <- BSAFD$dfTNoise);
      SigmaId = (1:length(colnames(BSAFD$CodaList[[ii]])))[
        substr(colnames(BSAFD$CodaList[[ii]]), 1, nchar("Sigma")) ==
        "Sigma"
      ];
      if (dfTNoise <= 0) {
        Noise = matrix(rep(sqrt(BSAFD$CodaList[[ii]][start:end,SigmaId]), length(FakeX[,1])),
           length(start:end), length(FakeX[,1]));
        Noise = Noise * matrix(rnorm(dim(Noise)[1] * dim(Noise)[2],0,1),
          dim(Noise)[1], dim(Noise)[2]);
      } else {
         Noise = matrix(rep(sqrt(BSAFD$CodaList[[ii]][start:end,SigmaId]), length(FakeX[,1])),
           length(start:end), length(FakeX[,1]));
         Noise = Noise * matrix(rt(dim(Noise)[1] * dim(Noise)[2],dfTNoise),
          dim(Noise)[1], dim(Noise)[2]);     
      }
      NewFakeCodaOb = NewFakeCodaOb + Noise;
    }
    NewFakeCodaOb = as.mcmc(NewFakeCodaOb);
    FakeCodaList[[ii]] = NewFakeCodaOb;
  }  
  try(FakeCodaList <- as.mcmc(FakeCodaList));
  try(FakeCodaList <- as.mcmc.list(FakeCodaList));
  ReturnFakeSet = list(FakeListjk= FakeListjk, FakeSexVector=FakeSexVector,
    FakeX = FakeX, FakeCoda=FakeCodaList);
  return(ReturnFakeSet);
}

###################################################################################
## rTruncT: Random Draws from Truncated T distibution
##
##   A function for testing function "rTT" which will try to draw from
## the truncated t distribution using truncated dimension.
##
rTruncT <- function(n=1, df = 1, L = NULL, U = NULL, sig=NULL, mu=NULL) {
  if (is.null(L) && is.null(U)) {
    print("rTruncT: Please supply actual bounds either upper or lower!");
    return(-1);
  }
  if (length(L) > 1 && length(U) > 1 && length(L) != length(U)) {
    print("rTruncT: Please supply L and U of equivalent lengths!");
  } else if (length(L) > 1 && length(U) == 1) {
    n = length(L);  U = rep(U[1], n);
  } else if (length(U) > 1 && length(L) == 1) {
    n = length(U); L = rep(L[1], n);
  }
  
 sRet = rep(0, n); smu = rep(0, 100)
 .Call("rTT", sRet, df, mu, sig, L, U);
 return(sRet);
}

################################################################################
##  .AFunctionNewUnCenteredCodaList (Hidden Function of a BayesSpikeCL Object MBS)
##
##      Sometimes group coordinates belong to a fixed equation 
##                        sum(j in group k) a_j = 0
##
##  As a result we choose to work on a "centered version" which rotates into a 
##   set of n_k-1 which have the same prior variance as the original coefficients.
##  This will reconvert the CodaList (which is precentered) to an UnCenteredCodaList.
##   and relabel correctly.  OtherNameCodaList then stores the desired renames.   
.AFunctionNewUnCenteredCodaList <- function(MBS) {
  if (is.null(MBS) || is.null(MBS$tauEndList)) { 
    print("Error: In UnCenteredCodaList: the observed MBS is null or does not have tauEndList!");
    flush.console();
    return(-1);
  }
  if (is.null(MBS$CenteredColumns) || length(MBS$CenteredColumns) <= 0)  { 
     print("Error: In UnCenteredCodaList: the CenteredColumns is null or zero length!");
     return(-1);
  }
  if (is.null(MBS$DoRecord) || length(MBS$DoRecord) < 1 ||
    MBS$DoRecord[1] != 1) {
    print("Error: in UnCenteredCodaList: the DoRecord did not record coefficients!");
    flush.console();
    return(-1);
  }
  if (is.null(MBS$CodaList) || length(MBS$CodaList) <= 0) {
    print("Error in UncenteredCodaList: the CodaList is not recorded!");
    flush.console();
    return(-1);
  }
  if (MBS$Verbose >= 1) {
    print("AFunctionNewUnCenteredCodaList: begin Assembly"); flush.console();
  }
  try(MBS$CenteredColumns <- sort(MBS$CenteredColumns));
  CenteredIndices <- list();
  NewCenteredNames <- list();
  ABAll <- 1:NCOL(MBS$CodaList[[1]]);
  ColsNames <- colnames(MBS$CodaList[[1]]);
  if (!is.null(MBS$OtherNameCodaList) && length(MBS$OtherNameCodaList) == NCOL(MBS$CodaList[[1]]) ) {
    ColsNames <- MBS$OtherNameCodaList
  }
  if (length(ColsNames) <= 0) { ColsNames<-NULL; }
  for (ii in 1:length(MBS$CenteredColumns)) {
    if (MBS$CenteredColumns[ii] == 0) {
      if (MBS$FirstRandom >= 1 && MBS$FirstRandom < MBS$p) {
        CenteredIndices[[ii]] = MBS$FirstRandom:(MBS$tauEndList[ii]+1);
      } else {
        CenteredIndices[[ii]] = 1:(MBS$tauEndList[ii]+1);
      }
    } else {
      CenteredIndices[[ii]] = (MBS$tauEndList[ii-1]+2):(MBS$tauEndList[ii]+1); 
    }
    if (length(CenteredIndices[[ii]]) >= 1) {
      ABAll = ABAll[!(ABAll %in% CenteredIndices[[ii]])];
    }
    if (length(CenteredIndices[[ii]]) >= 1 && !is.null(ColsNames)) {
      ABS <- strsplit(ColsNames[CenteredIndices[[ii]]], ":");
      AB1 <- -1;   VD = -1;
      TT <- 0;
      try(TT <- as.numeric(ABS[[length(ABS)]][length(ABS[[length(ABS)]])]), silent=TRUE); 
      if (!is.na(TT) && !is.null(TT) && length(TT) >= 1) {
        if (TT > AB1) {
          VD = length(ABS);  AB1 = TT;
        }
        NewCenteredNames[[ii]] <- c(ColsNames[CenteredIndices[[ii]]],
          paste(paste(ABS[[VD]][1:(length(ABS[[VD]])-1)], collapse=":"),
            ":", AB1+1, sep="") );
      }
    }
  }
  if (MBS$Verbose >= 1) {
    print("AFunctionNewUnCenteredCodaList:Finished GatheringCenteredIndices"); flush.console();
  }
  if (length(CenteredIndices) <= 1) {
    InBetweenLength = 0;
    InBetweenList = 0;
  } else {
    InBetweenLength <- rep(0, length(CenteredIndices)-1);
    InBetweenList <- list();
    for (ii in 2:length(CenteredIndices)) {
      if (max(CenteredIndices[[ii-1]])+1 == min(CenteredIndices[[ii]]) ) {
        InBetweenList[[ii-1]] <- -1;
        InBetweenLength[ii-1] = 0;
      } else {
        InBetweenList[[ii-1]] <- (max(CenteredIndices[[ii-1]])+1):(min(CenteredIndices[[ii]])-1); 
        InBetweenLength[ii-1] <- length(InBetweenList);
      }
    }
  }
  if (MBS$Verbose >= 2) {
    print("AFunctionNewUnCenteredCodaList: Finished Gathering InBetweenList"); flush.console();
  }  
  
MakeAM <- function(n) {
 ## A <- (n-2)^2+(n-2);
 ## B <- 2*(n-2)/ sqrt(n-1);
 ## C <- 1 - 1/(n-1);
 ## k <- (-B+ sqrt(B^2+ 4 * A * C)) / ( 2 * A) 
 k <- (-1 + sqrt(n)) / (n-1)^(3/2) 
 AMM <- matrix(-k, n-1,n-1);
 diag(AMM) <- k * (n-2) + 1 / sqrt(n-1);
 AM <- rbind(AMM, rep(-1/sqrt(n-1), n-1));
 return(AM);
}
  MinCentered = min(unlist(CenteredIndices));
  MaxCentered = max(unlist(CenteredIndices));
  ABLeft = ABAll[ABAll < MinCentered];
  ABRight = ABAll[ABAll > MaxCentered];
  RC <- rep(0, length(CenteredIndices));
  for (ii in 1:length(CenteredIndices)) {
    RC[ii] <- length(CenteredIndices[[ii]]);
  }
  if (MBS$Verbose >= 1) {
    print("AFunctionNewUnCenteredCodaList:NewMConstruction"); flush.console();
  }
  NewL <- list();
  for (tt in 1:length(MBS$CodaList)) {
    AP = 0;
    AGC <- c();
    try(NewM <- matrix(0, NROW(MBS$CodaList[[tt]]),
      length(ABLeft) + length(ABRight) + sum(InBetweenLength) + sum(RC) + length(RC)));
    if (length(ABLeft) >= 0) {
      try(NewM[,1:length(ABLeft)] <- MBS$CodaList[[tt]][, ABLeft]);
      try(AP <- AP + length(ABLeft));
      if (!is.null(ColsNames)) {
        AGC <- ColsNames[ABLeft];
      }
    } else {
      AP <- 0;  AGC <- c();
    }
    for (ii in 1:length(CenteredIndices)) {
      try(AM <- MakeAM(1+length(CenteredIndices[[ii]])));
      try(NewM[,AP+1:(length(CenteredIndices[[ii]])+1)] <-
        MBS$CodaList[[tt]][,CenteredIndices[[ii]] ] %*% t(AM));
      try(AP <- AP + NROW(AM));
      if (!is.null(ColsNames)) {
        try(AGC <- c(AGC, NewCenteredNames[[ii]]));
      }
      if (length(InBetweenLength) >= ii && InBetweenLength[ii] >= 1) {
        try(NewM[, AP+ 1:InBetweenLength[ii]] <-
          MBS$CodaList[[tt]][,InBetweenList[[ii]]]);
        if (!is.null(ColsNames)) {
          try(AGC <- c(AGC, ColsNames[InBetweenList[[ii]]]));
        }
        try(AP <- AP + InBetweenLength[ii]);
      }
    }
    if (length(ABRight) >= 0) {
      try(NewM[,AP+1:length(ABRight)] <- MBS$CodaList[[tt]][, ABRight]);
      if (!is.null(ColsNames[[1]])) {
        try(AGC <- c(AGC, ColsNames[ABRight]));
      }
    }
    if (!is.null(ColsNames) && length(AGC) >= 0 && !is.null(AGC)) {
      try(colnames(NewM) <- AGC);
    }
    try(NewM <- as.mcmc(NewM));
    try(NewL[[tt]] <- NewM);
  }
  try(NewL <- as.mcmc.list(NewL));
  try(MBS$DeCenteredCodaList <- NewL);
  try(MBS$OldDeCenterCodaNames <- colnames(NewL[[1]]));
  return(NewL);

}
