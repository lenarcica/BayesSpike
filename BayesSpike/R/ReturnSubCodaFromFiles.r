## ReturnSubCodaFromFiles;
ReturnSubCodaFromFiles <- function(
    TBSR5 = NULL, TBSCpp = NULL, StartIter = 0, EndIter = 2000,
    Verbose = 0, SubSetCoords = NULL,DoAllTemp = -1, TBSRoo = NULL,DoLogitPostProb = FALSE, DoProb = -1, 
    DoPostBuff = -1, DoILoc = -1, ...) {
  if (!exists("TBSRoo")) { TBSRoo <- NULL; }
  if (!exists("TBSR5")) {  if (exists(".self")) {
    TBSR5 <- .self;
    } else if (exists("GibbsSS")) {
      TBSR5 <- GibbsSS$TBSR5;
    } else if (exists("TBSCpp")) { TBSR5 <- TBSCpp; }
  }
  if (!exists("TBSCpp")) { if (exists("GibbsSS")) { TBSCpp <- GibbsSS; 
    } else if (exists("TBSR5")) {
      TBSCpp <- TBSR5$ABayesSpikeCL;
    }
  }
  if (is.null(TBSR5) && !is.null(TBSRoo)) {
    TBSR5 = TBSRoo;
  }
  if (is.numeric(DoProb) && length(DoProb) >= 1 && DoProb[1] == -1) {
    if (TBSR5$p >= 50000) {
      DoProb <- FALSE;
    } else {
      DoProb <- TRUE;
    }
  } else if (is.numeric(DoProb) && length(DoProb) >= 1 && DoProb[1] >= .5) {
    DoProb <- TRUE;
  } else {
    DoProb <- FALSE;
  }
  if (is.numeric(DoProb) && length(DoPostBuff) >= 1 && DoPostBuff[1] == -1) {
    if (TBSR5$p >= 50000) {
      DoPostBuff <- FALSE;
    } else {
      DoPostBuff <- TRUE;
    }
  } else if (is.numeric(DoPostBuff) && length(DoPostBuff) >= 1 && DoPostBuff[1] >= .5) {
    DoPostBuff <- TRUE;
  } else {
    DoPostBuff <- FALSE;
  }
  if (is.numeric(DoILoc) && length(DoILoc) >= 1 && DoILoc[1] == -1) {
    if (TBSR5$p >= 50000) {
      DoILoc <- FALSE;
    } else {
      DoILoc <- TRUE;
    }
  } else if (is.numeric(DoILoc) && length(DoILoc) >= 1 && DoILoc[1] >= .5) {
    DoILoc <- TRUE;
  } else {
    DoILoc <- FALSE;
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
    print("ReturnSubCodaFromFiles: Verbose is Null"); flush.console(); Verbose = 3;
  }
  if (Verbose > 0) {
    print("ReturnSubCodaFromFiles: Starting"); flush.console();
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
    print("ReturnSubCodaFromFiles:  Error Please give non null TBSRoo");
    return(NULL);
  }
  if (!exists("SubSetCoords") || is.null(SubSetCoords)) {
    SubSetCoords = 1:TBSRoo$p;
  }
  if (Verbose > 0) {
    print(paste("ReturnSubCodaFromFiles: SubSetCoords before fun = (",
      min(SubSetCoords), ", ", max(SubSetCoords), ")", sep="")); flush.console(); Verbose = 3;
  }
  SubSetCoords = SubSetCoords[SubSetCoords > 0 & SubSetCoords <= TBSRoo$p]
  CSubSetCoords = SubSetCoords -1;
  
  if (is.null(SaveDir)) {
    if (Verbose > 0) {
      print("SaveDir is Null, so can't use it"); flush.console();
    }
  }
  if (is.null(TBSRoo)) {
    if (Verbose > 0) {
      print("TBSRoo is Null so we can't use it"); flush.console();
    }
  }  else if (is.null(TBSRoo$ListSaveFiles)) {
    if (Verbose > 0) {
      print("TBSRoo is not null, but it does have null ListSaveFiles"); flush.console();
    }
  }
  if (!is.null(SaveDir)) {
    Dir = SaveDir;
    if (Verbose > 1) {
      print(paste("We're going to pull files from Save Dir = ", SaveDir, sep=""));
      print(paste("ListSaveFiles = (", paste(unlist(TBSRoo$ListSaveFiles), collapse=", "),
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
    if (length(TBSR5$ListSaveFiles) >= 1) {
      ATO <- rep(0, length(IFiles));
      AFTT <- sort(unique(unlist(TBSR5$ListSaveFiles)));
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
      
      for (jj in 1:length(TBSR5$ListSaveFiles)) {
         ALT <- substr(IFiles, 1, nchar(NewList[jj]));
         AM <-(1:length(ALT))[ ];
         ATO[ATO == 0 & ALT == NewList[jj]] = jj;
      }
      IFiles = IFiles[ATO > 0];
    }
    DFiles = MyNFiles[substr(MyNFiles, nchar(MyNFiles)-4, nchar(MyNFiles)) ==
      "D.bin"];  
    iTFiles = MyNFiles[substr(MyNFiles, nchar(MyNFiles)-5, nchar(MyNFiles)) ==
      "iT.bin"];               
    dTFiles = MyNFiles[substr(MyNFiles, nchar(MyNFiles)-5, nchar(MyNFiles)) ==
      "dT.bin"]; 
    PFiles = MyNFiles[substr(MyNFiles, nchar(MyNFiles)-4, nchar(MyNFiles)) ==
      "P.bin"];
    PostBuffFiles =   MyNFiles[
      substr(MyNFiles, nchar(MyNFiles)-nchar("PostBuff.bin")+1, 
      nchar(MyNFiles)) == "PostBuff.bin"];
    
    ILocFiles = MyNFiles[substr(MyNFiles, nchar(MyNFiles)-nchar("iLoc.bin")+1,
      nchar(MyNFiles)) == "iLoc.bin"]; 
    ProbFiles = MyNFiles[substr(MyNFiles, nchar(MyNFiles)-nchar("Prob.bin")+1,
      nchar(MyNFiles)) == "Prob.bin"];
      
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
    
    MDILoc = match(substr(ILocFiles, 1, nchar(ILocFiles) - nchar("ILoc.bin")),
      substr(ProbFiles, 1, nchar(ProbFiles)-nchar("Prob.bin")));
    MDILoc[is.na(MDILoc)] = -10000;
    sMDILoc = sort(MDILoc, index=TRUE)
    if (length(sMDILoc$x[sMDILoc$x != -10000]) > 0) {
      ILocFiles = ILocFiles[sMDILoc$ix[sMDILoc$x != -10000]];
      ProbFiles = ProbFiles[sMDILoc$x[sMDILoc$x != -10000]];
    } else {
    }  
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
      
    MPI = match(substr(PFiles, 1,nchar(PFiles)-5), 
      substr(IFiles, 1,nchar(IFiles)-5));  
    MPI[is.na(MPI)] = -10000;
    sMPI = sort(MPI, index=TRUE);
    if (length(sMPI$x[sMPI$x != -10000]) > 0) {
      PFiles = PFiles[sMPI$ix[sMPI$x != 10000]]; 
    }
     
    MPBI = match(substr(PostBuffFiles, 1
      ,nchar(PostBuffFiles)-nchar("PostBuff.bin")), 
      substr(IFiles, 1,nchar(IFiles)-nchar("I.bin") )); 
    MPBI[is.na(MPBI)] = -10000;
    sMPBI = sort(MPBI, index=TRUE);
    if (length(sMPBI$x[sMPBI$x != -10000]) > 0) {
      PostBuffFiles = PostBuffFiles[sMPBI$ix[sMPBI$x != -10000]];     
    }
    
    if (Verbose > 2) {
    print(paste("ReturnSubCodaFiles: Load from Save Dir IFiles ",
      length(IFiles), sep=""));
    print(paste("   IFiles was: ",  paste(IFiles, collapse=", "), ""));
    flush.console();
    }
  }  else if (!is.null(TBSRoo$ListSaveFiles)) {
    if (Verbose > 1) {
      print(paste("We're going to work with TBSRoo$ListSaveFiles", sep=""));
      print(paste("ListSaveFiles = (", paste(unlist(TBSRoo$ListSaveFiles), collapse=", "),
         ")", sep="")); flush.console();
      flush.console();
    }
    Dir = TBSRoo$SaveDir;
    LSFiles = unlist(TBSRoo$ListSaveFiles);
    for (ii in 1:length(LSFiles)) {
      if (substr(LSFiles[ii], 1, nchar(Dir)) == Dir) {
        LSFiles[ii] = substr(LSFiles[ii], nchar(Dir)+1, nchar(LSFiles[ii]));
        AT = unlist(strsplit(LSFiles[ii], "/"))
        LSFiles[ii] = AT[length(AT)]
      }
    }
    DFiles =  paste(LSFiles, "D.bin", sep="");
    IFiles =  paste(LSFiles, "I.bin", sep="");
    dTFiles =  paste(LSFiles, "dT.bin", sep="");
    iTFiles =  paste(LSFiles, "iT.bin", sep="");
    PFiles =  paste(LSFiles, "P.bin", sep=""); 

    ILocFiles =  paste(LSFiles, "ILoc.bin", sep="");
    ProbFiles =  paste(LSFiles, "Prob.bin", sep=""); 
    PostBuffFiles =   paste(LSFiles, "PostBuff.bin", sep=""); 
    if (Verbose > 2) {
    print(paste("ReturnSubCodaFiles: Load from ListSaveFiles IFiles ",
      length(IFiles), sep=""));
    print(paste("   IFiles was: ",  paste(IFiles, collapse=", "), ""));
    flush.console();
    }
  } else {
    print(paste("UhRoh:  TBSRoo$ListSaveFiles = ",
      paste(TBSRoo$ListSaveFiles, collapse=", "), sep=""));   flush.console();
    print(paste("  And SaveDir = ", SaveDir, sep=""));
    print("That Leaves us pretty S out of luck!"); flush.console();
  }
  if (Verbose > 1 && length(IFiles) >= 1) {
    print(paste("The IFiles are: ", paste(IFiles, collapse=", "), sep=""));
    flush.console();
    print(paste("The DFiles are: ", paste(DFiles, collapse=", "), sep=""));
    flush.console();
  }
  if (!exists("StartIter") || is.null(StartIter)) {
    StartIter = 1;
  }
  if (!exists("EndIter") || is.null(EndIter) || (is.numeric(EndIter) && EndIter <= 0)) {
    EndIter = TBSR5$MaxGibbsIters;
  }
  CodaBeta = list();
  if (length(IFiles) > 0) {
  for (ii in 1:length(IFiles)) {  
    NewCoda = NULL;
    if (Verbose > 0) {
      print(paste("  Going to GiveCodaBetaSubSet for IFiles[", ii, "] = ",
        IFiles[ii], sep="")); flush.console();
    }
    IOnFile = paste(Dir, "//", IFiles[ii], sep=""); 
    DOnFile = paste(Dir, "//", DFiles[ii], sep="");
    try( NewCoda <- .Call("GiveCodaSubset", CSubSetCoords, IOnFile, 
      DOnFile, StartIter, EndIter, Verbose-1), silent=FALSE);
    if (is.null(NewCoda)) {
    
    } else {
      if (is.null(names(TBSRoo$Beta))) {
        try(colnames(NewCoda) <- paste("Beta:", SubSetCoords, sep=""));
      } else {
        try(colnames(NewCoda) <- names(TBSRoo$Beta)[SubSetCoords]);
      }
      NewCoda = as.mcmc(NewCoda);
      CodaBeta[[length(CodaBeta) + 1]] = NewCoda;
    }
  }
  CodaBeta = as.mcmc.list(CodaBeta);  
  }

  if (Verbose > 1 && length(iTFiles) >= 1) {
    print(paste("The iTFiles are: ", paste(iTFiles, collapse=", "), sep=""));
    flush.console();
    print(paste("The dTFiles are: ", paste(dTFiles, collapse=", "), sep=""));
    flush.console();
  }  
  CodaTau = list();
  TauCoords = 1:length(TBSRoo$tauEndList);
  CTauCoords = TauCoords -1;
  if (length(iTFiles) > 0) {
  for (ii in 1:length(iTFiles)) {  
    NewCoda = NULL;
    if (Verbose > 0) {
      print(paste("  Going to GiveCodaTauSubSet for IFiles[", ii, "] = ",
        iTFiles[ii], sep="")); flush.console();
    }
    iTOnFile = paste(Dir, "//", iTFiles[ii], sep=""); 
    dTOnFile = paste(Dir, "//", dTFiles[ii], sep="");
    try( NewCoda <- .Call("GiveCodaSubset", CTauCoords, iTOnFile, 
      dTOnFile, StartIter, EndIter, Verbose-1), silent=FALSE);
    if (is.null(NewCoda)) {
    
    } else {
      if (is.null(names(TBSRoo$sOnTau))) {
        try(colnames(NewCoda) <- paste("tau:", TauCoords, sep=""));    
      } else {
        try(colnames(NewCoda) <- names(TBSRoo$sOnTau)[TauCoords]);
      }
      NewCoda = as.mcmc(NewCoda);
      CodaTau[[length(CodaTau) + 1]] = NewCoda;
    }
  }
  CodaTau = as.mcmc.list(CodaTau);  
  }
  
  if (is.numeric(DoProb) && length(DoProb) == 1 && DoProb == -1) {
    if (.self$p >= 50000) {
      DoProb <- FALSE;
    } else {
      DoProb <- TRUE;
    }
  } else if (is.numeric(DoProb) && length(DoProb) >= 1 && DoProb[1] >= .5) {
    DoProb <- TRUE;
  } else {
    DoProb <- FALSE;
  }
  if (is.numeric(DoPostBuff) && length(DoPostBuff) == 1 && DoPostBuff == -1) {
    if (.self$p >= 50000) {
      DoPostBuff <- FALSE;
    } else {
      DoPostBuff <- TRUE;
    }
  } else if (is.numeric(DoPostBuff) && length(DoPostBuff) >= 1 && DoPostBuff[1] >= .5) {
    DoPostBuff <- TRUE;
  } else {
    DoPostBuff <- FALSE;
  }

  CodaP = list();  ##CodaProbBuff <- list();
  if (length(PFiles) > 0 && DoProb == TRUE) {
  print(paste("  ABayesSpikeStart.r::1471:GiveCodaSubSet, PFiles = ", PFiles, " TBSRoo$DoProb = ", TBSRoo$DoProb, sep=""));
  flush.console();
  ##try(CodaP <- TBSRoo$GenerateSubProbCodaList());
  if (TRUE) {
  for (ii in 1:length(PFiles)) {
    NewCoda = NULL;
    if (Verbose > 0) {
      print(paste("  Going to GiveCoda SubSet for PFiles[", ii, "] = ",
        PFiles[ii], sep="")); flush.console();
    } 
    MyOA = NULL; 
    try(MyOA <- file(description = paste(Dir, "//", PFiles[ii], sep=""), 
      open = "rb", blocking = TRUE,
     encoding = getOption("encoding"), raw = FALSE));
    if (is.null(MyOA)) {
    } else {
      MyF = NULL;
      PreL = 0;
        if (sum(TBSRoo$DoRecord[3:4]) == 2) {
          PreL = 2;
        } else if (TBSRoo$DoRecord[3] == 1) {
          PreL = 1;
        } else if (TBSRoo$DoRecord[4] == 1) {
          PreL = 1;   
        }
      if (TBSRoo$p <= 2000) {
      try(MyF <- readBin(MyOA, what="double", n =PreL * TBSRoo$MaxGibbsIters ));
      if (!is.null(MyF)) {
        if (sum(TBSRoo$DoRecord[3:4]) == 2) {
          MyF = t(matrix(MyF, 2, length(MyF) / 2));
          colnames(MyF) = c("Sigma", "PiA");
        } else if (TBSRoo$DoRecord[3] == 1) {
          MyF = matrix(MyF, length(MyF), 1);
          colnames(MyF) = c("Sigma");
        } else if (TBSRoo$DoRecord[4] == 1) {
          MyF = matrix(MyF, length(MyF), 1);
          colnames(MyF) = c("PiA");      
        }
        NewCoda = as.mcmc(MyF);
        CodaP[[length(CodaP)+1]] = NewCoda;
      }
      } else {
      if (PreL <= 0) { PreL <- 1; }
      ATLoad <- 0;  MaxData <- 2000000;  ItersAtATime <- ceiling(MaxData /
        PreL);
      NewCoda <- matrix(0,  TBSRoo$MaxGibbsIters, PreL);
      while(ATLoad < TBSRoo$MaxGibbsIters) {
      try(MyF <- readBin(MyOA, what="double", n = PreL * ItersAtATime ));
      if (!is.null(MyF)) {
        if (sum(TBSRoo$DoRecord[3:4]) == 2) {
          MyF = t(matrix(MyF, 2, length(MyF) / 2));
          try(NewCoda[ATLoad + 1:NROW(MyF),] <- MyF);
          try(colnames(NewCoda) <- c("Sigma", "PiA"));
        } else if (TBSRoo$DoRecord[3] == 1) {
          MyF = matrix(MyF, length(MyF), 1);
          try(NewCoda[ATLoad+1:NROW(MyF),] <- MyF)
          try(colnames(NewCoda) <- c("Sigma"));
        } else if (TBSRoo$DoRecord[4] == 1) {
          MyF = matrix(MyF, length(MyF), 1);
          try(NewCoda[ATLoad+1:NROW(MyF),] <- MyF);
          try(colnames(NewCoda) <- c("PiA"));      
        }
        ATLoad <- ATLoad + ItersAtATime; 
      }
      }
      NewCoda <- as.mcmc(NewCoda);
      try(close(MyOA));
    }
  }
  CodaP = as.mcmc(CodaP);
  }
  }}
  
  CodaProb = list();
  if (length(ProbFiles) > 0 && DoProb == TRUE) {
    if (TBSRoo$Verbose >= 1) {
       print(paste("ABayesSpikeStart.r::1545::ReturnSubCodaFromList, ",
         "Getting ProbCodaList, TBSRoo$DoProb = ", TBSRoo$DoProb, sep=""));
       flush.console();
    }
    try(CodaProb <- TBSRoo$GenerateSubProbCodaList());
    if (FALSE) {
    for (ii in 1:length(ProbFiles)) {
      NewCoda = NULL;
      if (Verbose > 0) {
        print(paste("  Going to GiveCoda SubSet for ProbFiles[", ii, "] = ",
          ProbFiles[ii], sep="")); flush.console();
      } 
      MyOA = NULL; 
      try(MyOA <- file(description = paste(Dir, "//", ProbFiles[ii], sep=""), 
        open = "rb", blocking = TRUE,
       encoding = getOption("encoding"), raw = FALSE));
      if (is.null(MyOA)) {
      } else {
        MyF = NULL;
        try(MyF <- readBin(MyOA, what="double"));
        if (!is.null(MyF)) {
          CodaProb[[length(CodaProb)+1]] = as.mcmc(MyF);
        }
        try(close(MyOA));
      }
    }
    CodaProb = as.mcmc(CodaProb);
    }
  }
 
  CodaILoc = list();
  ## Note CodaILoc is a list of locations in IFile and JFile where the next iteration Starts
  ##  It is only of length 1.
  if (length(ILocFiles) > 0 && DoILoc == TRUE) {
  if (TBSRoo$Verbose >= 1) {
    print(paste("ABayesSpikeStart.r:::1580:::ReturnSubCodaFromList(), ",
      "doing ILoc because DoILoc == ", TBSRoo$DoILoc, sep=""));
    flush.console();
  }
  for (ii in 1:length(ILocFiles)) {
    NewCoda = NULL;
    if (Verbose > 1) {
      print(paste("  Going to GiveCoda SubSet for ILocFiles[", ii, "] = ",
        ILocFiles[ii], sep="")); flush.console();
    } 
    MyOA = NULL; 
    try(MyOA <- file(description = paste(Dir, "//", ILocFiles[ii], sep=""), 
      open = "rb", blocking = TRUE,
     encoding = getOption("encoding"), raw = FALSE));
    if (is.null(MyOA)) {
    } else {
      MyF = NULL;
      try(MyF <- readBin(MyOA, what="integer"));
      if (!is.null(MyF)) {
        CodaILoc[[length(CodaILoc)+1]] = as.mcmc(MyF);
      }
      try(close(MyOA));
    }
  }
  CodaILoc = as.mcmc(CodaILoc);
  }
  CodaPostBuff = list();
  if (length(PostBuffFiles) > 0 && DoPostBuff == TRUE) {
    if (TBSRoo$Verbose >= 1) {
      print(paste("  ABayesSpikeStart.r:::1610:::ReturnSubCodaFromFiles(), ",
        "we DoPostBuff because TBSRoo$DoPostBuff = ", TBSRoo$DoPostBuff, sep=""));
      flush.console();
    }
    try(CodaPostBuff <- TBSRoo$GenerateSubProbCodaList());
    if (FALSE) {
    LengthPostBuff  = 0;
    if (is.null(TBSRoo$sOnTau) || length(TBSRoo$sOnTau) <= 0) {
      LengthPostBuff = TBSRoo$p;
    } else {
      LengthPostBuff = TBSRoo$FirstRandom-1 + length(TBSRoo$sOnTau);
    }
    for (ii in 1:length(PostBuffFiles)) {
       NewCoda = NULL;
       if (Verbose > 1) {
        print(paste("  ABayesSpikeStart.r::1624::ReturnSubCodaFromFiles() TBSRoo$DoPostBuff = ", 
          TBSRoo$DoPostBuff, " and Going to GiveCoda SubSet for PostBuffFiles[", ii, "] = ",
          PostBuffFiles[ii], sep="")); flush.console();
       } 
       MyOA = NULL; 
       try(MyOA <- file(description = paste(Dir, "//", PostBuffFiles[ii], sep=""), 
         open = "rb", blocking = TRUE,
         encoding = getOption("encoding"), raw = FALSE));
       if (is.null(MyOA)) {
       } else {
         MyF = NULL;
         try(MyF <- readBin(MyOA, what="double", n = LengthPostBuff * TBSRoo$MaxGibbsIters));
         LLMyF = length(MyF) / LengthPostBuff;
         try(MyF <- t(matrix(as.numeric(MyF), LengthPostBuff, LLMyF)))
         if (!is.null(TBSRoo$sOnTau) && length(TBSRoo$sOnTau) > 0) {
           if (TBSRoo$FirstRandom > 1) {
             if (is.null(names(TBSRoo$sBeta)) || length(names(TBSRoo$sBeta)) <= 0) {
               NB = paste("Beta:", 1:length(TBSRoo$Beta), sep="");
             } else {
               NB = names(TBSRoo$sBeta);
             }
             if (is.null(names(TBSRoo$sOnTau)) || length(names(TBSRoo$sOnTau)) <= 0) {
               NT = paste("Tau:", 1:length(TBSRoo$sOnTau), sep="");
             } else {
               NT  = names(TBSRoo$sOnTau);
             }
             try(colnames(MyF) <- paste("Prob:", c(NB[1:(TBSRoo$FirstRandom-1)],
               NT), sep="") );
           } else {
             try(colnames(MyF) <- names(TBSRoo$sOnTau));
           }
         }  else {
           try(colnames(MyF) <- names(TBSRoo$sBeta));
         }
         if (DoLogitPostProb == FALSE) {
           MyF[is.na(MyF)] = -9999;
           NewMyF = matrix(0, NROW(MyF), NCOL(MyF));
              NewMyF[MyF >= 10] <- 1.0;
              NewMyF[MyF <= -10] <- exp(MyF[MyF <= -10]);
              NewMyF[MyF < 10 & MyF > -10] <- exp(MyF[MyF < 10 & MyF > -10])/
               (1.0 + exp(MyF[MyF < 10 & MyF > -10]));
              colnames(NewMyF) <- colnames(MyF);
              NewMyF[MyF == -9999] <- NA;
              MyF <- NewMyF;    
         }
         if (!is.null(MyF)) {
           CodaPostBuff[[length(CodaPostBuff)+1]] = as.mcmc(MyF);
         }
         try(close(MyOA));
       }
    }
    }
    try(CodaPostBuff <- as.mcmc.list(CodaPostBuff));
  }
  if (!is.null(CodaPostBuff) && length(CodaPostBuff) > 1) {
    try(CodaPostBuff <- as.mcmc.list(CodaPostBuff));
  }
  
  if (Verbose > 2) {
    print(paste("ReturnSubCodaFromFiles.r:::ReturnSubCodaFiles()::573:: Finished length(Beta) is ",
      length(CodaBeta), sep=""));
    print(paste("   IFiles was: ",  paste(IFiles, collapse=", "), ""));
    flush.console();
  }
  return(list(CodaBeta=CodaBeta, CodaTau = CodaTau, CodaP = CodaP,
    CodaProb = CodaProb, CodaILoc = CodaILoc, CodaPostBuff = CodaPostBuff));
}
