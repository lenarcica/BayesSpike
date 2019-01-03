

####################################################################################
##  ZeroOutFunction()
##
##   When p is large, this chooses a sparse subset of coordinates to start 
##  the sampler chains based upon one draw of MBX$SampleNewTaus() and MBS$SampleBFixed()
##  which fill the ProbTau and ProbFixed Vectors respectively
##
ZeroOutFunction <- function(MBS) {

Verbose <- MBS$Verbose; 

print("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$")
print("$$ ZeroOutFunction"); flush.console();
print("$$ Hey: MBS$ is large before a sample we are going to blank up and sample small");
flush.console();
MBS$DoAddCoordsOnSetBeta <- 0;
RfFlag <- MBS$RefreshBeta(rep(0,MBS$p));
if (Verbose >= -1) {
  print("$$ ZeroOutFunction: We Refreshed Beta into zero."); flush.console();
}
MyTryProbs <- NULL;
MyTryProbTau <- NULL;  MyTryProbFixed <- NULL;
if (MBS$TBSR5$OnChainIter >= 2) {
  if (Verbose >= -1) {
    print("$$ ZeroFunction: About to check for MBS$ProbVector"); flush.console();
  }
  try(MyTryProbs <- MBS$ProbVector);
  if (!is.null(MyTryProbs)) {
    if (is.null(MBS$tauEndList)) {
      MyTryProbFixed <- log(MyTryProbs/(1.0-MyTryProbs));  
      MyTryProbFixed[(MyTryProbFixed > 0 & !is.finite(MyTryProbFixed)) | MyTryProbs == 1.0] <- 4;
      MyTryProbFixed[(MyTryProbFixed < 0 & !is.finite(MyTryProbFixed)) | MyTryProbs == 0.0] <- -4;
    } else if (MBS$FirstRandom > 1 && MBS$FirstRandom <= MBS$p) {
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
} else {
  if (Verbose >= -1) {
    print("$$  ZeroOutFunction: Setting onSigma Based upon SigmaPrior"); flush.console();
    print("$$ Conducting Test of Sigma and Sigma Prior");   flush.console();
    MBS$TestSigma();
  }
  if (MBS$LengthSigmaPrior >= 2 && MBS$SigmaPrior[2] > 0) {
    try(MBS$OnSigma <- 1.5*MBS$SigmaPrior[2])
  }
  if (Verbose >= -1) {
    print("$$ ZeroOutFunction: Successfully placed OnSigma"); flush.console();
  }
}
if (Verbose >= -1) {
  print("$$ ZeroOutFunction: Conduct BlankAllNewCoords, set all Beta to Zero."); flush.console();
}
try(MBS$BlankAllNewCoords());
if (Verbose >= -1) {
  print("$$ ZeroOutFunction: Finished BlankAllNewCoords. Now check for MBS$tauEndList"); flush.console();
}
  if (length(MBS$tauEndList) >= 1) {
    if (Verbose >= -1) {
      print("$$ ZeroOutFunction: MBS$tauEndList is length greater than 1 or more."); flush.console();
    }
    if (!is.null(MyTryProbTau)) {
      MyProb <- MyTryProbTau;
    } else {
      MBS$SampleTausOnly <- 1;
      MBS$SampleNewTaus();
      MyProb <- MBS$ProbTau;
    }
    if (MBS$FirstRandom >= MBS$p) {
      print(paste("$$ ZeroOutFunction: ERROR ERROR ERROR, FirstRandom = ", MBS$FirstRandom, sep=""));
      flush.console();
      print("ERRORERRORERRORERRORERRORERRORERROR"); flush.console();
      MyT = "1 = 2";
      eval(parse(text=MyT));
      return(-1);
    }
    GoTau  <- MBS$OnTau;
    print(paste("$$ ZeroOutFunction: We sample Taus and the number of nonzero are ", length(GoTau[GoTau>0]),
     "/", length(GoTau), sep="")); flush.console();
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
    if (KeepCount <= 5) { KeepCount = 5; } 
    if (KeepCount > length(MBS$tauEndList)) {
      KeepCount <- length(MBS$tauEndList);
    }    
    if (KeepCount < length(GoTau[GoTau>0])) {
      if (Verbose >= -1) {
        print("$$ 0010ResortBayesSpike.r: ZeroOutFunction: Conducting Taus version."); flush.console();
      }
      Expt <- GiveMeNewBalance(lX=MyProb/2.0, WantTot = KeepCount, StartMove=1.0, MaxTTs = 300);
      if (Verbose >= -1) {
        print("0010ResortBayesSpike.r:ZeroOutFunction() We returned from GiveMeNewBalance."); flush.console();
      }
      if (length(Expt[Expt > 0]) < KeepCount) {
        KeepCount <- length(Expt[Expt > 0]);
      }
      print(paste("$$ after calculating Expt we have its sum is ", round(sum(Expt),4),
        "/", length(Expt), sep="")); flush.console();
      Keepers <- sample(1:length(GoTau), prob=Expt, replace=FALSE, size=KeepCount);
      print(paste("$$ Now we reduce active set to length(Keepers) = ", length(Keepers), sep=""));
      flush.console();
      DoTau <- rep(0, length(GoTau));
      DoTau[Keepers] <- GoTau[Keepers];
      try(MBS$DoAddCoordsOnSetBeta <- 0);
      try(MBS$BlankAllNewCoords());
      print(paste("AssignTau: about to Assign ")); flush.console();
      try(MBS$AssignTau(DoTau));
      print(paste("$$ Reduced DoTau length(Keepers) = ", length(Keepers), " is added", sep=""));
      print(paste("$$ After AssignTau, AllNewCoords = ", MBS$AllNewCoords, sep="")); flush.console();
      flush.console();
    } else {
      print(paste("$$ KeepCount = ", KeepCount, " but length(GoTau) = ",
        length(GoTau[GoTau>0]), "/", length(GoTau), " so no change.", sep=""));
      flush.console(); 
      if (length(GoTau[GoTau > 0]) <= KeepCount) {
        MyG <- sort(MyProb, decreasing=TRUE, index=TRUE)$ix[1:KeepCount];
        DG <- GoTau[MyG];
        DG[DG == 0.0] <- MBS$tauPriorMean
        GoTau[MyG] <- DG;
      } else {
        Keepers <- 1:length(GoTau);
      }
      try(MBS$BlankAllNewCoords());
      print(paste("$$ Now Assigning Tau to GoTau")); flush.console();
      MBS$AssignTau(GoTau);
      print(paste("$$ After AssigningTau to GoTau length ", length(GoTau), " we ",
        " have AllNewCoords = ", MBS$AllNewCoords, sep=""));
      flush.console();
    }
  }
  print("$$ Well looking into the length of MBS$tauEndList"); flush.console();
  if (length(MBS$tauEndList) >= 1) {
    AllRandomCoords <- rep(1,1);
    AllRandomCoords[1] <- MBS$AllNewCoords[1];
  } else {
    AllRandomCoords <- NULL;
  }
  if (Verbose >= -1) {
    print("$$ We have Checked to see that AllRandomCoords has been set due to tauEndList Length"); 
    flush.console();
  }
  if((MBS$FirstRandom > 1 && MBS$FirstRandom <= MBS$p) || length(MBS$tauEndList) <= 0) {
    print("$$ ZeroOutFunction(): We definitely have a fixed effects somehow. "); flush.console();
      if (length(MBS$tauEndList) == 0) {  
        if (Verbose >= -1) {
          print("$$ ZeroOutFunction(): tauEndList has zero length so pF is all of MBS.\n"); flush.console(); 
        }
        pF = MBS$p;  
      } else { 
        pF <- MBS$FirstRandom-1;
      }
      if (Verbose >= -1) {
        print("$$ ZeroOutFunction(): Going to try MyTryProbFixed"); flush.console();
      }
      if (!is.null(MyTryProbFixed)) {
        if (MBS$Verbose >= -1) {
          print("$$ ZeroOutFunction(): Set BProbs to MyTryProbFixed"); flush.console();
        }
        BProbs <- MyTryProbFixed;
        if (MBS$Verbose >= -1) {
          print(paste("$$ ZeroOutFunction(): Probability of Historical fixed variables ",
            "has Historical MIPS summing to ", round(sum(dLogit(BProbs)),4), "/", pF, sep=""));
         flush.console();
        }
      } else {
        if (MBS$Verbose >= -1) {
           print("$$ ZeroOutFunction(): setting MBS negative -1 to start."); flush.console();
        }
        MBS$tt <- -1;
        if (MBS$Verbose >= -1) {
          print(paste("$$  ZeroOutFunction() : We're refreshing fixed effects and conducting SampleFixedB, tt = ", MBS$tt)); flush.console();
          print("$$ ZeroOutFunction(): Set SampleTausOnly"); flush.console();
        }
        MBS$SampleTausOnly <- 1; 
        if (MBS$Verbose >=-1) {
          print("$$ ZeroOutFunction() Now SampleFixedB"); flush.console();
        }
        MBS$SampleFixedB();
        print(paste("$$ Now to reset SampleTausOnly back to zero.")); flush.console();
        MBS$SampleTausOnly <- 0;
        print(paste("$$ Now to get ProbFixed from MBS")); flush.console();
        BProbs <- MBS$ProbFixed + 0.0;
        print(paste("$$ Sample TausOnly, we retrieved BProbs first time", sep="")); flush.console();
        NewBProbs <- rep(0, length(BProbs));
        for (jtj in 1:length(BProbs)) {NewBProbs[jtj] <- BProbs[jtj];}
        print(paste("$$ SampleTausOnly, we copied NewBProbs", sep="")); flush.console();
        BProbs <- NewBProbs;
        MBS$tt <- 0
        if (Verbose >= -1) {
          print(paste("$$ ZeroOutFunction() : MIPs taken from MBS SampleFixedB()  average logit Probs is: ", 
            round(sum(dLogit(BProbs)),4), " over a count of ", pF, " coeficients", sep=""));
          flush.console();
        }
      }     
      
      OnPiAF <- MBS$OnPiA[1];
      KeepCount <- round(25*pF * OnPiAF);
      if (KeepCount <= 0) { KeepCount = 1; }
      if (KeepCount > pF) { KeepCount <- pF; }
    
      if (MBS$Verbose >= -1) {
        print("$$ ZeroOutFunction Estimating GiveMeNewBalance, filling BProbs and then going to GiveMeNewBalance"); flush.console();
      }  
      lX <- rep(0, length(BProbs));
      for (jtj in 1:length(lX)) {
        lX[jtj] <- BProbs[jtj] / 2.0;
      }
      eBProbs <- GiveMeNewBalance(lX=lX, WantTot = KeepCount, StartMove=1.0, MaxTTs = 300, Verbose = MBS$Verbose);
      if (MBS$Verbose >= -1) {
        print("$$ ZeroOutFunction: GiveMeNewBalance Returned"); flush.console();
      }
      if (length(eBProbs[eBProbs > 0])  < KeepCount) {
        KeepCount <- length(eBProbs[eBProbs > 0]); 
      }                                                   
      if (Verbose >= -1) {
        print("$$ ZeroOutFunction() Sampling Keepers."); flush.console();
      }
      Keepers <- sample(1:pF, prob=eBProbs, replace=FALSE, size=KeepCount);
      DKeepers <- (1:pF)[!(( 1:pF)  %in% Keepers )]; 
      print(paste("$$  We have chosen to eliminate a count of Keepers ", length(DKeepers)));
      flush.console();
      if (MBS$FirstRandom <= 0 || is.null(MBS$tauEndList) || length(MBS$tauEndList) <= 0) {
        ABeta <- MBS$Beta;
        ABeta[MBS$Beta == 0.0 && MBS$BFixed >= 1] <- 1;
      } else if (MBS$FirstRandom > 1 && MBS$FirstRandom <= MBS$p) {
        ABeta <- c(rep(1, MBS$FirstRandom-1), rep(0, MBS$p- MBS$FirstRandom+1)) * MBS$Beta;
        ABeta[(1:length(MBS$BFixed))[MBS$Beta == 0.0 && MBS$BFixed >= 1]] <- 1;
      }
      
      if (length(DKeepers) > 0) {
        ABeta[DKeepers] <- 0.0;
      }
      try(MBS$DoAddCoordsOnSetBeta <- 0);
     
      if (MBS$Verbose >= -1) {
        print("$$ About to RefreshBeta"); flush.console();
      }
      try(MBS$RefreshBeta(ABeta));
      ##try(MBS$SetAllNewCoords(length(AllRandomCoords)+length(Keepers), length(Keepers)));
    }
    print("About to Count All New Coords() "); flush.console();
    try(MBS$CountAllNewCoords());
    print(paste("$$  Adding initially ", MBS$AllNewCoords, "/", MBS$p, " new coords, maximum allotment is ",
      MBS$MaximumAllocation, ", OnKappaMem = ",  MBS$OnKappaMem, sep=""));  flush.console();
    print(paste("$$ length(MBS$Beta[MBS$Beta != 0.0]) == ", 
      length(MBS$Beta[MBS$Beta != 0.0]), sep="")); flush.console();
    print(paste("$$  Resizing now to run one regression and hopefully eliminate noise variables. ", sep=""));
    flush.console();
    try(MBS$AddAllNewCoords());
    print(paste("$$ We have AllNewCoords added, OnKappaS = ", MBS$OnKappaS, 
      ", OnKappaMem = ", MBS$OnKappaMem, ", NumActive = ", MBS$NumActive, sep=""));   flush.console();
    print(paste("$$ Now we Refresh OrderedActive")); flush.console();
    try(AGo <- NULL);
    try(AGo <- MBS$RefreshOrderedActive(1));
    print(paste("$$  Refreshing XtResid. ")); flush.console();
    if (MBS$dfRobit > 0) {
      MBS$RobitReplace();
    } else if (MBS$dfTNoise > 0) {
      MBS$UpdateTNoise();
    } else {
      MBS$UpdateFreshXtResid();
    }
    try(print(paste("$$  UpdateFreshXtResid conducted, p is ", MBS$p, ", OnKappaS is ",
      MBS$OnKappaS, ", MBS$NumActive is ", MBS$NumActive, "Our Test is: ", sep=""))); flush.console();
    ATestCount <- NULL;
    BackVerbose <- MBS$Verbose;
    ##MBS$Verbose = 4;
    try(ATestCount <- MBS$TestXtResid());
    MBS$Verbose <- BackVerbose;
    try(print(paste("$$ ATestCount from Test For XtResid is ", ATestCount, sep=""))); flush.console();
    if (!is.null(ATestCount) && ATestCount > 0)  {
      print("$$ 0010ReSortBayesSpike.r:::WE Quit"); flush.console();
      print("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$");
      flush.console();
      eval(parse(text=SetGText("MBS", "globalenv()", S=1)));
      return(-1);
    }
    try(print(paste("$$ Now about to test AGo.", sep=""))); flush.console();
    if (is.null(AGo) || (is.numeric(AGo) && AGo < 0)) {
      print("$$ Resample AGo, we failed to RefreshOrderedActive, this is a bad result, no regressions"); flush.console();
      return(-1);
    }
    print(paste("$$ Conducting PrepareForRegression with NumActive = ", MBS$NumActive, sep="")); flush.console();
    if (MBS$NumActiveFixed > MBS$NumActive) {
      print(paste("$$ NumActive: Why is NumActiveFixed = ", MBS$NumActiveFixed, " or NumActive is ", MBS$NumActive, sep="")); flush.console();
    }
    if (MBS$NumActive > 0) {
      if (MBS$NumActive == 1) {
        print(paste("$$ Obviously we have only 1 member. OrderedActive is: ", MBS$OrderedActive, sep=""));
        try(print(paste("$$ Which Suggests XLC[", MBS$OrderedActive, "] is ", MBS$XLC[MBS$OrderedActive[1]+1], sep=""))); flush.console();
        OA1 <- MBS$OrderedActive[1]; XLC1 <- MBS$XLC[OA1+1];
        try(print(paste("$$ the length of pXtX[", MBS$OrderedActive[1],",",
           XLC1,"] is then \"", 
           MBS$XtX[MBS$OrderedActive[1]+1, XLC1+1], "\"", sep=""))); flush.console();
        try(print(paste("$$ XtY[", OA1, "] = ", MBS$XtY[OA1+1], " and cXtY[", XLC1, "] is ", MBS$cXtY[XLC1+1], sep=""))); flush.console();
        print(paste("$$ OnKappaMem is ", MBS$OnKappaMem, " and OnKappaS is ", MBS$OnKappaS, " and NumberNew is ",
          MBS$NumberNew, sep="")); flush.console();
        try(print(paste("$$ We have get_LengthNoShrinkFixed() returns: ", MBS$LengthNoShrinkFixed, sep=""))); flush.console();
        print(paste("$$ NoShrinkFixedPrior is length ", length(MBS$NoShrinkFixedPrior), " and it is (", paste(MBS$NoShrinkFixedPrior, collapse=", ", sep=""), ")", sep="")); flush.console();
        try(print(paste("$$ rIChol is length", length(MBS$rIChol), " and length SmallXtX is ", length(MBS$SmallXtX), sep=""))); flush.console();
        print("$$ WE will now PrepareForRegression");
        MBS$Verbose = 4;
      }

      try(MBS$PrepareForRegression()); 
      if (MBS$NumActive == 1) {
        try(MBS$Verbose <- backV);
      }
      print(paste("$$  Now sampling PropBeta which will have length ", MBS$NumActive, sep=""));
      try(MBS$SamplePropBeta());
      CNBeta <- MBS$PropBeta[abs(MBS$PropBeta) > 999999]
      if (any(is.nan(MBS$PropBeta)) || length(CNBeta) >= 1) {
        print(paste("$$ CNBeta has errors in PropBeta, please new work.  Triggering an error", sep="")); flush.console();
        -1 <- 2;
      }
      try(print("$$ About to FillsBetaFromPropBetaAndCompute() ")); flush.console();
      try(MBS$FillsBetaFromPropBetaAndCompute()); 
    }
    try(MBS$Verbose <- BackVerbose);
    try(print("$$ About to UpdateSigma() ")); flush.console();
    try(MBS$UpdateSigma()); 
    try(print(paste("$$ ZeroOut Sampler, predictability after this is var(Y) = ",
      round(var(MBS$Y),6), " and YResid Squared is ", round(MBS$YResidSq,6), 
      " so Sigma = ", MBS$OnSigma, sep=""))); flush.console();
    print(paste("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$", sep=""));
    flush.console();
    return(1);
}

dLogit <- function(lX) {
  NI <- exp(lX) / (1.0+exp(lX));
  NI[lX < -10] = exp(-abs(lX[lX < -10]));
  NI[lX > 10] = 1.0 - exp(-abs(lX[lX > 10]));
  return(NI);
}
dLL <- function(lX, al) {
  NI <- rep(0, length(lX));
  for (ii in 1:length(NI)) {
    aV <- lX[ii]  + al;
    if (aV < -10) { NI[ii] <- exp(aV);
    } else if (aV > 10) { NI[ii] = 1.0-exp(-aV); 
    } else { NI[ii] = exp(aV)/(1.0+exp(aV)); }
  }
  return(NI);
}
##  GiveMeNewBalance: Helper function
##
##  Let lX be a real vector
##    We are looking for a vector LB = exp(lX + alpha) / (1.0 + exp(lX+alpha))
##     for some alpha such that  sum(LB) = WantTot
##
GiveMeNewBalance <- function(lX, WantTot = 1, StartMove = 1, MaxTTs = 300, epsilon = .00001, Verbose = -3) {
  if (WantTot <= 0) {
    print(paste("Error:GiveMeNewBalance, does not work if WantTot=", WantTot, sep=""));
    return(-1);
  }
  if (Verbose >= -1) {
    print("GiveMeNewBalance: Initiating"); flush.console();
  }
  if (length(lX) < WantTot) {
    print(paste("Error:GiveMeNeBalance, there are only ", length(lX),
      "members of lX not good enough to make ", WantTot, sep=""));
  }
  alphaOld = 0;
  
  if (Verbose >= -1) {
    print("0010ReSortBayesSpike:GiveMeNewBalance: writing deblurred form of lX to newlX"); flush.console();
  }
  newlX <- rep(0, length(lX));
  for (ii in 1:length(lX)) {
    newlX[ii] <- lX[ii];
  }
  if (MBS$Verbose >= -1) {
    print("0010ReSortBayesSpike:GiveMeNewBalance: now to remove lx"); flush.console();
  }
  rm(lX);
  lX = newlX;
  if (Verbose >= 2) {
    print("0010ReSortBayesSpike:GiveMeNewBalance: Concluded doing copy. "); flush.console();
  }
  
  AG <- dLogit(lX);
  if (Verbose >= -1) {
    print("0010ReSortBayesSpike:GiveMeNewBalance: We finished AG"); flush.console();
  }
  if (length(AG[AG >= 1.0]) >= WantTot) {
    print(paste("0010ReSortBayesSpike:GiveMeNewBalance, cannot be corrected, there are already ",
      " length(AG[AG >= 1.0]) = ",
      length(AG[AG >= 1.0]), " non zero, sample a subset then already!",
      sep="")); 
      flush.console();
    return(AG);
  }
  if (Verbose >= 2) {
    print("0010ReSortBayesSpike:GiveMeNewBalance: All Passed Trial 1.0"); flush.console();
  }
  if (length(AG[AG <= 0.0]) >= length(lX) - WantTot) {
    print(paste("0010ReSortBayesSpike:GiveMeNewBalance, cannot be corrected, there are already ",
      "length(AG[AG <= 0.0]) = ",
      length(AG[AG <= 0.0]), " are zero, be happy its small!",
      sep="")); 
      flush.console();
    return(dLogit(lX));
  }
  if (Verbose >= -1) {
    print("0010ReSortBayesSpike:GiveMeNewBalance: All Passed Trial 2.0"); flush.console();
  }
  LBOld <- sum(dLogit(lX));
  if (LBOld > WantTot) { StartMove = -1.0 *abs(StartMove) }
  alphaNew = alphaOld + StartMove;
  LBNew <- sum(dLogit(lX+alphaNew));
  tt = 0;
  if (Verbose >= 2) {
    print("GiveMeNewBalance: Initiating Trial 3.0"); flush.console();
  }
  if (LBNew > WantTot && LBOld > WantTot) {
    if (Verbose >= 2) {
      print("GiveMeNewBalance: while Loop 3.0a"); flush.console();
    }
    while(LBNew > WantTot && tt < MaxTTs) {
      alphaNew = alphaNew - abs(StartMove);
      LBNew <- sum(dLogit(lX+alphaNew)); tt <- tt+1;
      if (tt >= MaxTTs) {
        print(paste("GiveMeNewBalance, cannot succeed, tt = ", tt, 
          " and WantTot = ", WantTot, sep=""));  flush.console();
        return(-1);
      }
    }
    if (LBNew > WantTot) {
      print(paste("Error: GiveMeNewBalance in 300 moves we could not make WantTot = ",
        WantTot, " for LBNew = ", LBNew, " for alphaNew = ", alphaNew, sep=""));
      flush.console();  return(-1);
    }
  } else if (LBNew < WantTot && LBOld < WantTot) {
    if (Verbose >= -1) {
      print("GiveMeNewBalance: While Loop 3.0b"); flush.console();
    }
    while(LBNew < WantTot && tt < MaxTTs) {
      alphaNew = alphaNew + abs(StartMove);
      LBNew <- sum(dLogit(lX+alphaNew)); tt <- tt+1;
      if (tt >= MaxTTs) {
        print(paste("GiveMeNewBalance, cannot succeed, tt = ", tt, 
          " and WantTot = ", WantTot, sep=""));  flush.console();
        return(-1);
      }
    }
    if (LBNew < WantTot) {
      print(paste("Error: GiveMeNewBalance in 300 moves we could not make large WantTot = ",
        WantTot, " for LBNew = ", LBNew, " for alphaNew = ", alphaNew, sep=""));
      print(paste(" We will return a smaller form. ")); flush.console();
      return(dLogit(lX+1.0));
      flush.console();  return(-1);
    }
  }
  tt = 0;
  if (Verbose >= -1) {
    print(paste("GiveMeNewBalance: Testing alphaNew = ", alphaNew, 
      " versus alphaOld = ", alphaOld, sep="")); flush.console();
  }
  if (alphaNew < alphaOld) {
    LBMid <- LBNew;  alphaMid <- alphaNew;
    alphaNew <- alphaOld;  LBNew <- LBOld;
    alphaOld <- alphaMid;  LBOld <- LBMid;
  }   else {
    try(LBMid <- (LBNew + LBOld) / 2);
  }
  if (Verbose >= -1) {
    print(paste("0010ResortBayesSpike:GiveMeNewBalance: Intitiating While Loop 5.0, MaxTTs =", MaxTTs, sep="")); flush.console();
    print(paste("0010ResortBayesSpike:GiveMeNewBalance: length (lX) is ", length(lX), sep="")); flush.console();
    try(print(paste("0010ResortBayesSpike:   Start alphaNew = ", alphaNew, " and alphaOld = ", alphaOld, sep=""))); flush.console();
    try(print(paste("0010ResortBayesSpike:   Start with LBNew = ", LBNew, " and LBOld is ", LBOld, " and LBMid is ", LBMid, sep=""))); flush.console();
    try(print(paste("0010ResortBayesSpike:   epsilon is ", epsilon, " and MaxTTs is ", MaxTTs, sep="")));
  }
  while(abs(LBNew- WantTot) > epsilon && tt < MaxTTs) {
    tt <- tt+1;
    if (Verbose >= 2) {
      print(paste("0010ResortBayesSpike tt = ", tt, "/", MaxTTs, " we still have LBNew=", LBNew, " and WantTot = ", WantTot, 
        ", fbdiff = ", abs(LBNew-WantTot), " against ", epsilon, sep="")); flush.console();
    }
    alphaMid <- (alphaNew+alphaOld)  / 2.0;
    LBMid <- sum(dLogit(lX+rep(alphaMid, length(lX))));  
    if (abs(LBMid-WantTot) <= epsilon) {
      if (Verbose >= -1) {
        print(paste("0010ResortBayesSpike:GiveMeNewBalance: found alphaMid, LBMid which is ", LBMid,
         " is close to WantTot is ", WantTot, " and epsilon is ", epsilon, sep=""));
        print(paste("We Should Return a result now. using alphaMid = ", alphaMid, sep="")); flush.console();
        flush.console();
      }
      AVL <- dLL(lX, alphaMid);
      if (Verbose >= -1) {
        print("0010ResortBayesSpike:GiveMeNewBalance, we have AVL retrieved"); flush.console();
      }
      return(AVL);
      ##return(dLogit(alphaMid+lX));
    } else if (LBMid > WantTot) {
      alphaNew <- alphaMid;  LBNew <- LBMid;
    } else {
      alphaOld <- alphaMid;  LBOld <- LBMid;
    } 
  }
  if (Verbose >= -1) {
    print(paste("0010ResortBayesSpike:GiveMeNewBalance: Completing with alphaMid and lX found and tt = ", tt, " and MaxTTs is ", MaxTTs, sep="")); flush.console();
  }
  return(dLogit(alphaMid+lX));  
}