AlternateWeightSummary <- function(WantCodas, WeightCodas, burnin = 2,  
  quantilesGet = c(.01,.025,.05,.25,.5,.75,.95, .975,.99)) {
  if (!exists("burnin")) { burnin = 3; }
  if (!exists("quantilesGet")) { quantilesGet = c(.01,.025,.05,.25,.5,.75,.95, .975,.99); }
  ##if (burnin <= 1) {
  ##  print("AlternateWeightCI: Issue you set burnin = 1, this normally doesn't work");
  ##  flush.console();
  ##}  
  StackWantCodas <- AStackCoda(WantCodas, burnin);
  StackWeightCodas <- AStackCoda(WeightCodas, burnin);
  if (NROW(StackWeightCodas) != NROW(StackWantCodas)) {
    print(paste("AlternateWeightCI, uh oh, WantCodas had row ", NROW(WantCodas),
      " but WeightCoda is row ", NROW(WeightCodas), " with lengths ",
      "(", length(StackWantCodas), ", ", length(StackWeightCodas), ") respectively, not good.",
      sep="")); flush.console();
    return(NULL);
  }  
  if (sd(as.vector(StackWeightCodas)) >= 1.5) {
    StackWeightCodas[StackWeightCodas > quantile(StackWeightCodas,.99)] <- quantile(StackWeightCodas, .99);
    StackWeightCodas <- StackWeightCodas / sd(as.vector(StackWeightCodas)) * 1.5;
  }
  ReWeight <- exp(StackWeightCodas-max(StackWeightCodas));
  ReWeight <- ReWeight/ sum(ReWeight);
  if (NCOL(StackWantCodas) == 1) {
  Means = sum(as.vector(ReWeight) * StackWantCodas);
  SdQ <- sqrt(sum( as.vector(ReWeight)*t(t(StackWantCodas) - Means)^2));
  } else {
  Means = colSums(as.vector(ReWeight) * StackWantCodas);
  SdQ <- sqrt(colSums( as.vector(ReWeight)*t(t(StackWantCodas) - Means)^2));    
  }
  if (NCOL(StackWantCodas) == 1) {
    StackWantCodas <- matrix(StackWantCodas, length(StackWantCodas),1);
  }
  MyQuant <- QuantileMe(ReWeight, StackWantCodas, QuantileLocs = sort(quantilesGet), Verbose = 0);
  MyQuant <- t(MyQuant); ##colnames(MyQuant) <- paste(round(quantilesGet*100,2),"%", sep="");   
  MyHPD <- WeightHPD(ReWeight, StackWantCodas, quantilesGet, Verbose = 0)
  if (any(quantilesGet == .5)) {
     CN <- sort(quantilesGet[quantilesGet < .5]);
     if (NCOL(StackWantCodas) == 1) {
      MyHPD <-  c(MyHPD[1:(NROW(MyHPD)/2)], MyQuant[as.numeric(colnames(MyQuant)) == .5],
         MyHPD[(NROW(MyHPD)/2+1):(NROW(MyHPD))]);  
      MyHPD <- matrix(MyHPD, length(MyHPD),1);
      try(rownames(MyHPD) <- c(CN, .5, 1-CN[length(CN):1]));  
      try(colnames(MyHPD) <- colnames(WantCodas[[1]]));    
     } else {
       MyHPD <-  rbind(MyHPD[1:(NROW(MyHPD)/2),], MyQuant[,as.numeric(colnames(MyQuant)) == .5],
         MyHPD[(NROW(MyHPD)/2+1):(NROW(MyHPD)),]);
       try(rownames(MyHPD) <- c(CN, .5, 1-CN[length(CN):1]));  
       try(colnames(MyHPD) <- colnames(WantCodas[[1]])); 
     }
  }
  MyHPD <- t(MyHPD);
  MMT <- match(quantilesGet, as.numeric(colnames(MyQuant)));
  MyQuant <- MyQuant[,MMT, drop=FALSE];
  MMH <- match(quantilesGet, as.numeric(colnames(MyHPD)));
  MyHPD <- MyHPD[,MMT, drop=FALSE];  
  colnames(MyQuant) <- paste(quantilesGet*100, "%", sep="");
  colnames(MyHPD) <- paste(quantilesGet*100, "%", sep="");
  return(list(Means=Means, Sds = SdQ, MyQuant=MyQuant, MyHPD= MyHPD));
}
TakeSummary <- function(WantCodas, burnin = 2,  
  quantilesGet = c(.01,.025,.05,.25,.5,.75,.95,.975,.99)) {
  if (!exists("burnin")) { burnin = 3; }
  if (!exists("quantilesGet")) { quantilesGet = c(.01,.025,.05,.25,.5,.75,.95,.975,.99); }
  if (burnin <= 1) {
    ##print("AlternateWeightCI: Issue you set burnin = 1, this normally doesn't work");
    ##flush.console();
  }  
  StackWantCodas <- AStackCoda(WantCodas, burnin);

  if (NCOL(StackWantCodas)) {
  Means = colMeans(StackWantCodas);
  SdQ <- sqrt(colSums(t(t(StackWantCodas) - Means)^2)) / (NCOL(StackWantCodas)-1);
  } else {
  Means = sum(StackWantCodas);
  SdQ <- sqrt(sum(t(t(StackWantCodas) - Means)^2)) / (NCOL(StackWantCodas)-1);
  }
  ASum <- summary(as.mcmc(StackWantCodas), sort(quantilesGet));
  MyQuant <- ASum[[2]];  
  if (NCOL(StackWantCodas) == 1) {
    MyQuant <- matrix(ASum[[2]],1, length(ASum[[2]]));
  }
  colnames(MyQuant) <- sort(quantilesGet);
  MyHPD <- AllHPDs(as.mcmc(StackWantCodas), sort(quantilesGet))
  if (any(quantilesGet == .5)) {
    CN <- sort(quantilesGet[quantilesGet < .5]);
    if (NCOL(StackWantCodas) == 1) {
      MyHPD <- c(MyHPD[1:(NROW(MyHPD)/2)], MyQuant[as.numeric(colnames(MyQuant)) == .5],
        MyHPD[(NROW(MyHPD)/2+1):(NROW(MyHPD))]);
      MyHPD <- matrix(MyHPD, length(MyHPD),1);
      try(rownames(MyHPD) <- c(CN, .5, 1-CN[length(CN):1]));   
    } else {
      MyHPD <- rbind(MyHPD[1:(NROW(MyHPD)/2),], MyQuant[,as.numeric(colnames(MyQuant)) == .5],
       MyHPD[(NROW(MyHPD)/2+1):(NROW(MyHPD)),]);       
      try(rownames(MyHPD) <- c(CN, .5, 1-CN[length(CN):1]));   
    }
  } else {
    
  }
  
  MyHPD <- t(MyHPD);
  colnames(MyHPD) <- colnames(MyQuant);
  MMT <- match(quantilesGet, as.numeric(colnames(MyQuant)));
  MyQuant <- MyQuant[,MMT, drop=FALSE];
  MMH <- match(quantilesGet, as.numeric(colnames(MyHPD)));
  MyHPD <- MyHPD[,MMT, drop=FALSE];  
  colnames(MyQuant) <- paste(quantilesGet*100, "%", sep="");
  colnames(MyHPD) <- paste(quantilesGet*100, "%", sep="");
  return(list(Means=Means, Sds = SdQ, MyQuant=MyQuant, MyHPD= MyHPD));
}
AllHPDs <- function(StackWantCodas, quantilesGet) {
  quantilesLoc <- quantilesGet[quantilesGet < .5];
  HPDmatrix <- matrix(0, 2*length(quantilesLoc), NCOL(StackWantCodas));
  StackWantCodas <- as.mcmc(StackWantCodas);
  for (jj in 1:length(quantilesLoc)) {
    hGet <- HPDinterval(StackWantCodas, 1.0 - 2*quantilesLoc[jj]);
    HPDmatrix[jj,] <- hGet[,1];
    HPDmatrix[NROW(HPDmatrix)-jj+1,] <- hGet[,2];
  }
  colnames(HPDmatrix) <- colnames(StackWantCodas);
  rownames(HPDmatrix) <-c(quantilesLoc, 1.0-quantilesLoc[length(quantilesLoc):1]);
  return(HPDmatrix);
}
TakeSummaryZero <- function(WantCodas, burnin = 2,  
  quantilesGet = c(.01,.025,.05,.25,.5,.75,.95,.975,.99)) {
  if (!is.list(WantCodas) && length(dim(WantCodas)) == 2) {
    WeightCodas <- matrix(1.0, NROW(WantCodas), NCOL(WantCodas));
  } else {
    WeightCodas <- list();
    for (ii in 1:length(WantCodas)) {
      WeightCodas[[ii]] <- matrix(1.0, NROW(WantCodas[[ii]]), 1);  
    }
  }
  return(AlternateWeightZeroSummary(WantCodas, WeightCodas, burnin = 2,  
  quantilesGet = quantilesGet))  
}

AlternateWeightZeroSummary <- function(WantCodas, WeightCodas, burnin = 2,  
  quantilesGet = c(.01,.025,.05,.25,.5,.75,.95,.975,.99)) {
  if (!exists("burnin")) { burnin = 3; }
  if (!exists("quantilesGet")) { quantilesGet = c(.01,.025,.05,.25,.5,.75,.95, .975,.99); }
  if (burnin <= 1) {
    ##print("AlternateWeightCI: Issue you set burnin = 1, this normally doesn't work");
    ##flush.console();
  }  
  StackWantCodas <- AStackCoda(WantCodas, burnin);
  StackWeightCodas <- AStackCoda(WeightCodas, burnin);
  if (NROW(StackWeightCodas) != NROW(StackWantCodas)) {
    print(paste("AlternateWeightCI, uh oh, WantCodas had row ", NROW(StackWantCodas),
      " but WeightCoda is row ", NROW(StackWeightCodas), " with lengths ",
      "(", length(StackWantCodas), ", ", length(StackWeightCodas), ") respectively, not good.",
      sep="")); flush.console();
    return(NULL);
  }  
  if (sd(as.vector(StackWeightCodas)) >= 1.5) {
    StackWeightCodas[StackWeightCodas > quantile(StackWeightCodas,.99)] <- quantile(StackWeightCodas, .99);
    StackWeightCodas <- StackWeightCodas / sd(as.vector(StackWeightCodas)) * 1.5;
  }
  ReWeight <- exp(StackWeightCodas-max(StackWeightCodas));
  ReWeight <- ReWeight/ sum(ReWeight);
  if (NCOL(StackWantCodas) == 1) {
    Means = sum(as.vector(ReWeight) * StackWantCodas);
    SdQ <- sqrt(sum( as.vector(ReWeight)*t(t(StackWantCodas) - Means)^2));
  } else {
    Means = colSums(as.vector(ReWeight) * StackWantCodas);
    SdQ <- sqrt(colSums( as.vector(ReWeight)*t(t(StackWantCodas) - Means)^2));    
  }
  if (NCOL(StackWantCodas) == 1) {
    StackWantCodas <- matrix(StackWantCodas, length(StackWantCodas),1);
  }
  MyQuant <- QuantileMeZero(ReWeight, StackWantCodas, QuantileLocs = sort(quantilesGet), Verbose = 0)   
  MyQuant <- t(MyQuant); 
  MyHPD <- WeightHPDZero(ReWeight, StackWantCodas, quantilesGet, Verbose = 0)
  if (any(quantilesGet == .5)) {
     CN <- sort(quantilesGet[quantilesGet < .5]);
     if (NCOL(StackWantCodas) == 1) {
      MyHPD <-  c(MyHPD[1:(NROW(MyHPD)/2)], MyQuant[as.numeric(colnames(MyQuant)) == .5],
         MyHPD[(NROW(MyHPD)/2+1):(NROW(MyHPD))]);  
      MyHPD <- matrix(MyHPD, length(MyHPD),1);
      try(rownames(MyHPD) <- c(CN, .5, 1-CN[length(CN):1]));   
     } else {
       MyHPD <-  rbind(MyHPD[1:(NROW(MyHPD)/2),], MyQuant[,as.numeric(colnames(MyQuant)) == .5],
         MyHPD[(NROW(MyHPD)/2+1):(NROW(MyHPD)),]);
       try(rownames(MyHPD) <- c(CN, .5, 1-CN[length(CN):1]));   
     }
  }
  MyHPD <- t(MyHPD);
  MMT <- match(quantilesGet, as.numeric(colnames(MyQuant)));
  MyQuant <- MyQuant[,MMT, drop=FALSE];
  MMH <- match(quantilesGet, as.numeric(colnames(MyHPD)));
  MyHPD <- MyHPD[,MMT, drop=FALSE];  
  colnames(MyQuant) <- paste(quantilesGet*100, "%", sep="");
  colnames(MyHPD) <- paste(quantilesGet*100, "%", sep="");
  return(list(Means=Means, Sds = SdQ, MyQuant=MyQuant, MyHPD= MyHPD));
}

AStackCoda <- function(ACoda, MyBurn, KeepCols = NULL) {
  if (!is.list(ACoda) && length(dim(ACoda)) == 2) {
    return(as.mcmc(ACoda[MyBurn:NROW(ACoda),]));
  }
  if (!is.list(ACoda) && length(dim(ACoda)) == 1) {
    AM <- matrix(ACoda[MyBurn:NROW(ACoda)], length(MyBurn:NROW(ACoda)),1);
    try(rownames(AM) <- names(ACoda)[MyBurn:NROW(ACoda)]);
    return(as.mcmc(AM));
  }
  if (NCOL(ACoda[[1]]) == 1) {
    NeedCoda <- matrix(0, (NROW(ACoda[[1]])-MyBurn+1)*length(ACoda), 1);  colnames(NeedCoda) <- colnames(ACoda[[1]]);
    colnames(NeedCoda) <- colnames(ACoda[[1]]);
    AS <- 0;
    for (iii in 1:length(ACoda)) {
      ASL <- NROW(ACoda[[1]])-MyBurn+1;
      NeedCoda[AS+1:ASL,1] <-ACoda[[iii]][MyBurn:NROW(ACoda[[iii]]),1];
      AS <- AS + ASL;
    }
    library(coda);
    NeedCoda <- as.mcmc(NeedCoda);    return(NeedCoda);
  } else if (is.null(KeepCols)) {
  NeedCoda <- matrix(0, (NROW(ACoda[[1]])-MyBurn+1)*length(ACoda), NCOL(ACoda[[1]]));
  colnames(NeedCoda) <- colnames(ACoda[[1]]);
  AS <- 0;
  for (iii in 1:length(ACoda)) {
    ASL <- NROW(ACoda[[1]])-MyBurn+1;
    NeedCoda[AS+1:ASL,] <-ACoda[[iii]][MyBurn:NROW(ACoda[[iii]]),];
    AS <- AS + ASL;
  }
  library(coda);
  NeedCoda <- as.mcmc(NeedCoda);    return(NeedCoda);
  } else {
     NeedCoda <- matrix(0, (NROW(ACoda[[1]])-MyBurn+1)*length(ACoda),length(KeepCols));
     colnames(NeedCoda) <- colnames(ACoda[[1]])[KeepCols];
     AS <- 0;
     for (iii in 1:length(ACoda)) {
       ASL <- NROW(ACoda[[1]])-MyBurn+1;
       NeedCoda[AS+1:ASL,] <-ACoda[[iii]][MyBurn:NROW(ACoda[[iii]]),KeepCols];
       AS <- AS + ASL;
    }
    library(coda);
    NeedCoda <- as.mcmc(NeedCoda);
  }
  return(NeedCoda);
}

WeightHPD <- function(Weights, StackWantCodas, QuantileLocs, Verbose = 0, QuantileCoda = NULL) {
  QuantileLocs <- QuantileLocs[QuantileLocs >= 0 & QuantileLocs < .5];
  QuantileLocs <- sort(QuantileLocs);
  QuantileMatrix <- matrix(0, length(QuantileLocs)*2, NCOL(StackWantCodas));


  rownames(QuantileMatrix) <- c(QuantileLocs, 1.0-sort(QuantileLocs, decreasing=TRUE));

  Weights <- Weights / sum(Weights);
  colnames(QuantileMatrix) <- colnames(StackWantCodas);
  for (ii in 1:NCOL(StackWantCodas)) {
    MySt <- sort.int(StackWantCodas[,ii], index=TRUE)$ix;
    CumSumWeights <- cumsum(Weights[MySt]);
    SortX <-  StackWantCodas[MySt, ii];
    if (Verbose >= 1 && ii %% 10^(5-min(round(Verbose),5)) == 0) {
      print(paste(" -- Weighted Quantile on column ", ii, ", ", NCOL(StackWantCodas), sep="")); flush.console();
    }
    for (jj in 1:length(QuantileLocs)) {
       GetMinMax <- .Call("CallWeightedHPD", SortX, CumSumWeights, 1.0 - 2*QuantileLocs[jj]); 
       QuantileMatrix[jj,ii] <- SortX[GetMinMax[1]]; 
       QuantileMatrix[NROW(QuantileMatrix) - jj + 1,ii] <- SortX[GetMinMax[2]];
    }
  }
  return(QuantileMatrix);
}
WeightHPDZero <- function(Weights, StackWantCodas, QuantileLocs, Verbose = 0) {
  QuantileLocs <- QuantileLocs[QuantileLocs >= 0 & QuantileLocs < .5];
  QuantileLocs <- sort(QuantileLocs);
  QuantileMatrix <- matrix(0, length(QuantileLocs)*2, NCOL(StackWantCodas));


  rownames(QuantileMatrix) <- c(QuantileLocs, 1.0-sort(QuantileLocs, decreasing=TRUE));

  Weights <- Weights / sum(Weights);
  colnames(QuantileMatrix) <- colnames(StackWantCodas);
  for (ii in 1:NCOL(StackWantCodas)) {
    SKeep <- (1:NROW(StackWantCodas))[StackWantCodas[,ii] != 0.0];
    WeightKeep <- Weights[SKeep];  Xs <- StackWantCodas[SKeep,ii];
    MySt <- sort.int(Xs, index=TRUE)$ix;
    CumSumWeights <- cumsum(WeightKeep[MySt]) / sum(WeightKeep[MySt]);
    SortX <-  Xs[MySt];
    if (Verbose >= 1 && ii %% 10^(5-min(round(Verbose),5)) == 0) {
      print(paste(" -- Weighted Quantile on column ", ii, ", ", NCOL(StackWantCodas), sep="")); flush.console();
    }
    if (length(MySt) <= 1) {
      for (jj in 1:length(QuantileLocs)) {
         QuantileMatrix[jj,ii] <- min(SortX,0);
         QuantileMatrix[NROW(QuantileMatrix)-jj+1,ii] <-max(SortX,0);
      }
    } else {
      for (jj in 1:length(QuantileLocs)) {
         GetMinMax <- .Call("CallWeightedHPD", SortX, CumSumWeights, 1.0 - 2*QuantileLocs[jj]); 
         QuantileMatrix[jj,ii] <- SortX[GetMinMax[1]]; 
         QuantileMatrix[NROW(QuantileMatrix) - jj + 1,ii] <- SortX[GetMinMax[2]];
      }
    }
  }
  return(QuantileMatrix);
}

QuantileMe <- function(Weights, StackWantCodas, QuantileLocs, Verbose = 0) {
  QuantileMatrix <- matrix(0, length(QuantileLocs), NCOL(StackWantCodas));
  QuantileLocs <- sort(QuantileLocs);
  rownames(QuantileMatrix) <- QuantileLocs;
  colnames(QuantileMatrix) <- colnames(StackWantCodas);
  for (ii in 1:NCOL(StackWantCodas)) {
    MySt <- sort.int(StackWantCodas[,ii], index=TRUE)$ix;
    CumSumWeights <- cumsum(Weights[MySt]);
    tt <- 1
    if (Verbose >= 1 && ii %% 10^(5-min(round(Verbose),5)) == 0) {
      print(paste(" -- Weighted Quantile on column ", ii, ", ", NCOL(StackWantCodas), sep="")); flush.console();
    }
    for (jj in 1:length(QuantileLocs)) {
       MyMin <- sort(abs(CumSumWeights[tt:length(CumSumWeights)]-QuantileLocs[jj]),
         index=TRUE)$ix[1]
       QuantileMatrix[jj,ii] <- StackWantCodas[MySt[MyMin+tt-1],ii]; 
       tt <- MyMin;   
    }
  }
  return(QuantileMatrix)
}

QuantileMeZero <- function(Weights, StackWantCodas, QuantileLocs, Verbose = 0) {
  QuantileMatrix <- matrix(0, length(QuantileLocs), NCOL(StackWantCodas));
  QuantileLocs <- sort(QuantileLocs);
  rownames(QuantileMatrix) <- QuantileLocs;
  colnames(QuantileMatrix) <- colnames(StackWantCodas);
  for (ii in 1:NCOL(StackWantCodas)) {
    WantI <- (1:NROW(StackWantCodas))[StackWantCodas[,ii] != 0.0];
    MyX <- StackWantCodas[WantI,ii]; MyWeights <- Weights[WantI];
    MySt <- sort.int(MyX, index=TRUE)$ix;
    CumSumWeights <- cumsum(MyWeights[MySt]) / sum(MyWeights[MySt]);
    tt <- 1
    if (Verbose >= 1 && ii %% 10^(5-min(round(Verbose),5)) == 0) {
      print(paste(" -- Weighted Quantile on column ", ii, ", ", NCOL(StackWantCodas), sep="")); flush.console();
    }
    if (length(MyX) <= 1) {
      for (jj in 1:length(QuantileLocs)) {
        if (QuantileLocs[jj] <= .5) {
          QuantileMatrix[jj,ii] <- min(MyX,0); 
        } else {
          QuantileMatrix[jj,ii] <- max(MyX,0);
        }
      }
    } else {
    for (jj in 1:length(QuantileLocs)) {
       MyMin <- sort(abs(CumSumWeights[tt:length(CumSumWeights)]-QuantileLocs[jj]),
         index=TRUE)$ix[1]
       QuantileMatrix[jj,ii] <- MyX[MySt[MyMin+tt-1]]; 
       tt <- MyMin;   
    }
    }
  }
  return(QuantileMatrix)
}