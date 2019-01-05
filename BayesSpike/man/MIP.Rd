\name{MIP}
\alias{MIP}
\alias{mip}
\alias{postprobmip}
\alias{ModelInclusion}
\alias{ModelInclusionProbabilities}
\alias{PosteriorMIP}
\alias{posteriormip}
\alias{RunProbVector}
\alias{StartRunProbVector}
\alias{TotalRunProbVector}
\alias{DoShortMIP}
\alias{runprobvector}
\alias{Runprobvector}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{ "Model Inclusion Probabilities", Average MIP, BayesSpikeCpp \code{\link[Rcpp]{Rcpp}} Module Class}
\description{
  The posterior probability of inclusion averages, totalled through benefit from \code{\link{PostProbCodaBuffer}} .
  
  For all Fixed-Effects Beta/Tau parameters this is an average value from 0.0 to 1.0 expressing the posterior
  probability that this effect, or group of effects belongs within the buffer.
 
}
\usage{
BayesSpikeRegression <- function(X=NULL,Y=NULL, BetaStart = NULL,
  IndexFirstRandomEffect = 1, tauEndList=NULL,
  tauStartValues = NULL, PiAStart = .5, dfTauStart = -1, MStart = -1,
  dfTNoise = 0, Verbose = 0,
  DependenciesTau = NULL, DependenciesFixed=NULL,
  DoMax = 1, 
  PiAPrior = c(-1,-1), HowSample = 3,
  NumSpliceSampleReps = 5, SpliceMove= 1, CauchyEpsilon = .00001,
  MyPS = BayesSpike:::PriorStructure(.5,1),  tauFixed = 40,
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
  NoNoiseBetaStart = FALSE, AFD = NULL, NoShrinkVariables = NULL,
  NoShrinkTau =NULL, StartRunProbVector = -100,
  PriorProbFixed=NULL, PriorProbTau = NULL, dfRobit = -1, TemperatureList = NULL,
  burnin = NULL, EEMergeEvery = NULL, DoEEProbSort = 0, 
  EEProbSortWidth = -1, PrintIter = 200, ...)
  
  MBS$MIP
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{MIP}{Get this field from the BayesSpikeCL object after running the chains }
\item{RunProbVector}{Running total used to create MIP without loading draws from disk.  Accessible from MBS}
\item{StartRunProbVector}{Which iteration of chain to start saving to RunPRobVector}
\item{TotInProbVector}{Running total of how many draws populate the RunProbVector}
\item{DoShortMIP}{Flag whether the regression is performing a short MIP or gathering draws from saved to disk chains.}
}
\details{
 These buffers are updated during the chain using Posterior Model Inclusion, as opposed to just the average of (0,1's) for
 events when a coordinate is included or excluded.  Thus most parameters will have non-zero MIP, but they might still be low.
 
 Often performance from the sampler is good enough to simply run enough iterations to populate a stable MIP element, and then 
 run regression and average Bayesian mean only for parameters above a given desired MIP threshold.  If a Beta(1,p) prior is used
 then MIP larger than (1/p) after the regression could at least be argued to increase in inclusion probabilities.  MIP greater than \code{\link{OnPiA}} could be considered valable. 
}
\value{
    \item{MIP}{Get this vector from BayesSpikeCL Object}
}
\references{ Diallel Paper, pending BayesSpike submission }
\author{ Alan Lenarcic, Will Valdar }
\note{ 
  "TBSR5"  or \link{BayesSpikeR5} is a class using \code{\link[methods]{ReferenceClasses}} 
  but "TBSCpp" or \link{BayesSpikeCL} is a class using \link[Rcpp]{Rcpp} which 
  is the primary return object.  
}
\seealso{ \code{\link{BayesSpikeRegression}}, \code{\link{BayesSpikeFromAutoDiallel}},\code{link{Codajj}}, 
  \code{\link{BayesDiallel}}, \code{\link{coda}},
  \code{\link{Rcpp}}, \code{\link{R.oo}}, \code{\link{methods::ReferenceClasses}}
  }
\examples{
## Load in R packages
library(BayesSpike)
## Load in R packages
library(BayesSpike);
try(require(lattice, quietly=TRUE, warn.conflicts=FALSE), silent=TRUE);
try(require(coda, quietly=TRUE, warn.conflicts=FALSE), silent=TRUE);
try(require(R.methodsS3, quietly=TRUE, warn.conflicts=FALSE), silent=TRUE);
try(require(R.oo, quietly=TRUE, warn.conflicts=FALSE), silent=TRUE);
try(require(corpcor, quietly=TRUE, warn.conflicts=FALSE), silent=TRUE);
require(BayesDiallel, quietly = TRUE, warn.conflicts=FALSE);
require(R.oo, warn.conflicts=FALSE, quietly=TRUE);

 try(library(AlanDirectories), silent=TRUE)
 try(AlanDirectores:::SetSaveHomes(),silent=TRUE);
 try(SetSaveHomes());
 SaveDir = BSSaveDir;  ## Save Dir set by default, or set to ones own accessible harddrive space.
  
  NewWrite = 1; dfRobit = 7; dfRobit = -1;
  nop = 4; NoTau = 0;
  if (NoTau == 0) {
    FirstRandom =5 + nop+1; tauEndList = FirstRandom + c(7,14,22,28);
    CFirstRandom = FirstRandom -1;  CtauEndList = tauEndList -1;
    PiAPrior = c(5, nop, 10,10); PiAStart = c(5/(nop+5), 10/20); 
  } else {
    FirstRandom = -1;  tauEndList = NULL; CFirstRandom = -1; CtauEndList = NULL;
    PiAPrior = c(15, nop+10);  PiAStart = c(15 / (nop+25));
  }
  if (dfRobit >= 0) {
    n = 200;
    Sigma = 6;
  } else {
    n= 60;
    Sigma = 6;
  }
  dfTNoise = -1; ##dfTNoise = 30;  #dfTNoise = 0;
  Beta = 2*c(15,-15,0,15,15,rep(0,nop),-7,4,-9,12,-3,5,-4,.3, 0,0,0,0,0,0,0,4,-5,-6,7,-4,-2,1,.75,0,0,0,0,0,0);
  Beta = Beta / 5;

  SaveDir = paste(SaveDir, "//n", n, "p", length(Beta),  "dfR", dfRobit,
    "tEL", paste(tauEndList, collapse=""),
     sep=""); dir.create(SaveDir, showWarnings=FALSE);
     
  Verbose = 0;  InitKKs = 10;
  pTot = length(Beta);
  X = matrix(rnorm(n * pTot), n, pTot);
  EY = X %*% Beta;
  if (!is.null(dfRobit) && dfRobit == 0) {
    Y = EY + rnorm(n);  YY = Y;
    Z = Y; Z[Y > 0] = 1;  Z[Y < 0] = 0;
  } else if (is.null(dfRobit) || dfRobit < 0) {
    Y = EY + rnorm(n);  YY = Y;
  } else {
    Y = EY  + rt(n, dfRobit);
    Z = Y;  Z[Y > 0] = 1;  Z[Y < 0] = 0;
  }
  if (length(X[1,]) <= 400) {
    XtY = t(X) %*% Y;
    XtX = t(X) %*% X;
    XtResid = XtY - XtX %*% Beta;
    YResidSq = sum((Y-X %*% Beta)^2);  OutYResidSq = 0 * YResidSq;
  }
  
  if (is.null(dfRobit) || dfRobit < 0) {
    dfRobit = NULL;  YY = Y
  }
  
  ##BetaStart = MyEM2Lasso$ReturnBetas;
BetaStart = Beta;

PriorProbTau = c(.2,.8,.01,.99);
PriorProbFixed = rep(.5, FirstRandom-1);
t1 = proc.time();
TemperatureList = c(8.0, 4.0,2.0,1.0,1.0,1.0,1.0);
##TemperatureList = c(1.0,1.0)
TemperatureList = NULL;
ABS = BayesSpikeRegression(Y=YY,X=X, BetaStart = BetaStart,
  IndexFirstRandomEffect = FirstRandom, tauEndList=tauEndList,
  tauStartValues = NULL, PiAStart = PiAStart, dfTauStart = -1, MStart = -1,
  dfTNoise = dfTNoise, Verbose = 2,
  DependenciesTau = NULL, DependenciesFixed=NULL,
  DoMax = 1, dfRobit = dfRobit,
  PiAPrior = PiAPrior, HowSample = 3,
  NumSpliceSampleReps = 5, SpliceMove= 1, CauchyEpsilon = .00001,
  MyPS = PriorStructureR5$new(.5,1),  tauFixed = 800,
  TypeFixedPrior = 1, DoRNGState = 1,
  SigmaPrior = c(1,Sigma^2), InitKKs = 10, DoRecord = c(1,1,1,1,0, 1,1),
  NumSamples = 200, MaxIters = 200, MaxGibbsIters = 500,
  MaximizeMeCauchyTotalIters = 100,
  ttStart = 0,  NumChains = 4, FileSave = "MyS",
  SaveDir = SaveDir,  FileName = "ABSChain",
  NewWrite = 1, tauPriordf = 1, tauPriorMean = 1, 
  CodaTableNames = NULL, Run = TRUE,
  EarlyEndStep = -1, EarlyEndtt = -1, RStart = 0,
  StartRunProbVector = 100, PriorProbTau = PriorProbTau, 
  PriorProbFixed=PriorProbFixed,
  TemperatureList = TemperatureList, DoEEProbSort = 1, EEMergeEvery = 10,
  EEProbSortWidth = .5, DoLogitPostProb = FALSE);

 ## The following Digs variables from Saved Files.  It is similar to using
 ##  "CodaList", if that was not set up for this analysis.  
 ##  "CodaBeta" is a list of Beta coefficients
 ##  "CodaTau" is a list of Tau coefficients
 ##  "CodaP" is a list of Posterior likelihood probabilities of draws
 ##  "CodaPostBuff" is a list of Posterior Probability of activation of all variables
 ##  If SubSetCoords is set as a list of integers  that is a subset
 ##   of the coordinates of variables used in this sample, this will only
 ##   pull from files Beta coordinates of interest.
 ARSCF <- ReturnSubCodaFromFiles(
    TBSR5 = ABS$TBSR5, StartIter = 0, EndIter = 500,
    Verbose = 0, SubSetCoords = NULL, DoLogitPostProb = FALSE,
    DrawBeta = TRUE) 

 ## This will draw Posterior probability of Coda Likelihood.  
 plot(ABS$PostProbCodaList);
 
 ABS$PostProbBufferFile;  ## list last Buffer File saved to, if saved.
 
   
try(require(lattice, quietly=TRUE, warn.conflicts=FALSE), silent=TRUE);
try(require(coda, quietly=TRUE, warn.conflicts=FALSE), silent=TRUE);
try(require(R.methodsS3, quietly=TRUE, warn.conflicts=FALSE), silent=TRUE);
try(require(R.oo, quietly=TRUE, warn.conflicts=FALSE), silent=TRUE);
try(require(corpcor, quietly=TRUE, warn.conflicts=FALSE), silent=TRUE);
require(BayesDiallel, quietly = TRUE, warn.conflicts=FALSE);
require(R.oo, warn.conflicts=FALSE, quietly=TRUE);

###############################################################################
## First Example Piximus Diallel Mouse Data.
data(PiximusData);
##  Data is in "Piximus.Data"
##  Models are in "Examples.Piximus.Models"
##  tau.Prior is in "Examples.Piximus.tau.Prior.Info"
####################################
###### Run the Algorithm->
AFD = DiallelAnalyzer(data = Piximus.Data, father.strain="father.strain.name",
   mother.strain="mother.strain.name", phenotype="MouseWeight", is.female="is.female",
   sep="", FixedEffects = NULL,
   RandomEffects = NULL, 
   Models=Example.Piximus.Models, sigmasq.start = 1,  numChains = 5,  lengthChains=2500,
   burnin = 1,
   DIC.Only=FALSE,  tauPriorFile = Example.Piximus.tau.Prior.Info,
   SaveAFDFile = "SaveAFDBackUpPiximusData.RData");

## If SaveAFDFile saved a copy of AFD
if (FALSE) {
  library(R.oo);  library(BayesDiallel);
  load("SaveAFDBackUpPiximusData.RData");
}
AFD$ListElements;  ## List of Elements available to AFD
AFD$ListFunctions; ## List of Functions available to AFD  

AO = AFD$AllDiallelObs[[1]];
AO$ListElements  ## List of Elements available to AFD

## Chains are the first version, they are preshrunk towards group means.  You should
##   often use these by default;
summary(AO$cent.chains);   ## summary of the chains;
plot(AO$cent.chains);      ## Use Coda MCMC to plot the many member of chains consider "ask"/"wanted"
plot.hpd(AO$cent.chains);  ## Plot Hpd intervals for members of "chains" object.  Consider "ask"/"wanted"
summary(AO$raw.chains);   ## summary of the chains;
plot(AO$raw.chains);      ## Use Coda MCMC to plot the many member of chains consider "ask"/"wanted"
plot.hpd(AO$raw.chains);  ## Plot Hpd intervals for members of "chains" object.  Consider "ask"/"wanted"
  
## raw.chains are the "unshrunk version"
## Get a summary of the Coda Chains for first Diallel Model
summary(AFD$AllDiallelObs[[1]]$raw.chains)
plot(AFD$AllDiallelObs[[1]]$raw.chains);
plot.hpd(AFD$AllDiallelObs[[1]]$raw.chains);

##########################################################################
##  BayesSpike Selection uses Prior distributions to select a valid model from
##   subset of other models.  "tauPriorFile" can have an extra row "Prob"
##   with probabilities, if one wants to select certain groups with different
##   prior probabilities.  Default "Prob" will be .5
##
##  See library(BayesSpike) for more details on this selection process.
##
## Experimental, might not work if you don't have directory write access
if (length(list.files("C:")) > 0) {
  try(dir.create("C:/Stat"));  try(dir.create("C:/Stat/BayesSpikeSaves"));
  BSSaveDir = "C:/Stat/BayesSpikeSaves"
} else if (length(list.files("~/Documents")) > 0) {
  try(dir.create("~/Documents/Stat"));  try(dir.create("~/Documents/Stat/BayesSpikeSaves"));
  BSSaveDir = "~/Documents/Stat/BayesSpikeSaves"
} else if (length(list.files("~/")) > 0) {
 try(dir.create("~/Stat"));  try(dir.create("~/Stat/BayesSpikeSaves"));
  BSSaveDir = "~/Stat/BayesSpikeSaves"
}
AFD2 = DiallelAnalyzer(data = Piximus.Data, father.strain="father.strain.name",
   mother.strain="mother.strain.name", phenotype="MouseWeight", is.female="is.female",
   sep="", FixedEffects = NULL,
   RandomEffects = NULL, 
   Models=Example.Piximus.Models, sigmasq.start = 1,  numChains = 5,  lengthChains=2500,
   burnin = 1,
   DIC.Only=FALSE,  tauPriorFile = Example.Piximus.tau.Prior.Info,
   DoBayesSpike = TRUE, BSSaveDir = BSSaveDir);
   
## Get a summary of the Coda Chains for first Diallel Model
summary(AFD2$BSAFD$CodaList)
plot(AFD2$BSAFD$CodaList);
plot.hpd(AFD2$BSAFD$CodaList);

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }
\keyword{ ts }% __ONLY ONE__ keyword per line
