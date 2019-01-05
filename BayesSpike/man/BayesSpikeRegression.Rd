\name{BayesSpikeRegression}
\alias{BayesSpikeRegression}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ BayesSpikeRegression algorithm }
\description{
  The BayesSpike package for analyzing selection regression datasets via exclusionary spike sampling,
  BayesSpikeRegression is chief wrapper function.
  
  Warning, for successful use of \code{\link{BayesSpikeRegression}} the
  directory input \code{\link{SaveDir}} must be set to a user-accessible hard-drive
  location, see \code{\link{SaveDir}} for info on role of that data.
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
  EEProbSortWidth = -1, PrintIter = 200,DoLogitPostProb = FALSE,  
  DoLongCI = FALSE, CpBuffLongCI = -1, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{\eqn{\mathbf{X} \in \mathbb{R}^{n \times p}}{X in R[n * p]}
    Covariates, Independent Variables for both fixed and grouped paramters.  }
  \item{Y}{\eqn{\vec{Y} \in \mathbb{R}^{n}}{Y-vector in R^n} Dependent Variable. }
  \item{BetaStart}{Starting \eqn{\vec{\boldsymbol{\beta}}} draw.  }
  \item{IndexFirstRandomEffect}{Index in (1:NROW(X)) of the first "random effect", also \code{FirstRandom} }
  \item{\code{\link{tauEndList}}}{Ending Indices of each group of random effects (random effects groups must be contiguous)}
  \item{\code{\link{tauPriordf}}}{Degrees of freedom for \eqn{\tau_{j}^2}{tau_j^2} draws.}
  \item{\code{\link{tauPriorMean}}}{Mean value for prior of \eqn{\tau_{j}^2}{tau_j^2}.}
  \item{tauStartValues}{starting tau Values}
  \item{dfTauStart}{Degrees of freedom default to give to activated random effects groups inverse-gamma tau prior distribution}
  \item{mStart}{Mean default to give to activated random effects groups inverse-gamma tau prior distribution }
  \item{Verbose}{Toggle to deliver messages during algorithm run}
  \item{\code{\link{PiAPrior}}}{A Prior for PiA}
  \item{\code{\link{PriorProbTau}}}{Prior probability of inclusion optional extended user input beyond OnPiA for groups}
  \item{\code{\link{PriorProbFixed}}}{Prior probability of inclusion optional extended user input beyond OnPiA for solitary, ungrouped, coefficients}
  \item{\code{\link{SigmaPrior}}}{A prior for noise parameter OnSigma}
  \item{\code{\link{OnPiA}}}{Current \eqn{\pi_{\mbox{\tiny{$\mathcal{A}$}}}}{piA} parameter.}
  \item{\code{\link{OnSigma}}}{Current \eqn{\sigma^2}{sigma^2} noise paramter.}
  \item{\code{\link{dfRobit}}}{Set this to do Robit regression (Y entry must be 0-1)}
  \item{\code{\link{TemperatureList}}}{TemperatureList must be a set of descending temperatures (positive numbers greter than 1).
    Activating TemperatureList activates \link{EE-sampling} feature.  Only the lowest temperature Gibbs chain
    is valuable, but it features frequent mixing draws from higher temperature values.}
  \item{\code{\link{TypeFixedPrior}}}{Priors for Fixed parameters can take on multiple coding schemes.  0 is default, though
  exponential priors and uniform priors are possible as well.  See \code{\link{TypeFixedPrior}} for
  description of options.}
  \item{\code{\link{tauFixed}}}{Bayed upon integer flag \code{\link{TypeFixedPrior}},
    this double/vector/matrix should provide prior information on variances or bounds of the
    fixed paramters}
  \item{NumChains}{Number of chains to run.}
  \item{\code{\link{SaveDir}}}{Save information directory, see \code{\link{SaveDir}} for 
  description of BayesSpike saving buffer.}
  \item{\code{\link{FileName}}}{ Choose a unique header name for MCMC compressed files and RData}
  \item{\code{\link{FileSave}}}{ Additional File Save information.  }
  \item{CauchyEpsilon}{A Cauchy convergence parameter, useful in the \code{\link{SampleNewTaus}} operation}
  \item{MaxIters}{A Cauchy max iters parameter, useful in the \code{\link{SampleNewTaus}} operation}
  \item{EarlyEndStep}{Integer, if greater than zero sets a stopping test function, stops Gibbs sampler on a given step (samplign function) at time EarlyEndtt.}
  \item{EarlyEndtt}{Integer, if greater than zero setsiteration from gibbs sampler to stop.  Stop will occur
   after running EarlyEndStep.}
  \item{DoLongCI}{
    To save Long Intervals, this must be true in \code{\link{BayesSpikeRegression}} call.
    The "Long" Buffer is a storage of gibbs samples that are all non-zero for all p Beta coefficients
    This can only be saved on a file (Not stored in Coda List even if DoRecord[1] == 1).  This draws a nonzero
    Beta_j even if sampler itself drew a zero Beta_j.  This is used for robust confidence intervals for
    coefficients near zero.  These confidence intervals are overly wide, though they benefit from shrinkage on the other coefficients.
    See \code{\link{SubCodaList}} for more details, as well as \code{\link{BetaSymmetricLongConfidenceMatrix}}
    }
  \item{CpBuffLongCI}{Integer, How long should Long (No shrinkage to Zero) Buffer should be of size}
  \item{\code{\link{WritePostProbBuffer}}}{Whether to store Posterior Probabilities in Coda Chains}
  \item{burnin}{For calculations with MCMC chains, including drawing buffers, what initial position in chain to use as post-mixing time}
  \item{\code{\link{EEMergeEvery}}}{ How often to perform a EE merge step if necessary. }
  
}
\details{
 Analysis object is a "Rcpp modules" object or imbedded C++ class with multiple
   access vectors
   
   
  Note that recentering function is usually optional on data.  A matrix as described in
  BayesSpike Paper Section 3.1 is found at \code{\link{MakeAM}}
}
\value{
    Object analysis object is a "\code{\link[Rcpp]{Rcpp}} module" of time , with various properties and methods
     \item{CodaList}{Coda mcmc table list of output if generated}
     \item{\code{\link{OnBeta}}}{Current Beta Vector}
   \item{\link{PostProbBuffer}}{Buffer of information holding posterior probability information.}
     \item{\link{PostProbCodaList}}{Final Coda List of information holding posterior probability information.}
     \item{\code{\link{OnTau}}}{Current Draw of Tau random effects}
     \item{XtX}{Current Weighted part matrix of correlations}
     \item{XtY}{Current Weighted part vector of correlations}
     \item{XtResid}{Current XtResid vector or 
       \eqn{ \mathbf{X}^{T}( \vec{Y} - \mathbf{X} \boldsymbol{\beta}) }{ t(X) (Y-X*Beta)}
      useful in deciding posterior probabilities.  }
  \item{\code{\link{InitiateTime}}}{Time when BayesSpikeCL class was initiated}
  \item{\code{\link{SecondInitiateTime}}}{Time when GibbsSampler was started}
  \item{\code{\link{CompleteTime}}}{Time when Gibbs Sampler finished.}
  \item{\code{IamBayesSpike}}{Indicator that object is of BayesSpikeClass}
}
\references{ Diallel Paper, pending BayesSpike submission }
\author{ Alan Lenarcic, Will Valdar }
\note{ 
  "TBSRoo"  is a class using R.oo which we would like to deprecate in favor of 
  an R5 class
}
\seealso{ 
  \code{\link{BayesSpikeRegression}}, 
  \code{\link{RobitRegression}},
  \code{\link{MT}},
  \code{\link{OnPiA}},
  \code{\link{OnSigma}},
  \code{\link{PostProbCodaBuffer}},
  \code{\link{EESampling}},
  \code{\link{BayesSpikeFromAutoDiallel}}, 
  \code{\link[BayesDiallel]{BayesDiallel}}, 
  \code{\link[BayesDiallel]{DiallelAnalyzer}},
  \code{\link[coda]{mcmc}},
  \code{\link[Rcpp]{Rcpp}}, \code{\link[R.oo]{R.oo}}, 
  \code{\link[methods]{ReferenceClasses}},
  \code{\link{MakeAM}},
  \code{\link{Save}}, \code{\link{ReLoad}}
 }
\examples{
## Load in R packages
library(BayesSpike);
try(require(lattice, quietly=TRUE, warn.conflicts=FALSE), silent=TRUE);
try(require(coda, quietly=TRUE, warn.conflicts=FALSE), silent=TRUE);
try(require(R.methodsS3, quietly=TRUE, warn.conflicts=FALSE), silent=TRUE);
try(require(R.oo, quietly=TRUE, warn.conflicts=FALSE), silent=TRUE);
try(require(corpcor, quietly=TRUE, warn.conflicts=FALSE), silent=TRUE);
require(BayesDiallel, quietly = TRUE, warn.conflicts=FALSE);
require(R.oo, warn.conflicts=FALSE, quietly=TRUE);

################################################################################
## Not Usuable in sample package
##  try(library(AlanDirectories), silent=TRUE)
##  try(AlanDirectores:::SetSaveHomes(),silent=TRUE);
##  try(SetSaveHomes());
##  SaveDir = BSSaveDir;  ## Save Dir set by default, or set to ones own accessible harddrive space.
 ###############################################################################
 
 ###############################################################################
 ### Simulate Data 
  NewWrite = 1; dfRobit = 7; dfRobit = -1;
  nop = 4; NoTau = 0;
  if (NoTau == 0) {
    ## In this version the Random effects are simulated by "Non-Random" means
    FirstRandom =5 + nop+1; tauEndList = FirstRandom + c(7,14,22,28);
    CFirstRandom = FirstRandom -1;  CtauEndList = tauEndList -1;
    PiAPrior = c(5, nop, 10,10); PiAStart = c(5/(nop+5), 10/20); 
  } else {
    ##
    FirstRandom = -1;  tauEndList = NULL; CFirstRandom = -1; CtauEndList = NULL;
    PiAPrior = c(15, nop+10);  PiAStart = c(15 / (nop+25));
  }
  ## Set Noise Size at something relative to signal, Robit-ness, Beta
  if (dfRobit >= 0) {
    n = 500;
    Sigma = 3;
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
  }  else {
    YY = Z;
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
MyProcTime = tSeq(proc.time()[1]);
BSSaveDir  = paste(BSSaveDir, MyProcTime, sep="");  ## add a time stamp to analysis
  
AFD2 = DiallelAnalyzer(data = Piximus.Data, father.strain="father.strain.name",
   mother.strain="mother.strain.name", phenotype="MouseWeight", is.female="is.female",
   sep="", FixedEffects = NULL,
   RandomEffects = NULL, 
   Models=Example.Piximus.Models[1:2], sigmasq.start = 1,  numChains = 5,  lengthChains=2500,
   burnin = 1,
   DIC.Only=FALSE,  tauPriorFile = Example.Piximus.tau.Prior.Info,
   DoBayesSpike = TRUE, BSSaveDir = BSSaveDir);
   
## Get a summary of the Coda Chains for first Diallel Model
summary(AFD2$BS$CodaList)
plot(AFD2$BS$CodaList);
plot.hpd(AFD2$BS$CodaList);

SaveIn <- paste(BSAFD$SaveDir, "//AttemptSave//", sep="");
dir.create(SaveIn, showWarnings=FALSE);
try(BSAFD$save());

##try(save(BSAFD, paste(SaveIn, "//SaveBayesSpikeCLL.RData", sep="")));
##ARP <- BSAFD$TBSR5;
##myF <-  paste(SaveIn, "//TBSR5TrySave.RData", sep="")
##try(save(ARP, file = myF));


try(
  load(paste(BSAFD$TBSR5$SaveDir, "//AttemptSaveDir//PartSaveListTBSR5.RData"));
## This should attach a list object called "MyList" to the global environment
## To reload complete MyList object into a BayesSpike Object run  function
##                                 "ReLoad()"
##  
library("BayesSpike");
NewBayesSpikeLoadIn <- ReLoad(MyList);

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }
\keyword{ ts }% __ONLY ONE__ keyword per line
