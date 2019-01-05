\name{ReturnSubCodaFromFiles}
\alias{ReturnSubCodaFromFiles}
\alias{ReturnSubcodafromfiles}
\alias{ReturnSubCodafromfiles}
\alias{ReturnSubCodaFromfiles}
\alias{GiveCounts}
\alias{GiveCountFrequency}
\alias{GiveCodaSubset}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Load in Coda from compressed files}
\description{
   
   Standalone \code{\link{ReturnSubCodaFromFiles}()} is used to read sub-sized
   CodaChains from files if \code{\link{SaveDir}} was set.
}
\usage{
 ReturnSubCodaFromFiles <- function(
    TBSR5 = NULL, TBSCpp = NULL, StartIter = 0, EndIter = 2000,
    Verbose = 0, SubSetCoords = NULL,TBSRoo = NULL,DoLogitPostProb = FALSE, 
    DrawBeta = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{TBSCpp, TBSR5}{\code{\link{BayesSpikeCpp}} or \code{\link{BayesSpikeR5}}
   return objects from \code{\link{BayesSpikeRegression}} which has save file information.
   }
  \item{StartIter}{First iteration of chain to start reading in draws (a burnin).}
  \item{EndIter}{Last Iteration to consider.}
  \item{Verbose}{Printback information during read process (useful for debugging.) }
  \item{SubSetCoords}{A sub set of the total coodinates \eqn{j}{j} 
  \eqn{\in \{1, ..., p\} }{in 1,...,p} which saves memory by only loading in 
  coordinates of interest.}
  \item{DrawBeta}{Logical Flag, if FALSE, do not pull Beta}
  \item{DoLogitPostProb}{When calculating posterior probabilities, if FALSE, returns probabilities from 0 to 1,
    else returns probabilities on Logit Scale}

}
\details{
   \code{BayesSpike} saves compressed binary files to the harddrive to represent
   draws from the Gibbs sampling distribution.  Since a Gibbs sampling chain for 
   large p is memory-intensive, this allows considerable savings on memory.
   
   BayesSpike does save information to a buffer, so only when the buffer is filled is a chain of 
   multiple samples appended to the current binary stream.  Binary files are compressed such that
   only active variables are listed.  One file, the "I" file, lists the integer identity of
   active variables for every gibbs sample run, this is used to translate the values 
   file, the "D" file that records the actual Beta drawn values.  Tau parameters
   get their own files as well.  A file of Rao-blackwellized posterior probability of activation
   parameters is also saved as a file.
   
   Note that a model might include 100K plus parameters, and one cannot load a complete 10,000 iteration by 100K parameter chain
     easily into memory.  Usual use cases only want to select the top 10 to 100 parameters that had large MIP probability
     in the chain and this code will retrieve a code matrix only with those parameters.
     
   Note that Beta parameters will be a different chain than Tau parameters.  
   There are also posterior model probabilities, posterior inclusion parameter chains, and draws for Y,w in logistic/Probit samples in
     other stored buffer chains if necessary.
   
   \code{GiveCodaSubset()} is the .Call C function used.
}
\value{
    BayesSpikeRegression serves as a wrapper function to write and read these files.
    Setting \code{NewWrite=0} will by default set \code{\link{BayesSpike}} to attempt to read currently written
    files rather than start new chains.  If the files are only half written, then BayesSpike
    will read in current chains and attempt to append until the chains are finished.
    If chains are fully finished, \code{\link{BayesSpikeRegression}} will simply read in the files.
    
       Standalone \code{\link{ReturnSubCodaFromFiles}()} is used to read sub-sized
   CodaChains from files.
   
   The elements returned by this function are
     \item{CodaBeta}{A Coda \link[coda]{mcmc.list} object of the Beta draws}
     \item{CodaTau}{A Coda \link[coda]{mcmc.list} object of Tau draws}
     \item{CodaP}{A list of total Posterior probabilities for all draws.}
     \item{CodaPostBuff}{A Coda \link[coda]{mcmc.list} object of Posterior Probabilities of inclusion for all draws.}
     
   Alternatively, one can access \code{\link{SubCodaList}}, \code{\link{TauCodalist}},
     \code{\link{PiACodaList}}, \code{\link{SigmaCodaList}}, \code{\link{ProbCodaList}},
     \code{\link{WeightCodaList}}, \code{\link{YCodaList}} from a \code{\link{BayesSpikeCpp}}
     or \code{\link{BayesSpikeR5}} object to produce a similar effect.
}
\references{ Diallel Paper, pending BayesSpike submission }
\author{ Alan Lenarcic, Will Valdar }
\note{ 
  "TBSRoo"  is a class using R.oo which we would like to deprecate in favor of 
  an R5 class, "\code{\link{TBSR5}}"  or "\code{\link{BayesSpikeR5}}" is  the R5 storage class.
}
\seealso{ \code{\link{BayesSpikeRegression}}, \code{\link{BayesSpikeFromAutoDiallel}}, 
  \code{\link{BayesDiallel}}, \code{\link{coda}},
  \code{\link{Rcpp}}, \code{\link{R.oo}}, \code{\link{methods::ReferenceClasses}},
  \code{\link{SaveDir}}, \code{\link{CodaList}}, \code{\link{SubCodaList}}
  }
\examples{
## Load in R packages
try(require(lattice, quietly=TRUE, warn.conflicts=FALSE), silent=TRUE);
try(require(coda, quietly=TRUE, warn.conflicts=FALSE), silent=TRUE);
try(require(R.methodsS3, quietly=TRUE, warn.conflicts=FALSE), silent=TRUE);
try(require(R.oo, quietly=TRUE, warn.conflicts=FALSE), silent=TRUE);
try(require(corpcor, quietly=TRUE, warn.conflicts=FALSE), silent=TRUE);
require(BayesDiallel, quietly = TRUE, warn.conflicts=FALSE);
require(R.oo, warn.conflicts=FALSE, quietly=TRUE);
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
