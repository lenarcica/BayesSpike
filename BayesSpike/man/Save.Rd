\name{Save}
\alias{Save}
\alias{ssave}
\alias{save}
\alias{SSave}
\alias{Ssave}
\alias{reload}
\alias{Reload}
\alias{ReLoad}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ BayesSpike package }
\description{
  The BayesSpike package for analyzing selection regression datasets via exclusionary spike sampling
}
\usage{
  BS$save();
  
  load(paste(BSAFD$TBSR5$SaveDir, "//FinalOutputSave//PartSaveListTBSR5.RData"));
  NewBSObject <- ReLoad(MyList)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{BS}{BayesSpike Cpp Class}
  \item{RData File}{  BayesSpike automatically saves to its "SaveDir" directory
    a RData list into extention directory "AttemptSaveDir" into the element
    "PartSaveListTBSR5.RData"}
  \item{MyList}{This is the name of the list object saved by 
    PartSaveListTBSR5 Beta}
  \item{ReLoad}{Function that loads in MyList back to previous BayesSpike Object state}
}
\details{
  Uses object created for BayesSpike Analysis.

 Analysis object is a "Rcpp modules" object or imbedded C++ class with multiple
   access vectors
}
\value{
    Object analysis object is a "Rcpp module", with various properties and methods
     \item{CodaList}{Coda mcmc table list of output if generated}
     \item{Beta}{Current Beta Vector}
}
\references{ Diallel Paper, pending BayesSpike submission }
\author{ Alan Lenarcic, Will Valdar }
\note{ 
  "TBSRoo"  is a class using R.oo which we would like to deprecate in favor of 
  an R5 class
}
\seealso{ \code{\link{BayesSpikeRegression}}, \code{\link{BayesSpikeFromAutoDiallel}}, 
  \code{\link{BayesDiallel}}, \code{\link{coda}},
  \code{\link{Rcpp}}, 
  \code{\link{methods::ReferenceClasses}}, \code{\link{BaysSpike}}, \code{\link{Rcpp::Rcpp}}
  }
\examples{
library(BayesSpike)


NewWrite = 1;
NoTau = 0;
dfRobit = 7; nop = 4; dfTNoise = -1;
if (NoTau == 0) {
  FirstRandom =5 + nop+1; tauEndList = FirstRandom + c(7,14,22,28);
  CFirstRandom = FirstRandom -1;  CtauEndList = tauEndList -1;
  PiAPrior = c(5, nop, 10,10); PiAStart = c(5/(nop+5), 10/20); 
} else {
  FirstRandom = -1;  tauEndList = NULL; CFirstRandom = -1; CtauEndList = NULL;
  PiAPrior = c(15, nop+10);  PiAStart = c(15 / (nop+25));
}
## A group Beta
n = 200;
Sigma = 6;
Beta = 2*c(15,-15,0,15,15,rep(0,nop),-7,4,-9,12,-3,5,-4,.3, 0,0,0,0,0,0,0,4,-5,-6,7,-4,-2,1,.75,0,0,0,0,0,0);
Beta = Beta / 5;

## Generate Data
Verbose = 0;  InitKKs = 10;
pTot = length(Beta);
X = matrix(rnorm(n * pTot), n, pTot);
EY = X %*% Beta;
Y = EY  + rt(n, dfRobit);
Z = Y;  Z[Y > 0] = 1;  Z[Y < 0] = 0;
if (length(X[1,]) <= 400) {
XtY = t(X) %*% Y;
XtX = t(X) %*% X;
XtResid = XtY - XtX %*% Beta;
YResidSq = sum((Y-X %*% Beta)^2);  OutYResidSq = 0 * YResidSq;
}
##BetaStart = MyEM2Lasso$ReturnBetas;
BetaStart = Beta;


SaveDir = paste(BSSaveDir, "//n", n, "p", length(Beta), "tEL", paste(tauEndList, collapse=""),
  sep=""); dir.create(SaveDir, showWarnings=FALSE);

PriorProbTau = c(.2,.8,.01,.99);
PriorProbFixed = rep(.5, FirstRandom-1);
t1 = proc.time();
TemperatureList = c(8.0, 4.0,2.0,1.0,1.0,1.0,1.0);
##TemperatureList = c(1.0,1.0)
TemperatureList = NULL;
ABS = BayesSpikeRegression(Y=Y,X=X, BetaStart = BetaStart,
  IndexFirstRandomEffect = FirstRandom, tauEndList=tauEndList,
  tauStartValues = NULL, PiAStart = PiAStart, dfTauStart = -1, MStart = -1,
  dfTNoise = dfTNoise, Verbose = 1,
  DependenciesTau = NULL, DependenciesFixed=NULL,
  DoMax = 1, dfRobit = -1,
  PiAPrior = PiAPrior, HowSample = 3,
  NumSpliceSampleReps = 5, SpliceMove= 1, CauchyEpsilon = .00001,
  MyPS = PriorStructure(.5,1),  tauFixed = 800,
  TypeFixedPrior = 1, DoRNGState = 1,
  SigmaPrior = c(1,Sigma^2), InitKKs = 10, DoRecord = c(1,1,1,1,0, 1,1),
  NumSamples = 200, MaxIters = 200, MaxGibbsIters = 200,
  MaximizeMeCauchyTotalIters = 100,
  ttStart = 0,  NumChains = 4, FileSave = "MyS",
  SaveDir = SaveDir,  FileName = "ABSChain",
  NewWrite = 1, tauPriordf = 1, tauPriorMean = 1, 
  CodaTableNames = NULL, Run = TRUE,
  EarlyEndStep = -1, EarlyEndtt = -1, RStart = 0,
  StartRunProbVector = 100, PriorProbTau = PriorProbTau, 
  PriorProbFixed=PriorProbFixed,
  TemperatureList = NULL, DoEEProbSort = 1, EEMergeEvery = 10,
  EEProbSortWidth = .5);
ABS$save();

## Load in R packages
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
   mother.strain="mother.strain.name", 
   phenotype="MouseWeight", is.female="is.female",
   sep="", FixedEffects = NULL,
   RandomEffects = NULL, 
   Models=Example.Piximus.Models[1:2], sigmasq.start = 1, 
   numChains = 5,  lengthChains=500,
   burnin = 1,
   DIC.Only=FALSE,  tauPriorFile = Example.Piximus.tau.Prior.Info,
   DoBayesSpike = TRUE, BSSaveDir = BSSaveDir);
   
## Get a summary of the Coda Chains for first Diallel Model
summary(AFD2$BSAFD$CodaList)
plot(AFD2$BSAFD$CodaList);
plot.hpd(AFD2$BSAFD$CodaList);


SaveIn <- paste(BSAFD$SaveDir, "//FinalOutputSave//", sep="");
dir.create(SaveIn, showWarnings=FALSE);
try(BSAFD$save());

##try(save(BSAFD, paste(SaveIn, "//SaveBayesSpikeCLL.RData", sep="")));
##ARP <- BSAFD$TBSR5;
##myF <-  paste(SaveIn, "//TBSR5TrySave.RData", sep="")
##try(save(ARP, file = myF));


try(
  load(paste(BSAFD$TBSR5$SaveDir, "//FinalOutputSave//PartSaveListTBSR5.RData"));
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
