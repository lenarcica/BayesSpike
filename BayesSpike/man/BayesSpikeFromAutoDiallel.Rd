\name{BayesSpikeFromAutoDiallel}
\alias{BayesSpikeFromAutoDiallel}
\alias{AFDBayesSpike}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ BayesSpikeFromAutoDiallel, converting BayesDiallel, "FullDiallelAnalyze" object to use of BayesSpike algorithm }
\description{
  The BayesSpike package for analyzing selection regression datasets via exclusionary spike sampling,
  BayesSpikeRegression is chief wrapper function.  The Pakcage 
   \code{\link[BayesDiallel]{BayesDiallel}}l is an analysis method for 
  data generated from a "Diallel" or "F1-of-inbred-strains-single-phenotype" experiment.
  This function \code{BayesSpikeFromAutoDiallel} takes an AFD and converts it so that a BayesSpike Regression
  is run on that data.  
  
   See \code{\link[BayesDiallel]{BayesSpikeDiallel}}
  for information within package \code{\link[BayesDiallel]{BayesDiallel}}
   on how to used \code{\link[BayesDiallel]{FullDiallelAnalyze}} objects
   and \code{\link[BayesDiallel]{DiallelAnalyzer}} wrapper function to their fullest.
}
\usage{
       
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
  SigmaPrior = c(2,1), InitKKs = 5, DoRecord = rep(0,7),
  NumSamples = 1000, MaxIters = 1000, MaxGibbsIters = 1000,
  MaximizeMeCauchyTotalIters = 100,
  ttStart = 0, 
  EarlyEndStep, EarlyEndtt, FileSave = "MyS",
  SaveDir = NULL,  FileName = "AutoDiallelChain",
  NewWrite = 1, tauPriorMean = -1, tauPriordf = -1,
  WatchCoordinates = NULL, Run = FALSE,
  PriorProbTau = NULL, PriorProbFixed = NULL) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  The key input is \code{AFD} which is a type \code{\link[BayesDialle]{FullDiallelAnalyze}}
   object garnered after running \code{\link[BayesDiallel]{DiallelAnalyzer}}.  Using this
   object, \code{\link{BayesSpikeFromAutoDiallel}} culls information from \code{AFD}
   so that it can be added to \code{\link{BayesSpikeRegression}}.  
   Entry Elements to \code{\link{BayesSpikeFromAutoDiallel}} are thereafter similar
   to  \code{\link{BayesSpikeRegression}}, accepting that if information can be found
   in \code{AFD}, it will be used first.
   
  \item{SaveDir}{A Directory to which harddrive storage can be made for algorithm} 
  \item{X}{Covariates, Independent Variable}
  \item{Y}{Dependent Variable}
  \item{BetaStart}{Starting Beta}
  \item{IndexFirstRandomEffect}{Index in (1:NROW(X)) of the first "random effect"}
  \item{tauEndList}{Ending Indices of each group of random effects (random effects groups must be contiguous)}
  \item{tauStartValues}{starting tau Values}
  \item{PiAStart}{Value for Prior probability of Activation}
  \item{dfTauStart}{Degrees of freedom default to give to activated random effects groups inverse-gamma tau prior distribution}
  \item{mStart}{Mean default to give to activated random effects groups inverse-gamma tau prior distribution }
  \item{Verbose}{Toggle to deliver messages during algorithm run}
  \item{PiAPrior}{A Prior for PiA}
}
\details{
 Analysis object is a \code{\link[Rcpp]{Rcpp}} "Rcpp modules" object or imbedded C++ class with multiple
   access vectors of type \code{\link{BayesSpike}}.  
   
 This will return  an object of type  \code{\link{BayesSpikeCpp}} which performed the analysis.
 This object inside it will have a partial \code{\link[BayesDiallel]{FullDiallelAnalyze}}
 object which contains relevant copy information. 
 
 It is recommended however to simply run wrapper function \code{\link[BayesDiallel]{DiallelAnalyzer}}
 with \code{\link[BayesDiallel]{DoBayesSpike}} turned on.  This will run both selection and non-selection
 analyses, returning a single \code{\link[BayesDiallel]{FullDiallelAnalyze}} object
 to which the analysis \code{\link{BayesSpikeCpp}} analysis is attached.
 
   Because Diallels are small, and parameter lists are likely in the hundreds but 
   not in the thousands, saving and running \code{\link{CodaList}}
  is acceptable, and hence most analyses of the output will be based upon
  default large Coda-package \code{\link[coda]{mcmc}} objects.  Functions to consider
  are 
  such as \code{\link[BayesDiallel]{PlotHPD}} and \code{Summary}.  
  
The BayesSpikeCL member "BSAFD" or "BS" will be added to the AutoDiallelObject.

This object has Coda Chains in the form of \code{CAFD$BS$CodaList};

To get final estimated MIP Probabilities use  \code{CAFD$BS$MIP}

To get Posterior Probability chains use \code{CAFD$BS$ProbCodaList}

Documentation for other BayesSpike Features is in development.

BayesSpike Requires Harddrive storage.  Set "SaveDir" to a directory to
 which Coda chains can be stored.
}
\value{
This will return a copy of \code{AFD} (let's call it \code{CAFD}), 
  which has an element \code{CAFD$BS} which is the \code{\link{BayesSpikeCpp}}.
  
}
\references{ Diallel Paper, pending BayesSpike submission }
\author{ Alan Lenarcic, Will Valdar }
\note{ 
  "\code{\link{TBSR5}}"  is a class using R5  \code{\link[methods]{ReferenceClasses}} object.  
  A \code{BayesSpikeCpp} is a \code{\link[Rcpp]{Rcpp}}
  object.  Objects for package \code{\link[BayesDialle]{BayesDiallel}} are largely
  in \code{\link[R.oo]{R.oo}} format.  
}
\seealso{ \code{\link{BayesSpikeRegression}}, \code{\link{BayesSpikeFromAutoDiallel}}, 
  \code{\link[BayesDiallel]{BayesDiallel}}, 
  \code{\link[BayesDiallel]{DiallelAnalyzer}},
  \code{\link[coda]{mcmc}},
  \code{\link[Rcpp]{Rcpp}}, \code{\link[R.oo]{R.oo}}, 
  \code{\link[methods]{ReferenceClasses}},
  \code{\link{Save}}, \code{\link{ReLoad}}
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
   SaveAFDFile = "SaveAFDBackUpPiximusData.RData", LogTransform=TRUE);

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
  try(dir.create("~/Documents/Stat", showWarnings=FALSE));  
  try(dir.create("~/Documents/Stat/BayesSpikeSaves", showWarnings=FALSE));
  BSSaveDir = "~/Documents/Stat/BayesSpikeSaves"
} else if (length(list.files("~/")) > 0) {
 try(dir.create("~/Stat"));  try(dir.create("~/Stat/BayesSpikeSaves"));
  BSSaveDir = "~/Stat/BayesSpikeSaves"
}

MyProcTime = tSeq(proc.time()[1]);
BSSaveDir  = paste(BSSaveDir, MyProcTime, sep="");  ## add a time stamp to analysis

AFD2 = DiallelAnalyzer(data = Piximus.Data, father.strain="father.strain.name",
   mother.strain="mother.strain.name", phenotype="MouseWeight", 
   is.female="is.female",
   sep="", FixedEffects = NULL,
   RandomEffects = NULL, 
   Models=Example.Piximus.Models[1:2], sigmasq.start = 1,  
   numChains = 5,  lengthChains=800,
   burnin = 1,
   DIC.Only=FALSE,  tauPriorFile = Example.Piximus.tau.Prior.Info,
   DoBayesSpike = TRUE, BSSaveDir = BSSaveDir, LogTransform=TRUE);
   
BSAFD <- AFD2$BSAFD;
BSAFD$MIP;

plot(BSAFD$CodaList, ask=TRUE);
   
## Get a summary of the Coda Chains for first Diallel Model
summary(AFD2$BSAFD$CodaList)
plot(AFD2$BSAFD$CodaList);
plot.hpd(AFD2$BSAFD$CodaList);

#### Save the Work to Hard=drive
SaveIn <- paste(BSAFD$SaveDir, "//FinalOutputSave/", sep="");
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
