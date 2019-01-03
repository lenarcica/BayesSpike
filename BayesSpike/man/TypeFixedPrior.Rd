\name{TypeFixedPrior}
\alias{TypeFixedPrior}
\alias{typefixedprior}
\alias{tauFixed}
\alias{taufixed}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Prior for "fixed", individual Beta coefficients}
\description{
  This is a discussion on how to set prior information for fixed coefficients 
  (coefficients that do not belong to a grouping prior).
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
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{TemperatureList}{A List of Temperatures (such as c(10,5,2,1), necessary to activate EE-sampling.
    For typical regression needs, this list should be a list descending of postive, >= 1 real numbers.  
    If the last number is not 1.0, a 1.0 will be added to the end of the list in default.  

    }
  \item{EEMergeEvery}{Number of iterations of Gibbs sampler between every "merge step" (draw from higher temperature energy level)}
  \item{DoEEProbSort}{A flag whether to sort Energy level samples at end of simulation at a temperature level.}
  \item{EEProbSortWidth}{A width of energy to consider as an "equi-energy slice" for merge draws.  Although, the merge step will 
  consider all draws within EEProbSortWidth of the current draw, it will probabilistically
     weight these draws, locally, based upon differential likelihood of the current and proposed sample}
  \item{TemperatureDecreasingRate}{A number <= 1.0 to represent how many fewer draws from the next higher temperature should be taken.}
}
\details{
 Equi-Energy Sampling (Kou and Wong 2004) is a technique for MCMC under highly multimodal surfaces.
 
 It is simpler to Geyer (1998)'s parrellel tempering, in that MCMC draws are taken from a higher "temperature"
  of the true posterior under consideration.  In this case, that is { f(Theta|Y) }^(1/T) where T is a positive number.
Higher temperatures flatten out the posterior, producing chains that more readily cross multi modal divides at high temperatures.

However EE-sampling differs from Parrallel Samplings Parrallel chains by sampling from energy distributions in order.

First draws from the highest temperature are taken. 
 Then draws from the next temperature are taken, with occasional "merge" draws
 or draws from the previous samples in the draw from the highest temperature are taken. 
  One identifies quality "merge" draws from the previous posterior  by considering the
  current "energy" or posterior density value of the current point.  Then only draws within an equivalent energy band
  from the previous draws are considered.  These draws, if their posterior value is the same as the current value
  are probabalistically identical to the current draw and may serve as a new
  seed or jump-point for the rest of the chain of the Gibbs sampler.
  
EE sampling then repeats this routine until it reaches the "1" or original temperature.
Only draws from the 1 temperature constitute true draws from the posterior.

To perform EE-sampling "BSSaveDir" must be set, as it involves hard-drive storage of
previous draws.

A warning for selection problems, higher temperatures necessarily include less sparse
draws from the posterior.  The more active variables in a draw, the slower the routine.
EE-sampling typically uses a pyramid drawing schme, fewest draws from the highest temperatues,
a few more at each lower temperature, and most draws from the bottom temperature.  
Setting "TempertureDecreasingRate" to a number < 1, such as .6, 
achieves this pyramid structure.  


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
  \code{\link{Rcpp}}, \code{\link{R.oo}}, \code{\link{methods::ReferenceClasses}}
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
