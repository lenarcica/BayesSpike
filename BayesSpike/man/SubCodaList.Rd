\name{SubCodaList}
\alias{SubcodaList}
\alias{.TauCodaList}
\alias{TauCodaList}
\alias{PiACodaList}
\alias{SigmaCodaList}
\alias{SigCodaList}
\alias{SubsetCoords}
\alias{SubSetCoords}
\alias{.SubSetCoords}
\alias{SubSetTau}
\alias{.SubSetTau}
\alias{.SubsetTau}
\alias{subsettau}
\alias{SigmaCodaList}
\alias{.SigmaCodaList}
\alias{WeightCodaList}
\alias{ProbCodaList}
\alias{SubSetCoords}
\alias{PostProbCodaList}
\alias{PosteriorProbabilityCodaList}
\alias{PosteriorCodaList}
\alias{SubCodaLongList}
\alias{subcodalonglist}
\alias{longlist}
\alias{LongList}
\alias{YCodaList}
\alias{SubSetCoords}
\alias{subsetcoords}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Subset, disk-loaded Coda MCMC chains.}
\description{
   Package \link[coda]{mcmc} formatted chains for Bayesian posterior inference and plotting. 
   
   These are loaded from saved hard disk compact files.
   
   SubCodaList: refers to Beta parameter draws
   
   TauCodaList: refers to Tau parameter draws
   
   PiACodaList: refers to PiA parameter draws
   
   SigmaCodaList: refers to Sigma parameter draws
   
   YCodaList: refers to Y draws (in logistic regression)
   
   WeighCodaList: refers to Y weight draws (in t regression and logistic regression)
   
   ProbCodaList: refers to Beta activation posterior draw probabilities.
   
   
   By attaching "SubSetCoords" to a BayesSpikeObject, one can get a "SubCodaList"
   to be retrieved and returned featuring only the coordinates requested.  There is also a \code{SubSetTau}
   that can be set to get a subset of parameters introduced to \code{TauCodaList}.
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
  
  # BayesSpikeOB$SubSetCoords <- c(1,3,4,...);
  # BayesSpikeOB$SubSetTau <- c(1,4,5,3,...);
  #summary(BayesSpikeOB$SubCodaList)
  #summary(BayesSpikeOB$TauCodaList)
}
%- maybe also 'usage' for other objects documented here.
\arguments{

  \item{SubSetCoords}{Set this parameter to \code{\link{BayesSpikeCpp}} or
    \code{\link{BayesSpikeR5}} object to a set of ordered integers between 1 and p.
    SubCodaList will report only the Beta draws for these coordinates}
  \item{SubSetTau}{Set this parameter to  \code{\link{BayesSpikeCpp}} or
    \code{\link{BayesSpikeR5}} object to a set of ordered integers between 1 and Number of Random effects groups.
    SubCodaList will report only the Tau draws for these coordinates }
  \item{WritePostProbBuffer}{A logical flag.  If set to FALSE or 0, this will not write information
   to harddrive reporting Posterior Probabilities.  Thus ProbCodaList will be
   inaccessible post algorithm.}
  \item{WriteWeightBuffer}{A lgoical Flag. If set to FALSE or 0, this will not write information
   to harddrive reporting sample weights.  Thus WeightCodaList will be
   inaccessible post algorithm.}
  \item{WriteWeightBuffer}{A logical Flag. If set to FALSE or 0, this will not write information
   to harddrive reporting sample weights.  Thus WeightCodaList will be
   inaccessible post algorithm.}
  \item{WriteYBuffer}{A logical Flag. If set to FALSE or 0, this will not write information
   to harddrive reporting Y draws (in logistic regression).  Thus YCodaList will be
   inaccessible post algorithm.}
  \item{DoLongCI}{To save Long Intervals, this must be true in \code{\link{BayesSpikeRegression}} call.
    The "Long" Buffer is a storage of gibbs samples that are all non-zero for all p Beta coefficients
    This can only be saved on a file (Not stored in Coda List even if DoRecord[1] == 1).  This draws a nonzero
    Beta_j even if sampler itself drew a zero Beta_j.  This is used for robust confidence intervals for
    coefficients near zero.  These confidence intervals are overly wide, though they benefit from shrinkage on the other coefficients.
    See \code{\link{SubCodaList}} for more details, as well as \code{\link{BetaSymmetricLongQuantiles}}
    }
  \item{CpBuffLongCI}{Integer, input to \code{\link{BayesSpikeRegression}}, How long should Long (No shrinkage to Zero) Buffer should be of size}


}
\details{
   All sub coda elements, like \code{SubCodaList}, are accessible elements of \code{\link{BayesSpikeCpp}} object as
   well as \code{\link{BayesSpikeR5}} objects.  They are in \link[coda]{mcmc.list} format.
 
   By default, \code{\link{BayesSpikeRegression}()} does not record
   Gibbs sampler chains in RAM memory at all, instead writing to file buffers using
   \code{\link{SaveDir}}.  If directories are provided, this will write all draw information.
   
   Accessing these members causes the object to load these in from harddrive.
   
   These are important in calculating \code{\link{AutomaticQuantiles}}.
   
   These are accessible as a method from the \code{\link{BayesSpikeObject}}.  So "BayesSpikeObject$SubCodaList" will activate the retrieval
   function and attach the SubCodaList as a member coda-matrix to the the BaysSpike Object.
}
\value{

   \item{\code{SubCodaList}}{Is a mcmc.list of Beta Vector Draws (potentially smaller than
    all p draws if \code{SubSetCoords} is set).}
   \item{\code{TauCodaList}}{Is a mcmc.list of Tau Vector Draws (potentially smaller than
    all RandomEffects draws if \code{SubSetTau} is set).}
  \item{\code{PiACodaList}}{ Is mcmc.list of \eqn{\boldsymbol{\pi}^2_{
    \mbox{\tiny{$\mathcal{A}$}}} }{piA^2}  Draws}
  \item{\code{SigmaCodaList}}{ Is mcmc.list of \eqn{\sigma^2}{sigma^2}  Draws}
  \item{\code{YCodaList}}{ Is mcmc.list of Y draws.  Only valid for Logistic regression draws of Y}
  \item{\code{WeightCodaList}}{Is mcmc.list of weight draws.  Useful in logistic regression and t regression}
  \item{\code{ProbCodaList}}{Posterior Probability Draws}
  \item{\code{SubCodaLongList}}{Using Save Variables in SubCodaVector, that draws SubCodaList, 
    this is a list of "unshrunk" draws, that is, whenver Betaj was drawn zero, it is redrawn as a nonzero variable
    such that confidence intervals are robust  near zero if sparseness hypothesis fails.
    Must have set \code{DoLongCI} in \code{\link{BayesSpikeRegression}} to true.}
  
    
   This information can rapidly clog memory, so saving to compressed files with 
   \code{\link{SaveDir}} is considered a valuable alternative when \eqn{p}{p}
   is large.  In this case, the option of using \code{\link{ReturnSubCodaFromFiles}} must be used, which requires a user to choose a subset
   of paramters to reload.  

    \code{\link{BayesSpikeRegression}} serves as a wrapper function to write and read these files.
    Setting \code{NewWrite=0} will by default set \code{\link{BayesSpike}} to attempt to read currently written
    files rather than start new chains.  If the files are only half written, then BayesSpike
    will read in current chains and attempt to append until the chains are finished.
    If chains are fully finished, \code{\link{BayesSpikeRegression}} will simply read in the files.
}
\references{ Diallel Paper, pending BayesSpike submission }
\author{ Alan Lenarcic, Will Valdar }
\note{ 
  "TBSRoo"  is a class using R.oo which we would like to deprecate in favor of 
  an R5 class, "\code{\link{TBSR5}}"  or "\code{\link{BayesSpikeR5}}" is  the R5 storage class.
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
