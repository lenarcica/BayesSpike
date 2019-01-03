\name{CodaList}
\alias{codalist}
\alias{CodaList}
\alias{codalist}
\alias{CodaTable}
\alias{codaTable}
\alias{DoRecord}
\alias{dorecord}
\alias{Dorecord}
\alias{doRecord}
\alias{codaList}
\alias{Coda}
\alias{coda}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Coda MCMC chains.}
\description{
   Package \link[coda]{mcmc} formatted chains for Bayesian posterior inference and plotting. 
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
  \item{CodaList}{A \link[coda]{mcmc.list} list object with multiple gibbs sampler information.}
  \item{DoRecord}{Length 7 integer flag vector deciding what content to record in mcmc chains}
  \item{CodaTable}{Current \link[coda]{mcmc} chain.}


}
\details{
   \code{CodaList} is an accessible element of \code{\link{BayesSpikeCpp}} object.  It is in \link[coda]{mcmc.list} format.
 
   However, by default, \code{\link{BayesSpikeRegression}()} does not record
   Gibbs sampler chains in RAM memory at all, instead writing to file buffers using
   \code{\link{SaveDir}}.  
}
\value{
To store samples in RAM, the vector \code{DoRecord}
   which is length 7, must have values set to 1 (as opposed to 0) as logical flags.  For the seven locations in DoRecord toggle recording of te following information:
   
    \item{1}{Record \code{\link{OnBeta}} or \eqn{\boldsymbol{\beta}}{beta} Draws, even 
     for zero valued \eqn{\beta_{j}}{beta_j}.}
    \item{2}{Record \code{\link{OnTau}} or \eqn{\boldsymbol{\tau}^2}{tau^2} draws, even
    for zero valued \eqn{\tau^2_{j}}{tau^2_j}.}
    \item{3}{Record \code{\link{OnSigma}} or \eqn{\sigma^2}{sigma^2}
     noise dispersion paramter.}
     \item{4}{Record \code{\link{OnPiA}} or \eqn{\pi_{\mbox{\tiny{$\mathcal{A}$}}} }{piA} prior activation draws.}
    \item{5}{Record Rao-Blackwellized posterior probabilities for activation or 
    \eqn{P( \beta_{j} \neq 0)}{P(beta_j neq 0)} for all fixed paramters.}
   \item{6}{Record Rao-Blackwellized  posterior probabilities for activation
   or  \eqn{P( \tau^2_{k} \neq 0)}{P(tau^2_k neq 0)} for all random effects grouped paramters.}
   \item{7}{Record weights \eqn{\lambda_{i}}{lambda_i} for t-dispersion paramters
    if \code{\link{dfTNoise}} is > 0.}
    
   This information can rapidly clog memory, so saving to compressed files with 
   \code{\link{SaveDir}} is considered a valuable alternative when \eqn{p}{p}
   is large.  In this case, the option of using \code{\link{ReturnSubCodaFromFiles}}
    (See \code{\link{SubCodaList}}) must be used, which requires a user to choose a subset
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
  
  Set \code{\link{SubCodaList}} For how harddrive stored Coda objects are loaded post algoirthm.
}
\seealso{ \code{\link{BayesSpikeRegression}}, \code{\link{BayesSpikeFromAutoDiallel}}, 
  \code{\link{BayesDiallel}}, \code{\link{coda}},
  \code{\link{Rcpp}}, \code{\link{R.oo}}, \code{\link{methods::ReferenceClasses}},
  \code{\link{SubCodaList}}
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
