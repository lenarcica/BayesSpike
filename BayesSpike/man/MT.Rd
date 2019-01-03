\name{MT}
\alias{CurrentTauPriors}
\alias{TauOniiPrior}
\alias{MTOnii}
\alias{oniiMT}
\alias{mt}
\alias{TauOfFContainer}
\alias{tauoffcontainer}
\alias{Tauoffcontainer}
\alias{Mt}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{MT the loaded properties of Onii loaded in group of Tauii parameters.}
\description{
 MT object MT contains the critical data and parameters called upon to update
 a group paramter.
 
 The key innovation of BayesSpikeRegression is an efficient Gibbs sampler
 method for sampling selective groups of paramters indexed \eqn{\mathbb{J}(k)}{J(k)}
 for sub \eqn{\boldsymbol{\beta}}{beta} vector
 \eqn{\boldsymbol{\beta}_{\mathbb{J}(k)}}{beta_J(k)}.  
 
 The dispersion paramter \eqn{\tau^2_{k}}{tau^2_k} supposes that a-priori
 
 \eqn{\beta_{j} \sim N(0, \tau^2_{k}) \mbox{ for all $j \in \mathbb{J}(k)$}.}
 {beta_j ~ N(0, tau^2_k) for all j in J(k).}
 
 MT is a C structure that contains the data key iformation for calling 
   SampleANewTau(Onii)
 
 This shows all of the prior information useful to this part regression. 
 
 Please see the equation (18) featured in section 3.2 of the paper to understand
 the critical "D" and "R" pointers which are used by functions.
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
  \item{MT}{
     See \code{BayesSpikeCpp.h} and \code{BayesSpikeSliceSampling.cpp} for usage
  }
  \item{PriorOftau}{A Prior for the tau parameter}
  \item{staubarnu}{Mean of tau by default prior (Inv ChiSquared)}
  \item{staupriordf}{Degrees of freedom of tau prior};
  \item{sPriorXScaling}{SEXP with nu and df};   
  \item{D}{Diagonal Eigen values of X_g^T X_g for covariates in group g} 
  \item{R}{Residual correlation between diagonalized X_g coordinates and Y}
  \item{Onii}{Which iteration of Gibbs sampler (used in priting state)}
  \item{logOddsPi}{Information on "Odds" of PiA_Groups for activation default of a group}
  \item{Temperature}{Info on the Temperature of Sampler}
  \item{invTemperature}{Multiplicative inverse of that Temperature}
  \item{q3X2}{Elements q3X2, f3X2, f2X2,  f2X4, X2, X4, q4X4, xM2, xM3, xM4, lM3xM3, lB, lM4xM4, lM2xM2, lA, lLB42, lUB32, Df3, Df4,
    are all declared in first paragraph Section 5.3.  These are calculated quantities.
    "2" is the target density.
    "4" is the de bound inverse chisquared, and "3" is the lower bound. }
  \item{rUnif1}{Space to put a random uniform draw}

}
\details{
   At its heart, the major innovation of \code{BayesSpike} is the ability to do
   Group-Lasso type grouped paramter selection with a Bayesian prior of activation for a group of coefficients.
   
   If \eqn{\tau^2_{k}=0}{tau^2_k = 0} then all coefficients within
   \eqn{\mathbb{J}(k)}{J(k)} will be zero, conversely if \eqn{\tau^2_{k} > 0}{tau^2_k > 0}
   then all coefficients within \eqn{\mathbb{J}(k)}{J(k)}  are non-zero with prior
   sampling shrinkage distribution \eqn{N(0, \tau^2_{k})}{N(0,tau^2_k)}.  
   
   \code{\link{SampleANewTau}()} is  the \code{\link{BayesSpikeCpp}} function that samples
   a new tau.  
   
   For efficiency's sake, every \eqn{\mathbb{J}(k)}{J(k)} group is a consecutive list of
    indices.  
}
\value{
    Object analysis object is a "Rcpp module", with various properties and methods
  
     \item{TauOnii}{First member, the Onii integer that MT is configured to}
     \item{NumJ}{Number of elements to group Onii}
     \item{D}{D diagonal eigenvalues of the group Onii XtX}
     \item{R}{Residual Gamma *b XtResid eigenvalues of groupOnii XtX}
     \item{NumJ}{Number of elements to group Onii}
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
