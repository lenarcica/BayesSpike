\name{OnSigma}
\alias{Sigma}
\alias{sigma}
\alias{SigmaPrior}
\alias{sigmaPrior}
\alias{sigmaprior}
\alias{dfTNoise}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Noise distribution}
\description{
 Noise parameter OnSigma represents \eqn{\sigma^2}{sigma^2}  for
 a \eqn{\varepsilon_{i} \sim N(0, \sigma^2)}{epsilon_i ~ N(0, sigma^2)}
 or \eqn{\varpeilon_{i} \sim \sigma t_{\mbox{\tiny{d.f.}}} }{ epsilon_i ~ sigma * t(d.f.) }
 noise distributions.  When Robit information is used, \eqn{\sigma^2}{sigma^2}
 can be treated like a information magnification parameter.   Parameter \code{dfTNoise} sets the degrees of freedom for this T-noise,
 with setting it zero or negative supports gaussian distribution.
 
 The prior is entered in the form of vector \code{SigmaPrior}.  \code{SigmaPrior}
 should be length 2 \eqn{(\nu_{\sigma}, \bar{\sigma}^2)}{(nu_sigma, barsigma^2)}
   to program an inverse gamma 
 \eqn{\sigma^2 \sim   \bar{\sigma}^2 / \mbox{Gamma}[\nu_{\sigma}] }{sigma^2 sim barsigma^2 / Gamma[nu_sigma]}
 prior distribution for the noise.  \eqn{(1,1)} is the default prior for this noise distribution.
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
  \item{OnSigma}{Initial value to start \eqn{\sigma^2}{sigma^2}}
  \item{SigmaPrior}{Length 2 Vector for Inverse gamma prior for \eqn{\sigma^2}{sigma^2}}
  \item{dfTNoise}{Set degrees of freedom for T noise.}

}
\details{
  Noise is modelled as default gaussian 
   \eqn{\varepsilon_{i}\sim N(0, \sigma^2)
   }{epsilon_i~N(0,sigma^2)},
   though a t-weighting method 
   \eqn{
     \varepsilon_{i} \sim \sigma \times t_{\mbox{\tiny{d.f.}}} 
   }{epsilon_i ~ sigma*t(d.f.)}
   is possible, or a robit distribution
   \eqn{
     Y_{i} \sim \mbox{Bernoulli}(p(X))
   }{Y_i ~ Bernoulli(p(X))}
   which has 
   \eqn{ 
     p(X) = \Phi_{\mbox{\tiny{d.f.}}}^{-1}( \mathbf{X} \boldsymbol{\beta})
   }
   {p(X) = Phi_d.f.^-1(XBeta)}
   where 
   \eqn{
     \Phi_{\mbox{\tiny{d.f.}}}^{-1}
   }{Phi_d.f.^-1} 
   is the inverse cdf of 
   a t distribution.  (If d.f. is set to 7, this is very logit regression like,
    if d.f. is set to a large number like infity, this is a probit regression 
    with the normal cdf.)
    
  d.f. is set by setting \code{dfTNoise}. 
                        
   Hyper-Prior sampling for \code{dfTNoise}
  is currently not supported and must be considered as a user robustness choice.
  t-noise sampling is slower than Gaussian sampling since every sample requires a reweighting
  of \eqn{
   \wildetilde{\mathbf{X}^{T} Y} = 
   \sum_{i} X_{ij} \lambda_{i} Y_{i}
  }
  {
    tildeX^TY =
    sum_i X_ij lambda_i Y_i
  } and   
  \eqn{
    \widetilde{\mathbf{X}^{T} \mathbf{X}}=
    \sum_{i}  \lambda_{i} X_{i,}^{T}  X_{i,}
  }
  { tildeX^X=sum_i lambda_i X_{i,}^{T}  X_{i,} } 
    
and can take considerable computation time if sparsity is low (Only relevant vectors are reweighted if sparsity is high).
    This is also a significant computational cost for t-sampling. 
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
