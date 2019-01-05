\name{AM}
\alias{AM}
\alias{MakeAM}
\alias{MakeAm}
\alias{Makeam}
\alias{am}
\alias{MatrixM}
\alias{matrixm}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Matrix "AM" is a transform from n linearly dependent columns to n dependent columns}
\description{
Matrix AM was described in Section "3.1" of the BayesSpike Paper.

This matrix can be used to transform any linearly independent group of n columns 
assigned to a group parameters into
ideally a set of n-1 set that removes an additive dependence constraint (i.e. \eqn{a+b+c+d+e...=0}{a+b+c+d+e..==0})

It is used in \code{\link{BayesSpikeFromAutoDiallel}} so that selection is known
to operate on restricted columns.
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
  \enumerate{
   \item{n}{MakeAM(n) takes this many dimensions for function}
   \item{UseDetFAC}{\code{\link{BayesSpikeFromAutoDiallel}}() takes this as an optional parameter, usually set to zero, so as to specify no prior correction 
   while using the M matrix to make constrained parameters.}
  }

  MakeAM is a function that takes one input "n".
  
  It constructs a transformation matrix.  
  
  This matrix can be transposed to reconstruct
  data.  Note that \eqn{M^T{}M=n/(n-1)\times{}I(n-1)}{M^T*M==n/(n-1)*I(n-1)}, and thus for vectors
    \eqn{\beta\in R^n}{\beta in R^n} such that \eqn{\sum\beta=0}{sum(\beta)==0} centered then  
    \eqn{\frac{n-1}{n}M{}M^T\beta-\beta = 0}{(n-1)/n*M*M^T*\beta-\beta==0}.

}
\details{
   This is only used automatically with a Diallel structure, and in practice, a very
   particular user-friendly code would have to be written if groups were all of very different
   size and requirements within a given group sampler.
   
   To be pedantic, if we were fully certain of our prior opinions, albeit only on the
original X scale,   if we  manipulate of columns of X, to be particularly accurate, would require
   a slight adjustment of the "OnPiA" probability of inclusion.  
   
   This is because the transformation "M" adds a multiplicative factor of det|M| to a calculation.
    "\deqn{\int f(X \beta) p(\beta) d\beta = \int f(X \mathbf{M} \beta') p(\beta') \mathbf{M} d\beta' }{%
int f(X \beta) Prior(\beta) dbeta = int f(XM \beta') Prior(\beta') M d\beta'}"

   That said it isn't clear that, once one has demanded an \eqn{\mathbf{M}\mathbf{M}^{T}}{"M * t(M)"} correlation structure on Beta,
   that one's relative a opinion on the prior inclusion of the "tau^2_k" isn't intended for the
   transformed unconfounded parameters, instead of the confounded ones.  When groups are all confounded similarly and of similar size
   the relative effect to including or excluding a slight modification of \code{\link{PriorProbTau}} is probably more trouble than its worth. Given that 
   the confoundedness of X and the constraint condition  \eqn{\sum_{j\inJ(k)} \beta==0}{sum_{j in J(k)}\beta==0} mean that 
   \eqn{d\beta}{d\beta} can only  be integrated over the first \eqn{||J(k)|-1}{|J(k)|-1} coordinated.  The "UseDetFAC" feature is usually de-activated in the code for
   the paper.
   
   An effect to this is implemented in \code{\link{0009BayesSpikeAutoDiallel}} file, but it is generally a users responsibility to decide what sort
   of confounding their data has and which prior scale to choose for their data.  The modifications are performed using the inputs is effected on the \code{\link{PriorProbTau}}
   input to BayesSpike.  This is a relative change to the effect that \code{\link{OnPiA}} has on the overall group inclusion, shrinking it slightly on the Odds scale towards zero.
}


\value{
    Object analysis object is a "Rcpp module", with various properties and method 
}                                                                        
\references{ Diallel Paper, pending BayesSpike submission }
\author{ Alan Lenarcic, Will Valdar }
\note{ 
  "TBSRoo"  is a class using R.oo which we would like to deprecate in favor of 
  an R5 class, "\code{\link{TBSR5}}"  or "\code{\link{BayesSpikeR5}}" is  the R5 storage class.
}
\seealso{ \code{\link{BayesSpikeRegression}}, \code{\link{BayesSpikeFromAutoDiallel}}, \code{\link{PriorProbTau}},
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

  
