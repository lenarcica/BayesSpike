\name{EESampling}
\alias{EE-Sampling}
\alias{Eesampling}
\alias{EESampling}
\alias{eesampling}
\alias{EEsampling}
\alias{TemperatureList}
\alias{EEMergeEvery}
\alias{DoEEProbSort}
\alias{EEProbSortWidth}
\alias{BayesSpikeMergeAndPriors}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Equi-Energy Sampling for Bayes Spike}
\description{
  This is a slight discussion of the process of EE-sampling as it is 
used to accelerate mixing of this algorithm.  
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
  \item{BayesSpikeMergeAndPriors}{Code used to setup whether the Bayes Object should do EE sampling and which rate it should take higher temperature draws with}
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

##  Test Code To see how well the functions are running on server.

## Library TwoSimR5 contains a methodology to double 
library(TwoSimR5); library(BayesSpike);
n <- 100; sigma = 1.0; p = 1000; k = 6;

 DefaultCorrelationXmatrix <- function(kLen) {
      CorrelationXmatrix1(kLen, CovRho);
  }


SMS <- SimMeData(n = n, p = p/2, k = k, sigma = sigma,
     ElseGetX = NULL, 
     SigmaNoise = -999, GenerateBetaVec = DefaultGenerateBetaVec, 
     CorrelationXmatrix = DefaultCorrelationXmatrix, tNoiseDF = 0,
     LogitRegression = FALSE, LogitFactor = 1, Beta0 = 0, ExperimentName="", jobii= 0, 
     WorkingRow = 0, AlreadyLocked=TRUE, MeanCenterXX = TRUE,
     UniqueProcessIdentifier = "", ISample = 1,
     DontSave=TRUE, InOrder=FALSE, dfRobit = -1);
SMS <- DoubleSimulation(SMS, AlreadyLocked=FALSE, DontSave=TRUE);     
     
     
tauSqA = 1; CountSamp = 3200;
  NumChains = 1; PriorStrength = SMS$PriorStrength; SigmaPriorStrength = SMS$SigmaPriorStrength;
  CutOff = .1; DoCI = DefaultDoCI; DoMedian = TRUE; AutoInfo="Auto"; NumTemp=2;TestAtGoal=TRUE;
  EEMergeEvery = 5;        burnin = 200;
eval(parse(text=IntroGenerateBayesSpikeText()));

 BetaStart <- SMS$BetasPartReal;   
 TemperatureList <- c(1.5,1.3,1.1,1);
 MBS <- BayesSpike:::BayesSpikeRegression(X=SMS$X, Y=SMS$Y, 
      tauASq = tauSqA, PiAStart=SMS$puse,InitKKs=6, MaxGibbsIters = CountSamp,
      NumChains=NumChains,
      SigmaSq = SMS$SigmaSqNoise, Verbose = 1, tauEndList = NULL, PiAPrior = PiAPrior,
      SigmaPrior = SigmaPrior, DoSave=TRUE, DoSaveTBSR5=FALSE, DoRecord = c(0,0,0,0,0,0,0),
      SaveDir = BSSaveDir, IndexFirstRandomEffect = -1,
      StartRunProbVector = 100,
      DoLogitPreProb = DoLogitPreProb, DoLogitPostProb = DoLogitPostProb,
      dfRobit = dfRobit, PreRunMaxGibbsIters = PreRunMaxGibbsIters, 
      AlterWeightFlag = AlterWeightFlag, burnin=burnin, TemperatureList=TemperatureList,
      EEMergeEvery = EEMergeEvery,
      BetaStart=BetaStart, Run=TRUE);
  plot(MBS$MIP)
  length(MBS$MIP[MBS$MIP > .05]);
  MBS$MIP[SMS$BetasReal != 0]
  

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }
\keyword{ ts }% __ONLY ONE__ keyword per line
