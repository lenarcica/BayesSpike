\name{DoLongCI}
\alias{DoLongCI}
\alias{DolongCI}
\alias{dolongCI}
\alias{dolongci}
\alias{longci}
\alias{LongCI}
\alias{LongCi}
\alias{Credibility Intervals}
\alias{CI}
\alias{CredibilityIntervals}
\alias{Long Credibility Intervals}
\alias{credibility}
\alias{credibility intervals}
\alias{KeepPosteriorQuantiles}
\alias{keepposteriorquantiles}
\alias{GetConfidenceIntervals}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Credibility Intervals}
\description{
  Credibility Intervals are the main purpose of a Gibbs sampler.  HPD or symmetric quantile intervals are offered
  for BayesSpike, with HPD generally recommended.  Even more prefereable "Long Credibility Intervals" for 
  (\code{DoLongCI=TRUE}) setting will have the closest to correct coverage of true Beta_j (even Beta_j that are non-zero but close to zero)
  and attempt to have a "fair" interval for all coordinates  (Coordinates, even if the point estimate is zero and have significant 
  posterior probability of zero, will get a non-zero credibility interval demonstrating what the space of information
  suggested for this parameter in a possibly non-zero setting includes.)

  "DoLongCI=TRUE setting:"
 To save Long Intervals, this must be true in \code{\link{BayesSpikeRegression}} call.
    The "Long" Buffer is a storage of gibbs samples that are all non-zero for all p Beta coefficients
    This can only be saved on a file (Not stored in Coda List even if DoRecord[1] == 1).  This draws a nonzero
    Beta_j even if sampler itself drew a zero Beta_j.  This is used for robust confidence intervals for
    coefficients near zero.  These confidence intervals are overly wide, though they benefit from shrinkage on the other coefficients.
    See \code{\link{SubCodaList}} for more details, as well as \code{\link{BetaSymmetricLongConfidenceMatrix}}.
    
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
\item{
   \code{\link{CodaList}}}{\code{BayesSpikeOB$CodaList} will create an R list 
   object with Coda tables for each parameter set.  Confidence intervals can
   be calculated on this Coda Lsit by choice
   }
   
   
   To get BayesSpike to calculate a set of quantiles (say 5%, 2.5% or 1% credibility intervals)
   assign a vector of alpha values, say
   \code{BayesSpikeOB$TBSR5$KeepPosteriorQuantiles <- GetConfidenceIntervals}
   
   where in this case the alpha values are
   \code{GetConfidenceIntervals = c(0.500,0.250, 0.750, 0.050, 0.950, 0.025, 0.975, 0.010, 0.990, 0.005, 0.995)}
   Other than the first entry (set to give the median draw from each posterior for the draws), the others are two sided
   probability for symmetric 50\% (.25, .75), 10\% (.05, .95, 5\% (.025, .975), 2\% (.01, .99),
   1\% (.005,.995) two sided tails.  HPD intervals will detect the intended posterior value and give two sided invervals
   (even if non symmetric they should total to the target probability.) 
   
   
   Once KeepPosteriorQuantiles is set, if one accesses \code{BayesSpikeOB$TBSR5$(MemberMethod)}
   one will generate a calculation of critical CI of interest.
   
   \item{\code{BS$TBSR5$BetaSymmetricQuantiles}}{
     Symmetric Quantile (aka from (alpha to 1-alpa of posterior marginal quantity)
     around posterior.
   }
   \item{\code{BS$TBSR5$BetaHPDQuantiles}}{HPD quantile taken on the posterior draws.
     Note this could be only zeroes.  Due to the fact that this 
     could be completely to left or right of zero line, these HPD draws show there is technically something wrong 
     with an MCMC based  HPD calculation
     aren't really the "highest posterior density", which should always infinite posterior around zero.
   }
    \item{\code{BS$TBSR5$BetaSymmetricUnshrinkQuantiles}}{
     Symmetric Quantile (aka from (alpha to 1-alpa of posterior marginal quantity)
     around posterior.   Here we throw out all zero draws.
   }
   \item{\code{BS$TBSR5$BetaHPDUnshrinkQuantiles}}{HPD quantile taken on the posterior draws.
     This is based upon throwing out all zero draws for each coefficient.
   }
         \item{\code{BS$TBSR5$BetaSymmetricLongQuantiles}}{
     Symmetric Quantile (aka from (alpha to 1-alpa of posterior marginal quantity)
     around posterior.   These draws are possible if "DoLongCI" was set true.
     These are pulled from conditional posterior if Beta_j was non-zero each draw,
     condition on the sparse other coordinates.
   }
   \item{\code{BS$TBSR5$BetaHPDLongQuantiles}}{HPD quantile taken on the posterior draws.
     This is based upon throwing out all zero draws for each coefficient.
     hese draws are possible if "DoLongCI" was set true.
     These are pulled from conditional posterior if Beta_j was non-zero each draw,
     condition on the sparse other coordinates.
   }

}
\details{
For large models (large p), to get credibility on every coordinate, the data must be slowly read in from the saved file buffers, especially the LongCI buffers
if all coordinates are getting saved.  In this case we load in a CODA table of a subset of coordinates, calculate Credibility intervals and move on.
  
}
\value{
    Object analysis object is a "Rcpp module", with various properties and methods
     \item{\code{\link{CodaList}}}{
     Coda mcmc table list of output if generated.  OnPiA may be an element of this list if it is set through Beta Priors
     }
}
\references{ Diallel Paper, pending BayesSpike submission }
\author{ Alan Lenarcic, Will Valdar }
\note{ 
  "TBSRoo"  is a class using R.oo which we would like to deprecate in favor of 
  an R5 class, "\code{\link{TBSR5}}"  or "\code{\link{BayesSpikeR5}}" is  the R5 storage class.
}
\seealso{ \code{\link{BayesSpikeRegression}}, \code{\link{BayesSpikeFromAutoDiallel}}, \code{\link{OnBetea}}, \code{\link{OnPiA}},
  \code{\link{tauEndList}}, \code{\link{Y}},
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


GetConfidenceIntervals = c(0.500,0.250, 0.750, 0.050, 0.950, 0.025, 0.975, 0.010, 0.990, 0.005, 0.995)
## Load in R packages
library(BayesSpike);
try(require(lattice, quietly=TRUE, warn.conflicts=FALSE), silent=TRUE);
try(require(coda, quietly=TRUE, warn.conflicts=FALSE), silent=TRUE);
try(require(R.methodsS3, quietly=TRUE, warn.conflicts=FALSE), silent=TRUE);
try(require(R.oo, quietly=TRUE, warn.conflicts=FALSE), silent=TRUE);
try(require(corpcor, quietly=TRUE, warn.conflicts=FALSE), silent=TRUE);
require(BayesDiallel, quietly = TRUE, warn.conflicts=FALSE);
require(R.oo, warn.conflicts=FALSE, quietly=TRUE);

################################################################################
## Not Usuable in sample package
##  try(library(AlanDirectories), silent=TRUE)
##  try(AlanDirectores:::SetSaveHomes(),silent=TRUE);
##  try(SetSaveHomes());
##  SaveDir = BSSaveDir;  ## Save Dir set by default, or set to ones own accessible harddrive space.
 ###############################################################################
 
 ###############################################################################
 ### Simulate Data 
  NewWrite = 1; dfRobit = 7; dfRobit = -1;
  nop = 4; NoTau = 0;
  if (NoTau == 0) {
    ## In this version the Random effects are simulated by "Non-Random" means
    FirstRandom =5 + nop+1; tauEndList = FirstRandom + c(7,14,22,28);
    CFirstRandom = FirstRandom -1;  CtauEndList = tauEndList -1;
    PiAPrior = c(5, nop, 10,10); PiAStart = c(5/(nop+5), 10/20); 
  } else {
    ##
    FirstRandom = -1;  tauEndList = NULL; CFirstRandom = -1; CtauEndList = NULL;
    PiAPrior = c(15, nop+10);  PiAStart = c(15 / (nop+25));
  }
  ## Set Noise Size at something relative to signal, Robit-ness, Beta
  if (dfRobit >= 0) {
    n = 500;
    Sigma = 3;
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
  }  else {
    YY = Z;
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
  EEProbSortWidth = .5, DoLogitPostProb = FALSE,
  DoLongCI = TRUE);
  
  GetConfidenceIntervals = c(0.500,0.250, 0.750, 0.050, 0.950, 0.025, 0.975, 0.010, 0.990, 0.005, 0.995);
  ABS$TBSR5$KeepPosteriorQuantiles <- GetConfidenceIntervals;

  # Get credibility intervals caclulated on all ways and values.
   CIEst <- list();
	 try(CIEst[[1]] <- ABS$TBSR5$BetaSymmetricQuantiles);
	 try(CIEst[[2]] <- ABS$TBSR5$BetaSymmetricUnshrinkQuantiles);
	 try(CIEst[[3]] <- ABS$TBSR5$BetaHPDQuantiles);
	 try(CIEst[[4]] <- ABS$TBSR5$BetaHPDUnshrinkQuantiles);  
	 try(CIEst[[5]] <- ABS$TBSR5$BetaSymmetricLongQuantiles);
	 try(CIEst[[6]] <- ABS$TBSR5$BetaHPDLongQuantiles); 
   try(CIEst <- CleanCIEst(CIEst))   	   
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }
\keyword{ ts }% __ONLY ONE__ keyword per line
