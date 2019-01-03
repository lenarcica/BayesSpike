\name{FullDiallel}
\alias{FullDiallel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Create HappyR Object }
\description{
  Create "FullDiallel" either Analysis or Simulation object
}
\usage{
FullDiallel(numj=defaultnumj,  InListjk = NULL, InSexVector=NULL,
  InFixedEffects = NULL, nameFixedEffects = NULL,
  InRandomEffects = NULL, nameRandomEffects = NULL,
  InRandomPrintEffects=NULL, dfTNoise = NULL,
  Indf = NULL, Inm = NULL, IntauFixed=NULL, 
  Indftau = NULL,Inmtau = NULL,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{numj}{ Number of total Strains in Diallel}
  \item{Indf}{ A "RConstantsList", object representing prior inverse gamma degrees of freedom for parameters. 
    default NULL sets all values at .5}
  \item{Inm}{A "RConstantsList",object representing prior inverse gamma scaling for parameters. 
    default NULL sets all values at .5 }
  \item{dfTNoise}{If dfTNoise >0, then a outlier procedure based upon t-distributed noise is used}
  \item{InSexVector}{A vector declaring 0 if mouse is female, 1 if male}
  \item{InListjk}{A two column list, first column, strain of mother, second column, strain of father}
  \item{IntauFixed}{In Starting parameters for the fixed coefficients}
  \item{verbose}{ Decision to print actions to console }
  \item{Indftau}{Hierarchical prior information}
  \item{Inmtau}{Hierarchical prior information}
}
\details{
  All Analyses will use "FullDiallelAnalyze", an inherited class of Full Diallel
}
\value{
    Object  
     \item{FullDiallel Object}{An S3 Class designed to test multiple diallel models on the same dataset}
}
\references{ Diallel Paper }
\author{ Alan Lenarcic, Will Valdar }
\note{ 
  "FullDiallelAnalyze"  is a class that will analyze a dataset
}
\seealso{ FullDiallelAnalyze, FullDiallelSimulate, DiallelOb }
\examples{
## Load in R packages
try(require(lattice, quietly=TRUE, warn.conflicts=FALSE), silent=TRUE);
try(require(coda, quietly=TRUE, warn.conflicts=FALSE), silent=TRUE);
try(require(R.methodsS3, quietly=TRUE, warn.conflicts=FALSE), silent=TRUE);
try(require(R.oo, quietly=TRUE, warn.conflicts=FALSE), silent=TRUE);
try(require(corpcor, quietly=TRUE, warn.conflicts=FALSE), silent=TRUE);
require(RDiallel, quietly = TRUE, warn.conflicts=FALSE);
require(R.oo, warn.conflicts=FALSE, quietly=TRUE);
   
## Set the model shape to simulate
  ModelVector = c(1,1,2,1,3,1);
  ajModel = ModelVector[1];
  MotherModel = ModelVector[2];
  CrossModel = ModelVector[3];
  BetaHybridLevel = ModelVector[4];
  SexModel = ModelVector[5];
  BetaHybridTimesSex = ModelVector[6];
  IsFixed = 1; IsRandom = 1;
  NumReps = 2 ## Imagine dataset of size 2 times full diallel
  numj = 8 ## 8 strains of mice?
  
  
print("Sampling NonDiallel Random and Fake Effects");
InFixedEffects=NULL; InRandomEffects =NULL;
if (IsFixed || IsRandom) {
   WantLen = (numj *(numj+1)/2 ) *NumReps ;
   InRandomEffects = sample(size=WantLen, 1:5, replace=TRUE);
   InFixedEffects = rnorm(WantLen);
   if (IsFixed == 0) {InFixedEffects=NULL;}
   if (IsRandom == 0) {InRandomEffects=NULL;}
   mRandomEffects = 20;
   dfRandomEffects = 5; tauFixedEffects = 20;
   NumRepetitions=1;
}  else {
   NumRepetitions = NumReps;
}

## Simulate a data set given model of choice
print("Making first Simulation");
### Allocate AFDS Diallel Simulation Object
AFDS <- FullDiallelSimulate(numj = numj, InFixedEffects = InFixedEffects,
 InRandomEffects = InRandomEffects, mRandomEffects = 20,
 dfRandomEffects = 5, tauFixedEffects = 20, YesSex = TRUE, NumRepetitions = NumRepetitions);

## Sample the parameters for Diallel Simulation Object
defaultSimulateParam(AFDS);
if (is.null(InRandomEffects) && is.null(InFixedEffects)) {
  SetupFilledSimulatedListjk(AFDS, NumRepetitions = NumRepetitions);
}

## Set up the model for Diallel Simulation Object
SetupSimulateMethod(AFDS, ajModel = 1, MotherModel = MotherModel, CrossModel =
  CrossModel, SexModel = SexModel, BetaHybridLevel = BetaHybridLevel, 
  BetaHybridTimesSex = BetaHybridTimesSex);
## SimValues now simulates a random Beta vector, tau vector, and observed data Y      
SimValues(AFDS, Sigma = 1);

#############################################################################
###  Begin Analysis
###
### From Simulated Diallel in AFDS, we pull observed data to put in an 
###  Analysis Diallel AFD

Y = AFDS$Y; ListjkInput = AFDS$Listjk;
SexInputVec = AFDS$SexVector

## Pull Out the RandomEffects
ARandomEffects = rep(0, length(AFDS$RandomEffects[,1]));
for (ii in 1:length(AFDS$RandomEffects[1,])) {
  ARandomEffects[AFDS$RandomEffects[,ii] == 1] = ii;
}
if (is.null(InRandomEffects)) { ARandomEfects = NULL; }
##TargetLength=TARGET;
SampleBeta = AFDS$ADiallelOb$Beta;
FileDir = AllDir2;
quantile = .95;
INTIN = 0;

## These are the true simulation parameters from AFDS we would like experiment to learn  
TrueParametersList = c(AFDS$ADiallelOb$Beta, AFDS$ADiallelOb$tau,
    AFDS$ADiallelOb$Sigma)

#######################################################3
## Setup Analysis
##
## Allocate AFD : FullDiallelAnalyze Object
AFD  <- FullDiallelAnalyze(Y = Y, 
      SexVector = SexInputVec,
      Listjk = ListjkInput, dfStart = .5, mStart = .5,
      InFixedEffects = AFDS$FixedEffects,
      InRandomEffects = ARandomEffects);

## Assemble a List of models to be considered by AFD
##  This is a k * 6 matrix;
##
## In our case, let k = 1, and supply AFD the "true" simulation model 
ExampleModelsList = matrix(c(ajModel, MotherModel, CrossModel,
      SexModel, BetaHybridLevel, BetaHybridTimesSex), 1, 6);
      Sigma = 1;  
      names(Sigma) = "SigmaAll";

## Setup Analysis Methods for AFD, including allocating for MCMC chains
## "Sigma" is starting Sigma anticipated noise value
SetupAnalysisMethods(AFD, ExampleModelsList, Sigma = Sigma,
      numChains = 3, lengthChains = 2000);

#######################################################
## This will run the MCMC sampler on Model of choice
t1 = proc.time();
  RunChains(AFD);
t2 = proc.time();                  
print(t2-t1); ## Print Run time

ADO = AFD$AllDiallelObs[[1]];  ## This is the DiallelOb corresponding to our model

plot(ADO$CodaChains);  ## Plot the chains to understand
summary(ADO$CodaChains);  ## Get Summary information on a parameters
                                          
plot.hpd(ADO$CodaChains); ## Using Library_diallel_helper, plot chains.

plot(crosscorr(ADO$CodaChains)) ## Look at cross correlations of coda chains


CalcDICC(ADO, AFD) ## DIC value for this Coda Chain  

## This creates a summary table of applicable elements
MySum = summary.table(ADO, burnin = 400, AFD=AFD);

#MyPostSumGivesValues for all the diallel squares of every gender
MyPostSum = PosteriorPredSummary(ADO, AFD, burnin = 400, AFD=AFD, keep = TRUE);

SimmedGrid = GiveGrid(AFDS);
par(mfrow=c(1,2));
draw.diallel(data =SimmedGrid[MyPostSum$"Female/Male" == 1,], phenotype = "Mean.Simulated")
draw.diallel(data = MyPostSum[MyPostSum$"Female/Male" == 1,], phenotype = "Mean Posterior")
  
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }
\keyword{ ts }% __ONLY ONE__ keyword per line
