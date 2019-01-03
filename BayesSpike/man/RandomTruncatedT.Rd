\name{RandomTruncatedT}
\alias{rTT}
\alias{RTT}
\alias{RandomTruncT}
\alias{RandomTruncatedT}
\alias{RandomTruncNorm}
\alias{rTruncNorm}
\alias{rTruncT}
\alias{rtrunct}
\alias{rtruncnorm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Random Truncated Distribution Draws}
\description{
 These functions are important to a a Probit/Robit sampler, as draws from truncated distribution
 become the latent continuous measurements. 
 
 Some analysts with Binary 0/1 "Z" data will want to fit a Robit regression to the data.
 Generating draws from a truncated long tailed T distribution exactly is challenging, and for our purposes
 we only need an approximate solution that at least captures features of the long tail.  By slice sampling
 from the truncated t density we expect to get a reasonable approximation of an exact draw.  
 
 During the draw a "Z_i" is replaced by a new draw "Y_i" which is in the continuous Real Numbers and represents
 the truncated T draw (on right side of the regression line are the Z_i=1, on left are Z_i=0).
}
\usage{
rTruncT(n=1, df = 1, L = NULL, U = NULL, sig=NULL, mu=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
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
  \item{dfRobit}{Set this to do Robit regression (Y entry must be 0-1)}
}
\details{
 Analysis object is a "Rcpp modules" object or imbedded C++ class with multiple
   access vectors
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
library(BayesSpike);

 n = 20000;  snu = 7; ssig = 1; sL = rep(-100,n);  sU = rep(-99,n);
 sRet = rep(0, n); smu = rep(0, 100)
 t1 = proc.time();
 .Call("rTT", sRet, snu, smu, ssig, sL, sU)
 t2 = proc.time();
 
 rTruncT(df=snu, mu = smu, sig = ssig, L = sL, U = sU);
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }
\keyword{ ts }% __ONLY ONE__ keyword per line
