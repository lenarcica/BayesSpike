\name{tNoise}
\alias{tNoise}
\alias{tnoise}
\alias{dfTNoise}
\alias{dftnoise}
\alias{DFTnoise}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{ T Distributed Noise Truncated Distribution Draws}
\description{
 Although Y_i would be on the real line, and a \code{\link{dfRobit}} is unnecessary, instead weights for each sample i
 are drawn each iteration based upon current \eqn{X\beta}{X\beta} estimates.  

}
\usage{
rTruncT(n=1, df = 1, L = NULL, U = NULL, sig=NULL, mu=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{Covariates, Independent Variable}
  \item{Y}{Dependent Variable}
  \item{dfNoise}{Degrees of freedom for T Noise error}
}
\details{
 Analysis object is a "Rcpp modules" object or imbedded C++ class with multiple
   access vectors
   
  To do the T Noise, a weight W_i is drawn for every sample in the data.  The equations
  for this sample are contained in Section 2.1 of the paper.
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
