\name{mRMRe.Data-class}
\Rdversion{1.1}
\docType{class}
\alias{mRMRe.Data-class}
\alias{compressFeatureIndices,mRMRe.Data-method}
\alias{compressFeatureMatrix,mRMRe.Data-method}
\alias{expandFeatureIndices,mRMRe.Data-method}
\alias{expandFeatureMatrix,mRMRe.Data-method}
\alias{featureCount,mRMRe.Data-method}
\alias{featureData,mRMRe.Data-method}
\alias{featureNames,mRMRe.Data-method}
\alias{mim,mRMRe.Data-method}
\alias{priors,mRMRe.Data-method}
\alias{priors<-,mRMRe.Data-method}
\alias{sampleCount,mRMRe.Data-method}
\alias{sampleNames,mRMRe.Data-method}
\alias{sampleStrata,mRMRe.Data-method}
\alias{sampleStrata<-,mRMRe.Data-method}
\alias{sampleWeights,mRMRe.Data-method}
\alias{sampleWeights<-,mRMRe.Data-method}
\alias{subsetData,mRMRe.Data-method}

\title{Class \code{"mRMRe.Data"}}
\description{
%%  ~~ A concise (1-5 lines) description of what the class is. ~~
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("mRMRe.Data", data, strata, weights, priors)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{sample_names}:}{Object of class \code{"character"} ~~ }
    \item{\code{feature_names}:}{Object of class \code{"character"} ~~ }
    \item{\code{feature_types}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{data}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{strata}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{weights}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{priors}:}{Object of class \code{"matrix"} ~~ }
  }
}
\section{Methods}{
  \describe{
    \item{compressFeatureIndices}{\code{signature(object = "mRMRe.Data")}: ... }
    \item{compressFeatureMatrix}{\code{signature(object = "mRMRe.Data")}: ... }
    \item{expandFeatureIndices}{\code{signature(object = "mRMRe.Data")}: ... }
    \item{expandFeatureMatrix}{\code{signature(object = "mRMRe.Data")}: ... }
    \item{featureCount}{\code{signature(object = "mRMRe.Data")}: ... }
    \item{featureData}{\code{signature(object = "mRMRe.Data")}: ... }
    \item{featureNames}{\code{signature(object = "mRMRe.Data")}: ... }
    \item{mim}{\code{signature(object = "mRMRe.Data")}: ... }
    \item{priors}{\code{signature(object = "mRMRe.Data")}: ... }
    \item{priors<-}{\code{signature(object = "mRMRe.Data")}: ... }
    \item{sampleCount}{\code{signature(object = "mRMRe.Data")}: ... }
    \item{sampleNames}{\code{signature(object = "mRMRe.Data")}: ... }
    \item{sampleStrata}{\code{signature(object = "mRMRe.Data")}: ... }
    \item{sampleStrata<-}{\code{signature(object = "mRMRe.Data")}: ... }
    \item{sampleWeights}{\code{signature(object = "mRMRe.Data")}: ... }
    \item{sampleWeights<-}{\code{signature(object = "mRMRe.Data")}: ... }
    \item{subsetData}{\code{signature(object = "mRMRe.Data")}: ... }
	 }
}
\references{
%%  ~~put references to the literature/web site here~~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%%  ~~objects to See Also as \code{\link{~~fun~~}}, ~~~
%%  ~~or \code{\linkS4class{CLASSNAME}} for links to other classes ~~~
}
\examples{
showClass("mRMRe.Data")
}
\keyword{classes}
