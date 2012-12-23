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
	This class is a wrapper for data sets.  mRMRe operates using this mRMRe.Data class as input in order to simplify parameter passing.  Instanciate an mMRMe.Data object with a data frame contaning your data and optionally with strata, weights and priors.  In addition to basic accession functions, the mim() function can be called on an instance of this class that will produce a mutual information matrix based on the demanded estimator, strata, sample weights and priors.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("mRMRe.Data", data, strata, weights, priors)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{sample_names}:}{Object of class \code{"character"} Names of the samples for data rows}
    \item{\code{feature_names}:}{Object of class \code{"character"} Names of the feature for data columns}
    \item{\code{feature_types}:}{Object of class \code{"numeric"} Class of each feature (column), accepted values :{numeric, ordered factor, Surv}}
    \item{\code{data}:}{Object of class \code{"matrix"} The actual data set containing features as columns and obsverations as rows}
    \item{\code{strata}:}{Object of class \code{"numeric"} Sample stratification}
    \item{\code{weights}:}{Object of class \code{"numeric"}  Sample weights}
    \item{\code{priors}:}{Object of class \code{"matrix"} Feature priors (or weight) }
  }
}
\section{Methods}{
  \describe{
    \item{compressFeatureIndices}{\code{signature(object = "mRMRe.Data")}: ... }
    \item{compressFeatureMatrix}{\code{signature(object = "mRMRe.Data")}: Compresses multidimensional }
    \item{expandFeatureIndices}{\code{signature(object = "mRMRe.Data")}: }
    \item{expandFeatureMatrix}{\code{signature(object = "mRMRe.Data")}: Expands multidimensional feature (of type Surv) in unidimensional columns to simplify data passing between R and C}
    \item{featureCount}{\code{signature(object = "mRMRe.Data")}: Returns the number of features in the data set}
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
