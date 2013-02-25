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
	mRMRe.Data is a wrapper for data sets.  Most if not all of the routines in the mRMRe package use mRMRe.Data as primary input.
	Such an object is instantiated with a data frame containing the sample sets and optionally, stratum, weight vectors and a prior
	matrix.  In addition to basic accession functions, we describe several methods which serve to manipulate the contents of the data
	set.
}
\section{Instantiation}{
	Objects are created via calls of the form \code{new("mRMRe.Data", data, strata, weights, priors)}.
	
	\code{"data"}: is expected to be a data frame with samples and features respectively organized as rows and columns.  The columns
	have to be of type :{numeric, ordered factor, Surv} and respectively interpreted as :{continuous, discrete and survival} variables.
	
	\code{"strata"}: is expected to be a vector of type :{ordered factor} with the strata associated to the samples provided
	in \code{"data"}.
	
	\code{"weights"}: is expected to be a vector of type :{numeric} with the weights associated to the samples provided
	in \code{"data"}.
	
	\code{"priors"}: is expected to be a matrix of type :{numeric} where \code{"priors[i, j]"}: denotes an forced association between
	features i and j in \code{"data"}.  Directionality is important.   
}
\section{Slots}{
  \describe{
    \item{\code{sample_names}:}{Object of class \code{"character"} Vector containing the sample names}
    \item{\code{feature_names}:}{Object of class \code{"character"} Vector containing the feature names}
    \item{\code{feature_types}:}{Object of class \code{"numeric"} Vector containing the internal representation of features/variables}
    \item{\code{data}:}{Object of class \code{"matrix"} Matrix containing the internal representation of the data set}
    \item{\code{strata}:}{Object of class \code{"numeric"} Vector containing the feature strata}
    \item{\code{weights}:}{Object of class \code{"numeric"} Vector containing sample weights}
    \item{\code{priors}:}{Object of class \code{"matrix"} Matrix containing the priors}
  }
}
\section{Methods}{
  \describe{
%    \item{compressFeatureIndices}{\code{signature(object = "mRMRe.Data")}: ... }
%    \item{compressFeatureMatrix}{\code{signature(object = "mRMRe.Data")}: Compresses multidimensional }
%    \item{expandFeatureIndices}{\code{signature(object = "mRMRe.Data")}: }
%    \item{expandFeatureMatrix}{\code{signature(object = "mRMRe.Data")}: Expands multidimensional feature (of type Surv) in unidimensional columns to simplify data passing between R and C}
    \item{featureCount}{\code{signature(object = "mRMRe.Data")}: Returns the number of features}
    \item{featureData}{\code{signature(object = "mRMRe.Data")}: Returns a data set corresponding to the data set}
    \item{featureNames}{\code{signature(object = "mRMRe.Data")}: Returns a vector containing the feature names}
    \item{mim}{\code{signature(object = "mRMRe.Data", prior_weight = 0,
    continuous_estimator = c("pearson", "spearman", "kendall", "frequency"), outX = TRUE, bootstrap_count = 0)}:
        Computes and returns the mutual information matrix.
   	    \code{"prior_weight"}: expects a value between 0 and 1 that denotes the relative importance of the values in the prior matrix
   	    with respect to the values computed via the estimators.
   	    \code{"continuous_estimator"}: expects an estimator :{pearson, spearman, kendall, frequency} that is to be used for estimating
   	    association between two continuous features
   	    \code{"outX"}: ...
   	    \code{"bootstrap_count"}: ...
   	    }
    \item{priors}{\code{signature(object = "mRMRe.Data")}: Returns a matrix containing the priors}
    \item{priors<-}{\code{signature(object = "mRMRe.Data", value)}: Sets the prior matrix}
    \item{sampleCount}{\code{signature(object = "mRMRe.Data")}: Returns the number of samples}
    \item{sampleNames}{\code{signature(object = "mRMRe.Data")}: Returns a vector containing sample names}
    \item{sampleStrata}{\code{signature(object = "mRMRe.Data")}: Returns a vector containing sample strata}
    \item{sampleStrata<-}{\code{signature(object = "mRMRe.Data", value)}: Sets the sample strata}
    \item{sampleWeights}{\code{signature(object = "mRMRe.Data")}: Returns a vector containing sample weights}
    \item{sampleWeights<-}{\code{signature(object = "mRMRe.Data")}: Sets the sample weights}
    \item{subsetData}{\code{signature(object = "mRMRe.Data", row_indices, column_indices)}: Returns another mRMRe.Data object but with only the specified row and column (sample and feature) indices}
	 }
}
\references{

}
\author{

}
\note{

}
\seealso{

}
\examples{
showClass("mRMRe.Data")
}
\keyword{classes}
