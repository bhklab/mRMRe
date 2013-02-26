\name{mRMRe.Data-class}
\Rdversion{1.1}
\docType{class}
\alias{mRMRe.Data-class}
%\alias{compressFeatureIndices,mRMRe.Data-method}
%\alias{compressFeatureMatrix,mRMRe.Data-method}
%\alias{expandFeatureIndices,mRMRe.Data-method}
%\alias{expandFeatureMatrix,mRMRe.Data-method}
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
	
	%% perhaps explain our algorithms %%
}
\section{Instantiation}{
	Objects are created via calls of the form \code{new("mRMRe.Data", data, strata, weights, priors)}.
	
	\code{data}: is expected to be a data frame with samples and features respectively organized as rows and columns.  The columns
	have to be of type :{numeric, ordered factor, Surv} and respectively interpreted as :{continuous, discrete and survival} variables.
	
	\code{strata}: is expected to be a vector of type :{ordered factor} with the strata associated to the samples provided
	in \code{data}.
	
	\code{weights}: is expected to be a vector of type :{numeric} with the weights associated to the samples provided
	in \code{data}.
	
	\code{priors}: is expected to be a matrix of type :{numeric} where \code{priors[i, j]}: denotes an forced association between
	features i and j in \code{data}.  The latter takes into consideration the directionality of the relationship and must be a value
	between 0 and 1.
}
\section{Mutual Information Matrix}{
	The \code{mim} method computes and returns a mutual information matrix.  A correlation between continuous features is estimated
	using an estimator specified in \code{continuous_estimator}; currently, :{pearson, spearman, kendall, frequency} are supported.
	The estimator for discrete features is Cramér's V and for all other combinations, concordance index.
	
	When \code{outX} is set to \code{TRUE}, ties are ignored when computing the concordance index and otherwise, these are considered.
	The correlations are first computed per strata and these are then combined by the inverse variance weight mean of the estimates
	using a \code{bootstrap_count} number of bootstraps if the former parameter is greater than 0, and by the relative weights of each
	strata otherwise.  The resulting correlation is then summated with the corresponding value in the priors matrix with the
	latter being weighed for a proportion \code{prior_weight} of a final, biased correlation.
}
\section{Slots}{
  \describe{
    \item{\code{sample_names}:}{Object of class \code{"character"} containing the sample names.}
    \item{\code{feature_names}:}{Object of class \code{"character"} containing the feature names.}
    \item{\code{feature_types}:}{Object of class \code{"numeric"} containing the internal representation of features/variables.}
    \item{\code{data}:}{Object of class \code{"matrix"} containing the internal representation of the data set.}
    \item{\code{strata}:}{Object of class \code{"numeric"} containing the feature strata.}
    \item{\code{weights}:}{Object of class \code{"numeric"} containing sample weights.}
    \item{\code{priors}:}{Object of class \code{"matrix"} containing the priors.}
  }
}
\section{Methods}{
  \describe{
%    \item{compressFeatureIndices}{\code{signature(object = "mRMRe.Data")}: ... }
%    \item{compressFeatureMatrix}{\code{signature(object = "mRMRe.Data")}: Compresses multidimensional }
%    \item{expandFeatureIndices}{\code{signature(object = "mRMRe.Data")}: }
%    \item{expandFeatureMatrix}{\code{signature(object = "mRMRe.Data")}: Expands multidimensional feature (of type Surv) in unidimensional columns to simplify data passing between R and C}
    \item{featureCount}{\code{signature(object = "mRMRe.Data")}: Returns the number of features.}
    \item{featureData}{\code{signature(object = "mRMRe.Data")}: Returns a data frame corresponding to the data set.}
    \item{featureNames}{\code{signature(object = "mRMRe.Data")}: Returns a vector containing the feature names.}
    \item{mim}{\code{signature(object = "mRMRe.Data", prior_weight = 0,
    continuous_estimator = c("pearson", "spearman", "kendall", "frequency"), outX = TRUE, bootstrap_count = 0)}:
        Computes and returns the mutual information matrix.}
    \item{priors}{\code{signature(object = "mRMRe.Data")}: Returns a matrix containing the priors.}
    \item{priors<-}{\code{signature(object = "mRMRe.Data", value)}: Sets the prior matrix.}
    \item{sampleCount}{\code{signature(object = "mRMRe.Data")}: Returns the number of samples.}
    \item{sampleNames}{\code{signature(object = "mRMRe.Data")}: Returns a vector containing sample names.}
    \item{sampleStrata}{\code{signature(object = "mRMRe.Data")}: Returns a vector containing sample strata.}
    \item{sampleStrata<-}{\code{signature(object = "mRMRe.Data", value)}: Sets the sample strata.}
    \item{sampleWeights}{\code{signature(object = "mRMRe.Data")}: Returns a vector containing sample weights.}
    \item{sampleWeights<-}{\code{signature(object = "mRMRe.Data")}: Sets the sample weights.}
    \item{subsetData}{\code{signature(object = "mRMRe.Data", row_indices, column_indices)}: Returns another data object containing only the specified samples and features (rows and columns, respectively.)}
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
