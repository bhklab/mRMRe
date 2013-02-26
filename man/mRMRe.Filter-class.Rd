\name{mRMRe.Filter-class}
\Rdversion{1.1}
\docType{class}
\alias{mRMRe.Filter-class}
\alias{causality,mRMRe.Filter-method}
\alias{featureCount,mRMRe.Filter-method}
\alias{featureNames,mRMRe.Filter-method}
\alias{mim,mRMRe.Filter-method}
\alias{sampleCount,mRMRe.Filter-method}
\alias{sampleNames,mRMRe.Filter-method}
\alias{solutions,mRMRe.Filter-method}
\alias{target,mRMRe.Filter-method}

\title{Class \code{"mRMRe.Filter"}}
\description{
	mRMRe.Filter is a wrapper for various feature selection (filter.)
	
	%% perhaps explain our algorithms %%
}
\section{Instantiation}{
	Objects are created via calls of the form \code{new("mRMRe.Filter", data, prior_weight, target_indices, levels,
	method, continuous_estimator, outX, bootstrap_count)}.
	
	\code{data}: is expected to be a \code{mRMRe.Data} object.
	
	\code{target_indices}: is expected to be a vector of type \code{numeric} containing the indices of the features
	that will serve as targets for the feature selections.
	
	\code{levels}: is expected to be a vector of type \code{integer} containing the number of children of each element
	at each level of the resulting filter tree.
	
	\code{method}: is expected to be either :{exhaustive, bootstrap}.  The former picks siblings in the tree according
	to the mRMR metric and the latter, ...
	
	Since a mutual information matrix must be computed in order for feature selection to take place, the remaining arguments
	are identical to those required by the \code{mim} method of the \code{mRMRe.Data} object.
}
\section{Slots}{
  \describe{
    \item{\code{filters}:}{Object of class \code{"list"} containing for each target a solutions matrix.}
    \item{\code{mi_matrix}:}{Object of class \code{"matrix"} containing the combined MIM resulting of feature selection for the relevant targets.}
    \item{\code{causality_list}:}{Object of class \code{"list"} containing for each target a vector of causality coefficients between the target and its predictors.}
    \item{\code{sample_names}:}{Object of class \code{"character"} containing the sample names.}
    \item{\code{feature_names}:}{Object of class \code{"character"} containing the feature names.}
    \item{\code{target_indices}:}{Object of class \code{"integer"} containing the target indices.}
    \item{\code{levels}:}{Object of class \code{"integer"} containing the desired topology of the tree.}
  }
}
\section{Methods}{
  \describe{
    \item{causality}{\code{signature(object = "mRMRe.Filter")}: ... }
    \item{featureCount}{\code{signature(object = "mRMRe.Filter")}: Returns the number of features.}
    \item{featureNames}{\code{signature(object = "mRMRe.Filter")}: Returns a vector containing the feature names.}
    \item{mim}{\code{signature(object = "mRMRe.Filter")}: Returns the potentially partial mutual information matrix used for feature selection.}
    \item{sampleCount}{\code{signature(object = "mRMRe.Filter")}: Returns the number of samples.}
    \item{sampleNames}{\code{signature(object = "mRMRe.Filter")}: Returns a vector containing sample names.}
    \item{solutions}{\code{signature(object = "mRMRe.Filter", mi_threshold = -Inf, causality_threshold = Inf)}:
    Returns a matrix in which each column represents a different solution (path from root of the tree to a leaf.)
    The parameters serve respectively as high- and low- pass filters. %% Not sure how it's called? %%
    }
    \item{target}{\code{signature(object = "mRMRe.Filter")}: Returns a vector containing the target indices.}
  }
}
\references{

}
\author{

}
\note{

}
\seealso{
\code{\link{mRMRe.Data-class}}
}
\examples{
showClass("mRMRe.Filter")
}
\keyword{classes}
