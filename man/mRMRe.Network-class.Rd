\name{mRMRe.Network-class}
\Rdversion{1.1}
\docType{class}
\alias{mRMRe.Network-class}
\alias{adjacencyMatrix,mRMRe.Network-method}
\alias{adjacencyMatrixSum,mRMRe.Network-method}
\alias{causality,mRMRe.Network-method}
\alias{featureNames,mRMRe.Network-method}
\alias{mim,mRMRe.Network-method}
\alias{sampleNames,mRMRe.Network-method}
\alias{solutions,mRMRe.Network-method}
\alias{visualize,mRMRe.Network-method}

\title{Class \code{"mRMRe.Network"}}
\description{
	mRMRe.Network is a wrapper for inferring a network of features based on mRMR feature selection.

	%% perhaps explain our algorithms %%
}
\section{Instantiation}{
	Objects are created via calls of the form
	\code{new("mRMRe.Network", data, prior_weight, target_indices, levels, layers, ..., mi_threshold, causality_threshold)}.
	
	\code{layers}: is expected to be an \code{integer} specifying the number of layers of network inference desired.  When multiple
	layers are desired, the elements of the solutions found in the last step of feature selection are used as the targets of the next
	step.
	
	Since networking involves filter processing, the remaining arguments are identical to those required by \code{solutions} method of
	the \code{mRMRe.Filter} object and \code{mim} method of the \code{mRMRe.Data} object.
}
\section{Slots}{
  \describe{
    \item{\code{topologies}:}{Object of class \code{"list"} ~~ }
    \item{\code{mi_matrix}:}{Object of class \code{"matrix"} containing the combined mutual information matrix of the network elements.}
    \item{\code{causality_list}:}{Object of class \code{"list"} containing for each target a vector of causality coefficients between the target and its predictors.}
    \item{\code{sample_names}:}{Object of class \code{"character"} containing the sample names.}
    \item{\code{feature_names}:}{Object of class \code{"character"} containing the feature names.}
    \item{\code{target_indices}:}{Object of class \code{"integer"} containing the target indices.}
  }
}
\section{Methods}{
  \describe{
    \item{adjacencyMatrix}{\code{signature(object = "mRMRe.Network")}: Returns a matrix describing the topology of the network.}
    \item{adjacencyMatrixSum}{\code{signature(object = "mRMRe.Network")}: ... }
    \item{causality}{\code{signature(object = "mRMRe.Network")}: Returns a list containing vectors containing causality coefficients between targets and predictors.}
    \item{featureNames}{\code{signature(object = "mRMRe.Network")}: Returns a vector containing the feature names.}
    \item{mim}{\code{signature(object = "mRMRe.Network")}: ... }
    \item{sampleNames}{\code{signature(object = "mRMRe.Network")}: Returns a vector containing sample names.}
    \item{solutions}{\code{signature(object = "mRMRe.Network")}: ... }
    \item{visualize}{\code{signature(object = "mRMRe.Network")}: ... }
	 }
}
\references{
}
\author{
}
\note{
}
\seealso{
\code{\link{mRMRe.Filter-class}}, \code{\link{mRMRe.Data-class}}
}
\examples{
showClass("mRMRe.Network")
}
\keyword{classes}
