\name{mRMRe.Filter-class}
\Rdversion{1.1}
\docType{class}
\alias{mRMRe.Filter-class}
\alias{mRMR.ensemble}
\alias{mRMR.classic}
\title{Class \code{"mRMRe.Filter"}}
\description{
	mRMRe.Filter is a wrapper for various variants of the maximum relevance minimum redundancy (mRMR) feature selection/filter.
   
   Note that \code{mRMR.classic} and \code{mRMR.ensemble} functions are wrappers to easily perform classical (single) and ensemble mRMR feature selection.
}
\section{Instantiation}{
	Objects are created via calls of the form \code{new("mRMRe.Filter", data, prior_weight,
	target_indices,	levels,	method, continuous_estimator, outX, bootstrap_count)}.
	
	\code{data}: is expected to be a \code{mRMRe.Data} object.
	
	\code{target_indices}: is expected to be a vector of type \code{integer} containing the indices of the features
	that will serve as targets for the feature selections.
	
	\code{levels}: is expected to be a vector of type \code{integer} containing the number of children of each element
	at each level of the resulting filter tree.
	
	\code{method}: is expected to be either \code{exhaustive} or \code{bootstrap}.  The former uses the whole dataset to pick siblings in the tree according
	to the mRMR metric, while the latter perform the classical mRMR feature selection on several bootrstap selections of the dataset.
   
   \code{continuous_estimator}: it specifies the estimators for correlation between two continuous variables; value is either \code{pearson}, \code{spearman}, \code{kendall}, \code{frequency}, 
   
   \code{outX}: set to \code{TRUE} (default value) to not count pairs of 
   observations tied on \code{x} as a relevant pair. This results in a 
   Goodman-Kruskal gamma type rank correlation.
   
   \code{bootstrap_count}: Number of bootstraps to statistically compare
   the mRMR scores of each solution.
 	
	Since a mutual information matrix must be computed in order for feature selection to take place, the remaining arguments
	are identical to those required by the \code{mim} method of the \code{mRMRe.Data} object.
}
\section{Slots}{
  \describe{
    \item{\code{filters}:}{Object of class \code{"list"} containing for each target a solutions matrix.}
    \item{\code{mi_matrix}:}{Object of class \code{"matrix"} containing the combined mutual information matrix of the relevant targets.}
    \item{\code{causality_list}:}{Object of class \code{"list"} containing for each target a vector of causality coefficients between the target and its predictors.}
    \item{\code{sample_names}:}{Object of class \code{"character"} containing the sample names.}
    \item{\code{feature_names}:}{Object of class \code{"character"} containing the feature names.}
    \item{\code{target_indices}:}{Object of class \code{"integer"} containing the target indices.}
    \item{\code{fixed_feature_count}:}{Object of class \code{"integer"} containing the number of fixed features.}
    \item{\code{levels}:}{Object of class \code{"integer"} containing the desired topology of the tree.}
    \item{\code{scores}:}{Object of class \code{"list"} containing the mRMR score of selected features, respective to filters.}
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
    \item{solutions}{\code{signature(object = "mRMRe.Filter", mi_threshold = -Inf,
    						causality_threshold = Inf)}:
    Returns a matrix in which each column represents a different solution (path from root of the tree to a leaf.)
    }
    \item{target}{\code{signature(object = "mRMRe.Filter")}: Returns a vector containing the target indices.}
  }
}
\references{
Ding, C. and Peng, H. (2005). "Minimum redundancy feature selection from microarray gene expression data". \emph{Journal of bioinformatics and computational biology}, \bold{3}(2):185--205.
}
\author{
Nicolas De Jay, Simon Papillon-Cavanagh, Benjamin Haibe-Kains
}
\seealso{
\code{\link[mRMRe]{mRMRe.Data-class}}
}
\examples{
showClass("mRMRe.Filter")

set.thread.count(2)

## load data
data(cgps)

## build an mRMRe.Data object
ge <- mRMR.data(data = data.frame(cgps.ge[ , 1:100, drop=FALSE]))

## perform a classic (single) mRMR to select the 10 genes the most correlated with 
## the first gene but the less correlated between each other
exect <- system.time(fs <- new("mRMRe.Filter", data = ge, target_indices = 1,
					levels = c(8, 1, 1, 1, 1)))
print(exect)

## print the index of the selected features for each distinct mRMR solutions
print(solutions(fs)[[1]])

## print the names of the selected features for each distinct mRMR solutions
print(apply(solutions(fs)[[1]], 2, function(x, y) { return(y[x]) }, y=featureNames(ge)))
}
\keyword{classes}
