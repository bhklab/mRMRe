\name{mim}
\Rdversion{1.1}
\docType{methods}
\alias{mim}
\alias{mim,mRMRe.Data-method}
\alias{mim,mRMRe.Network-method}
\alias{mim,mRMRe.Filter-method}
\title{
	Accessor function for the 'mim' information in a mRMRe.Data, mRMRe.Filter and mRMRe.Network object
}
\description{
	In both mRMRe.Filter and mRMRe.Network objects, a sparse mutual information matrix is computed for the mRMRe procedure and this lazy-evaluated matrix is returned.
	In the context of a a mRMRe.Data 'mim', the full pairwise mutual information matrix is computed and returned.
}
\usage{
\S4method{mim}{mRMRe.Data}(object, prior_weight, continuous_estimator, outX, bootstrap_count)
\S4method{mim}{mRMRe.Filter}(object, method)
\S4method{mim}{mRMRe.Network}(object)
}
\arguments{
  \item{object}{a \code{mRMRe.Data, mRMRe.Filter or mRMRe.Network} object.}
  \item{prior_weight}{a numeric value [0,1] of indicating the impact of priors (mRMRe.Data only).}
  \item{continuous_estimator}{an estimator of the mutual information between features: either "pearson", "spearman", "kendall", "frequency" (mRMRe.Data only).}
  \item{outX}{a boolean used in the concordance index estimator to keep or throw out ties (mRMRe.Data only).}
  \item{bootstrap_count}{an integer indicating the number of bootstrap resampling used in estimation (mRMRe.Data only). }
  \item{method}{either "mi" or "cor"; the latter will return the correlation coefficients (rho) while the former will return the mutual information (-0.5 * log(1 - (rho^2))). }
}
\author{
Nicolas De Jay, Simon Papillon-Cavanagh, Benjamin Haibe-Kains
}
\examples{
set.thread.count(2)
data(cgps)
feature_data <- mRMR.data(data =  data.frame(cgps.ge))

# Calculate the pairwise mutual information matrix
mim(feature_data)
filter <- mRMR.classic("mRMRe.Filter", data = feature_data, target_indices = 3:5,
						feature_count = 2)

# Obtain the sparse (lazy-evaluated) mutual information matrix.
mim(filter)
}
\keyword{methods}
