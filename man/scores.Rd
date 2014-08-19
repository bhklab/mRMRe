\name{scores}
\Rdversion{1.1}
\docType{methods}
\alias{scores}
\alias{scores,mRMRe.Data-method}
\alias{scores,mRMRe.Filter-method}
\alias{scores,mRMRe.Network-method}
\title{
	mRMR Scores as per the MI gain for each feature}
\description{
	The \code{scores} method returns the scores of individual features in respect to previously selected features as per standard
	mRMR procedure.  For each target, the score of a feature is defined as the mutual information between the target and this feature
	minus the average mutual information of previously selected features and this feature.
}
\usage{
\S4method{scores}{mRMRe.Data}(object, solutions)
\S4method{scores}{mRMRe.Filter}(object)
\S4method{scores}{mRMRe.Network}(object)
}
\arguments{
  \item{object}{a \code{mRMRe.Data, mRMRe.Filter or mRMRe.Network} object.}
  \item{solutions}{a set of solutions from mRMRe.Filter or mRMRe.Network to be used in computing the scores from a mRMRe.Data set.}
}
\author{
Nicolas De Jay, Simon Papillon-Cavanagh, Benjamin Haibe-Kains
}
\examples{
set.thread.count(2)
data(cgps)
feature_data <- mRMR.data(data =  data.frame(cgps.ge))

# Create an mRMR filter and obtain the indices of selected features
filter <- mRMR.classic("mRMRe.Filter", data = feature_data, target_indices = 3:5,
						feature_count = 2)
scores(filter)
}
\keyword{methods}
