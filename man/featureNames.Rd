\name{featureNames}
\Rdversion{1.1}
\docType{methods}
\alias{featureNames}
\alias{featureNames,mRMRe.Filter-method}
\alias{featureNames,mRMRe.Network-method}
\alias{featureNames,mRMRe.Data-method}
\title{
	Accessor function for the 'featureNames' information in a mRMRe.Data, mRMRe.Filter and mRMRe.Network object
}
\description{
	featureNames are the names of the features given as input to the mRMRe procedure.
}
\usage{
\S4method{featureNames}{mRMRe.Data}(object)
\S4method{featureNames}{mRMRe.Filter}(object)
\S4method{featureNames}{mRMRe.Network}(object)
}
\arguments{
  \item{object}{a \code{mRMRe.Data, mRMRe.Filter or mRMRe.Network} object.}
}
\author{
Nicolas De Jay, Simon Papillon-Cavanagh, Benjamin Haibe-Kains
}
\examples{
set.thread.count(2)
data(cgps)
feature_data <- mRMR.data(data = data.frame(cgps.ge))
featureNames(feature_data)
filter <- mRMR.classic("mRMRe.Filter", data = feature_data, target_indices = 3:5,
						feature_count = 2)
featureNames(filter)
}
\keyword{methods}
