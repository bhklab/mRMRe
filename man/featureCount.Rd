\name{featureCount}
\Rdversion{1.1}
\docType{methods}
\alias{featureCount}
\alias{featureCount,mRMRe.Data-method}
\alias{featureCount,mRMRe.Filter-method}
\alias{featureCount,mRMRe.Network-method}
\title{
	Accessor function for the 'featureCount' information in a mRMRe.Data, mRMRe.Filter and mRMRe.Network object.
}
\description{
	The feature count is simply the total number of feature considered in the mRMRe procedure.
}
\usage{
\S4method{featureCount}{mRMRe.Data}(object)
\S4method{featureCount}{mRMRe.Filter}(object)
\S4method{featureCount}{mRMRe.Network}(object)
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
feature_data <- mRMR.data(data =  data.frame(cgps.ge))
featureCount(feature_data)
filter <- mRMR.classic("mRMRe.Filter", data = feature_data, target_indices = 3:5,
						feature_count = 2)
featureCount(filter)
}
\keyword{methods}
