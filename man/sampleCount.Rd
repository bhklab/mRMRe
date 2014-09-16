\name{sampleCount}
\Rdversion{1.1}
\docType{methods}
\alias{sampleCount}
\alias{sampleCount,mRMRe.Data-method}
\alias{sampleCount,mRMRe.Filter-method}
\alias{sampleCount,mRMRe.Network-method}

\title{
	Accessor function for the 'sampleCount' information in a mRMRe.Data, mRMRe.Filter and mRMRe.Network object.
}
\description{
	The feature count is simply the total number of samples considered in the mRMRe procedure.
}
\usage{
\S4method{sampleCount}{mRMRe.Data}(object)
\S4method{sampleCount}{mRMRe.Filter}(object)
\S4method{sampleCount}{mRMRe.Network}(object)
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
sampleCount(feature_data)
filter <- mRMR.classic("mRMRe.Filter", data = feature_data, target_indices = 3:5,
						feature_count = 2)
sampleCount(filter)
}
\keyword{methods}
