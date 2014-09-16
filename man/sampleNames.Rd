\name{sampleNames}
\Rdversion{1.1}
\docType{methods}
\alias{sampleNames}
\alias{sampleNames,mRMRe.Filter-method}
\alias{sampleNames,mRMRe.Network-method}
\alias{sampleNames,mRMRe.Data-method}
\title{
	Accessor function for the 'sampleNames' information in a mRMRe.Data, mRMRe.Filter and mRMRe.Network object.
}
\description{
	sampleNames are the names of the samples given as input to the mRMRe procedure.
}
\usage{
\S4method{sampleNames}{mRMRe.Data}(object)
\S4method{sampleNames}{mRMRe.Filter}(object)
\S4method{sampleNames}{mRMRe.Network}(object)
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
sampleNames(feature_data)
filter <- mRMR.classic("mRMRe.Filter", data = feature_data, target_indices = 3:5,
						feature_count = 2)
sampleNames(filter)
}
\keyword{methods}
