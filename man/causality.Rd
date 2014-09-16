\name{causality}
\Rdversion{1.1}
\docType{methods}
\alias{causality}
\alias{causality,mRMRe.Filter-method}
\alias{causality,mRMRe.Network-method}
\title{
	Accessor function for the 'causality' information in a mRMRe.Filter and mRMRe.Network object.
}
\description{
	The causality data is compute using the co-information lattice algorithm on each V-structure (feature, target, feature).  Given that this procedure
	is computed for each pair of features, the minimum result is kept. A negative score indicates putative causality of the feature to the target.
}
\usage{
\S4method{causality}{mRMRe.Filter}(object)
\S4method{causality}{mRMRe.Network}(object)
}
\arguments{
  \item{object}{a \code{mRMRe.Filter or mRMRe.Network} object.}
}
\author{
Nicolas De Jay, Simon Papillon-Cavanagh, Benjamin Haibe-Kains
}
\examples{
set.thread.count(2)
data(cgps)
feature_data <- mRMR.data(data =  data.frame(cgps.ge))
filter <- mRMR.classic("mRMRe.Filter", data = feature_data, target_indices = 3:5,
						feature_count = 2)
causality(filter)
}
\keyword{methods}