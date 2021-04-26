\name{solutions}
\Rdversion{1.1}
\docType{methods}
\alias{solutions}
\alias{solutions,mRMRe.Filter-method}
\alias{solutions,mRMRe.Network-method}
\title{
Basic result of the mRMR procedure}
\description{
The 'solutions' method allows one to access the set of selected features resulting of the mRMR algorithm.  More generally,
the set of feature are identified by their indices in the inputed feature set (1 being the first feature (column)).  At the
network level, 'solutions' consists of the topology of the network, identifying which features is connected to others.
}
\usage{
\S4method{solutions}{mRMRe.Filter}(object, mi_threshold, causality_threshold, with_fixed_features)
\S4method{solutions}{mRMRe.Network}(object)
}
\arguments{
  \item{object}{a \code{mRMRe.Filter or mRMRe.Network} object.}
  \item{mi_threshold}{a numeric value used in filtering the features based on their mRMR scores,
  		features that do not pass the threshold will be set at NA.}
  \item{causality_threshold}{a numeric value used in filtering the features based on their causality scores,
  		features that do not pass the threshold will be set at NA}
  \item{with_fixed_features}{a boolean indicating if fixed features are used in the computation, default TRUE}
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
solutions(filter)

# Build an mRMR-based network and obtain feature connections (topology)
network <- new("mRMRe.Network", data = feature_data, target_indices = c(1, 2),
			levels = c(2, 1), layers = 1)
solutions(network)
}
\keyword{methods}

