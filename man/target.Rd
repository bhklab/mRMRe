\name{target}
\Rdversion{1.1}
\docType{methods}
\alias{target}
\alias{target,mRMRe.Filter-method}
\alias{target,mRMRe.Network-method}
\title{
TODO
}
\description{
TODO
}
\usage{
\S4method{target}{mRMRe.Filter}(object)
\S4method{target}{mRMRe.Network}(object)
}
\arguments{
  \item{object}{a \code{mRMRe.Filter or mRMRe.Network} object.}
}
\author{
Nicolas De Jay, Simon Papillon-Cavanagh, Benjamin Haibe-Kains
}
\examples{
data(cgps)
feature_data <- mRMR.data(data =  data.frame(cgps.ge))

# Create an mRMR filter and obtain the targets of that filter
filter <- mRMR.classic("mRMRe.Filter", data = feature_data, target_indices = 3:5, feature_count = 2)
target(filter)

# Build an mRMR-based network and obtain targets (seeds) of the network
network <- new("mRMRe.Network", data = feature_data, target_indices = c(1, 2), levels = c(2, 1), layers = 1)
target(network)
}
\keyword{methods}