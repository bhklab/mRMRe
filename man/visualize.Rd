\name{visualize}
\Rdversion{1.1}
\docType{methods}
\alias{visualize}
\alias{visualize,mRMRe.Network-method}

\title{
mRMRe Network display}
\description{
The 'visualize' methods allows the visual display of an inferred mRMRe.Network topology.
}
\usage{
\S4method{visualize}{mRMRe.Network}(object)
}
\arguments{
  \item{object}{a \code{mRMRe.Network} object.}
}

\author{
Nicolas De Jay, Simon Papillon-Cavanagh, Benjamin Haibe-Kains
}
\examples{
set.thread.count(2)
data(cgps)
feature_data <- mRMR.data(data =  data.frame(cgps.ge))

# Build an mRMR-based network and display it
network <- new("mRMRe.Network", data = feature_data, target_indices = c(1),
				levels = c(3, 1), layers = 2)
visualize(network)
}
\keyword{methods}
