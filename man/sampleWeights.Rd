\name{sampleWeights}
\Rdversion{1.1}
\docType{methods}
\alias{sampleWeights}
\alias{sampleWeights<-}
\alias{sampleWeights,mRMRe.Data-method}
\alias{sampleWeights<-,mRMRe.Data-method}
\title{
	Accessor function for the 'sampleWeights' information in a mRMRe.Data object
}
\description{
TODO
}
\usage{
\S4method{sampleWeights}{mRMRe.Data}(object)
\S4method{sampleWeights}{mRMRe.Data}(object) <- value
}
\arguments{
  \item{object}{a \code{mRMRe.Data} object.}
  \item{value}{a numeric vector containing the biases of each sample in the mutual information computation.}
}
\author{
Nicolas De Jay, Simon Papillon-Cavanagh, Benjamin Haibe-Kains
}
\examples{
set.thread.count(2)
data(cgps)
feature_data <- mRMR.data(data = data.frame(cgps.ge))

# Uniform weight (default)
sampleWeights(feature_data)

# Random weighting  
sampleWeights(feature_data) <- runif(sampleCount(feature_data))
# Show result
sampleWeights(feature_data)
}
\keyword{methods}