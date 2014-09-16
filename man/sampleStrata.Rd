\name{sampleStrata}
\Rdversion{1.1}
\docType{methods}
\alias{sampleStrata}
\alias{sampleStrata<-}
\alias{sampleStrata,mRMRe.Data-method}
\alias{sampleStrata<-,mRMRe.Data-method}

\title{
	Accessor function for the 'sampleStrata' information in a mRMRe.Data object
}
\description{
	The sampleStrata vector consists of a sampling stratification that will be used in computing mutual information
	between features.  If known batch effects or sample stratification is present between samples, identify such subsets
	using this.
}
\usage{
\S4method{sampleStrata}{mRMRe.Data}(object)
\S4method{sampleStrata}{mRMRe.Data}(object) <- value
}
\arguments{
  \item{object}{a \code{mRMRe.Data} object.}
  \item{value}{a factor vector identifying the stratification of samples.}
}
\author{
Nicolas De Jay, Simon Papillon-Cavanagh, Benjamin Haibe-Kains
}
\examples{
set.thread.count(2)
data(cgps)
feature_data <- mRMR.data(data = data.frame(cgps.ge))

# No stratification (default)
sampleStrata(feature_data)

# Random stratification 
sampleStrata(feature_data) <- as.factor(sample(c(0,1), 
  sampleCount(feature_data), replace=TRUE))
# Show result
sampleStrata(feature_data)
}
\keyword{methods}
