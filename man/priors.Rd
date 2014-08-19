\name{priors}
\Rdversion{1.1}
\docType{methods}
\alias{priors}
\alias{priors<-}
\alias{priors,mRMRe.Data-method}
\alias{priors<-,mRMRe.Data-method}

\title{
	Accessor function for the 'priors' information in a mRMRe.Data object
}
\description{
	The priors matrix consists of a prior bias to be used in computation to mutual information between features.
}
\usage{
\S4method{priors}{mRMRe.Data}(object)
\S4method{priors}{mRMRe.Data}(object) <- value
}
\arguments{
  \item{object}{a \code{mRMRe.Data} object.}
  \item{value}{a numeric matrix containing values from 0 to 1 (or NA), one per pairwise feature bias.}
}
\author{
Nicolas De Jay, Simon Papillon-Cavanagh, Benjamin Haibe-Kains
}
\examples{
set.thread.count(2)
data(cgps)
feature_data <- mRMR.data(data = data.frame(cgps.ge))
priors(feature_data)
}
\keyword{methods}