\name{featureNames}
\Rdversion{1.1}
\docType{methods}
\alias{featureNames}
\alias{featureNames,mRMRe.Filter-method}
\alias{featureNames,mRMRe.Network-method}
\alias{featureNames,mRMRe.Data-method}
\title{
	Accessor function for the 'featureNames' information in a mRMRe.Data, mRMRe.Filter and mRMRe.Network object
}
\description{
	featureNames are the names of the features given as input to the mRMRe procedure.
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
data(featureNames)
## maybe str(featureNames) ; plot(featureNames) ...
}
\keyword{datasets}
