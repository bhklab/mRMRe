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
data(sampleCount)
## maybe str(sampleCount) ; plot(sampleCount) ...
}
\keyword{datasets}
