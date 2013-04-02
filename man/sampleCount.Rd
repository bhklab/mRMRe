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
  \item{sampleCount}{a \code{mRMRe.Data} object.}
  \item{sampleCount}{a \code{mRMRe.Filter} object.}
  \item{sampleCount}{a \code{mRMRe.Network} object.}
}

\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(sampleCount)
## maybe str(sampleCount) ; plot(sampleCount) ...
}
\keyword{datasets}
