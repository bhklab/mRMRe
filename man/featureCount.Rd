\name{featureCount}
\Rdversion{1.1}
\docType{methods}
\alias{featureCount}
\alias{featureCount,mMRRe.Data-method}
\alias{featureCount,mRMRe.Filter-method}
\alias{featureCount,mRMRe.Network-method}
\title{
	Accessor function for the 'featureCount' information in a mRMRe.Data, mRMRe.Filter and mRMRe.Network object.
}
\description{
	The feature count is simply the total number of feature considered in the mRMRe procedure.
}
\usage{
\S4method{featureCount}{mRMRe.Data}(object)
\S4method{featureCount}{mRMRe.Filter}(object)
\S4method{featureCount}{mRMRe.Network}(object)
}
\arguments{
  \item{object}{a \code{mRMRe.Data} object.}
  \item{object}{a \code{mRMRe.Filter} object.}
  \item{object}{a \code{mRMRe.Network} object.}
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
data(featureCount)
## maybe str(featureCount) ; plot(featureCount) ...
}
\keyword{datasets}
