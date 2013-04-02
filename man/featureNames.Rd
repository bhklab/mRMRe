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
data(featureNames)
## maybe str(featureNames) ; plot(featureNames) ...
}
\keyword{datasets}
