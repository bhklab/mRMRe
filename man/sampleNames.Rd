\name{sampleNames}
\Rdversion{1.1}
\docType{methods}
\alias{sampleNames}
\alias{sampleNames,mRMRe.Filter-method}
\alias{sampleNames,mRMRe.Network-method}
\alias{sampleNames,mRMRe.Data-method}
\title{
	Accessor function for the 'sampleNames' information in a mRMRe.Data, mRMRe.Filter and mRMRe.Network object.
}
\description{
	sampleNames are the names of the samples given as input to the mRMRe procedure.
}
\usage{
\S4method{sampleNames}{mRMRe.Data}(object)
\S4method{sampleNames}{mRMRe.Filter}(object)
\S4method{sampleNames}{mRMRe.Network}(object)
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
data(sampleNames)
## maybe str(sampleNames) ; plot(sampleNames) ...
}
\keyword{datasets}
