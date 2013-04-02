\name{priors}
\Rdversion{1.1}
\docType{methods}
\alias{priors}
\alias{priors,mMRRe.Data-method}
\alias{priors<-,mMRRe.Data,numeric-method}

\title{
	Accessor function for the 'priors' information in a mRMRe.Data object
}
\description{
	The priors matrix consists of a prior bias to be used in computation to mutual information between features.
}
\usage{
\S4method{priors}{mRMRe.Data}(object)
\S4method{priors}{mRMRe.Data, methods}(object) <- value
}
\arguments{
  \item{object}{a \code{mRMRe.Data} object.}
  \item{value}{a numeric matrix containing values from 0 to 1 (or NA), one per pairwise feature bias.}
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
data(priors)
## maybe str(priors) ; plot(priors) ...
}
\keyword{datasets}
