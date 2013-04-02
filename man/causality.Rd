\name{causality}
\Rdversion{1.1}
\docType{methods}
\alias{causality}
\alias{causality,mRMRe.Filter-method}
\alias{causality,mRMRe.Network-method}
\title{
	Accessor function for the 'causality' information in a mRMRe.Filter and mRMRe.Network object.
}
\description{
	The causality data is compute using the co-information lattice algorithm on each V-structure (feature, target, feature).  Given that this procedure
	is computed for each pair of features, the minimum result is kept. A negative score indicates putative causality of the feature to the target.
}
\usage{
\S4method{causality}{mRMRe.Filter}(object)
\S4method{causality}{mRMRe.Network}(object)
}
\arguments{
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
data(causality)
## maybe str(causality) ; plot(causality) ...
}
\keyword{datasets}
