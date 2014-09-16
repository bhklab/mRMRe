\name{featureData}
\Rdversion{1.1}
\docType{methods}
\alias{featureData}
\alias{featureData,mRMRe.Data-method}
\title{
	Accessor function for the 'featureData' information in a mRMRe.Data object
}
\description{
	the featureData consists of the numerical value of each feature for each sample considered
}
\usage{
\S4method{featureData}{mRMRe.Data}(object)
}
\arguments{
  \item{object}{a \code{mRMRe.Data} object.}
}
\author{
Nicolas De Jay, Simon Papillon-Cavanagh, Benjamin Haibe-Kains
}
\examples{
set.thread.count(2)
data(cgps)
feature_data <- mRMR.data(data =  data.frame(cgps.ge))
featureData(feature_data)
}
\keyword{methods}
