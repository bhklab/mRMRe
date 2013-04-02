\name{subsetData}
\Rdversion{1.1}
\docType{methods}
\alias{subsetData}
\alias{subsetData,mRMRe.Data-method}

\title{
	Returns a mRMRe.Data object using a subset of the current mRMRe.Data object.
}
\description{
	This method is used to extract a subset of the current mRMRe.Data object.
}
\usage{
\S4method{sampleStrata}{mRMRe.Data}(object, row_indices, column_indices)
}
\arguments{
  \item{object}{a \code{mRMRe.Data} object.}
  \item{row_indices}{An integer vector of the rows to be included in the subset.}
  \item{column_indices}{An integer vector of the columns to be included in the subset.}
}
\author{
Nicolas De Jay, Simon Papillon-Cavanagh, Benjamin Haibe-Kains
}
\examples{
data(subsetData)
## maybe str(subsetData) ; plot(subsetData) ...
}
\keyword{datasets}
