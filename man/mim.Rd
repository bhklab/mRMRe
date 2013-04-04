\name{mim}
\Rdversion{1.1}
\docType{methods}
\alias{mim}
\alias{mim,mRMRe.Data-method}
\alias{mim,mRMRe.Network-method}
\alias{mim,mRMRe.Filter-method}
\title{
	Accessor function for the 'mim' information in a mRMRe.Data, mRMRe.Filter and mRMRe.Network object
}
\description{
	In both mRMRe.Filter and mRMRe.Network objects, a sparse mutual information matrix is computed for the mRMRe procedure and this lazy-evaluated matrix is returned.
	In the context of a a mRMRe.Data 'mim', the full pairwise mutual information matrix is computed and returned.
}
\usage{
\S4method{mim}{mRMRe.Data}(object)
\S4method{mim}{mRMRe.Filter}(object)
\S4method{mim}{mRMRe.Network}(object)
}
\arguments{
  \item{object}{a \code{mRMRe.Data, mRMRe.Filter or mRMRe.Network} object.}
}
\author{
Nicolas De Jay, Simon Papillon-Cavanagh, Benjamin Haibe-Kains
}
\examples{
data(cgps)
feature_data <- mRMR.data(data =  data.frame(cgps.ge))

# Calculate the pairwise mutual information matrix
mim(feature_data)
filter <- mRMR.classic("mRMRe.Filter", data = feature_data, target_indices = 3:5, feature_count = 2)

# Obtain the sparse (lazy-evaluated) mutual information matrix.
mim(filter)
}
\keyword{methods}
