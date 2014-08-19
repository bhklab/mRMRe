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
\S4method{subsetData}{mRMRe.Data}(object, row_indices, column_indices)
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
set.thread.count(2)
data(cgps)
feature_data <- mRMR.data(data = data.frame(cgps.ge))
# Subset the same dimensions, equivalent to making a copy
feature_data_copy <- subsetData(feature_data, 
    row_indices=sampleCount(feature_data),
		column_indices=featureCount(feature_data))

# Use only half of the samples
feature_data_samples <- subsetData(feature_data, row_indices=sampleCount(feature_data)/2)

# Use only half of the features
feature_data_features <- subsetData(feature_data,
    column_indices=featureCount(feature_data))
}
\keyword{methods}
