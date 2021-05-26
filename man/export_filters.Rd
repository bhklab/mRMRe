\name{export_filters}
\alias{export_filters}
\title{
Export filters
}
\description{
Export the filters
}
\arguments{
   \item{childrenCountPerLevel}{...}
   \item{dataMatrix}{...}
   \item{priorsMatrix}{...}
   \item{priorsWeight}{...}
   \item{sampleStrata}{...}
   \item{sampleWeights}{...}
   \item{featureTypes}{...}
   \item{sampleCount}{...}
   \item{featureCount}{...}
   \item{sampleStratumCount}{...}
   \item{targetFeatureIndices}{...}
   \item{continuousEstimator}{...}
   \item{outX}{...}
   \item{bootstrapCount}{...}
   \item{miMatrix}{...}
}
\value{
Exhaustively computes the minimum redudancy maximim relevance features from
the mutual information matrix, returing a list of solutions where each item
is a numeric index of selected features.
}
\author{
Nicolas De Jay, Simon Papillon-Cavanagh, Benjamin Haibe-Kains
}
\keyword{ methods }
