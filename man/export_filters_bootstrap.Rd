\name{export_filters_bootstrap}
\alias{export_filters_bootstrap}
\title{
Export filters bootstrap
}
\description{
Export the filters
}
\arguments{
   \item{solutionCount}{...}
   \item{solutionLength}{...}
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
Bootstraps and estimate of the minimum redundancy maximim relevance features
from the mutual information, returning a list where each item is a numeric
vector of selected feature indices.
}
\author{
Nicolas De Jay, Simon Papillon-Cavanagh, Benjamin Haibe-Kains
}
\keyword{ methods }
