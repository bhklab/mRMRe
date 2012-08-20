\name{build.mim}
\alias{build.mim}
\title{Constructs a mutual information matrix.}
\description{}
\usage{build.mim(data, strata, weights, feature_types, uses_ranks, outX, bootstrap_count = 0)}
\arguments{
  \item{data}{A data frame with rows and columns respectively corresponding to samples and features.}
  \item{strata}{A vector containing the stratum associated with each sample. Defaults to all samples being in the same stratum.}
  \item{weights}{A vector containing the weight associated with each sample. Defaults to all samples having the same weight.}
  \item{uses_ranks}{If set to TRUE, the Spearman correlation is used between continuous features. If set to FALSE, the Pearson correlation is used. Defaults to TRUE.}
  \item{outX}{Boolean value on whether ties should be counted. Defaults to TRUE.}
  \item{bootstrap_count}{If a positive non-zero integer is given,
  samples of each stratum will be recombined by the generic inverse-variance
  weighted average based on sample bootstrapping. Defaults to 0.}
}
\details{}
\value{A symmetric matrix whose rows and columns correspond to features. Each (i, j) in the matrix is the mutual information between the features i and j.}
\references{}
\author{Nicolas De Jay, Simon Papillon-Cavanagh, Benjamin Haibe-Kains}
\note{}
\seealso{}
\examples{
data(cgps)

# The variables are all of continuous type
build.mim(data.frame(cgps_ic50, cgps_ge)) # Uses Spearman as correlation estimator
build.mim(data.frame(cgps_ic50, cgps_ge), uses_ranks=FALSE) # Uses Pearson as correlation estimator
}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
