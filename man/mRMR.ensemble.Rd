\name{mRMR.ensemble}
\alias{mRMR.ensemble}
\title{
Performs an ensemble mRMR feature selection
}
\description{
This function performs an ensemble mRMR feature selection, that is multiple mRMR selections in parallel.
}
\usage{
mRMR.ensemble(data, target_index, levels, strata, weights, uses_ranks, outX, bootstrap_count)
}
\arguments{
  \item{data}{
  A data frame with rows and columns respectively corresponding to samples and features.
   Only columns of the following types are supported: "numeric", "ordered_factor" and "Surv".
}
  \item{target_index}{
  Index of the target feature (column) from the data matrix.
}
  \item{levels}{
  A vector containing the number of features to be selected at each turn of ensemble mRMR feature selection.
}
  \item{strata}{
  A vector containing the stratum associated with each sample. Defaults to all samples being in the same stratum.
}
  \item{weights}{
  A vector containing the weight associated with each sample. Defaults to all samples having the same weight.
}
  \item{uses_ranks}{
  If set to TRUE, the Spearman correlation is used between continuous features. If set to FALSE, the Pearson correlation is used. Defaults to TRUE.
}
  \item{outX}{
  Boolean value on whether ties should be counted. Defaults to TRUE.
}
  \item{bootstrap_count}{
  If a positive non-zero integer is given,
  samples of each stratum will be recombined by the generic inverse-variance
  weighted average based on sample bootstrapping. Defaults to 0.
}
}
%\details{
%}
\value{
Returns a mRMReObject object with the following attributes:
  \item{paths}{Index matrix in which rows correspond to solutions and columns to the features selected for each solution.}
  \item{scores}{Score matrix associated with ensemble mRMR feature selection.}
  \item{mim}{Mutual information matrix used for ensemble mRMR feature selection.}
}
%\references{
%}
\author{
Nicolas De Jay, Simon Papillon-Cavanagh, Benjamin Haibe-Kains
}
%\note{
%}
\seealso{build.mim, mRMR.classic}
\examples{

data(cgps)

data <- data.frame(target=cgps_ic50, cgps_ge)
classic <- mRMR.ensemble(data, 1, rep.int(1, 30)) # For mRMR.classic-like results
ensemble <- mRMR.ensemble(data, 1, 10:1)
}
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
