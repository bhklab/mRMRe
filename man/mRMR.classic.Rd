\name{mRMR.classic}
\alias{mRMR.classic}
\title{Performs an mRMR feature selection}
\description{This function performs an mRMR feature selection.}
\usage{mRMR.classic(data, target_index, feature_count, strata, weights, uses_ranks, outX, bootstrap_count)}
\arguments{
  \item{data}{A data frame with rows and columns respectively corresponding to samples and features. Only columns of the following types are supported: "numeric", "ordered_factor" and "Surv".}
  \item{target_index}{Index of the target feature (column) from the data matrix.}
  \item{feature_count}{The number of features to be selected by mRMR feature selection.}
  \item{strata}{A vector containing the stratum associated with each sample. Defaults to all samples being in the same stratum.}
  \item{weights}{A vector containing the weight associated with each sample. Defaults to all samples having the same weight.}
  \item{uses_ranks}{If set to TRUE, the Spearman correlation is used between numeric features. If set to FALSE, the Pearson correlation is used. Defaults to TRUE.}
  \item{outX}{Boolean value on whether ties should be counted. Defaults to TRUE.}
  \item{bootstrap_count}{If a positive non-zero integer is given, samples of each stratum will be recombined by the generic inverse-variance weighted average based on sample bootstrapping. Defaults to 0.}
}
\details{
The resulting mutual information matrix is computed using Pearson or Spearman's rho between numeric features (depending on the uses_ranks argument);
using Cramer's V between ordered factors; and Concordance index for any other combination.
}
\value{
Returns a mRMReObject object with the following attributes:
  \item{paths}{Index vector in which elements correspond to the features selected.}
  \item{scores}{Score vector associated with mRMR feature selection.}
  \item{mim}{Mutual information matrix used for mRMR feature selection.}
}
%\references{}
\author{Nicolas De Jay, Simon Papillon-Cavanagh, Benjamin Haibe-Kains}
%\note{}
\seealso{mRMR.classic, mRMR.ensemble}
\examples{
data(cgps)

data <- data.frame(target=cgps_ic50, cgps_ge)
mRMR.classic(data, 1, 30)
}
%\keyword{ ~kwd1 }
