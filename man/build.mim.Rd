\name{build.mim}
\alias{build.mim}
\title{Constructs a mutual information matrix (MIM)}
\description{This function constructs a mutual information matrix (MIM).}
\usage{build.mim(data, priors, prior_weight, strata, weights, uses_ranks, outX, bootstrap_count, as_mi)}
\arguments{
  \item{data}{A data frame with rows and columns respectively corresponding to samples and features.
  Only columns of the following types are supported: "numeric", "ordered_factor" and "Surv".}
  \item{priors}{A non-symmetric matrix of prior causal relations between features where values can be [0, 1]; 
  0 meaning i is not causal of j and 1 indicates that i is a cause of j. Defaults to having no priors.}
  \item{prior_weight}{A value [0, 1] representing the weight the importance of priors for biasing the estimation of mutual information. Defaults to 0.}
  \item{strata}{A vector of factors specifying the stratum associated with each sample. Defaults to all samples being in the same stratum.}
  \item{weights}{A vector containing the weight associated with each sample. Defaults to all samples having the same weight.}
  \item{uses_ranks}{If set to TRUE, the Spearman correlation is used between numeric features. If set to FALSE, the Pearson correlation is used. Defaults to TRUE.}
  \item{outX}{Boolean value on whether ties should be counted. Defaults to TRUE.}
  \item{bootstrap_count}{If a positive non-zero integer is given,
  samples of each stratum will be recombined by the generic inverse-variance
  weighted average based on sample bootstrapping. Defaults to 0.}
  \item{as_mi}{If set to TRUE, the function returns a MI matrix. If set to FALSE, the function returns a correlation matrix. Defaults to TRUE.}
}
%\details{}
\value{A symmetric matrix whose rows and columns correspond to features. Each (i, j) in the matrix is the mutual information between the features i and j.}
%\references{}
\author{Nicolas De Jay, Simon Papillon-Cavanagh, Benjamin Haibe-Kains}
%\note{}
\seealso{mRMR.classic, mRMR.ensemble}
\examples{
## Test on a dummy dataset
set.seed(0)
library(survival)
dd <- data.frame("surv1"=Surv(runif(100), sample(0:1, 100, replace=TRUE)),
	"cont1"=runif(100),
	"cat1"=factor(sample(1:5, 100, replace=TRUE), ordered=TRUE),
	"surv2"=Surv(runif(100), sample(0:1, 100, replace=TRUE)),
	"cont2"=runif(100),
	"cont3"=runif(100),
	"surv3"=Surv(runif(100),
	sample(0:1, 100, replace=TRUE)),
	"cat2"=factor(sample(1:5, 100, replace=TRUE), ordered=TRUE))
	library(mRMRe)
	data <- new("mRMRe.Data", data=dd)
	filter <- new("mRMRe.Filter", data=data, target_index=1, levels=c(4))
	#mRMR.network(data = dd, target_indices = c(1, 2), solution_count = 2, feature_count = 2)
		#mRMR.ensemble(data = dd, target_index = 1, solution_count = 2, feature_count = 2)
	
print(head(dd))
build.mim(data=dd)

## Test on the 'cgps' dataset
data(cgps)
## The variables are all of continuous type
mim <- build.mim(data.frame(cgps_ic50, cgps_ge)) # Uses Spearman as correlation estimator
print(head(mim))
mim <- build.mim(data.frame(cgps_ic50, cgps_ge), uses_ranks=FALSE) # Uses Pearson as correlation estimator
print(head(mim))
}
%\keyword{mim}
