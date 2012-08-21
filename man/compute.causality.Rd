\name{compute.causality}
\alias{compute.causality}
\title{Computes the causality}
\description{This function computes causality scores based on mutual information theory.}
\usage{compute.causality(data, target_index, mim, solutions, estimator)}
\arguments{
  \item{data}{
  	An mRMReObject
  	OR
	A data frame with rows and columns respectively corresponding to samples and features.
    Only columns of the following types are supported: "numeric", "ordered_factor" and "Surv".
    
    If mRMReObject is used, the other arguments need not be used.
  }
  \item{target_index}{Index of the target feature (column) from the data matrix.}
  \item{mim}{A mutual information matrix that will be used in computing causality.}
  \item{solutions}{A matrix where each row represents an mRMR solution.}
  \item{estimator}{String containing the estimator to be used for computing the correlation between elements. Must be one of the following: "pearson", "spearman" or "kendall".}
}
%\details{}
\value{Returns a matrix of causality coefficients for every pair of features considered in solutions.}
\references{
Bell AJ: The co-information lattice. In Proceedings of the Fifth International Workshop on Independent Component Analysis and Blind Signal Separation: ICA 2003. edited by Amari S et al.2003.
}
\author{Nicolas De Jay, Simon Papillon-Cavanagh, Benjamin Haibe-Kains}
%\note{}
\seealso{mRMR.classic, mRMR.ensemble}
\examples{
data(cgps)

data <- data.frame(target=cgps_ic50, cgps_ge)
ensemble <- mRMR.ensemble(data, 1, c(10, 5, rep(1, 5)))

compute.causality(ensemble)
compute.causality(data=data, target_index=1, mim=NULL, solutions=ensemble$paths, estimator="spearman")
}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
