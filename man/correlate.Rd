\name{correlate}
\alias{correlate}
\title{
Computes the degree of association between two variables
}
\description{
This function computes the degree of association between two variables (correlation)
}
\usage{
correlate(x, y, strata, weights, method, outX, bootstrap_count)
}
\arguments{
  \item{x, y}{
Two vector containing sample data.
}
  \item{strata}{
A vector containing the stratum associated with each sample. Defaults to all samples being in the same stratum.
}
  \item{weights}{
A vector containing the weight associated with each sample. Defaults to all samples having the same weight.
}
  \item{method}{
A string containing the estimator to be used for computing the correlation between x and y.
Must be one of the following methods: "cramer", "pearson", "spearman", "cindex".
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
%\details{}
\value{
A value depending on the estimator. "cramer" returns Cramer's V,
"pearson" returns Pearson's r, "spearman" returns Spearman's r, "cindex" returns the concordance index.
}
%\references{}
\author{
Nicolas De Jay, Simon Papillon-Cavanagh, Benjamin Haibe-Kains
}
%\note{}
\seealso{build.mim}
\examples{
data(cgps)
# Compute c-index between feature 1 and 2
correlate(cgps_ge[,1],cgps_ge[,2], method="cindex")

# Compute Cramer's V
x <- sample(c(0,1), 100, replace=TRUE)
y <- sample(c(0,1), 100, replace=TRUE)
correlate(x,y, method="cramer")

# Compute pearson coefficient with random strata and sample weights
# between feature 1 and 2
strata <- sample(c(0,1,2), nrow(cgps_ge), replace=TRUE)
weights <- runif(nrow(cgps_ge))
correlate(cgps_ge[,1],cgps_ge[,2], strata=strata, weights=weights, method="pearson", bootstrap_count=1000)
}
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line