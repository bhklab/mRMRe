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
\details{
}
\value{
A value depending on the estimator. "cramer" returns Cramer's V,
"pearson" returns Pearson's r, "spearman" returns Spearman's r, "cindex" returns the concordance index.
}
\references{
}
\author{
Nicolas De Jay, Simon Papillon-Cavanagh, Benjamin Haibe-Kains
}
\note{
}
\seealso{
}
\examples{
}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line