\name{correlate}
\alias{correlate}

\title{
Correlate
}
\description{
	Correlate is a function that computes the measure of association between two variables.
}
\usage{
correlate(X, Y, method = "pearson", strata, weights, outX = TRUE, bootstrap_count = 0, alpha = 0.05, alternative = c("two.sided", "less", "greater"))
}
\arguments{
  \item{X}{Vector of type {numeric, ordered factor, Surv}.}
  \item{Y}{Vector of type {numeric, ordered factor, Surv} of same length as \code{X}.}
  \item{method}{One of the following values: {pearson, spearman, kendall, frequency, cramersv, cindex}.}
  \item{strata}{Vector of type {ordered factor} corresponding to the sample strata.}
  \item{weights}{Vector of type {numeric} corresponding to the sample weights.}
  \item{outX}{If set to \code{TRUE}, ignore ties. Otherwise, take them into account.}
  \item{bootstrap_count}{If set to 0, combine strata by their relative weights. Otherwise, use the inverse variance weight mean of the estimates.}
  \item{alpha}{...}
  \item{alternative}{...}
}
\details{
	Since a mutual information matrix must be computed in order for feature selection to take place, the remaining arguments
	are very similar to those required by the \code{mim} method of the \code{mRMRe.Data} object.
}
%\value{
%
%}
\references{
}
\author{
}
\note{
}
\seealso{
\code{\link{mRMRe.Data-class}}
}
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
% How does the below work ??
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
