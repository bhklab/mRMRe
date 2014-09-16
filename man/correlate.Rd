\name{correlate}
\alias{correlate}

\title{
Function to compute various correlation measures between two variables
}
\description{
	Correlate is a function that cestimates correlation between two variables, which can be either continuous, categorical (ordered factor) or censored (survival data).
}
\usage{
correlate(X, Y, method = c("pearson", "spearman", "kendall", "frequency", "cramersv",
		"cindex"),strata, weights, outX = TRUE, bootstrap_count = 0, alpha = 0.05,
 		alternative = c("two.sided", "less", "greater"))
}
\arguments{
  \item{X}{Vector of type {numeric, ordered factor, or Surv}.}
  \item{Y}{Vector of type {numeric, ordered factor, or Surv} of same length as \code{X}.}
  \item{method}{One of the following values: {pearson, spearman, kendall, frequency, cramersv, or cindex}.}
  \item{strata}{Vector of type {factor} corresponding to the sample strata.}
  \item{weights}{Vector of type {numeric} corresponding to the sample weights.}
  \item{outX}{For cindex, if set to \code{TRUE}, ignore ties; otherwise, take them into account when computing the concordance index.}
  \item{bootstrap_count}{If set to \code{0}, analytical standard error for the correlation estimate in each strata is used to compute the meta-estimate (inverse-variance weighting avarega); otherwise a number of bootstraps are used to computes standard errors.}
  \item{alpha}{The probability of Type I error that is, rejecting a null hypothesis when it is in fact true}
  \item{alternative}{a character string specifying the alternative hypothesis,must be one of {two.sided (default), greater or less}.  You can specify just the initial letter.}
}
\details{
   The \code{correlate} function could be used to measure correlation between any types of variables:
   \describe{
   \item{numeric vs. numeric}{Pearson, Spearman, Kendall or concordance index}
   \item{numeric vs. ordered factor}{concordance index (Somers' Dxy)}
   \item{numeric vs. survival data}{concordance index (Somers' Dxy)}
   \item{ordered factor vs. ordered factor}{Carmer's V}
   \item{ordered factor vs. survival data}{concordance index (Somers' Dxy)}
   \item{survival data vs. survival data}{concordance index (Somers' Dxy)}
   }
 
	Part of the code underlying \code{correlate} is also used in \code{mim} method of the \code{mRMRe.Data} object because correlations are used to build the mutual information matrix in order for feature selection to take place. This is why these two functions have many argiuments in common.
}
\value{
   \item{estimate}{point estimate}
   \item{se}{standard error}
   \item{lower}{lower confidence bound}
   \item{upper}{upper confidence bound}
   \item{p}{p-value}
   \item{n}{sample size}
}

\author{
Nicolas De Jay, Simon Papillon-Cavanagh, Benjamin Haibe-Kains
}

\seealso{
\code{\link[mRMRe]{mRMRe.Data-class}}
}
\examples{
set.thread.count(2)

## load data
data(cgps)
   
## spearman correlation coefficent between the first gene and Camptothecin IC50
correlate(X=cgps.ge[ ,1], Y=cgps.ic50, method="spearman")
   
## concordance index between the first gene and Camptothecin IC50
correlate(X=cgps.ge[ ,1], Y=cgps.ic50, method="cindex")
}

\keyword{ univar }
