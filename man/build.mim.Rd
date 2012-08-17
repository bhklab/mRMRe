\name{build.mim}
\alias{build.mim}
\title{Constructs a mutual information matrix.}
\description{%%  ~~ A concise (1-5 lines) description of what the function does. ~~}
\usage{build.mim(data, strata = rep.int(0, nrow(data)), weights = rep.int(1, nrow(data)), feature_types = rep.int(0, ncol(data)), uses_ranks = TRUE, outX = TRUE, bootstrap_count = 0)}
\arguments{
  \item{data}{A data frame with rows and columns respectively corresponding to samples and features.}
  \item{strata}{A vector containing the stratum associated with each sample. Optional argument.}
  \item{weights}{A vector containing the weight associated with each sample. Optional argument.}
  \item{feature_types}{A vector specifying the type of each feature. A value of 0 is for a continuous feature; 1, for a discrete feature; and 2, for survival data.}
  \item{uses_ranks}{If set to TRUE, the Spearman correlation is used between continuous features. If set to FALSE, the Pearson correlation is used.}
  \item{outX}{Boolean value on whether ties should be counted.}
  \item{bootstrap_count}{%%     ~~Describe \code{bootstrap_count} here~~}
}
\details{%%  ~~ If necessary, more details than the description above ~~}
\value{A symmetric matrix whose rows and columns correspond to features. Each (i, j) in the matrix is the mutual information between the features i and j.}
\references{%% ~put references to the literature/web site here ~}
\author{%%  ~~who you are~~}
\note{%%  ~~further notes~~}
\seealso{%% ~~objects to See Also as \code{\link{help}}, ~~~}
\examples{}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
