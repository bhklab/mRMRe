\name{compute.causality}
\alias{compute.causality}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Computes the causality
}
\description{
This function computes causality scores based on mutual infomration theory.
}
\usage{
compute.causality(data, target_index, mim, solutions, estimator = c("pearson", "spearman", "kendall"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
  	An mRMReObject
  	OR
	A data frame with rows and columns respectively corresponding to samples and features.
    Only columns of the following types are supported: "numeric", "ordered_factor" and "Surv".
    
    If mRMReObject is used, the other arguments need not to be used.
}
  \item{target_index}{
	Index of the target feature (column) from the data matrix.
}
  \item{mim}{
	(optional) A mutual information matrix that will be used in computing causality.
}
  \item{solutions}{
	A matrix where each row represents an mRMR solution
}
  \item{ estimator=c("pearson", "spearman", "kendall")}{
	A string containing the estimator to be used for computing the correlation between elements.
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
Returns a matrix of causality coefficients for every pair of features considered in solutions.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Nicolas De Jay, Simon Papillon-Cavanagh, Benjamin Haibe-Kains
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{mRMR.classic, mRMR.ensemble}
\examples{
data(cgps)

data <- data.frame(target=as.numeric(cgps_ic50), as.numeric(cgps_ge))
ensemble <- mRMR.ensemble(data, 1, 10:1)
compute.causality(ensemble)
compute.causality(data=data, target_index=1, mim=NULL, solutions=ensemble$paths, estimator="spearman")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line