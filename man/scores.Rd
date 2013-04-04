\name{scores}
\Rdversion{1.1}
\docType{methods}
\alias{scores}
\alias{scores,mRMRe.Filter-method}
\title{
	mRMR Scores as per the MI gain for each feature}
\description{
	The 'scores' method returns the scores of individual features in respect to previously selected features as per standard
	mRMR procedure.  For each target, the score of a feature is defined as the mutual information between the target and this feature
	minus the average mutual information of previously selected features and this feature.
}
\usage{
\S4method{scores}{mRMRe.Filter}(object)
}
\arguments{
  \item{object}{a \code{mRMRe.Filter} object.}
}
\author{
Nicolas De Jay, Simon Papillon-Cavanagh, Benjamin Haibe-Kains
}
\examples{
data(cgps)
feature_data <- mRMR.data(data =  data.frame(cgps.ge))

# Create an mRMR filter and obtain the indices of selected features
filter <- mRMR.classic("mRMRe.Filter", data = feature_data, target_indices = 3:5,
						feature_count = 2)
scores(filter)
}
\keyword{methods}
