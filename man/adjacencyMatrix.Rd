\name{adjacencyMatrix}
\Rdversion{1.1}
\docType{methods}
\alias{adjacencyMatrix}
\alias{adjacencyMatrixSum}
\alias{adjacencyMatrix,mRMRe.Network-method}
\alias{adjacencyMatrixSum,mRMRe.Network-method}
\title{
	Accessor function for the 'adjacencyMatrix' information in a mRMRe.Network object.
	
}
\description{
	The adjency matrix is a directed matrix of 0's and 1's indicating if there is a link between features.
}
\usage{
\S4method{adjacencyMatrix}{mRMRe.Network}(object)
}
\arguments{
  \item{object}{a \code{mRMRe.Network} object.}
}
\value{
A square numeric matrix with feature names as row and column names. Entry \code{[i, j]} is 1 if feature \code{i} was selected for target \code{j} (parents in rows, targets in columns); for networks built without causality the matrix is also filled symmetrically.
}
\author{
Nicolas De Jay, Simon Papillon-Cavanagh, Benjamin Haibe-Kains
}
\examples{
set.thread.count(2)
data(cgps)
feature_data <- mRMR.data(data =  data.frame(cgps.ge))

# Build an mRMR-based network and display adjacency matrix (topology)
network <- new("mRMRe.Network", data = feature_data, target_indices = c(1, 2),
			levels = c(2, 1), layers = 1)
adjacencyMatrix(network)
}
\keyword{methods}
