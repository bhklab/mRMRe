\name{get.thread.count}
\alias{get.thread.count}
\title{
	openMP Thread Count
}
\description{
This methods allows you to retrieve the number of cores currently accessible to openMP}
\usage{
get.thread.count()
}
\value{
Return the number of cores accessible to openMP for C level parallelization.
}
\author{
Nicolas De Jay, Simon Papillon-Cavanagh, Benjamin Haibe-Kains
}
\examples{
get.thread.count() }
\keyword{ methods }