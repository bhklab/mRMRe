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
An integer: the number of cores accessible to openMP for C level parallelization. It is 0 when the package was built without OpenMP support (e.g. the macOS CRAN binaries), in which case the C code runs single-threaded.
}
\author{
Nicolas De Jay, Simon Papillon-Cavanagh, Benjamin Haibe-Kains
}
\examples{
get.thread.count() }
\keyword{ methods }