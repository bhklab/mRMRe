\name{set.thread.count}
\alias{set.thread.count}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
openMP Thread Count}
\description{
This methods allows you to set the number of cores currently accessible to openMP
}
\usage{
set.thread.count(thread_count)
}
\arguments{
  \item{thread_count}{number of OPENMP threads to be used}
}
\author{
Nicolas De Jay, Simon Papillon-Cavanagh, Benjamin Haibe-Kains
}
\value{
No return, sets the number of cores available to openMP for C level parallelization.
}
\examples{
# Access to number of available threads
threads <- get.thread.count()
# Force a single threaded openMP job
set.thread.count(1)

# Revert back to all accessible threads
set.thread.count(threads)
}
\keyword{methods}