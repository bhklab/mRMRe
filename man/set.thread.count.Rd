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
  \item{thread_count}{a single whole number >= 1: the number of OPENMP threads to request. Anything else (\code{NA}, zero, negative, non-integer or longer than one) is an error.}
}
\author{
Nicolas De Jay, Simon Papillon-Cavanagh, Benjamin Haibe-Kains
}
\value{
Invisibly, the thread count in effect: \code{thread_count} clamped to between 1 and the number of processors available to openMP. It is 0 (and nothing is changed) when the package was built without OpenMP support. An error is raised if \code{thread_count} is not a single whole number >= 1.
}
\examples{
# Access to number of available threads
threads <- get.thread.count()
# Force a single threaded openMP job
set.thread.count(1)

# Revert back to all accessible threads
if (threads >= 1) set.thread.count(threads)
}
\keyword{methods}