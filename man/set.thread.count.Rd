\name{set.thread.count}
\alias{set.thread.count}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
TODO}
\description{
TODO}
\usage{
set.thread.count(thread_count)
}
\arguments{
  \item{thread_count}{number of OPENMP threads to be used}
}
\author{
Nicolas De Jay, Simon Papillon-Cavanagh, Benjamin Haibe-Kains
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (thread_count) 
{
    thread_count <- as.integer(thread_count)
    .Call(mRMRe:::.C_set_thread_count, thread_count)
    return(thread_count)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
