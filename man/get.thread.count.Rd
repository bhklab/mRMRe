\name{get.thread.count}
\alias{get.thread.count}
\title{
TODO
}
\description{
TODO}
\usage{
get.thread.count()
}

\author{
Nicolas De Jay, Simon Papillon-Cavanagh, Benjamin Haibe-Kains
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function () 
{
    thread_count <- vector(mode = "integer", length = 1)
    .Call(mRMRe:::.C_get_thread_count, thread_count)
    return(thread_count)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
