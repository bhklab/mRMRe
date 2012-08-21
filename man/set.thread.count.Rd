\name{set.thread.count}
\alias{set.thread.count}
\title{Sets the maximum number of threads that can be launched}
\description{This function sets the maximum number of threads that can be launched.}
\usage{set.thread.count(thread_count)}
\arguments{
  \item{thread_count}{Maximum number of threads.}
}
%\details{}
\value{None}
%\references{}
\author{Nicolas De Jay, Simon Papillon-Cavanagh, Benjamin Haibe-Kains}
%\note{}
%\seealso{}
\examples{
data(cgps)

set.thread.count(4)
mim <- build.mim(data.frame(target=cgps_ic50, cgps_ge))
head(mim)

}
%\keyword{ ~kwd1 }
