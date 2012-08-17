\name{mRMR.classic}
\alias{mRMR.classic}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Performs an mRMR feature selection.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
mRMR.classic(data, target, feature_count, strata = rep.int(0, nrow(data)), weights = rep.int(1, nrow(data)), feature_types = rep.int(0, ncol(data)), uses_ranks = TRUE, outX = TRUE, bootstrap_count = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
%%     ~~Describe \code{data} here~~
}
  \item{target}{
%%     ~~Describe \code{target} here~~
}
  \item{feature_count}{
%%     ~~Describe \code{feature_count} here~~
}
  \item{strata}{
%%     ~~Describe \code{strata} here~~
}
  \item{weights}{
%%     ~~Describe \code{weights} here~~
}
  \item{feature_types}{
%%     ~~Describe \code{feature_types} here~~
}
  \item{uses_ranks}{
%%     ~~Describe \code{uses_ranks} here~~
}
  \item{outX}{
%%     ~~Describe \code{outX} here~~
}
  \item{bootstrap_count}{
%%     ~~Describe \code{bootstrap_count} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
