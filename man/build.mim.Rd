\name{build.mim}
\alias{build.mim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
build.mim(data, strata = rep.int(0, nrow(data)), weights = rep.int(1, nrow(data)), feature_types = rep.int(0, ncol(data)), uses_ranks = TRUE, outX = TRUE, bootstrap_count = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
%%     ~~Describe \code{data} here~~
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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (data, strata = rep.int(0, nrow(data)), weights = rep.int(1, 
    nrow(data)), feature_types = rep.int(0, ncol(data)), uses_ranks = TRUE, 
    outX = TRUE, bootstrap_count = 0) 
{
    data <- as.matrix(data)
    mi_matrix <- .Call(C_build_mim, as.vector(data), as.vector(strata), 
        as.vector(weights), as.vector(feature_types), as.integer(nrow(data)), 
        as.integer(ncol(data)), as.integer(length(unique(strata))), 
        as.integer(uses_ranks), as.integer(outX), as.integer(bootstrap_count))
    mi_matrix <- matrix(mi_matrix, nrow = ncol(data), ncol = ncol(data))
    rownames(mi_matrix) <- colnames(data)
    colnames(mi_matrix) <- colnames(data)
    return(mi_matrix)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
