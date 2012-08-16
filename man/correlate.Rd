\name{correlate}
\alias{correlate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Computes the degree of association between two variables.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
correlate(x, y, strata = rep.int(0, length(x)), weights = rep.int(1, length(x)), method = c("cramer", "pearson", "spearman", "cindex"), outX = TRUE, bootstrap_count = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{y}{
%%     ~~Describe \code{y} here~~
}
  \item{strata}{
%%     ~~Describe \code{strata} here~~
}
  \item{weights}{
%%     ~~Describe \code{weights} here~~
}
  \item{method}{
%%     ~~Describe \code{method} here~~
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
function (x, y, strata = rep.int(0, length(x)), weights = rep.int(1, 
    length(x)), method = c("cramer", "pearson", "spearman", "cindex"), 
    outX = TRUE, bootstrap_count = 0) 
{
    x <- as.vector(x)
    y <- as.vector(y)
    strata <- as.vector(strata)
    weights <- as.vector(weights)
    method <- match.arg(method)
    stratum_count <- as.integer(length(unique(strata)))
    bootstrap_count <- as.integer(bootstrap_count)
    value <- NA
    if (method == "cramer") 
        value <- .Call(C_compute_cramers_v, x, y, weights, strata, 
            stratum_count, bootstrap_count)
    else if (method == "pearson") 
        value <- .Call(C_compute_pearson_correlation, x, y, weights, 
            strata, stratum_count, bootstrap_count)
    else if (method == "spearman") 
        value <- .Call(C_compute_spearman_correlation, x, y, 
            weights, strata, stratum_count, bootstrap_count)
    else if (method == "cindex") 
        value <- .Call(C_compute_concordance_index, x, y, weights, 
            strata, stratum_count, outX)
    return(value)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
