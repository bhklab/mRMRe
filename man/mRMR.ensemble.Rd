\name{mRMR.ensemble}
\alias{mRMR.ensemble}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
mRMR.ensemble(data, target, levels, strata = rep.int(0, nrow(data)), weights = rep.int(1, nrow(data)), feature_types = rep.int(0, ncol(data)), uses_ranks = TRUE, outX = TRUE, bootstrap_count = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
%%     ~~Describe \code{data} here~~
}
  \item{target}{
%%     ~~Describe \code{target} here~~
}
  \item{levels}{
%%     ~~Describe \code{levels} here~~
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
function (data, target, levels, strata = rep.int(0, nrow(data)), 
    weights = rep.int(1, nrow(data)), feature_types = rep.int(0, 
        ncol(data)), uses_ranks = TRUE, outX = TRUE, bootstrap_count = 0) 
{
    if (is.data.frame(data)) 
        stop("data must be of type data frame")
    levels <- as.vector(levels)
    data <- cbind(target, data)
    tree <- .Call(C_build_mRMR_tree_from_data, levels, as.vector(data), 
        as.vector(strata), as.vector(weights), as.vector(feature_types), 
        nrow(data), ncol(data), as.integer(length(unique(strata))), 
        0, as.integer(uses_ranks), as.integer(outX), as.integer(bootstrap_count))
    tree$mim <- matrix(tree$mim, ncol = sqrt(length(tree$mim)), 
        nrow = sqrt(length(tree$mim)))
    rownames(tree$mim) <- colnames(data)
    colnames(tree$mim) <- colnames(data)
    wrap <- function(i) t(matrix(i[length(i):1], nrow = length(levels), 
        ncol = length(i)/length(levels)))
    object <- list(paths = wrap(tree$paths) + 1, scores = wrap(tree$scores), 
        tree$mim)
    class(object) <- "mRMReObject"
    return(object)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
