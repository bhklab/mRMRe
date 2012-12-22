\name{mRMRe.Filter-class}
\Rdversion{1.1}
\docType{class}
\alias{mRMRe.Filter-class}
\alias{causality,mRMRe.Filter-method}
\alias{featureNames,mRMRe.Filter-method}
\alias{mim,mRMRe.Filter-method}
\alias{sampleNames,mRMRe.Filter-method}
\alias{solutions,mRMRe.Filter-method}
\alias{target,mRMRe.Filter-method}

\title{Class \code{"mRMRe.Filter"}}
\description{
%%  ~~ A concise (1-5 lines) description of what the class is. ~~
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("mRMRe.Filter", data, prior_weight, target_indices, levels, continuous_estimator, outX, bootstrap_count)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{filters}:}{Object of class \code{"list"} ~~ }
    \item{\code{mi_matrix}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{causality_list}:}{Object of class \code{"list"} ~~ }
    \item{\code{sample_names}:}{Object of class \code{"character"} ~~ }
    \item{\code{feature_names}:}{Object of class \code{"character"} ~~ }
    \item{\code{target_indices}:}{Object of class \code{"integer"} ~~ }
    \item{\code{levels}:}{Object of class \code{"integer"} ~~ }
  }
}
\section{Methods}{
  \describe{
    \item{causality}{\code{signature(object = "mRMRe.Filter")}: ... }
    \item{featureNames}{\code{signature(object = "mRMRe.Filter")}: ... }
    \item{mim}{\code{signature(object = "mRMRe.Filter")}: ... }
    \item{sampleNames}{\code{signature(object = "mRMRe.Filter")}: ... }
    \item{solutions}{\code{signature(object = "mRMRe.Filter")}: ... }
    \item{target}{\code{signature(object = "mRMRe.Filter")}: ... }
	 }
}
\references{
%%  ~~put references to the literature/web site here~~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%%  ~~objects to See Also as \code{\link{~~fun~~}}, ~~~
%%  ~~or \code{\linkS4class{CLASSNAME}} for links to other classes ~~~
}
\examples{
showClass("mRMRe.Filter")
}
\keyword{classes}
