\name{mRMRe.Network-class}
\Rdversion{1.1}
\docType{class}
\alias{mRMRe.Network-class}
\alias{adjacencyMatrix,mRMRe.Network-method}
\alias{adjacencyMatrixSum,mRMRe.Network-method}
\alias{causality,mRMRe.Network-method}
\alias{featureNames,mRMRe.Network-method}
\alias{mim,mRMRe.Network-method}
\alias{sampleNames,mRMRe.Network-method}
\alias{solutions,mRMRe.Network-method}
\alias{visualize,mRMRe.Network-method}

\title{Class \code{"mRMRe.Network"}}
\description{
%%  ~~ A concise (1-5 lines) description of what the class is. ~~
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("mRMRe.Network", data, prior_weight, target_indices, levels, layers, ..., mi_threshold, causality_threshold)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{topologies}:}{Object of class \code{"list"} ~~ }
    \item{\code{mi_matrix}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{causality_list}:}{Object of class \code{"list"} ~~ }
    \item{\code{sample_names}:}{Object of class \code{"character"} ~~ }
    \item{\code{feature_names}:}{Object of class \code{"character"} ~~ }
    \item{\code{target_indices}:}{Object of class \code{"integer"} ~~ }
  }
}
\section{Methods}{
  \describe{
    \item{adjacencyMatrix}{\code{signature(object = "mRMRe.Network")}: ... }
    \item{adjacencyMatrixSum}{\code{signature(object = "mRMRe.Network")}: ... }
    \item{causality}{\code{signature(object = "mRMRe.Network")}: ... }
    \item{featureNames}{\code{signature(object = "mRMRe.Network")}: ... }
    \item{mim}{\code{signature(object = "mRMRe.Network")}: ... }
    \item{sampleNames}{\code{signature(object = "mRMRe.Network")}: ... }
    \item{solutions}{\code{signature(object = "mRMRe.Network")}: ... }
    \item{visualize}{\code{signature(object = "mRMRe.Network")}: ... }
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
showClass("mRMRe.Network")
}
\keyword{classes}
