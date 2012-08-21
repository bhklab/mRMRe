\name{cgps}
\alias{cgps_annot}
\alias{cgps_ge}
\alias{cgps_ic50}
\docType{data}
\title{Part of the large pharmacogenomic dataset published by Carnett et al. in the Cancer Genome Project (CGP)}
\description{This dataset contains gene expression of 200 cancer cell lines for which sensitivity (IC50) to Irinotecan was measured.}
\usage{data(cgps)}
\format{
The \code{\link[mRMRe]{cgps}} dataset is composed of three objects
 \describe{
    \item{cgps_annot}{Dataframe containing gene annotations}
    \item{cgps_ge}{Matrix containing expressions of 1000 genes; cell lines in rows, genes in columns}
    \item{cgps_ic50}{Drug sensitivity measurements (IC50) for Irinotecan}
}
}
\details{Irinotecan is a drug mainly used in colorectal cancer.}
\source{
  http://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-783
  http://www.nature.com/nature/journal/v483/n7391/extref/nature11005-s2.zip
}
\references{Garnett MJ et al. "Systematic identification of genomic markers of drug sensitivity in cancer cells", \emph{Nature}, \bold{483}:570-575, 2012.}
\examples{
  data(cgps)

  message("Gene expression data:")
  print(cgps_ge[1:3, 1:3])

  message("Gene annotations:")
  print(head(cgps_annot))

  message("Drug sensitivity (IC50) values:")
  print(head(cgps_ic50))
}
\keyword{datasets}
