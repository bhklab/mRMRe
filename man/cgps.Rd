\name{cgps}
\alias{cgps.annot}
\alias{cgps.ge}
\alias{cgps.ic50}
\docType{data}
\title{Part of the large pharmacogenomic dataset published by Garnett et al. within the Cancer Genome Project (CGP)}
\description{This dataset contains gene expression of 200 cancer cell lines for which sensitivity (IC50) to Camptothecin was measured (release 2).}
\usage{data(cgps)}
\format{
The \code{\link[mRMRe]{cgps}} dataset is composed of three objects
 \describe{
    \item{cgps.annot}{Dataframe containing gene annotations}
    \item{cgps.ge}{Matrix containing expressions of 1000 genes; cell lines in rows, genes in columns}
    \item{cgps.ic50}{Drug sensitivity measurements (IC50) for Camptothecin}
}
}
\details{Camptothecin is a drug mainly used in colorectal cancer.}
\source{
  http://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-783
    
  http://www.nature.com/nature/journal/v483/n7391/extref/nature11005-s2.zip
}
\references{Garnett MJ et al. "Systematic identification of genomic markers of drug sensitivity in cancer cells", \emph{Nature}, \bold{483}:570-575, 2012.}
\examples{
set.thread.count(2)
data(cgps)

message("Gene expression data:")
print(cgps.ge[1:3, 1:3])

message("Gene annotations:")
print(head(cgps.annot))

message("Drug sensitivity (IC50) values:")
print(head(cgps.ic50))
}
\keyword{datasets}
