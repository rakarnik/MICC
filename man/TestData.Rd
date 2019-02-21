\name{TestData}
\alias{TestData}
\docType{data}
\title{
	TestData: A PET cluster file to test MICC model
}
\description{
	This data is a standard input PET cluster file for MICC model
}
\usage{data(TestData)}
\format{
  A data frame with 5613 observations on the following 9 variables.
  \describe{
    \item{\code{chr.}}{Chromosome label with levels \code{chr1} \code{chr10} \code{chr11} \code{chr12} \code{chr13} \code{chr14} \code{chr15} \code{chr16} \code{chr17} \code{chr18} \code{chr19} \code{chr2} \code{chr20} \code{chr21} \code{chr22} \code{chr3} \code{chr4} \code{chr5} \code{chr6} \code{chr7} \code{chr8} \code{chr9} \code{chrX} \code{chrY}}
    \item{\code{start}}{Starting position of left anchor regions }
    \item{\code{end}}{Ending position of left anchor regions}
    \item{\code{chr.}}{Chromosome label with levels \code{chr1} \code{chr10} \code{chr11} \code{chr12} \code{chr13} \code{chr14} \code{chr15} \code{chr16} \code{chr17} \code{chr18} \code{chr19} \code{chr2} \code{chr20} \code{chr21} \code{chr22} \code{chr3} \code{chr4} \code{chr5} \code{chr6} \code{chr7} \code{chr8} \code{chr9} \code{chrX} \code{chrY}}
    \item{\code{start}}{Starting position of right anchor regions}
    \item{\code{end}}{Ending position of right anchor regions}
    \item{\code{cAB}}{PET-count between two anchor regions}
    \item{\code{cA}}{Total PET-count in the left anchor region}
    \item{\code{cB}}{Total PET-count in the right anchor region}
  }
}

\examples{
data(TestData)
}
\keyword{datasets}
