\name{InputMatrixFormatted}
\alias{InputMatrixFormatted}
\title{InputMatrixFormatted: make the input data meet the format requirements}
\description{
	This function is used to make the input PET clusters meet the format requirements of MICC model
}
\usage{
InputMatrixFormatted(data)
}
\arguments{
  \item{data}{
		PET clusters. The input PET clusters should have 9 columns as:
		chromosome; start; end; chromosome; start; end; PET-count between two anchor regions; total PET-count in left anchor region; total PET-count in right anchor region
	}
}
\value{
	\item{cAB }{PET-count between two anchor regions}
	\item{cA }{Total PET-count in left anchor region}
	\item{cB }{Total PET-count in right anchor region}
	\item{distance }{genomic distance between two anchor regions, unit: kb; for inter-chromosomal PET clusters, distance="inf"}
}
\author{
	Chao He
}
\seealso{
\code{\link{MICC_1.0-package}}
}
\examples{
library(MICC)

## Import data
data(TestData)

# format the data
data_formatted <- InputMatrixFormatted(TestData)

## The function is currently defined as
function (data) 
{
    x <- data
    cAB <- x[, 7]
    cA <- x[, 8]
    cB <- x[, 9]
    intra <- as.character(x[, 1]) == as.character(x[, 4])
    inter <- as.character(x[, 1]) != as.character(x[, 4])
    distance <- (x[, 5] + x[, 6] - x[, 3] - x[, 2])/2/1000
    distance[inter] <- Inf
    data_formatted <- list(cAB = cAB, cA = cA, cB = cB, distance = distance)
  }
}

\keyword{ format }