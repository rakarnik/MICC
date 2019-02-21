\name{FDRcompute}
\alias{FDRcompute}
\title{
	FDRcompute: compute false discovery date by comparing posterior probability from real data and simulated data
}
\description{
	This function is used to estimate false discovery date of MICC model by comparing posterior probability from real data and simulated data
}

\usage{
FDRcompute(data_formatted, params, PostProb)
}

\arguments{
  \item{data_formatted}{
		Formatted data matrix by "InputMatrixFormatted" function
	}
  \item{params}{
		Parameters after fitting MICC model. See "MICCMainLearn" for more details.
	}
  \item{PostProb}{
		Posterior probability that a PET cluster is true interaction clusters
	}
}

\value{
	A vector of estimated false discovery date.
}

\author{Chao He}
\seealso{
\code{\link{InputMatrixFormatted}},
\code{\link{MICC_1.0-package}}.
}
\examples{
library(MICC)

## Import data
data(TestData)

# format the data
data_formatted <- InputMatrixFormatted(TestData)

# train the model
Par <- MICCMainLearn( data_formatted, reltol=1e-5, step=200 )

# compute FDR
params <- Par$params
PostProb <- Par$PostProb
fdr <- FDRcompute( data_formatted, params, PostProb[,1] )

## The function is currently defined as
function (data_formatted, params, PostProb) 
{
    fdr <- FDRestimate(data, params, PostProb)
    fdr
  }
}
\keyword{ method }
