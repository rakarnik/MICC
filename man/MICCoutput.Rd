\name{MICCoutput}
\alias{MICCoutput}
\title{
	MICCoutput: main function to implement MICC model and get output
}
\description{
	This function is used to implement MICC model and get the output files
}
\usage{
MICCoutput(data, outfilename, params.init = NULL, reltol = 1e-05, abstol = 0.001, step = 200, restart = 5, MinConfident = 5)
}

\arguments{
  \item{data}{
		Input matrix of PET clusters
	}
  \item{outfilename}{
		Output filename
	}
  \item{params.init}{
		Initialized paramters, see "MICCMainLearn" function for more details
	}
  \item{reltol}{
		Relative tolerance, default value: 1e-5 
	}
  \item{abstol}{
		Absolute tolerance, default value: 1e-5 
	}
  \item{step}{
		Max number of steps before convergence, default value: 200
	}
  \item{restart}{
		Times to restart before convergence, default value: 5
	}
  \item{MinConfident}{
		Minimal number of PET-count to classify true interaction PET clusters when initializing the paramters, default value: 5 
	}
}

\author{Chao He}
\seealso{
\code{\link{InputMatrixFormatted}},
\code{\link{MICCMainLearn}},
\code{\link{MICC_1.0-package}}.
}
\examples{
library(MICC)

## Import data
data(TestData)

# implement the model
MICCoutput( TestData, "./TestData.txt" )

## The function is currently defined as
function (data, outfilename, params.init = NULL, reltol = 1e-05, 
    abstol = 0.001, step = 200, restart = 5, MinConfident = 5) 
{
    data_formatted <- InputMatrixFormatted(data)
    Par <- MICCMainLearn(data_formatted, params.init = params.init, 
        reltol = reltol, abstol = abstol, step = step, restart = restart, 
        MinConfident = MinConfident)
    params <- Par$params
    PostProb <- Par$PostProb
    fdr <- FDRcompute(data, params, PostProb[, 1])
    output.colnames <- c("chr.", "start", "end", "chr.", "start", 
        "end", "cAB", "cA", "cB", "-log10(1-PostProb)", "fdr")
    y <- cbind(x, -log10(1 - PostProb[, 1]), fdr)
    colnames(y) <- output.colnames
    write.table(y, file = outfilename, sep = "\t", row.names = F, 
        quote = F)
  }
}

\keyword{ model }

