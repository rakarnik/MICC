\name{MICCMainLearn}
\alias{MICCMainLearn}
\title{MICCMainLearn: main function to fit the model}
\description{
	This function is used to learn the model parameters from formatted PET clusters
}
\usage{
MICCMainLearn(data_formatted, params.init = NULL, reltol = 1e-05, abstol = 0.001, step = 200, restart = 5, MinConfident = 5)
}
\arguments{
  \item{data_formatted}{
		Formatted data matrix by "InputMatrixFormatted" function
	}
  \item{params.init}{
		Initialized paramters, see value section for more details
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

\value{
	\item{params }{ A list of model parameters
		\itemize{
			\item{theta12.par} {\cr 
				A vector of parameters for zeta distribution to describe true interaction PET clusters and random collision PET clusters. \cr 
				Let (par1,par2,par3,par3,theta0) = theta12.par, and distance denote genomic distance between two anchor regions of a PET cluster, then \cr 
				zeta paramter for true interaction PET clusters is \cr 
					theta1 = ( par1*distance + par3*par2*1000 ) / ( par2*1000 + distance ) + par4 / (10*distance) \cr 
				zeta paramter for true interaction PET clusters is \cr 
					theta2 = ( par1*distance + par3*par2*1000 ) / ( par2*1000 + distance ) + par4 / (10*distance) + theta0
			}
			\item{lambda.par} { \cr 
				A vector of parameters for prior probability to describe random ligation PET clusters. \cr 
				Let (par1,par2,lambda0) = lambda.par, and distance denote genomic distance between two anchor regions of a PET cluster, then \cr 
				lambda <- lambda0 * exp( par1*log(distance) + par2 ) / ( 1 + log(distance) )
			}
			\item{mu.par} { \cr 
				A vector of parameters for prior probability to describe random collision PET clusters. \cr 
				Let (par1,par2) = mu.par, and (cA, cB) denote total PET-count in the two anchor regions, then \cr 
				mu <- exp( par1*log(cA) + par2 ) * exp( par1*log(cB) + par2 ) / ((1+exp( par1*log(cA) + par2 ) )*(1+exp( par1*log(cB) + par2 )))
			}
		}
	}	
	\item{LogLik }{Logliklihood after fitting the model}
	\item{PostProb }{A matrix with 3 columns to describe the posterior probability of each PET clusters for each state, respectively}
	\item{CONVERGED!}{ The training is absolutely converged! }
	\item{RELCONVERGED!}{ The training is relatively converged! }
	\item{NOT CONVERGED!}{ The training is neither absolutely nor relatively converged! }

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

## The function is currently defined as
function (data_formatted, params.init = NULL, reltol = 1e-05, 
    abstol = 0.001, step = 200, restart = 5, MinConfident = 5) 
{
    Par <- EMIter(data_formatted, params.init = params.init, 
        reltol = reltol, abstol = abstol, step = step, restart = restart, 
        MinConfident = MinConfident)
    Par
  }
}

\keyword{ method }
