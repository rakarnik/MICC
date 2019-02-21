\name{MICC_1.0-package}
\alias{MICC_1.0-package}
\alias{MICC_1.0}
\docType{package}
\title{
MICC: an R package to identify chromatin interactions from ChIA-PET data
}
\description{
This package combined Bayesian prior and a three-component mixture model to systematically remove random ligation and random collison noise. It could detect chromatin interactions at a higher sensitivity from ChIA-PET sequencing data
}
\details{
\tabular{ll}{
Package: \tab MICC_1.0\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2015-03-18\cr
License: LGPL (>=2)
}
InputMatrixFormatted: process input PET clusters into MICC main model input\cr
MICCMainLearn: fit model parameters\cr
FDRcompute: compute FDR\cr
MICCfoutput: implement MICC model and get output\cr
}
\author{
Chao He\cr
Maintainer: hechaobnu@gmail.com
}
\references{
Dempster, P. et al. (1977). Maximum likelihood from incomplete data via the EM algorithm. Journal of the Royal statistical Society, 39, 1¨C38.\cr
Jessen, B. and Winter, A. (1935).Distribution functions and the Riemann zeta function. Transactions of the American Mathematical Society, 38, 48-88.\cr
Li, G. et al. (2010). ChIA-PET tool for comprehensive chromatin interaction analysis with paired-end tag sequencing. Genome Biology, 11, 1¨C13.\cr
Newman, E. (2005). Power laws, pareto distributions and zipf¡¯s law. Contemporary physics, 46, 323¨C351.\cr
}
\keyword{ MICC }
\seealso{
\code{\link{InputMatrixFormatted}},
\code{\link{MICCMainLearn}},
\code{\link{MICCoutput}}.
}
\examples{
library(MICC)

## Import data
data(TestData)

# implement the model
MICCoutput( TestData, "./TestData.txt" )

}
