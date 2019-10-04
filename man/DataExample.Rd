\name{DataExample}
\alias{A}
\alias{b}
\alias{DataExample}
\docType{data}
\title{
Simulated data for linear, logistic and Poisson Regression
}
\description{Simple simulated data, used to demonstrate the feature of glmADMM
}
\usage{
data(GaussianExample)
data(BinomialExample)
data(PoissonExample)
}
\format{
  Data objects used to demonstrate the feature of glmADMM
}

\examples{
data(GaussianExample)
glmADMM(A,b,family="gaussian")
}
\keyword{datasets}
