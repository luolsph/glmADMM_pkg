\name{plotglmADMM}
\alias{plotglmADMM}
\title{
plot the solution path from a "glmADMM" object
}
\description{
This function plot the solution of coefficients at each \code{lambda} value}
\usage{
plotglmADMM(fit)
}
\arguments{
  \item{fit}{the output of the \code{glmADMM} function}
}

\details{Note that we do not plot \code{intercept} if \code{intercept=TRUE}, because \code{intercept} is not involved in the shrinkage}

\examples{
data(GaussianExample)
fit=glmADMM(A,b,family="gaussian")
plotglmADMM(fit)
}

\keyword{plot}

% \item{\dots}{Other graphical parameters to plot}
