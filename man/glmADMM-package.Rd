\name{glmADMM-package}
\alias{glmADMM-package}
\docType{package}
\title{
Lasso regularization path for some generalized linear models
}
\description{
This package fits lasso paths for linear regression, logistic and Poisson regression. It can also incorporate an equality constraint on the coefficients. In addition, this package can solve logistic and Poisson regression via inexact ADMM which is more efficient.
}
\details{
\tabular{ll}{
Package: \tab glmADMM\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2015-07-13\cr
License: \tab GPL-2\cr
}
Accepts \code{A,b} data for regression models, and produces the regularization path over a default of 100 tuning parameter \code{lambda} values. Note that \code{C} and \code{d} may need to be specified only when \code{equality=TRUE}.
Only 2 function to use:
\code{glmADMM}\cr
\code{plotglmADMM}
}
\author{Hui Jiang, Lan Luo\cr
Maintainer: Lan Luo <luolsph@umich.edu>
}
\references{
Boyd S., Parikh N., Chu E., Peleato B. and Eckstein J. \emph{Distributed optimization and statistical learning via the alternating direction method of multipliers}, \url{https://web.stanford.edu/~boyd/papers/pdf/admm_distr_stats.pdf}\cr \emph{Foundations and Trends in Machine Learning, 3(1)(1-122), 2011}\cr
Chang T., Hong M. and Wang X. \emph{Multi-agent distributed optimization via inexact consensus admm}, \url{http://arxiv.org/pdf/1402.6065v2.pdf}\cr \emph{IEEE TRANSACTIONS
ON SIGNAL PROCESSING, 63(482-497), 2015}\cr
Ivanoff S., Picard F. and Rivoirard V. \emph{Adaptive Lasso and group-Lasso for functional Poisson regression}, \url{http://arxiv.org/pdf/1412.6966v2.pdf}\cr \emph{arXiv:1412.6966 [stat.ME] 26 Dec 2014}\cr
}
\keyword{ package }

