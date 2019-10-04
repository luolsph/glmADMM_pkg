\name{glmADMM}
\alias{glmADMM}
\title{fit a GLM with lasso regularization and equality constraint if specified}
\description{
  Fit a sequence of generalized linear model by ADMM (Alternating Direction Method of Multipliers). The sequence is implied by \code{lambda} which is the regularization parameter. Fits linear, logistic and Poisson regression models. Users can add an equality constraint to the coefficients. In order to solve logistic and Poisson regression efficiently, users can also use inexact-ADMM to get the solution.
}
\usage{
glmADMM(A, b, lambda_seq, family, rho = 1, intercept = TRUE, equality=FALSE, C, d, inexact=FALSE)
}

\arguments{
  \item{A}{input feature matrix, of dimension nobs x nvars; each row is an observation vector.}
  \item{b}{the response vector. Quantitative for \code{family="gaussian"},
  or non-negative counts for \code{family="poisson"}. For
  \code{family="binomial"}, it should be a factor with two levels indicated by either \code{1} or \code{-1}.}
  \item{lambda_seq}{A user supplied \code{lambda} sequence. By default, the program compute its own \code{lambda} sequence based on the data and run the regression based on a decreasing sequence of \code{lambda} values. The default number of \code{lambda} is \code{100} and we use a warm-start in runing this sequence of regressions.}
  \item{family}{Response type, including "gaussian", "binomial" and "poisson".}
  \item{rho}{the augmented Lagrangian parameter in ADMM, usually set to 1.}
  \item{intercept}{Should intercept(s) be fitted (default=TRUE) or igonored (FALSE)}
  \item{equality}{Should there be an equality constraint to the coefficients (TRUE) or no constraint (default=FALSE)}
  \item{C}{a designed matrix that add equality constraint to the coefficients, dimension of \code{1} row and number of columns matches the number of observations (number of rows in \code{A}). Note that if an intercept is fitted, adding a zero to the left of the first column, meaning do not constrain the intercept. If \code{C} is not supplied and \code{equality=TRUE}, the default design is to sum up all the coefficients.}
  \item{d}{specify the value on the right-hand side of the equality constraint. If not specified, \code{d=0}.}
  \item{inexact}{Should the program solve \code{x} via inexact-ADMM (TRUE) or not (default=FALSE).}
  }
\value{
\item{solution}{solution of \code{x}, each column represents the solution at a \code{lambda} value.}
\item{time}{Program running time.}
\item{iter}{iteration number at each \code{lambda}.}
\item{lambda}{The actual sequence of \code{lambda} values used.}
\item{intercept}{record if an intercept is fitted (intercept=1) or not (intercept=0)}
}
\references{
Boyd S., Parikh N., Chu E., Peleato B. and Eckstein J. \emph{Distributed optimization and statistical learning via the alternating direction method of multipliers}, \url{https://web.stanford.edu/~boyd/papers/pdf/admm_distr_stats.pdf}\cr \emph{Foundations and Trends in Machine Learning, 3(1)(1-122), 2011}\cr
Chang T., Hong M. and Wang X. \emph{Multi-agent distributed optimization via inexact consensus admm}, \url{http://arxiv.org/pdf/1402.6065v2.pdf}\cr \emph{IEEE TRANSACTIONS
ON SIGNAL PROCESSING, 63(482-497), 2015}\cr
Ivanoff S., Picard F. and Rivoirard V. \emph{Adaptive Lasso and group-Lasso for functional Poisson regression}, \url{http://arxiv.org/pdf/1412.6966v2.pdf}\cr \emph{arXiv:1412.6966 [stat.ME] 26 Dec 2014}\cr
}
\author{Hui Jiang, Lan Luo\cr Maintainer: Lan Luo \email{luolsph@umich.edu}}

\examples{
# Gaussian
data(GaussianExample)
fit1=glmADMM(A,b,family="gaussian")

# Binomial
data(BinomialExample)
fit2=glmADMM(A,b,family="binomial",inexact=TRUE)

# Poisson
data(PoissonExample)
fit3=glmADMM(A,b,family="poisson",inexact=TRUE)

}


