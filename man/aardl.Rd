\name{aardl}
\alias{aardl}
\title{Augmented ARDL Cointegration Analysis}
\usage{
aardl(formula, data, type = "aardl",
      decompose = NULL, max_lag = 4L, max_k = 3,
      ic = "bic", reps = 999L, level = 95,
      case = 3L, bootstrap = "bvz")
}
\arguments{
\item{formula}{A formula of the form \code{y ~ x1 + x2 + ...}.}

\item{data}{A data frame containing all variables in \code{formula}.}

\item{type}{Character string specifying the model variant. One of
\code{"aardl"} (default), \code{"baardl"}, \code{"faardl"},
\code{"fbaardl"}, \code{"nardl"}, \code{"fanardl"}, \code{"banardl"},
\code{"fbanardl"}.}

\item{decompose}{Character vector naming variables to decompose into
positive and negative partial sums (required for NARDL variants).}

\item{max_lag}{Positive integer. Maximum lag order for selection.
Default \code{4}.}

\item{max_k}{Numeric. Maximum Fourier frequency for selection in Fourier
variants. Default \code{3}.}

\item{ic}{Character. Information criterion: \code{"bic"} (default) or
\code{"aic"}.}

\item{reps}{Positive integer. Number of bootstrap replications for
bootstrap variants. Default \code{999}.}

\item{level}{Numeric (0, 100). Confidence level. Default \code{95}.}

\item{case}{Integer 1--5. PSS deterministic case. Default \code{3}.}

\item{bootstrap}{Character. Bootstrap method: \code{"bvz"} (default) or
\code{"mcnown"}.}
}
\value{
An object of class \code{"aardl"} with components:
\describe{
  \item{coefficients}{EC-form coefficients (ADJ, LR terms).}
  \item{vcov}{Variance-covariance matrix.}
  \item{F_pss}{F-statistic: joint significance of all lagged levels.}
  \item{t_pss}{t-statistic on the lagged dependent variable.}
  \item{F_ind}{F-statistic: lagged independent levels (Sam 2019).}
  \item{p_F_pss}{p-value of \code{F_pss}.}
  \item{p_t_pss}{Two-sided p-value of \code{t_pss}.}
  \item{p_F_ind}{p-value of \code{F_ind}.}
  \item{coint_status}{Cointegration conclusion.}
  \item{r2}{R-squared.}
  \item{r2_adj}{Adjusted R-squared.}
  \item{aic}{AIC.}
  \item{bic}{BIC.}
  \item{nobs}{Number of observations.}
  \item{opt_p}{Selected AR lag.}
  \item{opt_q}{Selected lags for each regressor.}
  \item{kstar}{Optimal Fourier frequency (Fourier models only).}
  \item{type}{Model type.}
  \item{depvar}{Dependent variable name.}
  \item{indepvars}{Independent variable names.}
  \item{ic}{Information criterion used.}
}
}
\description{
Estimates an Augmented Autoregressive Distributed Lag (A-ARDL) model and
tests for cointegration using the three-test framework of Sam, McNown and
Goh (2019). Eight model variants are supported.
}
\examples{
\donttest{
set.seed(42)
n  <- 80
x1 <- cumsum(rnorm(n))
x2 <- cumsum(rnorm(n))
y  <- 0.5 * x1 + 0.3 * x2 + rnorm(n, sd = 0.5)
df <- data.frame(y = y, x1 = x1, x2 = x2)
result <- aardl(y ~ x1 + x2, data = df, max_lag = 2, ic = "aic")
print(result)
}
}
\references{
Sam, C. Y., McNown, R. and Goh, S. K. (2019). An augmented autoregressive
distributed lag bounds test for cointegration. \emph{Economics Letters},
174, 47--50. \doi{10.1016/j.econlet.2018.12.007}

McNown, R., Sam, C. Y. and Goh, S. K. (2018). Bootstrapping the
Autoregressive Distributed Lag Test for Cointegration.
\emph{Applied Economics}, 50(13), 1509--1521.
\doi{10.1080/00036846.2017.1366150}

Shin, Y., Yu, B. and Greenwood-Nimmo, M. (2014). Modelling Asymmetric
Cointegration and Dynamic Multipliers in a Nonlinear ARDL Framework.
In Horrace, W. C. and Sickles, R. C. (Eds.), \emph{Festschrift in Honor
of Peter Schmidt}. Springer. \doi{10.1007/978-1-4899-8008-3_9}

Pesaran, M. H., Shin, Y. and Smith, R. J. (2001). Bounds testing
approaches to the analysis of level relationships.
\emph{Journal of Applied Econometrics}, 16(3), 289--326.
\doi{10.1002/jae.616}
}
