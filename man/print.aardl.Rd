\name{print.aardl}
\alias{print.aardl}
\alias{summary.aardl}
\title{Print and Summary Methods for aardl Objects}
\usage{
\method{print}{aardl}(x, digits = 4L, ...)

\method{summary}{aardl}(object, ...)
}
\arguments{
\item{x}{An object of class \code{"aardl"}.}

\item{digits}{Integer. Significant digits.}

\item{...}{Further arguments (ignored).}

\item{object}{An object of class \code{"aardl"}.}
}
\value{
Invisibly returns the input object.
}
\description{
Print and summary methods for objects of class \code{"aardl"} returned
by \code{\link{aardl}}.
}
\examples{
\donttest{
set.seed(1)
n  <- 60
df <- data.frame(y = cumsum(rnorm(n)), x1 = cumsum(rnorm(n)))
res <- aardl(y ~ x1, data = df, max_lag = 2)
print(res)
summary(res)
}
}
