# aardl

> **This repository is superseded and no longer maintained.**
> At the request of the CRAN team, this package was merged into the CRAN package
> [ardlverse](https://cran.r-project.org/package=ardlverse). The function `aardl()` is maintained there,
> with corrections that are not in this repository. The code here is an older version
> and should not be used for new work.
>
> ```r
> install.packages("ardlverse")
> ```

**Augmented ARDL Cointegration Analysis** for R

## Overview

`aardl` implements the Augmented ARDL (A-ARDL) cointegration framework
proposed by Sam, McNown and Goh (2019). The package resolves the degenerate
cases of the standard ARDL bounds test through an augmented three-test
framework (F-overall, t-DV, F-independent).

## Model Variants

| `type`      | Description                                      |
|-------------|--------------------------------------------------|
| `"aardl"`   | Augmented ARDL (default)                         |
| `"baardl"`  | Bootstrap Augmented ARDL                         |
| `"faardl"`  | Fourier Augmented ARDL                           |
| `"fbaardl"` | Fourier Bootstrap Augmented ARDL                 |
| `"nardl"`   | Augmented NARDL (nonlinear)                      |
| `"fanardl"` | Fourier Augmented NARDL                          |
| `"banardl"` | Bootstrap Augmented NARDL                        |
| `"fbanardl"`| Fourier Bootstrap Augmented NARDL                |

## Installation

```r
# Old version, GitHub only (see the notice at the top of this page)
# install.packages("remotes")
remotes::install_github("muhammedalkhalaf/aardl")
```

## Usage

```r
library(aardl)

set.seed(42)
n  <- 80
x1 <- cumsum(rnorm(n))
x2 <- cumsum(rnorm(n))
y  <- 0.5 * x1 + 0.3 * x2 + rnorm(n, sd = 0.5)
df <- data.frame(y = y, x1 = x1, x2 = x2)

result <- aardl(y ~ x1 + x2, data = df, max_lag = 3, ic = "bic")
print(result)
```

## References

Sam, C. Y., McNown, R. and Goh, S. K. (2019). Economics Letters, 174, 47–50.

Pesaran, M. H., Shin, Y. and Smith, R. J. (2001). Journal of Applied
Econometrics, 16(3), 289–326.

## License

GPL-3
