# aardl

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
install.packages("aardl")
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
