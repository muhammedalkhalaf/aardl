#' Print method for aardl objects
#'
#' @param x An object of class \code{"aardl"}.
#' @param digits Integer. Significant digits.
#' @param ... Further arguments (ignored).
#' @return Invisibly returns \code{x}.
#' @export
print.aardl <- function(x, digits = 4L, ...) {
  type_labels <- c(
    aardl    = "Augmented ARDL (A-ARDL)",
    baardl   = "Bootstrap Augmented ARDL (BA-ARDL)",
    faardl   = "Fourier Augmented ARDL (FA-ARDL)",
    fbaardl  = "Fourier Bootstrap Augmented ARDL (FBA-ARDL)",
    nardl    = "Augmented NARDL (A-NARDL)",
    fanardl  = "Fourier Augmented NARDL (FA-NARDL)",
    banardl  = "Bootstrap Augmented NARDL (BA-NARDL)",
    fbanardl = "Fourier Bootstrap Augmented NARDL (FBA-NARDL)"
  )

  message(strrep("-", 65))
  message(type_labels[x$type])
  message(strrep("-", 65))
  message(sprintf("Dependent variable  : %s", x$depvar))
  message(sprintf("Independent vars    : %s",
                  paste(x$indepvars, collapse = ", ")))
  message(sprintf("N (used)            : %d", x$nobs))
  message(sprintf("Optimal p           : %d", x$opt_p))
  for (nm in names(x$opt_q)) {
    message(sprintf("  q[%s]             : %d", nm, x$opt_q[[nm]]))
  }
  if (x$kstar > 0) {
    message(sprintf("Fourier k*          : %.2f", x$kstar))
  }
  message(sprintf("R-squared           : %.6f", x$r2))
  message(sprintf("Adj. R-squared      : %.6f", x$r2_adj))
  message(sprintf("AIC                 : %.4f", x$aic))
  message(sprintf("BIC                 : %.4f", x$bic))
  message(strrep("-", 65))
  message("Augmented 3-Test Cointegration Framework (Sam et al. 2019)")
  message(strrep("-", 65))
  message(sprintf("  F_overall (F_pss)  : %.*f   p = %.*f  %s",
                  digits, x$F_pss,
                  digits, x$p_F_pss,
                  .aardl_stars(x$p_F_pss)))
  message(sprintf("  t_DV (t_pss)       : %.*f   p = %.*f  %s",
                  digits, x$t_pss,
                  digits, x$p_t_pss,
                  .aardl_stars(x$p_t_pss)))
  message(sprintf("  F_ind              : %.*f   p = %.*f  %s",
                  digits, x$F_ind,
                  digits, x$p_F_ind,
                  .aardl_stars(x$p_F_ind)))
  message(strrep("-", 65))
  coint_msg <- switch(x$coint_status,
    cointegrated    = ">>> COINTEGRATION: All 3 tests significant",
    degenerate_1    = ">>> Degenerate case 1: No cointegration",
    degenerate_2    = ">>> Degenerate case 2: Dependent may be I(0)",
    no_cointegration = ">>> NO cointegration"
  )
  message(coint_msg)
  message(strrep("-", 65))
  message("Long-run / EC-form coefficients:")
  cfs <- x$coefficients
  for (nm in names(cfs)) {
    message(sprintf("  %-20s : %.*f", nm, digits, cfs[[nm]]))
  }
  message(strrep("-", 65))
  message("*** p<0.01, ** p<0.05, * p<0.10")
  invisible(x)
}


#' Summary method for aardl objects
#'
#' @param object An object of class \code{"aardl"}.
#' @param ... Further arguments passed to \code{print.aardl}.
#' @return Invisibly returns \code{object}.
#' @export
summary.aardl <- function(object, ...) {
  print(object, ...)
}


#' @keywords internal
.aardl_stars <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.01) "***" else if (p < 0.05) "**" else if (p < 0.10) "*" else ""
}
