#' Trimming
#'
#' The function trims the weights for IPTW methods.
#'
#' @param x A numeric vector
#' @param trim_perc A numeric vector with length 2 indicating trimming
#' percentile.
#'
#' @return A numeric vector
#'
trunc_fun <- function(x, trim_perc = 0.05) {
    pmin(stats::quantile(x, trim_perc[2]),
         pmax(stats::quantile(x, trim_perc[1]), x))
}
