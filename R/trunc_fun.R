#' Function to truncate weight
#'
#' @param x weight
#' @param trim_perc the percentile at which the inverse probability of treatment weights should be trimmed.
#'
#' @return Truncated weights
#' @export
#'
#' @examples
#' trunc_fun(rnorm(1000,0,1),0.05)
trunc_fun <- function(x, trim_perc = 0.05) {
  if (length(trim_perc) == 1 && trim_perc < 0.5) {
    pmax(stats::quantile(x, trim_perc),x)
  } else if (length(trim_perc) == 1 && trim_perc >= 0.5) {
    pmin(stats::quantile(x, trim_perc), x)
  } else if (length(trim_perc) == 2) {
    pmin(stats::quantile(x, trim_perc[2]), pmax(stats::quantile(x, trim_perc[1]), x))
  }
}
