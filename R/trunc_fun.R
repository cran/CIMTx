
trunc_fun <- function(x, trim_perc = 0.05) {
  if (length(trim_perc) == 1 && trim_perc < 0.5) {
    pmax(stats::quantile(x, trim_perc),x)
  } else if (length(trim_perc) == 1 && trim_perc >= 0.5) {
    pmin(stats::quantile(x, trim_perc), x)
  } else if (length(trim_perc) == 2) {
    pmin(stats::quantile(x, trim_perc[2]), pmax(stats::quantile(x, trim_perc[1]), x))
  }
}
