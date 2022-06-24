#' Posterior distribution summary
#'
#' @param rd_est Posterior distribution for rd estimates.
#' @param rr_est Posterior distribution for rr estimates.
#' @param or_est Posterior distribution for rr estimates.
#'
#' @return A list of causal estimands including risk difference (rd),
#' odds ratios (or) and relative risk (rr)
#' between different treatment groups.
posterior_summary <- function(rd_est, rr_est, or_est) {
  # Risk difference (rd)
  rd_mean <- mean(rd_est)
  rd_se <- stats::sd(rd_est)
  rd_lower <- stats::quantile(rd_est, probs = 0.025, na.rm = T)
  rd_upper <- stats::quantile(rd_est, probs = 0.975, na.rm = T)

  # Relative risk (rr)
  rr_mean <- mean(rr_est)
  rr_se <- stats::sd(rr_est)
  rr_lower <- stats::quantile(rr_est, probs = 0.025, na.rm = T)
  rr_upper <- stats::quantile(rr_est, probs = 0.975, na.rm = T)

  # Odds ratio (or)
  or_mean <- mean(or_est)
  or_se <- stats::sd(or_est)
  or_lower <- stats::quantile(or_est, probs = 0.025, na.rm = T)
  or_upper <- stats::quantile(or_est, probs = 0.975, na.rm = T)

  # summarize results
  rd <- c(rd_mean, rd_se, rd_lower, rd_upper)
  rr <- c(rr_mean, rr_se, rr_lower, rr_upper)
  or <- c(or_mean, or_se, or_lower, or_upper)

  res <- rbind(rd, rr, or)
  rownames(res) <- c("RD", "RR", "OR")
  colnames(res) <- c("EST", "SE", "LOWER", "UPPER")
  return(res)
}
