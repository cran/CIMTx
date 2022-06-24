#' Summarize a CIMTx_IPTW object
#'
#' @param object a \code{CIMTx_IPTW} object
#' @param ... further arguments passed to or from other methods.
#'
#' @return a data frame with ATT/ATE effect estimates in terms of
#' RD, RR and OR.
#' @importFrom stringr str_detect str_sub
#' @export
#'
#' @references
#' Hadley Wickham (2019).
#' \emph{stringr: Simple, Consistent Wrappers for Common String Operations}.
#' R package version 1.4.0.
#' URL:\url{https://CRAN.R-project.org/package=stringr}
#'
#' @examples
#' lp_w_all <-
#'   c(
#'     ".4*x1 + .1*x2  - .1*x4 + .1*x5", # w = 1
#'     ".2 * x1 + .2 * x2  - .2 * x4 - .3 * x5"
#'   ) # w = 2
#' nlp_w_all <-
#'   c(
#'     "-.5*x1*x4  - .1*x2*x5", # w = 1
#'     "-.3*x1*x4 + .2*x2*x5"
#'   ) # w = 2
#' lp_y_all <- rep(".2*x1 + .3*x2 - .1*x3 - .1*x4 - .2*x5", 3)
#' nlp_y_all <- rep(".7*x1*x1  - .1*x2*x3", 3)
#' X_all <- c(
#'   "rnorm(0, 0.5)", # x1
#'   "rbeta(2, .4)", # x2
#'   "runif(0, 0.5)", # x3
#'   "rweibull(1,2)", # x4
#'   "rbinom(1, .4)" # x5
#' )
#'
#' set.seed(111111)
#' data <- data_sim(
#'   sample_size = 300,
#'   n_trt = 3,
#'   x = X_all,
#'   lp_y = lp_y_all,
#'   nlp_y = nlp_y_all,
#'   align = FALSE,
#'   lp_w = lp_w_all,
#'   nlp_w = nlp_w_all,
#'   tau = c(-1.5, 0, 1.5),
#'   delta = c(0.5, 0.5),
#'   psi = 1
#' )
#' iptw_multi_res <- ce_estimate(
#'   y = data$y, x = data$covariates,
#'   w = data$w, method = "IPTW-Multinomial", estimand = "ATE"
#' )
#' summary(iptw_multi_res)
summary.CIMTx_IPTW <- function(object, ...) {
  if (object$estimand == "ATE") {
    object <- object[stringr::str_detect(names(object), "ATE")]
    n_trt <-
      length(unique(as.integer(stringr::str_sub(
        names(object), 5, 5
      )))) + 1
    result_summary <- NULL
    for (i in 1:(n_trt - 1)) {
      for (j in (i + 1):n_trt) {
        result_once <- object[[paste0("ATE", i, j)]]
        result_once <- round(result_once, digits = 2)
        colnames(result_once) <- paste0("ATE", i, j)
        result_summary <- cbind(result_summary, result_once)
      }
    }
    return(result_summary)
  }
  if (object$estimand == "ATT") {
    object <- object[stringr::str_detect(names(object), "ATT")]
    reference_trt <-
      as.integer(stringr::str_sub(names(object)[1], 4, 4))
    trt_indicator_no_reference <-
      unique(as.integer(stringr::str_sub(names(object), 5, 5)))
    result_summary <- NULL
    for (j in seq_len(length(trt_indicator_no_reference))) {
      result_once <-
        object[[paste0("ATT", reference_trt, trt_indicator_no_reference[j])]]
      result_once <- round(result_once, digits = 2)
      colnames(result_once) <-
        paste0("ATT", reference_trt, trt_indicator_no_reference[j])
      result_summary <- cbind(result_summary, result_once)
    }
    return(result_summary)
  }
}
