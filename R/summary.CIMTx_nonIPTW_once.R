#' Summarize a CIMTx_nonIPTW_once object
#'
#' @param object a \code{CIMTx_nonIPTW_once} object
#' @param ... further arguments passed to or from other methods.
#'
#' @return a data frame with ATT/ATE effect estimates in terms of
#' RD, RR and OR.
#' @importFrom dplyr select slice pull
#' @importFrom stringr str_detect str_sub
#' @export
#'
#' @references
#'
#' Hadley Wickham, Romain François, Lionel Henry and Kirill Müller (2021).
#' \emph{dplyr: A Grammar of Data Manipulation}.
#' R package version 1.0.7.
#' URL: \url{https://CRAN.R-project.org/package=dplyr}
#'
#' Hadley Wickham (2019).
#' \emph{stringr: Simple, Consistent Wrappers for Common String Operations}.
#' R package version 1.4.0.
#' URL:\url{https://CRAN.R-project.org/package=stringr}
#'
#' @examples
#' \donttest{
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
#' iptw_tmle_res <- ce_estimate(
#'   y = data$y, x = data$covariates,
#'   w = data$w, method = "TMLE", estimand = "ATE",
#'   sl_library = c("SL.glm", "SL.glmnet")
#' )
#' summary(iptw_tmle_res)
#' }
summary.CIMTx_nonIPTW_once <- function(object, ...) {
  if (object$method == "VM") {
    object <- object[stringr::str_detect(names(object), "ATT")]
    reference_trt <-
      as.integer(stringr::str_sub(names(object)[1], 4, 4))
    trt_indicator_no_reference <-
      unique(as.integer(stringr::str_sub(names(object), 5, 5)))
    result_final <- NULL
    for (j in seq_len(length(trt_indicator_no_reference))) {
      result_once <-
        object[[paste0("ATT", reference_trt, trt_indicator_no_reference[j])]]
      result_once <- round(result_once, digits = 2)
      colnames(result_once) <-
        paste0("ATT", reference_trt, trt_indicator_no_reference[j])
      result_final <- cbind(result_final, result_once)
    }
    return(result_final)
  }
  if (object$method == "TMLE") {
    n_trt <- object$n_trt
    w_result_one_repetition <- object$result_TMLE
    for (i in 1:n_trt) {
      assign(
        paste0("mu_", i, "_hat"),
        w_result_one_repetition %>%
          dplyr::select("EYt") %>%
          dplyr::slice(i) %>%
          dplyr::pull("EYt")
      )
    }
    for (i in 1:n_trt) {
      assign(
        paste0("lower_", i, "_hat"),
        w_result_one_repetition %>%
          dplyr::select("CI1") %>%
          dplyr::slice(i) %>%
          dplyr::pull("CI1")
      )
    }
    for (i in 1:n_trt) {
      assign(
        paste0("upper_", i, "_hat"),
        w_result_one_repetition %>%
          dplyr::select("CI2") %>%
          dplyr::slice(i) %>%
          dplyr::pull("CI2")
      )
    }
    result_final <- NULL
    for (i in 1:(n_trt - 1)) {
      result_once <- NULL
      for (j in (i + 1):n_trt) {
        assign(paste0("RD", i, j), eval(parse(text = paste0(
          "mu_", i, "_hat"
        ))) - eval(parse(text = paste0(
          "mu_", j, "_hat"
        ))))
        assign(paste0("RR", i, j), eval(parse(text = paste0(
          "mu_", i, "_hat"
        ))) / eval(parse(text = paste0(
          "mu_", j, "_hat"
        ))))
        assign(paste0("OR", i, j), (eval(parse(
          text = paste0("mu_", i, "_hat")
        )) / (1 - eval(
          parse(text = paste0("mu_", i, "_hat"))
        ))) / (eval(parse(
          text = paste0("mu_", j, "_hat")
        )) / (1 - eval(
          parse(text = paste0("mu_", j, "_hat"))
        ))))
        result_once <-
          rbind(eval(parse(text = paste0(
            "round(RD", i, j, ",2)"
          ))), eval(parse(text = paste0(
            "round(RR", i, j, ",2)"
          ))), eval(parse(text = paste0(
            "round(OR", i, j, ",2)"
          ))))
        colnames(result_once) <- "EST"
        rownames(result_once) <- c("RD", "RR", "OR")
        colnames(result_once) <- paste0("ATE", i, j)
        result_final <- cbind(result_final, result_once)
      }
    }
    return(result_final)
  }
  if (object$method %in% c("RAMS-Multinomial", "RAMS-SL",
                           "RAMS-GBM")) {
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
}
