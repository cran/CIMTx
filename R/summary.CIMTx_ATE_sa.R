#' Summarize a CIMTx_ATE_sa object
#'
#' @param object a \code{CIMTx_ATE_sa} object obtained with
#' \code{\link{sa}} function.
#' @param ... further arguments passed to or from other methods.
#'
#' @return a data frame containing the estimation, standard error,
#' lower and upper 95\% CI for the causal estimand in terms of RD.
#' @importFrom stringr str_detect str_sub
#' @export
#'
#' @references
#' Hadley Wickham (2019).
#' \emph{stringr: Simple, Consistent Wrappers for Common String Operations}.
#' R package version 1.4.0.
#' URL:\url{https://CRAN.R-project.org/package=stringr}
#' @examples
#' \donttest{
#' lp_w_all <-
#'   c(
#'     ".4*x1 + .1*x2  - 1.1*x4 + 1.1*x5", # w = 1
#'     ".2 * x1 + .2 * x2  - 1.2 * x4 - 1.3 * x5"
#'   ) # w = 2
#' nlp_w_all <-
#'   c(
#'     "-.5*x1*x4  - .1*x2*x5", # w = 1
#'     "-.3*x1*x4 + .2*x2*x5"
#'   ) # w = 2
#' lp_y_all <- rep(".2*x1 + .3*x2 - .1*x3 - 1.1*x4 - 1.2*x5", 3)
#' nlp_y_all <- rep(".7*x1*x1  - .1*x2*x3", 3)
#' X_all <- c(
#'   "rnorm(0, 0.5)", # x1
#'   "rbeta(2, .4)", # x2
#'   "runif(0, 0.5)", # x3
#'   "rweibull(1,2)", # x4
#'   "rbinom(1, .4)" # x5
#' )
#' set.seed(1111)
#' data <- data_sim(
#'   sample_size = 100,
#'   n_trt = 3,
#'   x = X_all,
#'   lp_y = lp_y_all,
#'   nlp_y = nlp_y_all,
#'   align = FALSE,
#'   lp_w = lp_w_all,
#'   nlp_w = nlp_w_all,
#'   tau = c(0.5, -0.5, 0.5),
#'   delta = c(0.5, 0.5),
#'   psi = 2
#' )
#' c_grid <- c(
#'   "runif(-0.6, 0)", # c(1,2)
#'   "runif(0, 0.6)", # c(2,1)
#'   "runif(-0.6, 0)", # c(2,3)
#'   "seq(-0.6, 0, by = 0.3)", # c(1,3)
#'   "seq(0, 0.6, by = 0.3)", # c(3,1)
#'   "runif(0, 0.6)" # c(3,2)
#' )
#' sensitivity_analysis_parallel_ATE_result <-
#'   sa(
#'     m1 = 1,
#'     x = data$covariates,
#'     y = data$y,
#'     w = data$w,
#'     prior_c_function = c_grid,
#'     nCores = 1,
#'     estimand = "ATE",
#'   )
#' summary(sensitivity_analysis_parallel_ATE_result)
#' }
summary.CIMTx_ATE_sa <- function(object, ...) {
  object <- object[stringr::str_detect(names(object), "ATE")]
  n_trt <-
    length(unique(as.integer(stringr::str_sub(
      names(object), 8, 8
    )))) + 1
  result_final <- NULL
  counter <- 1
  for (k in 1:(n_trt - 1)) {
    for (m in (k + 1):n_trt) {
      assign(paste0("mean", k, m), mean(object[[paste0("ATE_RD", k, m)]]))
      assign(paste0("sd", k, m), stats::sd(object[[paste0("ATE_RD", k, m)]]))
      assign(paste0("lower", k, m), eval(parse(text = paste0("mean", k, m))) -
               1.96 * eval(parse(text = paste0("sd", k, m))))
      assign(paste0("upper", k, m), eval(parse(text = paste0("mean", k, m))) +
               1.96 * eval(parse(text = paste0("sd", k, m))))
      assign(paste0("RD", k, m), round(c(
        eval(parse(text = paste0("mean", k, m))),
        eval(parse(text = paste0("sd", k, m))),
        eval(parse(text = paste0("lower", k, m))),
        eval(parse(text = paste0("upper", k, m)))
      ), 2))
      result_final <-
        rbind(result_final, eval(parse(text = paste0("RD", k, m))))
      rownames(result_final)[[counter]] <- paste0("ATE_RD", k, m)
      counter <- counter + 1
    }
  }
  colnames(result_final) <- c("EST", "SE", "LOWER", "UPPER")
  return(result_final)
}
