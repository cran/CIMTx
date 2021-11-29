
#' Summarize a CIMTx_ATE_posterior object
#'
#' @param object a \code{CIMTx_ATE_posterior} object obtained with \code{\link{ce_estimate}} function.
#' @param ... further arguments passed to or from other methods.
#'
#' @return a list with w*(w-1)/2 elements for ATE effect. Each element of the list contains the estimation, standard error, lower and upper 95\% CI for RD/RR/OR.
#' @export
#'
#' @examples
#' library(CIMTx)
#'lp_w_all <-
#'  c(".4*x1 + .1*x2  - .1*x4 + .1*x5",    # w = 1
#'    ".2 * x1 + .2 * x2  - .2 * x4 - .3 * x5")  # w = 2
#'nlp_w_all <-
#'  c("-.5*x1*x4  - .1*x2*x5", # w = 1
#'    "-.3*x1*x4 + .2*x2*x5")# w = 2
#'lp_y_all <- rep(".2*x1 + .3*x2 - .1*x3 - .1*x4 - .2*x5", 3)
#'nlp_y_all <- rep(".7*x1*x1  - .1*x2*x3", 3)
#'X_all <- c(
#'  "rnorm(300, 0, 0.5)",# x1
#'  "rbeta(300, 2, .4)",   # x2
#'  "runif(300, 0, 0.5)",# x3
#'  "rweibull(300,1,2)",  # x4
#'  "rbinom(300, 1, .4)"# x5
#')

#'set.seed(111111)
#'data <- data_sim(
#'  sample_size = 300,
#'  n_trt = 3,
#'  X = X_all,
#'  lp_y = lp_y_all,
#'  nlp_y  = nlp_y_all,
#'  align = FALSE,
#'  lp_w = lp_w_all,
#'  nlp_w = nlp_w_all,
#'  tau = c(-1.5,0,1.5),
#'  delta = c(0.5,0.5),
#'  psi = 1
#')
#'ce_estimate_ra_ate_result <- ce_estimate(y = data$y, x = data$covariates
#', w = data$w, ndpost = 10,method = "RA", estimand = "ATE")
#'summary(ce_estimate_ra_ate_result)
summary.CIMTx_ATE_posterior <- function(object,...){
  object <- object[stringr::str_detect(names(object), "ATE")]
  n_trt <- length(unique(as.integer(stringr::str_sub(  names(object),8,8)))) + 1
  result_summary <- NULL
  for (i in 1:(n_trt-1)){
    for (j in (i + 1):n_trt){
      assign(paste0("ate",i,j), posterior_summary(object[[paste0("ATE_RD",i,j)]], object[[paste0("ATE_RR",i,j)]], object[[paste0("ATE_OR",i,j)]]))
      assign(paste0("ATE",i,j), list(round(eval(parse(text =(paste0("ate",i,j)))), digits = 2)))
      assign(paste0("ATE",i,j), stats::setNames(eval(parse(text =(paste0("ATE",i,j)))), paste0("ATE",i,j)))
      result_summary <- c(result_summary, (eval(parse(text =(paste0("ATE",i,j))))))
    }
  }
  return(result_summary)
}



