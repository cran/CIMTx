#' Summarize a CIMTx_ATT_sa object
#'
#' @param object a \code{CIMTx_ATT_sa} object obtained with \code{\link{sa}} function.
#' @param ... further arguments passed to or from other methods.
#'
#' @return a data frame containing the estimation, standard error, lower and upper 95\% CI for the causal estimand in terms of RD.
#' @export
#'
#' @examples
#' \donttest{
#'lp_w_all <-
#'  c(".4*x1 + .1*x2  - 1.1*x4 + 1.1*x5",    # w = 1
#'    ".2 * x1 + .2 * x2  - 1.2 * x4 - 1.3 * x5")  # w = 2
#'nlp_w_all <-
#'  c("-.5*x1*x4  - .1*x2*x5", # w = 1
#'    "-.3*x1*x4 + .2*x2*x5")# w = 2
#'lp_y_all <- rep(".2*x1 + .3*x2 - .1*x3 - 1.1*x4 - 1.2*x5", 3)
#'nlp_y_all <- rep(".7*x1*x1  - .1*x2*x3", 3)
#'X_all <- c(
#'  "rnorm(100, 0, 0.5)",# x1
#'  "rbeta(100, 2, .4)",   # x2
#'  "runif(100, 0, 0.5)",# x3
#'  "rweibull(100,1,2)",  # x4
#'  "rbinom(100, 1, .4)"# x5
#')
#'set.seed(1111)
#'data <- data_sim(
#'  sample_size = 100,
#'  n_trt = 3,
#'  X = X_all,
#'  lp_y = lp_y_all,
#'  nlp_y  = nlp_y_all,
#'  align = FALSE,
#'  lp_w = lp_w_all,
#'  nlp_w = nlp_w_all,
#'  tau = c(0.5,-0.5,0.5),
#'  delta = c(0.5,0.5),
#'  psi = 2
#')
#'c_grid <- c(
#'  "runif(-0.6, 0)",# c(1,2)
#'  "runif(0, 0.6)",# c(2,1)
#'  "runif(-0.6, 0)", # c(2,3)
#'  "seq(-0.6, 0, by = 0.3)", # c(1,3)
#'  "seq(0, 0.6, by = 0.3)", # c(3,1)
#'  "runif(0, 0.6)" # c(3,2)
#')
#'sensitivity_analysis_parallel_ATT_result <-
#'  sa(
#'    M1 = 1,
#'    x = data$covariates,
#'    y = data$y,
#'    w = data$w,
#'    prior_c_function = c_grid,
#'    nCores = 1,
#'    estimand = "ATE",
#'  )
#'  summary(sensitivity_analysis_parallel_ATT_result)
#'  }
summary.CIMTx_ATT_sa <- function(object,...){
  object <- object[stringr::str_detect(names(object), "ATT")]
  reference_trt <- as.integer(stringr::str_sub(  names(object)[1],7,7))
  n_trt <- length(unique(as.integer(stringr::str_sub(  names(object),8,8)))) + 1
  w_ind_no_reference <- unique(as.integer(stringr::str_sub(  names(object),8,8)))
result_final <- NULL
counter <- 1
for (k in 1:(n_trt-1)){
  assign(paste0("mean",reference_trt,w_ind_no_reference[k]), mean(object[[paste0("ATT_RD",reference_trt,w_ind_no_reference[k])]]))
  assign(paste0("sd",reference_trt,w_ind_no_reference[k]), stats::sd(object[[paste0("ATT_RD",reference_trt,w_ind_no_reference[k])]]))
  assign(paste0("lower",reference_trt,w_ind_no_reference[k]), eval(parse(text = paste0("mean",reference_trt,w_ind_no_reference[k])))- 1.96 * eval(parse(text = paste0("sd",reference_trt,w_ind_no_reference[k]))))
  assign(paste0("upper",reference_trt,w_ind_no_reference[k]), eval(parse(text = paste0("mean",reference_trt,w_ind_no_reference[k])))+ 1.96 * eval(parse(text = paste0("sd",reference_trt,w_ind_no_reference[k]))))
  assign(paste0("RD",reference_trt,w_ind_no_reference[k]), round(c(eval(parse(text = paste0("mean",reference_trt,w_ind_no_reference[k]))), eval(parse(text = paste0("sd",reference_trt,w_ind_no_reference[k]))), eval(parse(text = paste0("lower",reference_trt,w_ind_no_reference[k]))), eval(parse(text = paste0("upper",reference_trt,w_ind_no_reference[k])))),2))

  result_final <- rbind(result_final, eval(parse(text = paste0("RD",reference_trt,w_ind_no_reference[k]))))
  rownames(result_final)[[counter]] <- paste0("ATT_RD",reference_trt,w_ind_no_reference[k])
  counter <- counter +1
}
colnames(result_final) <- c("EST","SE","LOWER","UPPER")
return(result_final)
}
