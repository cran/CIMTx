#' Regression adjustment with multivariate spline of GPS (RAMS) for ATE estimation with bootstrapping
#'
#'This function implements bootstrapping for the RAMS method when estimand is ATE. Please use our main function ce_estimate.R.
#'
#' @param y numeric vector for the binary outcome
#' @param w numeric vector for the treatment indicator
#' @param x dataframe including the treatment indicator and the covariates
#' @param method a character string. Users can selected from the following methods including "RAMS-Multinomial", "RAMS-GBM", "RAMS-SL"
#' @param nboots a numeric value representing the number of bootstrap samples
#' @param ... Other paramters to be passed to twang::mnps()
#'
#' @return a list with w*(w-1)/2 elements for ATE effect. Each element of the list contains the estimation, standard error, lower and upper 95\% CI for RD/RR/OR
#' @export
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
#'"rnorm(300, 0, 0.5)",# x1
#'"rbeta(300, 2, .4)",   # x2
#'"runif(300, 0, 0.5)",# x3
#'"rweibull(300,1,2)",  # x4
#'"rbinom(300, 1, .4)"# x5
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
#'ce_estimate_rams_ate_boot(y = data$y, x = data$covariates ,
#'w = data$w, method = "RAMS-Multinomial",nboots = 1)
ce_estimate_rams_ate_boot <- function(y, w, x, method, nboots,...){

  trim_alpha <- parent.frame()$trim_alpha
  SL.library <- parent.frame()$SL.library
  n_trt <- length(unique(w))
  for (i in 1:n_trt){
    assign(paste0("RAMS_ate_result_",i, "_all"), NULL)
  }

  names_result <- NULL
  for (j in 1:nboots) {
    bootstrap_id <- sample(length(y), replace = T)
    y_boot <- y[bootstrap_id]
    trt_boot <- w[bootstrap_id]
    x_boot <- x[bootstrap_id,]
    RAMS_ate_result <- ce_estimate_rams_ate(y = y_boot,
                                x = x_boot,
                                w = trt_boot,
                                method = method, ...)
    names_result <- names(RAMS_ate_result)
    for (i in 1:n_trt){
      assign(paste0("RAMS_ate_result_",i, "_all"), cbind(eval(parse(text = paste0("RAMS_ate_result_",i, "_all"))), RAMS_ate_result[[i]]))
    }
    print(paste0("Finish bootstrapping ", j))
  }

  for (i in 1:n_trt){
    assign(paste0("RD_mean_",i), mean(eval(parse(text = paste0("RAMS_ate_result_", i, "_all")))[1,]))
    assign(paste0("RD_se_",i), stats::sd(eval(parse(text = paste0("RAMS_ate_result_", i, "_all")))[1,]))
    assign(paste0("RD_lower_",i), stats::quantile(eval(parse(text = paste0("RAMS_ate_result_", i, "_all")))[1,], probs=0.025, na.rm = T))
    assign(paste0("RD_upper_",i), stats::quantile(eval(parse(text = paste0("RAMS_ate_result_", i, "_all")))[1,], probs=0.975, na.rm = T))
    assign(paste0("RR_mean_",i), mean(eval(parse(text = paste0("RAMS_ate_result_", i, "_all")))[2,]))
    assign(paste0("RR_se_",i), stats::sd(eval(parse(text = paste0("RAMS_ate_result_", i, "_all")))[2,]))
    assign(paste0("RR_lower_",i), stats::quantile(eval(parse(text = paste0("RAMS_ate_result_", i, "_all")))[2,], probs=0.025, na.rm = T))
    assign(paste0("RR_upper_",i), stats::quantile(eval(parse(text = paste0("RAMS_ate_result_", i, "_all")))[2,], probs=0.975, na.rm = T))
    assign(paste0("OR_mean_",i), mean(eval(parse(text = paste0("RAMS_ate_result_", i, "_all")))[3,]))
    assign(paste0("OR_se_",i), stats::sd(eval(parse(text = paste0("RAMS_ate_result_", i, "_all")))[3,]))
    assign(paste0("OR_lower_",i), stats::quantile(eval(parse(text = paste0("RAMS_ate_result_", i, "_all")))[3,], probs=0.025, na.rm = T))
    assign(paste0("OR_upper_",i), stats::quantile(eval(parse(text = paste0("RAMS_ate_result_", i, "_all")))[3,], probs=0.975, na.rm = T))
    # summarize results
    assign(paste0("RD_",i), c(eval(parse(text = paste0("RD_mean_",i))),
                              eval(parse(text = paste0("RD_se_",i))),
                              eval(parse(text = paste0("RD_lower_",i))),
                              eval(parse(text = paste0("RD_upper_",i)))
    ))
    assign(paste0("RR_",i), c(eval(parse(text = paste0("RR_mean_",i))),
                              eval(parse(text = paste0("RR_se_",i))),
                              eval(parse(text = paste0("RR_lower_",i))),
                              eval(parse(text = paste0("RR_upper_",i)))
    ))
    assign(paste0("OR_",i), c(eval(parse(text = paste0("OR_mean_",i))),
                              eval(parse(text = paste0("OR_se_",i))),
                              eval(parse(text = paste0("OR_lower_",i))),
                              eval(parse(text = paste0("OR_upper_",i)))
    ))
    assign(paste0("res_",i), rbind(eval(parse(text = paste0("RD_",i))),
                                   eval(parse(text = paste0("RR_",i))),
                                   eval(parse(text = paste0("OR_",i)))
    ))

  }
  result_list <- NULL
  for (i in 1:n_trt){
    result_list <- c(result_list, list(round(eval(parse(text = paste0("res_",i))),2)))
  }
  for (i in 1:n_trt){
    colnames(result_list[[i]]) <- c("EST","SE","LOWER","UPPER")
    rownames(result_list[[i]]) <- c("RD", "RR", "OR")
  }

  names(result_list) <- names_result[1:n_trt]

  return(result_list)
}
