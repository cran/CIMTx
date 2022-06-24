#' Calculate the true c functions with 3 treatments and a binary predictor
#'
#' This function calculates the true confounding functions
#' with 3 treatments and a binary predictor for simulated data.
#'
#' @param x A matrix with one column for the binary predictor
#' with values 0 and 1
#' @param w A treatment indicator
#'
#' @return A matrix with 2 rows and 6 columns
#' @export
#'
#' @examples
#' set.seed(111)
#' data_SA <- data_sim(
#'   sample_size = 100,
#'   n_trt = 3,
#'   x = c(
#'     "rbinom(1, .5)", # x1:measured confounder
#'     "rbinom(1, .4)"
#'   ), # x2:unmeasured confounder
#'   lp_y = rep(".2*x1+2.3*x2", 3), # parallel response surfaces
#'   nlp_y = NULL,
#'   align = FALSE, # w model is not the same as the y model
#'   lp_w = c(
#'     "0.2 * x1 + 2.4 * x2", # w = 1
#'     "-0.3 * x1 - 2.8 * x2"
#'   ),
#'   nlp_w = NULL,
#'   tau = c(-2, 0, 2),
#'   delta = c(0, 0),
#'   psi = 1
#' )
#' x1 <- data_SA$covariates[, 1, drop = FALSE]
#' w <- data_SA$w
#' Y1 <- data_SA$Y_true[, 1]
#' Y2 <- data_SA$Y_true[, 2]
#' Y3 <- data_SA$Y_true[, 3]
#' true_c_fun <- true_c_fun_cal(x = x1, w = w)
true_c_fun_cal <-
  function(x, w) {
    # Calculate the true confounding functions within x = 1 and x= 0 stratum
    c_truth_list <- vector("list", length(unique(x[, 1])))
    x_unique <- unique(x[, 1])
    for (i in seq_len(length(unique(x[, 1])))) {
      assign(paste0("c_1_x_", x_unique[i]), eval(parse(text = (
        paste0("mean(Y1[w==1&x==", x_unique[i], "])")
      ))) - eval(parse(text = (
        paste0("mean(Y1[w==2&x==", x_unique[i], "])")
      )))) # This is c(1,2)
      assign(paste0("c_2_x_", x_unique[i]), eval(parse(text = (
        paste0("mean(Y2[w==2&x==", x_unique[i], "])")
      ))) - eval(parse(text = (
        paste0("mean(Y2[w==1&x==", x_unique[i], "])")
      )))) # This is c(2,1)
      assign(paste0("c_3_x_", x_unique[i]), eval(parse(text = (
        paste0("mean(Y2[w==2&x==", x_unique[i], "])")
      ))) - eval(parse(text = (
        paste0("mean(Y2[w==3&x==", x_unique[i], "])")
      )))) # This is c(2,3)
      assign(paste0("c_4_x_", x_unique[i]), eval(parse(text = (
        paste0("mean(Y1[w==1&x==", x_unique[i], "])")
      ))) - eval(parse(text = (
        paste0("mean(Y1[w==3&x==", x_unique[i], "])")
      )))) # This is c(1,3)
      assign(paste0("c_5_x_", x_unique[i]), eval(parse(text = (
        paste0("mean(Y3[w==3&x==", x_unique[i], "])")
      ))) - eval(parse(text = (
        paste0("mean(Y3[w==1&x==", x_unique[i], "])")
      )))) # This is c(3,1)
      assign(paste0("c_6_x_", x_unique[i]), eval(parse(text = (
        paste0("mean(Y3[w==3&x==", x_unique[i], "])")
      ))) - eval(parse(text = (
        paste0("mean(Y3[w==2&x==", x_unique[i], "])")
      )))) # This is c(3,2)
      assign(paste0("c_x_", x_unique[i]), eval(parse(text = (
        paste0(
          "cbind(c_1_x_",
          x_unique[i],
          ", c_2_x_",
          x_unique[i],
          ", c_3_x_",
          x_unique[i],
          ", c_4_x_",
          x_unique[i],
          ", c_5_x_",
          x_unique[i],
          ", c_6_x_",
          x_unique[i],
          ")"
        )
      ))))
      c_truth_list[[i]] <-
        eval(parse(text = paste0("c_x_", x_unique[i])))
    }
    do.call(rbind, c_truth_list)
  }
