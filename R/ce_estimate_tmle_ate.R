#' Targeted Maximum Likelihood (TMLE) for ATE estimation
#'
#'This function implements the TMLE method when estimand is ATE. Please use our main function ce_estimate.R.
#'
#' This function implements the TMLE method. Please use our main function causal_multi_treat.R.
#' @param y numeric vector for the binary outcome
#'
#' @param w numeric vector for the treatment indicator
#' @param x data frame containing the treatment indicator and covariates
#' @param SL.library a character vector of prediction algorithms. A list of functions included in the SuperLearner package can be found with listWrappers()
#' @param ... Other arguments
#'
#' @importFrom magrittr "%>%"
#' @import SuperLearner
#' @export
#' @return a list with w*(w-1)/2 elements for ATE effect. Each element of the list contains the estimation for RD/RR/OR
#' @examples
#' \donttest{
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
#'ce_estimate_tmle_ate(y = data$y, x = data$covariates ,
#'w = data$w, SL.library = c("SL.glm", "SL.mean"))
#' }

ce_estimate_tmle_ate <- function(y, w, x, SL.library,...){

  n_trt <- length(unique(w))
  Xmat <- cbind(w, x)
  for (i in 1:n_trt){
    Xmat <- Xmat %>%
      dplyr::as_tibble() %>%
      # dplyr::mutate(!!sym(new_col_name) := dplyr::case_when(w == i ~ 1,
      #                                                       TRUE ~ 0))%>%
      dplyr::mutate(dplyr::case_when(w == i ~ 1,
                                                            TRUE ~ 0))%>%
        as.data.frame()
    names(Xmat)[length(Xmat)] <- paste0("w",i)
  }

  K <- n_trt

  n <- dim(Xmat)[1]
  t_mat <- NULL
  for (i in 1:n_trt){
    t_mat_once <- Xmat %>%
      dplyr::select(paste0("w",i))
    t_mat <- dplyr::bind_cols(t_mat, t_mat_once)
  }

  W <- x
  #---------------------------------------------#
  ###Create Counterfactual Treatment Scenarios###
  #---------------------------------------------#
  for (i in 1:n_trt){
    assign(paste0("w", i, "_countfactual"), W)
  }

  for (j in 1:n_trt){
    new_col_name <- paste0("w",j)
    for (i in 1:n_trt){
      if (i == j) {
        names_w_countfactual <- names(eval(parse(text =(paste0("w", i, "_countfactual")))))
        assign(paste0("w", i, "_countfactual"), eval(parse(text = paste0("w", i, "_countfactual"))) %>% dplyr::mutate(1))
        assign(paste0("w", i, "_countfactual"), stats::setNames(eval(parse(text =(paste0("w", i, "_countfactual")))), c(names_w_countfactual, paste0("w", j))))
      } else {
        names_w_countfactual <- names(eval(parse(text =(paste0("w", i, "_countfactual")))))
        assign(paste0("w", i, "_countfactual"), eval(parse(text = paste0("w", i, "_countfactual"))) %>% dplyr::mutate(0))
        assign(paste0("w", i, "_countfactual"), stats::setNames(eval(parse(text =(paste0("w", i, "_countfactual")))), c(names_w_countfactual, paste0("w", j))))
      }
    }
  }

  w_countfactual_combined <- NULL
  for (i in 1:n_trt){
    w_countfactual_combined <- as.data.frame(rbind(w_countfactual_combined, eval(parse(text = paste0("w", i, "_countfactual")))))
  }



  ###Run Super Learner Once, Obtain Initial Predicted Values for All Counterfactual Settings###

  #Step 1: Estimating the outcome regression using super learner

  sl_fit <-
    SuperLearner::SuperLearner(
      Y = y,
      X = Xmat[,-1],
      newX = w_countfactual_combined,
      SL.library = SL.library,
      family = stats::binomial(),
      verbose = FALSE,...
    )

  q_0 <- rep(0, n * K)
  q_tvector <- cbind(q_0, sl_fit$SL.predict)
  q_tmat <-
    matrix(unlist(split(
      as.data.frame(q_tvector), rep(1:K, each = n)
    )), ncol = 2 * K)

  #-----------------------------------------------------------------------------------#
  ###Run TMLE to Calculate Point Estimates of each T=t with Custom Cluster-Based SEs###
  #-----------------------------------------------------------------------------------#
  #Steps 2-5 are performed in this code chunk#

  w_results <- matrix(NA, nrow = K, ncol = 4)
  w_results_row_names <- NULL
  for (i in 1:n_trt){
    w_results_row_names_once <- paste0("w",i)
    w_results_row_names <- c(w_results_row_names, w_results_row_names_once)
  }
  rownames(w_results) <- w_results_row_names

  colnames(w_results) <- c("EYt", "SE", "CI1", "CI2")
  start <- 1
  end <- 2
  gbound <- 0.025
  for(t in 1:K) {
    #Step 2: Super learner fit for P(T_k=t|W) specified with g.SL.library=SL.library in tmle call#
    #Steps 3-4: Performed in tmle call, target  EYt parameter using A=NULL and Delta=Tmat[,t]#
    fit <- tmle::tmle(
      Y = y,
      A = NULL,
      Delta = t_mat[, t],
      W = W,
      Q = q_tmat[, c(start, end)],
      g.SL.library = SL.library,
      family = "binomial",
      verbose = FALSE,...
    )
    #Step 5: The parameter estimates are stored in fit$estimates$EY1$psi#
    w_results[t, 1] <- fit$estimates$EY1$psi
    w_results[t, 2] <- fit$estimates$EY1$var.psi
    w_results[t, 3] <- fit$estimates$EY1$CI[1]
    w_results[t, 4] <- fit$estimates$EY1$CI[2]
    start <- start + 2
    end <- end + 2
  }

  w_result_one_repetition <- w_results %>%
    dplyr::as_tibble(rownames = "treatment")

  for (i in 1:n_trt){
    assign(paste0("mu_",i, "_hat"), w_result_one_repetition %>% dplyr::select("EYt") %>% dplyr::slice(i) %>% dplyr::pull("EYt"))
  }
  for (i in 1:n_trt){
    assign(paste0("lower_",i, "_hat"), w_result_one_repetition %>% dplyr::select("CI1") %>% dplyr::slice(i) %>% dplyr::pull("CI1"))
  }
  for (i in 1:n_trt){
    assign(paste0("upper_",i, "_hat"), w_result_one_repetition %>% dplyr::select("CI2") %>% dplyr::slice(i) %>% dplyr::pull("CI2"))
  }
  result_list <- NULL
  for (i in 1:(n_trt-1)){
    result_once <- NULL
    for (j in (i + 1):n_trt){
      assign(paste0("RD",i,j),eval(parse(text = paste0("mu_",i, "_hat"))) - eval(parse(text = paste0("mu_",j, "_hat"))))
      assign(paste0("RR",i,j),eval(parse(text = paste0("mu_",i, "_hat"))) / eval(parse(text = paste0("mu_",j, "_hat"))))
      assign(paste0("OR",i,j), (eval(parse(text = paste0("mu_",i, "_hat"))) /(1 - eval(parse(text = paste0("mu_",i, "_hat"))))) / (eval(parse(text = paste0("mu_",j, "_hat"))) /(1 - eval(parse(text = paste0("mu_",j, "_hat"))))))
      result_once <- rbind(eval(parse(text = paste0("round(RD",i,j,",2)"))), eval(parse(text = paste0("round(RR",i,j,",2)"))), eval(parse(text = paste0("round(OR",i,j,",2)"))))
      colnames(result_once) <- "EST"
      rownames(result_once) <- c("RD", "RR", "OR")
      result_once_list <- list(result_once)
      names(result_once_list) <- paste0("ATE",i,j)
      result_list <- c(result_list, result_once_list)
    }
  }
  return(result_list)
}
