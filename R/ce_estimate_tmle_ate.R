

ce_estimate_tmle_ate <- function(y, w, x, SL.library,...){
  if (any((SL.library %in% getNamespaceExports("SuperLearner")[grepl(pattern = "^[S]L", getNamespaceExports("SuperLearner"))]) == F)) stop("SL.library argument unrecgonized; please use listWrappers() in SuperLearner to find the list of supported values", call. = FALSE)
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
  result_final <- NULL
  for (i in 1:(n_trt-1)){
    result_once <- NULL
    for (j in (i + 1):n_trt){
      assign(paste0("RD",i,j),eval(parse(text = paste0("mu_",i, "_hat"))) - eval(parse(text = paste0("mu_",j, "_hat"))))
      assign(paste0("RR",i,j),eval(parse(text = paste0("mu_",i, "_hat"))) / eval(parse(text = paste0("mu_",j, "_hat"))))
      assign(paste0("OR",i,j), (eval(parse(text = paste0("mu_",i, "_hat"))) /(1 - eval(parse(text = paste0("mu_",i, "_hat"))))) / (eval(parse(text = paste0("mu_",j, "_hat"))) /(1 - eval(parse(text = paste0("mu_",j, "_hat"))))))
      result_once <- rbind(eval(parse(text = paste0("round(RD",i,j,",2)"))), eval(parse(text = paste0("round(RR",i,j,",2)"))), eval(parse(text = paste0("round(OR",i,j,",2)"))))
      colnames(result_once) <- "EST"
      rownames(result_once) <- c("RD", "RR", "OR")
      # result_once_list <- list(result_once)
      colnames(result_once) <- paste0("ATE",i,j)
      result_final <- cbind(result_final, result_once)
    }
  }
  return(result_final)
}
