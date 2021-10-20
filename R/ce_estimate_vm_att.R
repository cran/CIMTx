#' Vector Macthing (VM) for ATT estimation
#'
#'This function implements the VM method when estimand is ATT. Please use our main function ce_estimate.R.
#'
#' @param y numeric vector for the binary outcome
#' @param x dataframe including the treatment indicator and the covariates
#' @param w numeric vector for the treatment indicator
#' @param caliper is a numeric value denoting the caliper which should be used when matching on the logit of GPS within each cluster formed by K-means clustering. The caliper is in standardized units. For example, caliper = 0.25 means that all matches greater than 0.25 standard deviations of the logit of GPS are dropped. The default value is 0.25
#' @param n_cluster a numeric value denoting the number of clusters to form using K means clustering on the logit of GPS. The default value is 5
#' @param reference_trt reference_trt group for ATT
#'
#' @return a list with w-1 elements for ATT effect. Each element of the list contains the estimation, standard error, lower and upper 95\% CI for RD/RR/OR
#' @export
#' @examples
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
#'ce_estimate_vm_att(y = data$y, x = data$covariates,
#'w = data$w,reference_trt = 1, caliper = 0.25, n_cluster = 5)


ce_estimate_vm_att <-  function(y, x, w, reference_trt,caliper,n_cluster) {
  n_trt <- length(unique(w))
  if (n_trt > 3){
    stop("We do not recommend using VM for more than 3 treatments")
  }
  # estimate generalized propensity scores using multinomial logistic regression
  xwdata = cbind(w, x)
  # WARNING: reorder the dataset by "Treatment 1", "Treatment 2", "Treatment 3"
  # it will affect the identification of matched indices in row 249
  xwdata = xwdata[order(xwdata$w),]
  # ps model 1
  ps_fit = nnet::multinom(as.factor(w) ~ ., data = xwdata, trace = FALSE)
  probs_logit1 = data.frame(stats::fitted(ps_fit))
  colnames_probs_logit1 <- NULL
  for (i in 1:n_trt){
    colnames_probs_logit1_once <- paste0("p",i)
    colnames_probs_logit1 <- c(colnames_probs_logit1, colnames_probs_logit1_once)
  }
  colnames(probs_logit1) = colnames_probs_logit1
  xwdata = cbind(xwdata, probs_logit1)

  # Determine eligibility
  min_max_Ps <- NULL
  for (i in 1: n_trt){
    xwdata_summarise_once <- xwdata %>%
      dplyr::group_by(w) %>%
      dplyr::summarise(min(eval(parse(text = paste0("p",i)))), max(eval(parse(text = paste0("p",i))))) %>%
      dplyr::ungroup()
    names(xwdata_summarise_once)[c(2,3)] <- c(paste0("min",i), paste0("max",i))
    if (i == 1){
      min_max_Ps <- xwdata_summarise_once
    } else {
      min_max_Ps <- min_max_Ps %>% dplyr::inner_join(xwdata_summarise_once, by = "w")
    }
  }

  #min_max_Ps
  Eligible <- TRUE
  for (i in 1:n_trt){
    Eligible_once <- xwdata[[paste0("p", i)]] >= max(min_max_Ps[[paste0("min", i)]]) & xwdata[[paste0("p", i)]] <= min(min_max_Ps[[paste0("max", i)]])
    Eligible <- Eligible & Eligible_once
  }
  xwdata <- xwdata %>%
    dplyr::mutate(Eligible = Eligible)

  xwydata = cbind(y, xwdata)
  xwydata = dplyr::filter(xwydata, Eligible)

  # Calculate new propensity scores for eligible subjects
  ps_fit_E = nnet::multinom(w ~ ., data = xwydata[,-1], trace = FALSE)
  probs_logit1_E = stats::fitted(ps_fit_E)

  colnames_probs_logit1_E <- NULL
  for (i in 1:n_trt){
    colnames_probs_logit1_E_once <- paste0("p",i)
    colnames_probs_logit1_E <- c(colnames_probs_logit1_E, colnames_probs_logit1_E_once)
  }
  colnames(probs_logit1_E) = colnames_probs_logit1_E
  # colnames(probs_logit1_E) = c("p1", "p2", "p3")
  xwydata <- xwydata %>%
    dplyr :: select(-"p1", -"p2", -"p3")
  xwydata = cbind(xwydata, probs_logit1_E)

  for (i in 1:n_trt){
    assign(paste0("n",i), sum(xwydata$w == i))
  }

  ### Vector Matching for ATT for outcome 1 (comp_resp_obs)
  # Stratify logit(r(ti, Xi)) using K-means clustering
  # n_cluster <- 5

  for (i in 1:n_trt){
    temp_once <-  stats::kmeans(logit(xwydata[[paste0("p",i)]]), n_cluster)
    xwydata <- xwydata %>%
      dplyr::mutate(temp_once$cluster)
    names(xwydata)[length(names(xwydata))] <- paste0("Quint",i)
  }
  treat <- NULL
  colnames(xwydata)[1:2] <-  c("Y_obs","treat")
  trt_indicator = 1:n_trt
  trt_indicator_no_reference <- trt_indicator[trt_indicator!=reference_trt]

  for (i in 1:length((trt_indicator_no_reference))){
    assign(paste0("temp", reference_trt,trt_indicator_no_reference[i]), dplyr::filter(xwydata, treat %in% c(reference_trt,trt_indicator_no_reference[i])))
  }
  assign(paste0("temp", trt_indicator_no_reference[1],trt_indicator_no_reference[2]), dplyr::filter(xwydata, treat %in% trt_indicator_no_reference))

  # matching
  for (i in 1:length((trt_indicator_no_reference))){
    trt_indicator_left <- trt_indicator[!(trt_indicator %in% c(reference_trt, trt_indicator_no_reference[i]))]
    assign(paste0("match",reference_trt,trt_indicator_no_reference[i]),   Matching::Matchby(Y = eval(parse(text = paste0("temp", reference_trt,trt_indicator_no_reference[i])))[["Y_obs"]], Tr = eval(parse(text = paste0("temp", reference_trt,trt_indicator_no_reference[i])))[["treat"]] == reference_trt,X = logit(eval(parse(text = paste0("temp", reference_trt,trt_indicator_no_reference[i])))[[paste0("p",reference_trt)]]), by = eval(parse(text = paste0("temp", reference_trt,trt_indicator_no_reference[i])))[[paste0("Quint", trt_indicator_left)]],
  caliper = caliper,
  replace = T, estimand = "ATT", print.level = 0))
  }

  # Identify the matched subgroups
  rownames(xwydata) = 1:nrow(xwydata)
  xwydata$id = 1:nrow(xwydata)
  eligible_matching <- TRUE
  for (i in 1:length((trt_indicator_no_reference))){
    eligible_matching_once <- xwydata$id %in% eval(parse(text = paste0("match",reference_trt,trt_indicator_no_reference[i])))[["index.treated"]]
    eligible_matching <- eligible_matching & eligible_matching_once
  }
  if (sum(eligible_matching) == 0){
    stop("No eligible matches could be found.")
  }
  xwydata$both_1 <- eligible_matching
  temp = xwydata[xwydata$both_1 == "TRUE", ]


  for (i in 1:length((trt_indicator_no_reference))){
    if ( i == 1){
      assign(paste0("m",reference_trt,trt_indicator_no_reference[i]), cbind( eval(parse(text = paste0("match",reference_trt,trt_indicator_no_reference[i])))[["index.treated"]], eval(parse(text = paste0("match",reference_trt,trt_indicator_no_reference[i])))[["index.control"]]))
    } else {
      trt_indicator_left <- trt_indicator[!(trt_indicator %in% c(reference_trt, trt_indicator_no_reference[i]))]
      assign(paste0("m",reference_trt,trt_indicator_no_reference[i]), cbind( eval(parse(text = paste0("match",reference_trt,trt_indicator_no_reference[i])))[["index.treated"]], eval(parse(text = paste0("match",reference_trt,trt_indicator_no_reference[i])))[["index.control"]] + sum(xwydata[["treat"]] == trt_indicator_left)))
    }
  }
  # WARNING: reorder the dataset by "Treatment 1", "Treatment 2", "Treatment 3" before running the following lines
  for (i in 1:length((trt_indicator_no_reference))){
    assign(paste0("m",reference_trt,trt_indicator_no_reference[i]), eval(parse(text = paste0("m",reference_trt,trt_indicator_no_reference[i])))[eval(parse(text = paste0("m",reference_trt,trt_indicator_no_reference[i])))[,1] %in% rownames(temp),])
  }
  triplets <- NULL
  for (i in 1:length(trt_indicator_no_reference)){
    triplets_once <- eval(parse(text = paste0("m",reference_trt,trt_indicator_no_reference[i])))[order(eval(parse(text = paste0("m",reference_trt,trt_indicator_no_reference[i])))[,1]), ]
    triplets <- cbind(triplets,triplets_once)
  }
  triplets = as.matrix(triplets[,c(1, 2, 4)])
  n_trip = nrow(triplets)
  # Matching Estimator
  # For subjects receiving reference_trt treatment
  for (i in 1:3){
    assign(paste0("Y",i,"_imp"), xwydata$Y_obs[triplets[,i]])
    assign(paste0("y",i,"_hat"), mean(eval(parse(text = paste0("Y",i,"_imp")))))
  }

  for (i in 1:length((trt_indicator_no_reference))){
    assign(paste0("RD", reference_trt, trt_indicator_no_reference[i],"_est"), eval(parse(text = paste0("y",reference_trt,"_hat"))) - eval(parse(text = paste0("y",trt_indicator_no_reference[i],"_hat"))))
    assign(paste0("RR", reference_trt, trt_indicator_no_reference[i],"_est"), eval(parse(text = paste0("y",reference_trt,"_hat"))) / eval(parse(text = paste0("y",trt_indicator_no_reference[i],"_hat"))))
    assign(paste0("OR", reference_trt, trt_indicator_no_reference[i],"_est"), (eval(parse(text = paste0("y",reference_trt,"_hat"))) / ( 1- eval(parse(text = paste0("y",reference_trt,"_hat"))))) / (eval(parse(text = paste0("y",trt_indicator_no_reference[i],"_hat"))) / ( 1- eval(parse(text = paste0("y",trt_indicator_no_reference[i],"_hat"))))))
  }

  result <- NULL
  result_list <- NULL
  for (j in 1:length(trt_indicator_no_reference)){
    result_once <- rbind(eval(parse(text = paste0("round(RD",reference_trt,trt_indicator_no_reference[j],"_est,2)"))), eval(parse(text = paste0("round(RR",reference_trt,trt_indicator_no_reference[j],"_est,2)"))), eval(parse(text = paste0("round(OR",reference_trt,trt_indicator_no_reference[j],"_est,2)"))))
    colnames(result_once) <- "EST"
    rownames(result_once) <- c("RD", "RR", "OR")
    result_once_list <- list(result_once)
    names(result_once_list) <- paste0("ATT",reference_trt,trt_indicator_no_reference[j])
    result_list <- c(result_list, result_once_list)
  }
  return(c(result_list,number_matched = n_trip))
}


