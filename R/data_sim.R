#' Simulate data for binary outcome with multiple treatments
#'
#' This function simulate data for binary outcome with multiple treatments. Users can adjust the follwoing 7 design factors: (1) sample size, (2) ratio of units across treatment groups, (3) whether the treatment assignment model and the outcome generating model are linear or nonlinear, (4) whether the covariates that best predict the treatment also predict the outcome well, (5) whether the response surfaces are parallel across treatment groups, (6) outcome prevalence, and (7) degree of covariate overlap.
#'
#' @param sample_size total number of units.
#' @param n_trt the number of treatments
#' @param X a vector of characters representing covariates, with each covariate being generated from the standard probability distributions in the stats package
#' @param lp_y a vector of characters of length n_trt, representing the linear effects in the outcome generating model
#' @param nlp_y a vector of characters of length n_trt, representing the nonlinear effects in the outcome generating model
#' @param align logical,indicating whether the predictors in the treatment assignment model are the same as the predictors for the outcome generating model. The default is TRUE. If the argument is set to FALSE, users need to specify additional two arguments lp_w and nlp_w.
#' @param tau a numeric vector of length n_trt inducing different outcome event probabilities across treatment groups
#' @param delta a numeric vector of length n_trt-1 inducing different ratio of units across treatment groups.
#' @param psi a numeric value for the parameter psi in the treatment assignment model, governing the sparsity of covariate overlap.
#' @param lp_w is a vector of characters of length n_trt - 1, representing in the treatment assignment model
#' @param nlp_w is a vector of characters of length n_trt - 1, representing in the treatment assignment model
#'
#' @import dplyr
#'
#' @return list with 7 elements for simulated data. It contains
#' \item{covariates:}{X matrix}
#' \item{w:}{treatment indicators}
#' \item{y:}{observed binary outcomes}
#' \item{y_prev:}{outcome prevalence rates}
#' \item{ratio_of_units:}{the proportions of units in each treatment group}
#' \item{overlap_fig:}{the visualization of covariate overlap via boxplots of the distributions of true GPS}
#' \item{Y_true_matrix:}{simulated true outcome in each treatment group}
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
data_sim <-
  function(sample_size,
           n_trt,
           X,
           lp_y,
           nlp_y,
           align = TRUE,
           tau,
           delta,
           psi,
           lp_w,
           nlp_w) {

    if (align == TRUE) {
      lp_w <- lp_y
      nlp_w <- nlp_y
    }
    if (is.null(lp_y)) {
      lp_y <- rep(0, n_trt)
    }
    if (is.null(nlp_y)) {
      nlp_y <- rep(0, n_trt)
    }
    if (is.null(lp_w)) {
      lp_w <- rep(0, n_trt)
    }
    if (is.null(nlp_w)) {
      nlp_w <- rep(0, n_trt)
    }

  for (i in 1:length(X)){
    assign(paste0("x",i), eval(parse(text = X[i])))
  }
  X_matrix <- matrix(NA, nrow = sample_size, ncol = length(X))
  for (i in 1:dim(X_matrix)[2]){
    X_matrix[,i] <- eval(parse(text = paste0("x",i)))
  }
  treatment_exp_matrix_all <- NULL
  for (i in 1:(n_trt-1)){
    treatment_exp_matrix <- exp(eval(parse(text = paste0( psi, "*(",delta[i], "+", lp_w[i], "+", nlp_w[i], ")"))))
    treatment_exp_matrix_all <- cbind(treatment_exp_matrix_all, treatment_exp_matrix)
  }
  treatment_exp_matrix_all[is.infinite(treatment_exp_matrix_all)] <- 10^8
  probs <- sweep(treatment_exp_matrix_all, 1, (rowSums(treatment_exp_matrix_all)+1), '/')

  prob_last_all <- rep(NA, dim(probs)[1])
  for (i in 1:dim(probs)[1]){
    prob_last_all[i] <- 1-sum(probs[i,])
  }
  probs_all <- cbind(probs, prob_last_all)
  w <- rep(NA, dim(probs)[1])
  for (i in 1:dim(probs)[1]){
    w[i] <- sample(x = 1:n_trt, size = 1, prob = c(probs[i,], 1-sum(probs[i,])), replace = TRUE)
  }
  w
  expit = function(x) {exp(x)/(1+exp(x))}

  for (i in 1:n_trt){
    assign(paste0("Y",i,"_final"),NULL)
  }

  for (j in 1:n_trt){
    assign(paste0("y",j),  expit(eval(parse(text = paste0(tau[j], "+", lp_y[j], "+", nlp_y[j])))))
  }

  for (j in 1:n_trt){
    assign(paste0("Y",j,"_final"), stats::rbinom(n = sample_size, size = 1, prob = eval(parse(text = paste0("y",j)))))
  }

  Y_true_matrix <- matrix(NA, nrow = sample_size, ncol = n_trt)
  for (i in 1:n_trt){
    Y_true_matrix[,i] <- eval(parse(text = paste0("Y",i,"_final")))
  }
  Y_true_matrix_with_treatment <- cbind(Y_true_matrix, w)

  #observed outcomes
  Yobs = apply(Y_true_matrix_with_treatment, 1, function(x) x[1:n_trt][x[n_trt+1]])

  y_prev <- tibble(w = as.character(w), Yobs) %>%
    group_by(w) %>%
    summarise(y_prev = mean(Yobs)) %>%
    bind_rows(tibble(w = "Overall",
                     y_prev = mean(Yobs))) %>%
    dplyr::mutate(y_prev = round(y_prev, 2)) %>%
    as.data.frame()


    p_covariate_overlap <- covariate_overlap(treatment = w,
                            prob = probs_all)
  return(
    list(
      covariates = as.data.frame(X_matrix),
      w = w,
      y = Yobs,
      y_prev = y_prev,
      ratio_of_units = round(table(w) / length(w),2),
      overlap_fig = p_covariate_overlap,
      Y_true_matrix = Y_true_matrix
    )
  )
}
