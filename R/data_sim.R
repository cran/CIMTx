#' Simulate data for binary outcome with multiple treatments
#'
#' The function \code{data_sim} simulate data for binary outcome with
#' multiple treatments. Users can adjust the following 7 design factors:
#' (1) sample size, (2) ratio of units across treatment groups,
#' (3) whether the treatment assignment model and the outcome generating model
#' are linear or nonlinear, (4) whether the covariates that best predict
#' the treatment also predict the outcome well,
#' (5) whether the response surfaces are parallel across treatment groups,
#' (6) outcome prevalence, and (7) degree of covariate overlap.
#'
#' @param sample_size A numeric value indicating the total number of units.
#' @param n_trt A numeric value indicating the number of treatments.
#' The default is set to 3.
#' @param x A vector of characters representing covariates,
#' with each covariate being generated from the standard probability.
#' The default is set to "rnorm(0, 1)".
#' \code{\link[stats:Distributions]{distributions}} in the
#' \code{\link[stats:stats-package]{stats}} package.
#' @param lp_y A vector of characters of length \code{n_trt},
#' representing the linear effects in the outcome generating model.
#' The default is set to rep("x1", 3).
#' @param nlp_y A vector of characters of length \code{n_trt},
#' representing the nonlinear effects in the outcome generating model.
#' The default is set to NULL.
#' @param align A logical indicating whether the predictors in the
#' treatment assignment model are the same as the predictors for
#' the outcome generating model.
#' The default is \code{TRUE}. If the argument is set to \code{FALSE},
#' users need to specify additional two arguments \code{lp_w} and \code{nlp_w}.
#' @param tau A numeric vector of length \code{n_trt} inducing different
#' outcome event probabilities across treatment groups.
#' Higher values mean higher outcome event probability for the treatment group;
#' lower values mean lower outcome event probability for the treatment group.
#' The default is set to c(0, 0, 0), which corresponds to an
#' approximately equal outcome event probability across three treatment groups.
#' @param delta A numeric vector of length \code{n_trt}-1 inducing different
#' ratio of units across treatment groups.
#' Higher values mean higher proportion for the treatment group;
#' lower values mean lower proportion for the treatment group.
#' The default is set to c(0,0), which corresponds to an approximately
#' equal sample sizes across three treatment groups.
#' @param psi A numeric value for the parameter governing the sparsity
#' of covariate overlap. Higher values mean weaker covariate overlap;
#' lower values mean stronger covariate overlap. The default is set to 1,
#' which corresponds to a moderate covariate overlap.
#' @param lp_w is a vector of characters of length \code{n_trt} - 1,
#' representing in the treatment assignment model
#' @param nlp_w is a vector of characters of length \code{n_trt} - 1,
#' representing in the treatment assignment model
#'
#' @import dplyr
#'
#' @return A list with 7 elements for simulated data. It contains
#' \item{covariates:}{x matrix}
#' \item{w:}{treatment indicators}
#' \item{y:}{observed binary outcomes}
#' \item{y_prev:}{outcome prevalence rates}
#' \item{ratio_of_units:}{the proportions of units in each treatment group}
#' \item{overlap_fig:}{the visualization of covariate overlap
#' via boxplots of the distributions of true GPS}
#' \item{y_true:}{simulated true outcome in each treatment group}
#' @importFrom stringr str_locate str_sub str_locate
#' @importFrom dplyr mutate
#' @export
#'
#' @references
#'
#' Hadley Wickham (2019).
#' \emph{stringr: Simple, Consistent Wrappers for Common String Operations}.
#' R package version 1.4.0.
#' URL:\url{https://CRAN.R-project.org/package=stringr}
#'
#' Hadley Wickham, Romain François,
#' Lionel Henry and Kirill Müller (2021).
#' \emph{dplyr: A Grammar of Data Manipulation}.
#' R package version 1.0.7.
#' URL: \url{https://CRAN.R-project.org/package=dplyr}
#' @examples
#' library(CIMTx)
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
#'   "rbeta(2,0.4)", # x2
#'   "runif(0, 0.5)", # x3
#'   "rweibull(1,2)", # x4
#'   "rbinom(1,0.4)" # x5
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
data_sim <-
  function(sample_size,
           n_trt = 3,
           x = "rnorm(0, 1)",
           lp_y = rep("x1", 3),
           nlp_y = NULL,
           align = TRUE,
           tau = c(0, 0, 0),
           delta = c(0, 0),
           psi = 1,
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

    for (i in seq_len(length(x))) {
      str_locate_parenthesis <- stringr::str_locate(x[i], "\\(")
      assign(paste0("x", i), eval(parse(text = paste0(
        paste0(
          stringr::str_sub(x[i], 1, str_locate_parenthesis[1]),
          sample_size,
          ",",
          stringr::str_sub(x[i], str_locate_parenthesis[1] + 1)
        )
      ))))
    }
    x_matrix <- matrix(NA, nrow = sample_size, ncol = length(x))
    for (i in 1:dim(x_matrix)[2]) {
      x_matrix[, i] <- eval(parse(text = paste0("x", i)))
    }
    treatment_exp_matrix_all <- NULL
    for (i in 1:(n_trt - 1)) {
      treatment_exp_matrix <-
        exp(eval(parse(
          text = paste0(psi, "*(", delta[i], "+", lp_w[i], "+", nlp_w[i], ")")
        )))
      treatment_exp_matrix_all <-
        cbind(treatment_exp_matrix_all, treatment_exp_matrix)
    }
    treatment_exp_matrix_all[is.infinite(treatment_exp_matrix_all)] <-
      10 ^ 8
    probs <-
      sweep(treatment_exp_matrix_all, 1,
            (rowSums(treatment_exp_matrix_all) + 1), "/")

    prob_last_all <- rep(NA, dim(probs)[1])
    for (i in 1:dim(probs)[1]) {
      prob_last_all[i] <- 1 - sum(probs[i, ])
    }
    probs_all <- cbind(probs, prob_last_all)
    w <- rep(NA, dim(probs)[1])
    for (i in 1:dim(probs)[1]) {
      w[i] <-
        sample(
          x = 1:n_trt,
          size = 1,
          prob = c(probs[i, ], 1 - sum(probs[i, ])),
          replace = TRUE
        )
    }
    w


    for (i in 1:n_trt) {
      assign(paste0("Y", i, "_final"), NULL)
    }

    for (j in 1:n_trt) {
      assign(paste0("y", j), stats::plogis(eval(parse(
        text = paste0(tau[j], "+", lp_y[j], "+", nlp_y[j])
      ))))
    }

    for (j in 1:n_trt) {
      assign(paste0("Y", j, "_final"),
             stats::rbinom(
               n = sample_size,
               size = 1,
               prob = eval(parse(text = paste0("y", j)))
             ))
    }

    y_true <- matrix(NA, nrow = sample_size, ncol = n_trt)
    for (i in 1:n_trt) {
      y_true[, i] <- eval(parse(text = paste0("Y", i, "_final")))
    }
    y_true_with_treatment <- cbind(y_true, w)

    # observed outcomes
    y_obs <-
      apply(y_true_with_treatment, 1, function(x)
        x[1:n_trt][x[n_trt + 1]])

    y_prev <- tibble(w = as.character(w), y_obs) %>%
      group_by(w) %>%
      summarise(y_prev = mean(y_obs)) %>%
      bind_rows(tibble(w = "Overall",
                       y_prev = mean(y_obs))) %>%
      dplyr::mutate(y_prev = round(y_prev, 2)) %>%
      as.data.frame()


    p_covariate_overlap <- covariate_overlap(treatment = w,
                                             prob = probs_all)
    covariates <- as.data.frame(x_matrix)
    names(covariates) <- paste0("x", seq_len(length(x)))
    return(
      list(
        covariates = covariates,
        w = w,
        y = y_obs,
        y_prev = y_prev,
        ratio_of_units = round(table(w) / length(w), 2),
        overlap_fig = p_covariate_overlap,
        y_true = y_true
      )
    )
  }
