#' Causal inference with multiple treatments using RA for ATE effects
#'
#' The function \code{ce_estimate_ra_ate} implements
#' RA to estimate ATE effect with
#' multiple treatments using observational data.
#'
#' @param y A numeric vector (0, 1) representing a binary outcome.
#' @param x A dataframe, including all the covariates but not treatments.
#' @param w A numeric vector representing the treatment groups.
#' @param ndpost A numeric value indicating the number of posterior draws.
#'
#' @return A summary of the effect estimates can be obtained
#' with \code{summary} function. The output also
#' contains a list of the posterior samples of causal estimands.
#' @importFrom arm bayesglm sim invlogit
#' @importFrom dplyr mutate
#' @references
#'
#' Hadley Wickham, Romain François, Lionel Henry and Kirill Müller (2021).
#' \emph{dplyr: A Grammar of Data Manipulation}.
#' R package version 1.0.7.
#' URL: \url{https://CRAN.R-project.org/package=dplyr}
#'
#' Andrew Gelman and Yu-Sung Su (2020).
#' \emph{arm: Data Analysis Using Regression and
#' Multilevel/Hierarchical Models}.
#' R package version 1.11-2.
#' \emph{https://CRAN.R-project.org/package=arm}
ce_estimate_ra_ate <- function(y, x, w, ndpost) {
  n_trt <- length(unique(w)) # Number of unique treatments
  for (i in 1:n_trt) {
    # Number of individuals receiving each treatment
    assign(paste0("n", i), sum(w == i))
  }
  n <- length(w)
  xwdata <- cbind(w, x)
  xwydata <- cbind(y, xwdata)
  # Fit Bayesian logistic regression
  reg_mod <-
    arm::bayesglm(
      y ~ .,
      data = as.data.frame(xwydata),
      family = stats::binomial(link = "logit"),
      x = TRUE
    )
  mod_sims <- arm::sim(reg_mod, n.sims = ndpost)
  sim_beta <- as.matrix(stats::coef(mod_sims))
  x_tilde <- stats::model.matrix(reg_mod)
  for (i in 1:n_trt) {
    assign(paste0("x_tilde", i),
           as.data.frame(x_tilde) %>%
             dplyr::mutate(w = i))
  }
  # Predict potential outcomes for each treatment group
  for (i in 1:n_trt) {
    assign(paste0("p", i, "_tilde"),
           arm::invlogit(as.matrix(sim_beta %*% t(eval(
      parse(text = paste0("x_tilde", i))
    )))))
    assign(paste0("y", i, "_tilde"),
           matrix(stats::rbinom(ndpost * n, 1, eval(
             parse(text = paste0("p", i, "_tilde"))
           )), nrow = ndpost))
  }
  for (i in 1:(n_trt - 1)) {
    for (j in (i + 1):(n_trt)) {
      assign(paste0("RD", i, j, "_est"), NULL)
      assign(paste0("RR", i, j, "_est"), NULL)
      assign(paste0("OR", i, j, "_est"), NULL)
    }
  }

  for (i in 1:n_trt) {
    assign(paste0("y", i, "_pred"), rowMeans(eval(parse(
      text = paste0("y", i, "_tilde")
    ))))
  }
  # Estimate causal effects in terms of OR, RR and RD
  result <- NULL
  for (i in 1:(n_trt - 1)) {
    for (j in (i + 1):n_trt) {
      assign(paste0("RD", i, j, "_est"), list(eval(parse(
        text = (paste0("y", i, "_pred"))
      )) - eval(parse(
        text = (paste0("y", j, "_pred"))
      ))))
      assign(paste0("RR", i, j, "_est"), list(eval(parse(
        text = (paste0("y", i, "_pred"))
      )) / eval(parse(
        text = (paste0("y", j, "_pred"))
      ))))
      assign(paste0("OR", i, j, "_est"), list((eval(
        parse(text = (paste0(
          "y", i, "_pred"
        )))
      ) / (
        1 - eval(parse(text = (
          paste0("y", i, "_pred")
        )))
      )) / (eval(
        parse(text = (paste0(
          "y", j, "_pred"
        )))
      ) / (
        1 - eval(parse(text = (
          paste0("y", j, "_pred")
        )))
      ))))
      assign(paste0("RD", i, j, "_est"),
             stats::setNames(eval(parse(
               text = (paste0("RD", i, j, "_est"))
             )), paste0("ATE_RD", i, j)))
      assign(paste0("RR", i, j, "_est"),
             stats::setNames(eval(parse(
               text = (paste0("RR", i, j, "_est"))
             )), paste0("ATE_RR", i, j)))
      assign(paste0("OR", i, j, "_est"),
             stats::setNames(eval(parse(
               text = (paste0("OR", i, j, "_est"))
             )), paste0("ATE_OR", i, j)))
      result <-
        c(result, (eval(parse(
          text = (paste0("RD", i, j, "_est"))
        ))), (eval(parse(
          text = (paste0("RR", i, j, "_est"))
        ))), (eval(parse(
          text = (paste0("OR", i, j, "_est"))
        ))))
    }
  }
  result <- c(result, list(method = parent.frame()$method))
  class(result) <- "CIMTx_ATE_posterior"
  return(result)
}
