#' Causal inference with multiple treatments using BART for ATE effects
#'
#' The function \code{ce_estimate_bart_ate} implements
#' BART to estimate ATE effect with
#' multiple treatments using observational data.
#'
#' @param y A numeric vector (0, 1) representing a binary outcome.
#' @param x A dataframe, including all the covariates but not treatments.
#' @param w A numeric vector representing the treatment groups.
#' @param discard A logical indicating whether to use the discarding rules.
#' The default is \code{FALSE}.
#' @param ndpost A numeric value indicating the number of posterior draws.
#' @param ... Other parameters that can be passed through to functions.
#' @return A summary of the effect estimates can be obtained
#' with \code{summary} function. The output also
#' contains a list of the posterior samples of causal estimands. When
#' \code{discard = TRUE}, the output contains number of discarded
#' individuals.
#' @importFrom BART pbart
#' @importFrom dplyr mutate select filter
#' @references
#' Sparapani R, Spanbauer C, McCulloch R
#' Nonparametric Machine Learning and
#' Efficient Computation with Bayesian Additive Regression Trees:
#' The BART R Package.
#' \emph{Journal of Statistical Software}, \strong{97}(1), 1-66.
#'
#' Hadley Wickham, Romain François, Lionel Henry and Kirill Müller (2021).
#' \emph{dplyr: A Grammar of Data Manipulation}. R package version 1.0.7.
#' URL: \url{https://CRAN.R-project.org/package=dplyr}
ce_estimate_bart_ate <- function(y, x, w, discard = FALSE,
                                 ndpost = 1000, ...) {
  # Get the number of treatment groups
  n_trt <- length(unique(w))
  # Get the sample size for each treatment group
  for (i in 1:n_trt) {
    assign(paste0("n", i), sum(w == i))
  }
  xwdata <- cbind(w, x)
  # Fit BART model
  bart_mod <- BART::pbart(x.train = xwdata, y.train = y,
                          ndpost = ndpost, ...)
  # Predict potential outcomes for each treatment group w
  for (i in 1:n_trt) {
    assign(paste0("xp", i), xwdata)
    for (j in 1:(n_trt)) {
      assign(paste0("xp", j),
             as.data.frame(eval(parse(text = paste0("xp", i)))) %>%
        dplyr::mutate(w = j))
      assign(paste0("xp", i, j),
             as.data.frame(eval(parse(text = paste0("xp", i)))) %>%
        dplyr::filter(w == i) %>%
        dplyr::mutate(w = j))
      assign(paste0("bart_pred", i, j),
             stats::predict(bart_mod,
                            newdata = eval(parse(text = paste0("xp", i, j)))))
      assign(paste0("bart_pred", j),
             stats::predict(bart_mod,
                            newdata = eval(parse(text = paste0("xp", j)))))
      assign(paste0("pred_prop", i, j),
             eval(parse(text = paste0("bart_pred", i, j)))[["prob.test"]])
      assign(paste0("pred_prop", j),
             eval(parse(text = paste0("bart_pred", j)))[["prob.test"]])
    }
  }
  # Implement the discarding rules
  if (discard == FALSE) {
    for (i in 1:n_trt) {
      for (j in 1:(n_trt)) {
        assign(paste0("pred_prop", i, j),
               eval(parse(text = paste0("bart_pred", i, j)))[["prob.test"]])
      }
    }
    discard_all <- 0
  } else if (discard == TRUE) {
    for (i in 1:n_trt) {
      for (j in 1:(n_trt)) {
        assign(paste0("post.ind.sd", i, j),
               apply(eval(parse(text =
                                  paste0("pred_prop", i, j))), 2, stats::sd))
      }
      assign(paste0("threshold", i),
             max(eval(parse(text = paste0("post.ind.sd", i, i)))))
    }
    # Individuals with a large variability in the predicted potential outcomes
    # among each of the treatment groups are discarded
    for (i in 1:n_trt) {
      n_trt_no_i <- unique(w)[unique(w) != i]
      assign(paste0("eligible", i), TRUE)
      assign(paste0("criteria", i), TRUE)
      for (j in seq_len(n_trt_no_i)) {
        assign(paste0("criteria", i, n_trt_no_i[j]),
               eval(parse(text = paste0("post.ind.sd", i, n_trt_no_i[j]))) <=
                 eval(parse(text = paste0("threshold", i))))
        assign(paste0("eligible", i),
               eval(parse(text = paste0("criteria", i, n_trt_no_i[j]))) &
                 eval(parse(text = paste0("criteria", i))))
      }
      # Record the number of discarded individuals in each treatment group
      assign(paste0("n_", i, "_discard"),
             sum(eval(parse(text = paste0("eligible", i))) == FALSE))
    }
    discard_all <- NULL
    for (i in 1:n_trt) {
      discard_all <- c(discard_all,
                       eval(parse(text = paste0("n_", i, "_discard"))))
    }
    # Only use the undiscarded individuals for the final causal estimands
    for (i in 1:n_trt) {
      for (j in 1:(n_trt)) {
        assign(paste0("pred_prop", i, j),
               eval(parse(text = paste0("pred_prop", i, j))) %>%
          as.data.frame() %>%
          dplyr::select(which(eval(parse(text = paste0("eligible", i))))) %>%
          as.matrix())
      }
    }
  }
  # Set up the final matrix to save the RD, RR and OR
  for (i in 1:(n_trt - 1)) {
    for (j in (i + 1):n_trt) {
      assign(paste0("RD", i, j, "_est"), NULL)
      assign(paste0("RR", i, j, "_est"), NULL)
      assign(paste0("OR", i, j, "_est"), NULL)
    }
  }

  for (i in 1:n_trt) {
    assign(paste0("y", i, "_pred"), NULL)
    assign(paste0("y", i, "_pred"), matrix(stats::rbinom(
      dim(eval(parse(text = (
        paste0("pred_prop", i)
      ))))[1] * dim(eval(parse(text = (
        paste0("pred_prop", i)
      ))))[2], 1, eval(parse(text = (
        paste0("pred_prop", i)
      )))
    ), nrow = ndpost))
  }
  # Calculate the final RD, RR and OR using
  # the predictive posterior distributions from BART model
  for (i in 1:(n_trt - 1)) {
    for (j in (i + 1):n_trt) {
      assign(paste0("RD", i, j, "_est"), list(rowMeans(eval(parse(text = (
        paste0("y", i, "_pred")
      )))) - rowMeans(eval(parse(text = (
        paste0("y", j, "_pred")
      ))))))
      assign(paste0("RR", i, j, "_est"), list(rowMeans(eval(parse(text = (
        paste0("y", i, "_pred")
      )))) / rowMeans(eval(parse(text = (
        paste0("y", j, "_pred")
      ))))))
      assign(paste0("OR", i, j, "_est"), list(rowMeans((eval(parse(
        text = (paste0("y", i, "_pred"))
      ))) / (1 - rowMeans(eval(
        parse(text = (paste0(
          "y", i, "_pred"
        )))
      )))) / rowMeans((eval(parse(
        text = (paste0("y", j, "_pred"))
      ))) / (1 - rowMeans(eval(
        parse(text = (paste0(
          "y", j, "_pred"
        )))
      ))))))
    }
  }
  result <- NULL
  for (i in 1:(n_trt - 1)) {
    for (j in (i + 1):n_trt) {
      assign(paste0("RD", i, j, "_est"),
             stats::setNames(eval(parse(text = (paste0("RD", i, j, "_est")))),
                             paste0("ATE_RD", i, j)))
      assign(paste0("RR", i, j, "_est"),
             stats::setNames(eval(parse(text = (paste0("RR", i, j, "_est")))),
                             paste0("ATE_RR", i, j)))
      assign(paste0("OR", i, j, "_est"),
             stats::setNames(eval(parse(text = (paste0("OR", i, j, "_est")))),
                             paste0("ATE_OR", i, j)))
      result <- c(result,
                  (eval(parse(text = (paste0("RD", i, j, "_est"))))),
                  (eval(parse(text = (paste0("RR", i, j, "_est"))))),
                  (eval(parse(text = (paste0("OR", i, j, "_est"))))))
    }
  }
  result <- c(result, list(n_discard = discard_all),
              list(method = parent.frame()$method))
  class(result) <- "CIMTx_ATE_posterior"

  return(result)
}
