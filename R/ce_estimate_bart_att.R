#' Causal inference with multiple treatments using BART for ATT effects
#'
#' The function \code{ce_estimate_bart_att} implements
#' BART to estimate ATT effect with
#' multiple treatments using observational data.
#'
#' @param y A numeric vector (0, 1) representing a binary outcome.
#' @param x A dataframe, including all the covariates but not treatments.
#' @param w A numeric vector representing the treatment groups.
#' @param discard A logical indicating whether to use the discarding rules.
#' The default is \code{FALSE}.
#' @param ndpost A numeric value indicating the number of posterior draws.
#' @param reference_trt A numeric value indicating reference treatment group
#' for ATT effect.
#' @param ... Other parameters that can be passed through to functions.
#'
#' @return A summary of the effect estimates can be obtained
#' with \code{summary} function. The output also
#' contains a list of the posterior samples of causal estimands. When
#' \code{discard = TRUE}, the output contains number of discarded
#' individuals.
#' @importFrom BART pbart
#' @importFrom dplyr mutate select filter
#' @references
#'
#' Sparapani R, Spanbauer C, McCulloch R
#' Nonparametric Machine Learning and  Efficient Computation with
#' Bayesian Additive Regression Trees: The BART R Package.
#' \emph{Journal of Statistical Software}, \strong{97}(1), 1-66.
#'
#' Hadley Wickham, Romain François, Lionel Henry and Kirill Müller (2021).
#' \emph{dplyr: A Grammar of Data Manipulation}. R package version 1.0.7.
#' URL: \url{https://CRAN.R-project.org/package=dplyr}
ce_estimate_bart_att <- function(y, x, w, discard = FALSE, ndpost = 1000,
                                 reference_trt, ...) {
  # Get the number of treatment groups
  n_trt <- length(unique(w))
  trt_ind <- seq_len(n_trt)
  trt_ind_no_ref <- trt_ind[trt_ind != reference_trt]
  # Get the sample size for each treatment group
  for (i in seq_len(n_trt)) {
    assign(paste0("n", i), sum(w == i))
  }

  xwdata <- cbind(w, x)

  # Fit BART model
  bart_mod <- BART::pbart(x.train = xwdata, y.train = y, ndpost = ndpost, ...)

  assign(paste0("xp", reference_trt), xwdata[w == reference_trt, ])
  for (i in trt_ind[trt_ind != reference_trt]) {
    assign(paste0("xp", i), xwdata[w == reference_trt, ])
    assign(paste0("xp", i),
           as.data.frame(eval(parse(text = paste0("xp", i)))) %>%
      dplyr::mutate(w = i))
  }
  # Predict potential outcomes among those in reference_trt group
  for (j in seq_len(n_trt)) {
    assign(paste0("bart_pred", reference_trt, j),
           BART::pwbart(eval(parse(text = paste0("xp", j))),
                        bart_mod$treedraws, mu = mean(y)))
    assign(paste0("pred_prop", reference_trt, j),
           stats::pnorm(eval(parse(text =
                                     paste0("bart_pred", reference_trt, j)))))
  }
  # Implement the discarding rules
  if (discard == FALSE) {
    for (j in seq_len(n_trt)) {
      assign(paste0("pred_prop", reference_trt, j),
             stats::pnorm(eval(parse(text =
                                    paste0("bart_pred", reference_trt, j)))))
    }
    n_discard_att <- 0
  } else if (discard == TRUE) {
    for (j in seq_len(n_trt)) {
      assign(paste0("post.ind.sd", j),
             apply(eval(parse(text =
                                paste0("pred_prop", reference_trt, j))), 2,
                   stats::sd))
    }
    threshold <- max(eval(parse(text = paste0("post.ind.sd", reference_trt))))

    for (j in seq_len(length(trt_ind_no_ref))) {
      assign(paste0("eligible", trt_ind_no_ref[j]),
             (eval(parse(text =
                           paste0("post.ind.sd",
                                  trt_ind_no_ref[j]))) <= threshold))
    }
    # Individuals with a large variability in the predicted potential outcomes
    # in the reference group are discarded
    eligible <- rep(TRUE, dim(eval(parse(text =
                                           paste0("pred_prop",
                                                  reference_trt, 1))))[2])
    for (j in seq_len(length(trt_ind_no_ref))) {
      eligible <- (eval(parse(text =
                                paste0("eligible",
                                       trt_ind_no_ref[j]))) &
                     eligible)
    }
    # Record the number of discarded individuals
    # Only use the undiscarded individuals for the final causal estimands
    n_discard_att <- sum(eligible == FALSE)
    for (j in seq_len(n_trt)) {
      assign(paste0("pred_prop", reference_trt, j),
             (eval(parse(text = paste0("pred_prop", reference_trt, j))) %>%
        as.data.frame() %>%
        dplyr::select(which(eligible)) %>%
        as.matrix()))
    }
  }
  # Set up the final matrix to save the RD, RR and OR
  for (j in seq_len(length(trt_ind_no_ref))) {
    assign(paste0("RD",
                  reference_trt, trt_ind_no_ref[j], "_est"), NULL)
    assign(paste0("RR",
                  reference_trt, trt_ind_no_ref[j], "_est"), NULL)
    assign(paste0("OR",
                  reference_trt, trt_ind_no_ref[j], "_est"), NULL)
  }

  for (i in seq_len(n_trt)) {
    assign(paste0("y", i, "_pred"), NULL)
    assign(paste0("y", i, "_pred"), matrix(stats::rbinom(
      dim(eval(parse(text = (
        paste0("pred_prop", reference_trt, i)
      ))))[1] * dim(eval(parse(text = (
        paste0("pred_prop", reference_trt, i)
      ))))[2], 1, eval(parse(text = (
        paste0("pred_prop", reference_trt, i)
      )))
    ), nrow = ndpost))
  }
  # Calculate the final RD, RR and OR using the
  # predictive posterior distributions from BART model
  for (j in seq_len(length(trt_ind_no_ref))) {
    assign(paste0("RD",
                  reference_trt, trt_ind_no_ref[j], "_est"),
           list(rowMeans(eval(parse(text = (
      paste0("y", reference_trt, "_pred")
    )))) - rowMeans(eval(parse(text = (
      paste0("y", trt_ind_no_ref[j], "_pred")
    ))))))
    assign(paste0("RR", reference_trt, trt_ind_no_ref[j], "_est"),
           list(rowMeans(eval(parse(text = (
      paste0("y", reference_trt, "_pred")
    )))) / rowMeans(eval(parse(text = (
      paste0("y", trt_ind_no_ref[j], "_pred")
    ))))))
    assign(paste0("OR", reference_trt, trt_ind_no_ref[j], "_est"),
           list(rowMeans((eval(parse(
      text = (paste0("y", reference_trt, "_pred"))
    ))) / (1 - rowMeans(eval(
      parse(text = (paste0(
        "y", reference_trt, "_pred"
      )))
    )))) / rowMeans((eval(parse(
      text = (paste0("y", trt_ind_no_ref[j], "_pred"))
    ))) / (1 - rowMeans(eval(
      parse(text = (paste0(
        "y", trt_ind_no_ref[j], "_pred"
      )))
    ))))))
  }

  result <- NULL
  for (i in seq_len(n_trt - 1)) {
    for (j in seq_len(length(trt_ind_no_ref))) {
      assign(paste0("RD", reference_trt, trt_ind_no_ref[j], "_est"),
             stats::setNames(eval(
               parse(text =
                       (paste0("RD",
                               reference_trt, trt_ind_no_ref[j], "_est")))),
                             paste0("ATT_RD",
                                    reference_trt, trt_ind_no_ref[j])))
      assign(paste0("RR", reference_trt, trt_ind_no_ref[j], "_est"),
             stats::setNames(eval(
               parse(text = (paste0("RR",
                                    reference_trt,
                                    trt_ind_no_ref[j], "_est")))),
                             paste0("ATT_RR",
                                    reference_trt, trt_ind_no_ref[j])))
      assign(paste0("OR", reference_trt, trt_ind_no_ref[j], "_est"),
             stats::setNames(eval(
               parse(text = (paste0("OR",
                                    reference_trt,
                                    trt_ind_no_ref[j], "_est")))),
                             paste0("ATT_OR",
                                    reference_trt, trt_ind_no_ref[j])))
      result <- c(result,
                  (eval(parse(text =
                                (paste0("RD",
                                        reference_trt,
                                        trt_ind_no_ref[j], "_est"))))),
                  (eval(parse(text =
                                (paste0("RR",
                                        reference_trt,
                                        trt_ind_no_ref[j], "_est"))))),
                  (eval(parse(text =
                                (paste0("OR",
                                        reference_trt,
                                        trt_ind_no_ref[j], "_est"))))))
    }
  }
  result <- c(result, list(n_discard = n_discard_att),
              list(method = parent.frame()$method))
  class(result) <- "CIMTx_ATT_posterior"

  return(result)
}
