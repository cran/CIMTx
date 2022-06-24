#' Causal inference with multiple treatments using IPTW for ATT effects
#' (bootstrapping for CI)
#'
#' The function \code{ce_estimate_iptw_att_boot} implements
#' IPTW with bootstrapping to estimate ATT effect with
#' multiple treatments using observational data.
#'
#' @param y A numeric vector (0, 1) representing a binary outcome.
#' @param x A dataframe, including all the covariates but not treatments.
#' @param w A numeric vector representing the treatment groups.
#' @param reference_trt A numeric value indicating reference treatment group
#' for ATT effect.
#' @param method A character string. Users can selected from the
#' following methods including \code{"IPTW-Multinomial"},
#'  \code{"IPTW-GBM"}, \code{"IPTW-SL"}.
#' @param nboots A numeric value representing the number of bootstrap samples.
#' @param verbose_boot A logical value indicating whether to print
#' the progress of nonparametric bootstrap. The default is \code{TRUE}.
#' @param ... Other parameters that can be passed through to functions.
#'
#' @return A summary of the effect estimates can be obtained
#' with \code{summary} function.
#' @importFrom stringr str_sub
#' @references
#' Hadley Wickham (2019).
#' \emph{stringr: Simple, Consistent Wrappers for Common String Operations}.
#' R package version 1.4.0.
#' URL:\url{https://CRAN.R-project.org/package=stringr}
ce_estimate_iptw_att_boot <-
  function(y,
           x,
           w,
           reference_trt,
           method,
           nboots,
           verbose_boot,
           ...) {
    sl_library <- parent.frame()$sl_library
    # Get the number of treatment group
    n_trt <- length(unique(w))
    for (i in 1:(n_trt - 1)) {
      assign(paste0("iptw_multitrt_att_result_", i, "_all"), NULL)
    }
    names_result <- NULL
    # Start bootstrapping
    for (j in 1:nboots) {
      bootstrap_id <- sample(length(y), replace = T)
      y_boot <- y[bootstrap_id]
      trt_boot <- w[bootstrap_id]
      x_boot <- x[bootstrap_id, ]
      iptw_multitrt_att_result <- ce_estimate_iptw_att(
        y = y_boot,
        x = x_boot,
        w = trt_boot,
        reference_trt,
        method = method,
        ...
      )
      names_result <- names(iptw_multitrt_att_result)
      for (i in 1:(n_trt - 1)) {
        assign(
          paste0("iptw_multitrt_att_result_", i, "_all"),
          cbind(eval(parse(
            text = paste0("iptw_multitrt_att_result_", i, "_all")
          )), iptw_multitrt_att_result[[i]])
        )
      }
      if (verbose_boot == TRUE) {
        print(paste0("Finish bootstrapping ", j))
      }
    }
    # Save the results of bootstrapping
    result <- NULL
    for (i in 1:(n_trt - 1)) {
      assign(paste0("RD_", i), list(as.double(eval(
        parse(text = paste0(
          "iptw_multitrt_att_result_", i, "_all"
        ))
      )[1, ])))
      assign(paste0("RR_", i), list(as.double(eval(
        parse(text = paste0(
          "iptw_multitrt_att_result_", i, "_all"
        ))
      )[2, ])))
      assign(paste0("OR_", i), list(as.double(eval(
        parse(text = paste0(
          "iptw_multitrt_att_result_", i, "_all"
        ))
      )[3, ])))

      assign(paste0("RD_", i), stats::setNames(eval(parse(text = (
        paste0("RD_", i)
      ))), paste0(
        "ATT_RD", stringr::str_sub(names_result[i], 4, 5)
      )))
      assign(paste0("RR_", i), stats::setNames(eval(parse(text = (
        paste0("RR_", i)
      ))), paste0(
        "ATT_RR", stringr::str_sub(names_result[i], 4, 5)
      )))
      assign(paste0("OR_", i), stats::setNames(eval(parse(text = (
        paste0("OR_", i)
      ))), paste0(
        "ATT_OR", stringr::str_sub(names_result[i], 4, 5)
      )))
      result <-
        c(result, (eval(parse(text = (
          paste0("RD_", i)
        )))), (eval(parse(text = (
          paste0("RR_", i)
        )))), (eval(parse(text = (
          paste0("OR_", i)
        )))))
    }
    result <- c(result, list(method = parent.frame()$method))
    class(result) <- "CIMTx_ATT_posterior"
    return(result)
  }
