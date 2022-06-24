#' Print the ATE results for non-IPTW results
#'
#' The \code{\link{print}} method for class "CIMTx_ATE_posterior"
#'
#' @param x a \code{CIMTx_ATE_posterior} object obtained
#' from \code{\link{ce_estimate}} function
#' @param ... further arguments passed to or from other methods.
#'
#' @return The output from \code{\link{print}}
#' @export
#'
#' @examples
#' result <- list(method = "RA")
#' class(result) <- "CIMTx_ATE_posterior"
#' print(result)
print.CIMTx_ATE_posterior <- function(x, ...) {
  if (x$method %in% c("RA", "BART")) {
    cat(
      "This is ATE results from estimation method",
      x$method,
      "with confidence interval estimated by Bayesian posterior samples.",
      "For effect estimates, please use summary function."
    )
  } else {
    cat(
      "This is ATE results from estimation method",
      x$method,
      "with confidence interval estimated by bootstrapping.",
      " For effect estimates, please use summary function."
    )
  }
}
