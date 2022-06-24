#' Print the ATE/ATT results for non-IPTW results
#'
#' The \code{\link{print}} method for class "CIMTx_nonIPTW_once"
#'
#' @param x a \code{CIMTx_ATE_posterior} x obtained
#' from \code{\link{ce_estimate}} function
#' @param ... further arguments passed to or from other methods.
#'
#' @return The output from \code{\link{print}}
#' @export
#'
#' @examples
#' result <- list(method = "TMLE", estimand = "ATE")
#' class(result) <- "CIMTx_nonIPTW_once"
#' print(result)
print.CIMTx_nonIPTW_once <- function(x, ...) {
  cat(
    "This is",
    x$estimand,
    "results from estimation method",
    x$method,
    "using CIMTx. For effect estimates,",
    "please use summary function."
  )
}
