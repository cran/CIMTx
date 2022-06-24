#' Print the ATE/ATT results for IPTW results
#'
#' The \code{\link{print}} method for class "CIMTx_IPTW"
#'
#' @param x a \code{CIMTx_ATE_posterior} x obtained
#' @param ... further arguments passed to or from other methods.
#'
#' @return The output from \code{\link{print}}
#' @export
#'
#' @examples
#' result <- list(method = "IPTW-Multinomial", estimand = "ATE")
#' class(result) <- "CIMTx_IPTW"
#' print(result)
print.CIMTx_IPTW <- function(x, ...) {
  cat(
    "This is",
    x$estimand,
    "results from estimation method",
    x$method,
    "using CIMTx. For effect estimates,",
    "please use summary function.",
    "To visualize weight,",
    "please use plot function."
  )
}
