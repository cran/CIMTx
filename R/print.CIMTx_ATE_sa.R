#' Print the ATE results for from sensitivity analysis
#'
#' The \code{\link{print}} method for class "CIMTx_ATE_sa"
#'
#' @param x a \code{CIMTx_ATE_sa} object obtained
#' from \code{\link{sa}} function
#' @param ... further arguments passed to or from other methods.
#'
#' @return The output from \code{\link{print}}
#' @export
#'
#' @examples
#' result <- list(estimand = "ATE")
#' class(result) <- "CIMTx_ATE_sa"
#' print(result)
print.CIMTx_ATE_sa <- function(x, ...) {
  cat(
    "This is sa results using CIMTx pakcage",
    "For effect estimates, please use summary function."
  )
}
