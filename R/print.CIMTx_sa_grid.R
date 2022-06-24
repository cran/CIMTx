#' Print the ATT results for from sensitivity analysis
#'
#' The \code{\link{print}} method for class "CIMTx_sa_grid"
#'
#' @param x a \code{CIMTx_sa_grid} object obtained from
#' \code{\link{sa}} function
#' @param ... further arguments passed to or from other methods.
#'
#' @return The output from \code{\link{print}}
#' @export
#'
#' @examples
#' result <- list(estimand = "ATE")
#' class(result) <- "CIMTx_sa_grid"
#' print(result)
print.CIMTx_sa_grid <- function(x, ...) {
  cat(
    "This is sa results using CIMTx pakcage",
    "when the c functions involves a range of point mass priors",
    "To visualize effect estimates, please plot function."
  )
}
