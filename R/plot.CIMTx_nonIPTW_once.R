#' Plot for non-IPTW estimation methods for ATE effect
#'
#' @param x a \code{CIMTx_nonIPTW_once} object obtained
#' @param ... further arguments passed to or from other methods.
#'
#' @return an error message
#' @export
#'
#' @examples
#' \dontrun{
#' result <- list(method = "TMLE")
#' class(result) <- "CIMTx_nonIPTW_once"
#' plot(result)
#' }
plot.CIMTx_nonIPTW_once <- function(x, ...) {
  stop(
    "Plot function is not supported for estimation method using ",
    x$method,
    ". Plot function is only supported for estimation method using IPTW"
  )
}
