#' Plot for non-IPTW estimation methods for ATT effect
#'
#' @param x a \code{CIMTx_ATT_posterior} object obtained
#' @param ... further arguments passed to or from other methods.
#'
#' @return an error message
#' @export
#'
#' @examples
#' \dontrun{
#' result <- list(method = "RA")
#' class(result) <- "CIMTx_ATT_posterior"
#' plot(result)
#' }
plot.CIMTx_ATT_posterior <- function(x, ...) {
  if (x$method %in% c("IPTW-Multinomial", "IPTW-GBM", "IPTW-SL")) {
    stop(
      "Plot function is not supported for estimation method using ",
      x$method,
      " when boot is set to TRUE. Plot function is only supported for",
      " estimation method using IPTW when boot is set to FALSE."
    )
  } else {
    stop(
      "Plot function is not supported for estimation method using ",
      x$method,
      ". Plot function is only supported for estimation method using",
      " IPTW when boot is set to FALSE."
    )
  }
}
