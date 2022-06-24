#' Contour plot for the grid specification of sensitivity analysis
#'
#' This function make the countor plot after the grid specification
#' of sensitivity analysis.
#' The input of the function is from the output of the sa.R function.
#'
#' @param x Object from sa function
#' @param ate a character indicating the ATE effect to plot,
#' eg, "1,3" or "2,3"
#' @param att a character indicating the ATT effect to plot,
#' eg, "1,3" or "1,2"
#' @param ... further arguments passed to or from other methods.
#'
#' @return A ggplot figure
#' @importFrom stringr str_sub
#' @importFrom ggplot2 ggplot aes_string geom_contour
#' labs theme element_text theme_bw
#' @importFrom metR geom_text_contour
#' @export
#'
#' @references
#'
#' Hadley Wickham (2019).
#' \emph{stringr: Simple, Consistent Wrappers for Common String Operations}.
#' R package version 1.4.0.
#' URL:\url{https://CRAN.R-project.org/package=stringr}
#'
#' Elio Campitelli (2021).
#'  \emph{metR: Tools for Easier Analysis of Meteorological Fields}.
#'   R package version 0.11.0.
#'   URL:\url{https://github.com/eliocamp/metR}
#'
#' @examples
#' sa_object_example <- list(
#'   ATE13 = seq(0, 1, length.out = 25), grid_index = c(4, 5),
#'   c_functions = data.frame(
#'     c4 = rep(seq(-0.6, 0, 0.15), each = 5),
#'     c5 = rep(seq(0, 0.6, 0.15), 5)
#'   )
#' )
#' class(sa_object_example) <- "CIMTx_sa_grid"
#' plot(sa_object_example, ate = "1,3")
plot.CIMTx_sa_grid <- function(x, ate = NULL, att = NULL, ...) {
  if (!is.null(ate)) {
    estimand <- paste0("ATE", as.numeric(gsub("\\,", "", ate)))
    m <- stringr::str_sub(ate, 1, 1)
    n <- stringr::str_sub(ate, 3, 3)
    plot_data <- x$c_functions %>%
      as.data.frame() %>%
      select(paste0("c", x$grid_index)) %>%
      mutate(estimand = x[[estimand]])
    plot_result <- plot_data %>%
      ggplot2::ggplot(ggplot2::aes_string(
        x = names(plot_data)[1],
        y = names(plot_data)[2],
        z = names(plot_data)[3]
      )) +
      ggplot2::geom_contour() +
      metR::geom_text_contour() +
      ggplot2::labs(
        x = bquote(italic(c) ~ "(" ~ .(m) ~ "," ~ .(n) ~ ")"),
        y = bquote(italic(c) ~ "(" ~ .(n) ~ "," ~ .(m) ~ ")"),
        title = bquote(ATE[.(m) ~ "," ~ .(n)])
      ) +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
      ggplot2::theme_bw() +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
    return(plot_result)
  }
  if (!is.null(att)) {
    estimand <- paste0("ATT", as.numeric(gsub("\\,", "", att)))
    m <- stringr::str_sub(att, 1, 1)
    n <- stringr::str_sub(att, 3, 3)
    plot_data <- x$c_functions %>%
      as.data.frame() %>%
      select(paste0("c", x$grid_index)) %>%
      mutate(estimand = x[[estimand]])
    plot_result <- plot_data %>%
      ggplot2::ggplot(ggplot2::aes_string(
        x = names(plot_data)[1],
        y = names(plot_data)[2],
        z = names(plot_data)[3]
      )) +
      ggplot2::geom_contour() +
      metR::geom_text_contour() +
      ggplot2::labs(
        x = bquote(italic(c) ~ "(" ~ .(m) ~ "," ~ .(n) ~ ")"),
        y = bquote(italic(c) ~ "(" ~ .(n) ~ "," ~ .(m) ~ ")"),
        title = bquote(ATT[.(m) ~ "," ~ .(n)])
      ) +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
      ggplot2::theme_bw() +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
    return(plot_result)
  }
}
