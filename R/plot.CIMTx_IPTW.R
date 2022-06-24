#' Boxplot for weight distribution
#'
#' This function make the boxplot plot for the weights estimated by
#' different IPTW methods.
#' The inputs of the function are from the output of ce_estimate.R function
#' when the methods are "IPTW-Multinomial", "IPTW-GBM", "IPTW-SL".
#'
#' @param ... Objects from IPTW related methods
#'
#' @return A ggplot figure
#' @importFrom stringr str_detect
#' @importFrom dplyr mutate
#' @importFrom ggplot2 ggplot geom_boxplot scale_fill_manual
#' facet_wrap labs theme_bw theme element_blank
#' @importFrom cowplot plot_grid get_legend
#' @export
#'
#' @references
#' Hadley Wickham (2019).
#' \emph{stringr: Simple, Consistent Wrappers for Common String Operations}.
#' R package version 1.4.0.
#' URL:\url{https://CRAN.R-project.org/package=stringr}
#'
#' Hadley Wickham, Romain François, Lionel Henry and Kirill Müller (2021).
#' \emph{dplyr: A Grammar of Data Manipulation}.
#' R package version 1.0.7.
#' URL: \url{https://CRAN.R-project.org/package=dplyr}
#'
#' H. Wickham. \emph{ggplot2: Elegant Graphics for Data Analysis}.
#' Springer-Verlag New York, 2016.
#'
#' Claus O. Wilke (2020).
#' \emph{cowplot: Streamlined Plot Theme and Plot Annotations for 'ggplot2'}.
#' R package version 1.1.1.
#' URL:\url{https://CRAN.R-project.org/package=cowplot}
#' @examples
#' iptw_object_example <- list(weight = rnorm(1000, 1, 1), method = "IPTW-SL")
#' class(iptw_object_example) <- "CIMTx_IPTW"
#' plot(iptw_object_example)
plot.CIMTx_IPTW <- function(...) {
  iptw_object <- list(...)
  weight <- NULL
  method <- NULL
  weight_all <- NULL
  method_all <- NULL
  trim <- NULL
  for (i in seq_len(length(iptw_object))) {
    weight_all <- c(weight_all, iptw_object[[i]]$weight)
    method_all <-
      c(method_all,
        rep(iptw_object[[i]]$method,
            length(iptw_object[[i]]$weight)))
  }
  boxplot_figure_data <- data.frame(weight = weight_all,
                                    method = method_all)
  if (any(stringr::str_detect(method_all, "Trim")) == TRUE &&
      length(unique(method_all)) == 6) {
    boxplot_figure_data <- boxplot_figure_data %>%
      dplyr::mutate(trim = ifelse(
        stringr::str_detect(method_all, "Trim"),
        "(b) trimmed weights",
        "(a) untrimmed weights"
      ))

    boxplot_figure_data_reordered <- boxplot_figure_data %>%
      dplyr::mutate(method = factor(
        method,
        levels = c(
          "IPTW-Multinomial",
          "IPTW-SL",
          "IPTW-GBM",
          "IPTW-Multinomial-Trim",
          "IPTW-SL-Trim",
          "IPTW-GBM-Trim"
        )
      ))
    boxplot_figure_main <-
      ggplot2::ggplot(ggplot2::aes(x = method, y = weight, fill = method),
                      data = boxplot_figure_data_reordered) +
      ggplot2::geom_boxplot() +
      ggplot2::scale_fill_manual(values = c(
        "#CCCCCC",
        "#4D4D4D",
        "white",
        "#CCCCCC",
        "#4D4D4D",
        "white"
      )) +
      ggplot2::facet_wrap(~ trim, scales = "free_x") +
      ggplot2::labs(color = "", y = "Weights", x = "") +
      ggplot2::theme_bw() +
      ggplot2::theme(
        legend.position = "none",
        axis.title.x = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank()
      )
    boxplot_figure_data_for_legend <- boxplot_figure_data %>%
      dplyr::mutate(method = ifelse(
        trim == "(b) trimmed weights",
        stringr::str_sub(method, 1, -6),
        method
      )) %>%
      dplyr::mutate(method = factor(
        method,
        levels = c("IPTW-Multinomial", "IPTW-SL", "IPTW-GBM")
      ))

    boxplot_figure_legend <-
      ggplot2::ggplot(ggplot2::aes(x = method, y = weight, fill = method),
                      data = boxplot_figure_data_for_legend) +
      ggplot2::geom_boxplot() +
      ggplot2::scale_fill_manual(values = c(
        "#CCCCCC",
        "#4D4D4D",
        "white",
        "#CCCCCC",
        "#4D4D4D",
        "white"
      )) +
      ggplot2::facet_wrap(~ trim, scales = "free_x") +
      ggplot2::labs(color = "",
                    y = "Weights",
                    fill = "") +
      ggplot2::theme_bw() +
      ggplot2::theme(legend.position = "top")

    boxplot_figure <-
      cowplot::plot_grid(
        cowplot::get_legend(boxplot_figure_legend),
        boxplot_figure_main,
        rel_heights = c(0.1, 0.9),
        nrow = 2
      )
  } else {
    boxplot_figure <-
      ggplot2::ggplot(ggplot2::aes(x = method, y = weight, fill = method),
                      data = boxplot_figure_data) +
      ggplot2::geom_boxplot() +
      ggplot2::scale_fill_manual(values = c(
        "#CCCCCC",
        "#4D4D4D",
        "white",
        "#CCCCCC",
        "#4D4D4D",
        "white"
      )) +
      ggplot2::labs(color = "",
                    y = "Weights",
                    fill = "") +
      ggplot2::theme_bw() +
      ggplot2::theme(legend.position = "top")
  }

  return(boxplot_figure)
}
