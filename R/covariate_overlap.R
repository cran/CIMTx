#' Covariate overlap figure
#'
#' @param treatment A numeric vector representing the treatment groups.
#' @param prob Probability of being assigned to the treatment group.
#'
#' @return A ggplot figure.
#'
#' @importFrom ggplot2 ggplot aes geom_boxplot ggtitle xlab ylab
#' theme_bw ylim theme element_text
#' @importFrom cowplot ggdraw draw_label plot_grid
#' @references
#' H. Wickham.
#' \emph{ggplot2: Elegant Graphics for Data Analysis}.
#' Springer-Verlag New York, 2016.
#'
#' Claus O. Wilke (2020).
#' \emph{cowplot: Streamlined Plot Theme and Plot Annotations for 'ggplot2'}.
#' R package version 1.1.1.
#' URL:\url{https://CRAN.R-project.org/package=cowplot}
covariate_overlap <- function(treatment, prob) {
  treatment <- as.factor(treatment)
  df <- data.frame(cbind(treatment, prob))
  plot_list <- list()
  for (i in seq_len(length(unique(treatment)))) {
    # gather boxplot for each treatment
    plot_list[[i]] <- local({
      i <- i
      ggplot2::ggplot(ggplot2::aes(
        y = prob[, i],
        x = treatment,
        group = treatment
      ),
      data = df) +
        ggplot2::geom_boxplot() +
        ggplot2::ggtitle(bquote(
          italic(P) ~ "(" ~ italic(W) ~ "=" ~ .(i) ~ "|" ~ bolditalic(X) ~ ")"
        )) +
        ggplot2::xlab("Treatment") +
        ggplot2::ylab("") +
        ggplot2::theme_bw()
    }) +
      ggplot2::ylim(0, 1) +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  }

  # set title
  title <- cowplot::ggdraw() +
    cowplot::draw_label(
      "Covariate Overlap",
      fontface = "bold",
      x = 0.42,
      hjust = 0
    )

  # if the number of treatment equal or less than 3,
  # then plot them in a row together
  if (length(unique(treatment)) <= 3) {
    cowplot::plot_grid(
      title,
      cowplot::plot_grid(plotlist = plot_list, nrow = 1),
      ncol = 1,
      rel_heights = c(0.1, 1)
    )
  } else {
    cowplot::plot_grid(
      title,
      cowplot::plot_grid(plotlist = plot_list),
      ncol = 1,
      rel_heights = c(0.1, 1)
    )
  }
}
