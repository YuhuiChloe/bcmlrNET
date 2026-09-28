#' Plot credible intervals for changes in network-statistic coefficients
#' @param Beta Finite numeric array: samples by coefficients by chronological segments.
#' @param mfrow Optional c(number of rows, number of columns). NULL chooses a layout.
#' @param coef_labels One label per coefficient; defaults to sequential labels.
#' @param ci_level Equal-tailed posterior credible-interval level in (0, 1).
#' @param xlim Optional increasing numeric axis limits; zooms without dropping samples.
#' @param base_size Base font size.
#' @returns A ggplot object. Points are posterior means of next-minus-previous
#'   coefficient differences, with equal-tailed credible intervals.
#' @export
plot_beta_diff_CIs <- function(Beta, mfrow = NULL, coef_labels = NULL,
                               ci_level = 0.95, xlim = NULL, base_size = 12) {
  if (!is.numeric(Beta) || is.complex(Beta) || length(dim(Beta)) != 3L ||
      any(dim(Beta) < 1L) || any(!is.finite(Beta)) || dim(Beta)[3] < 2L) {
    stop("Beta must be a finite array: samples x coefficients x at least two segments.",
         call. = FALSE)
  }
  .bcmlr_number(ci_level, "ci_level")
  if (ci_level <= 0 || ci_level >= 1) stop("ci_level must lie in (0, 1).")
  .bcmlr_number(base_size, "base_size", lower = .Machine$double.eps)
  p <- dim(Beta)[2]
  J <- dim(Beta)[3]
  if (is.null(coef_labels)) coef_labels <- paste0("coef ", seq_len(p))
  .bcmlr_labels(coef_labels, p, "coef_labels")
  if (is.null(mfrow)) {
    nc <- min(3L, J - 1L)
    mfrow <- c(ceiling((J - 1L) / nc), nc)
  }
  if (!is.numeric(mfrow) || is.complex(mfrow) || length(mfrow) != 2L || any(!is.finite(mfrow)) ||
      any(mfrow < 1 | mfrow != floor(mfrow)) || prod(mfrow) < J - 1L) {
    stop("mfrow must provide enough panels for all coefficient differences.")
  }
  beta_diff <- Beta[, , -1, drop = FALSE] - Beta[, , -J, drop = FALSE]
  probs <- c((1 - ci_level) / 2, 1 - (1 - ci_level) / 2)
  means <- apply(beta_diff, c(2, 3), mean)
  lower <- apply(beta_diff, c(2, 3), stats::quantile, probs = probs[1], names = FALSE)
  upper <- apply(beta_diff, c(2, 3), stats::quantile, probs = probs[2], names = FALSE)
  df <- expand.grid(coef = seq_len(p), diff = seq_len(J - 1L), KEEP.OUT.ATTRS = FALSE)
  df$mean <- as.vector(means)
  df$lower <- as.vector(lower)
  df$upper <- as.vector(upper)
  df$coef <- factor(df$coef, levels = seq_len(p))
  facet_labels <- paste0("beta[", 2:J, "] - beta[", seq_len(J - 1L), "]")
  # Explicit levels retain chronological order when segment numbers exceed nine.
  df$diff_lab <- factor(facet_labels[df$diff], levels = facet_labels)
  if (is.null(xlim)) xlim <- range(df$lower, df$upper)
  if (!is.numeric(xlim) || is.complex(xlim) || length(xlim) != 2L ||
      any(!is.finite(xlim)) || xlim[1] > xlim[2]) stop("xlim must be finite ordered limits.")
  if (diff(xlim) == 0) xlim <- xlim + c(-1, 1)
  ggplot2::ggplot(df, ggplot2::aes(x = .data$mean, y = .data$coef)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed") +
    ggplot2::geom_errorbar(ggplot2::aes(xmin = .data$lower, xmax = .data$upper),
                           orientation = "y", width = 0.2) +
    ggplot2::geom_point(size = 1.6) +
    ggplot2::facet_wrap(~ diff_lab, nrow = mfrow[1], ncol = mfrow[2],
                         labeller = ggplot2::label_parsed) +
    ggplot2::scale_y_discrete(labels = coef_labels) +
    ggplot2::scale_x_continuous(breaks = scales::pretty_breaks(n = 6),
                               expand = ggplot2::expansion(mult = c(0.05, 0.05))) +
    ggplot2::coord_cartesian(xlim = xlim) +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme_bw(base_size = base_size, base_family = "serif") +
    ggplot2::theme(strip.text = ggplot2::element_text(face = "bold", size = base_size + 1),
                   panel.spacing.x = grid::unit(6, "pt"),
                   panel.spacing.y = grid::unit(6, "pt"),
                   plot.margin = ggplot2::margin(3, 3, 3, 3))
}
