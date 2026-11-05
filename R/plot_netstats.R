#' Plot network-statistic trajectories
#' @param S Finite numeric matrix/data frame with one row per observation.
#' @param network_stats Optional distinct column names to display; NULL uses all.
#' @param cp Optional valid changepoint indices (last observations of segments).
#' @param title Optional plot title.
#' @param shade_segments Shade alternating segments when cp is supplied?
#' @param outfile Optional file path for saving the ggplot.
#' @param width,height Plot dimensions in inches for saving.
#' @param dpi Raster resolution when saving.
#' @returns A ggplot object. All numeric columns are supported, not only the
#'   four statistics computed by [compute_netstats()].
#' @export
plot_netstats_paper <- function(S, network_stats = NULL, cp = NULL,
                                title = NULL, shade_segments = TRUE,
                                outfile = NULL, width = 7, height = 4.8, dpi = 600) {
  if (is.data.frame(S) && all(vapply(S, is.numeric, logical(1)))) {
    S <- as.matrix(S)
  }
  if (!is.matrix(S) || !is.numeric(S) || is.complex(S) ||
      any(dim(S) < 1L) || any(!is.finite(S))) {
    stop("S must be a nonempty finite numeric matrix or numeric data frame.")
  }
  Tt <- nrow(S)
  if (!is.logical(shade_segments) || length(shade_segments) != 1L ||
      is.na(shade_segments)) stop("shade_segments must be TRUE or FALSE.")
  if (is.null(colnames(S))) colnames(S) <- paste0("Statistic ", seq_len(ncol(S)))
  if (anyDuplicated(colnames(S))) stop("S must have distinct column names.")
  if (is.null(network_stats)) network_stats <- colnames(S)
  if (!is.character(network_stats) || !length(network_stats) || anyNA(network_stats)) {
    stop("network_stats must contain column names of S.")
  }
  network_stats[network_stats == "triangle"] <- "triangles"
  if (anyDuplicated(network_stats) || any(!network_stats %in% colnames(S))) {
    stop("network_stats must contain distinct column names of S.")
  }
  if (length(cp)) {
    if (!is.numeric(cp) || is.complex(cp) || any(!is.finite(cp)) ||
        any(cp != floor(cp)) || any(cp < 1 | cp >= Tt)) {
      stop("cp must contain integer indices in 1, ..., nrow(S) - 1.")
    }
    cp <- sort(unique(cp))
  }
  S <- S[, network_stats, drop = FALSE]
  nice <- c(edges = "Edges", density = "Density", triangles = "Triangles", isolates = "Isolates")
  labels <- network_stats
  known <- network_stats %in% names(nice)
  labels[known] <- nice[network_stats[known]]
  # Keep facet identities separate from display labels.
  df <- data.frame(time = rep(seq_len(Tt), ncol(S)),
                   statistic = factor(rep(network_stats, each = Tt), levels = network_stats),
                   value = as.vector(S))
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$time, y = .data$value))
  if (shade_segments && length(cp)) {
    cuts <- c(0, cp, Tt)
    shade <- data.frame(xmin = head(cuts, -1) + 0.5,
                        xmax = tail(cuts, -1) + 0.5, segment = seq_len(length(cuts) - 1L))
    shade <- shade[shade$segment %% 2L == 1L, , drop = FALSE]
    p <- p + ggplot2::geom_rect(data = shade,
                                ggplot2::aes(xmin = .data$xmin, xmax = .data$xmax),
                                ymin = -Inf, ymax = Inf, inherit.aes = FALSE,
                                fill = "grey94", color = NA)
  }
  p <- p + ggplot2::geom_point(size = 0.85, alpha = 0.55) +
    ggplot2::facet_wrap(~ statistic, ncol = 1, scales = "free_y", strip.position = "left",
                         labeller = ggplot2::as_labeller(stats::setNames(labels, network_stats))) +
    ggplot2::scale_x_continuous(breaks = scales::pretty_breaks(n = 8),
                               expand = ggplot2::expansion(mult = c(0.01, 0.01))) +
    ggplot2::scale_y_continuous(breaks = scales::pretty_breaks(n = 6),
                               labels = scales::label_number(big.mark = ","),
                               expand = ggplot2::expansion(mult = c(0.08, 0.12))) +
    ggplot2::labs(x = "Time", y = NULL, title = title) +
    ggplot2::theme_classic(base_size = 12, base_family = "serif") +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
                   strip.background = ggplot2::element_blank(), strip.placement = "outside",
                   strip.text.y.left = ggplot2::element_text(angle = 0, face = "bold"),
                   panel.spacing.y = grid::unit(10, "pt"),
                   plot.margin = ggplot2::margin(6, 10, 6, 6))
  if (length(cp)) p <- p + ggplot2::geom_vline(xintercept = cp, linetype = "dashed", linewidth = 0.55)
  if (!is.null(outfile)) ggplot2::ggsave(outfile, p, width = width, height = height, dpi = dpi)
  return(p)
}
