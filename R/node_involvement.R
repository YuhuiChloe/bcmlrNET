#=============================================#
# Plot node involvement by row-wise L2 norms ####
#=============================================#

#' Plot node involvement across changepoints
#'
#' @param heat_out Output of [plot_B_diff_heatmaps()] or [compute_B_diffs_from_out()],
#'   with a nonempty `B_diff` list retaining posterior draws.
#' @param node_labels Character vector of node names in adjacency-matrix order.
#' @param top_n Number of highest-ranked nodes selected at each changepoint;
#'   default eight. The plot displays their union.
#' @param summary Posterior median or mean of row-wise L2 norms computed within draws.
#' @returns Invisibly, a list with `l2_rank`, `active_nodes`, and `L2_mat`.
#'   Draws a heatmap on the current base-graphics device.
#' @details Uses [rank_nodes_by_Bdiff_L2()]. The original node ordering, colors,
#'   and graphics layout are retained. The left margin expands when needed to fit
#'   long node labels at the original text size.
#' @examples
#' d <- list(draws_upper = rbind(c(3, 4, 0), c(-3, -4, 0)), n = 3,
#'           label = "B2 - B1")
#' plot_node_involvement_L2(list(B_diff = list(d)), c("A", "B", "C"))
#' @export
plot_node_involvement_L2 = function(heat_out,
                                     node_labels,
                                     top_n = 8,
                                     summary = c("median", "mean")) {
  summary = match.arg(summary)
  K = length(heat_out$B_diff)
  node_names = unname(node_labels)

  # Sanity checks
  if (K < 1) stop("heat_out$B_diff is empty.")
  if (length(node_names) != heat_out$B_diff[[1]]$n) {
    stop("node_labels must have length equal to the number of nodes.")
  }
  if (length(top_n) != 1 || top_n < 1 || top_n != as.integer(top_n)) {
    stop("top_n must be a positive integer.")
  }
  top_n = min(as.integer(top_n), length(node_names))

  #......Posterior L2 norms for each changepoint
  l2_ls = lapply(seq_len(K), function(k) {
    rank_nodes_by_Bdiff_L2(
      B_diff = heat_out$B_diff[[k]],
      node_labels = node_names,
      summary = summary
    )
  })

  # Verify that node indices map directly to node_labels
  check_labels = data.frame(
    node = seq_along(node_names),
    expected = node_names,
    returned = l2_ls[[1]]$node_name[
      match(seq_along(node_names), l2_ls[[1]]$node)
    ]
  )
  stopifnot(identical(check_labels$expected, check_labels$returned))

  # Union of top nodes from each changepoint
  top_nodes = lapply(l2_ls, function(x) utils::head(x$node, top_n))
  active_nodes = sort(unique(unlist(top_nodes)))

  # L2 matrix for selected nodes across all changepoints
  L2_mat = sapply(seq_len(K), function(k) {
    x = l2_ls[[k]]
    x$l2_norm[match(active_nodes, x$node)]
  })
  if (K == 1) L2_mat = matrix(L2_mat, ncol = 1)

  rownames(L2_mat) = node_names[active_nodes]
  colnames(L2_mat) = paste0("B", 2:(K + 1), " - B", 1:K)

  # Order selected nodes by maximum L2 norm across changepoints
  ord = order(apply(L2_mat, 1, max), decreasing = TRUE)
  L2_mat = L2_mat[ord, , drop = FALSE]
  row_nodes = active_nodes[ord]
  rownames(L2_mat) = node_names[row_nodes]

  # Colors
  zmax = max(L2_mat)
  if (!is.finite(zmax) || zmax <= 0) zmax = 1

  cols = grDevices::colorRampPalette(
    c("white", "#DCEAF7", "#9ECAE1", "#4292C6", "#2171B5", "#084594")
  )(101)

  # Heatmap and shared color bar
  old_par = graphics::par(no.readonly = TRUE)
  on.exit({
    graphics::layout(1)
    graphics::par(old_par)
  }, add = TRUE)

  graphics::layout(
    matrix(c(1, 2), nrow = 1),
    widths = c(0.9, 0.1)
  )

  graphics::par(
    mar = c(5, 10, 2, 1),
    family = "serif",
    font.lab = 1,
    font.axis = 1,
    font.main = 2,
    cex.lab = 1.5,
    cex.main = 1.5
  )

  # Reserve enough space for long country names at the existing axis text size.
  plot_mai = graphics::par("mai")
  plot_mai[2] = max(plot_mai[2], max(graphics::strwidth(
    rownames(L2_mat), units = "inches", cex = 1.3, font = 1
  )) + 0.2)
  graphics::par(mai = plot_mai)

  graphics::image(
    x = seq_len(ncol(L2_mat)),
    y = seq_len(nrow(L2_mat)),
    z = t(L2_mat[nrow(L2_mat):1, , drop = FALSE]),
    zlim = c(0, zmax),
    col = cols,
    axes = FALSE,
    xlab = "Changepoint",
    ylab = "",
    main = bquote(
      bold("Posterior " * .(summary) * " of row-wise " * L[2] * " norms")
    )
  )

  graphics::axis(
    1,
    at = seq_len(ncol(L2_mat)),
    labels = seq_len(K),
    las = 1
  )

  graphics::axis(
    2,
    at = seq_len(nrow(L2_mat)),
    labels = rev(rownames(L2_mat)),
    las = 2,
    cex.axis = 1.3
  )

  graphics::box()

  #......Shared color bar
  graphics::par(
    mar = c(5, 0.5, 2, 2.5),
    family = "serif",
    font.axis = 1
  )

  graphics::plot.new()
  graphics::plot.window(
    xlim = c(0, 1),
    ylim = c(0, zmax)
  )

  brks = seq(
    0,
    zmax,
    length.out = length(cols) + 1
  )

  for (i in seq_along(cols)) {
    graphics::rect(
      0.15, brks[i],
      0.75, brks[i + 1],
      col = cols[i],
      border = NA
    )
  }

  graphics::axis(4, las = 1)
  graphics::box()

  invisible(list(
    l2_rank = l2_ls,
    active_nodes = row_nodes,
    L2_mat = L2_mat
  ))
}
