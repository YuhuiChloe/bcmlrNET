#' Plot node involvement across changepoints
#' @param heat_out Output of [compute_B_diffs_from_out()] or
#'   [plot_B_diff_heatmaps()] with sample-based B_diff objects.
#' @param node_labels Node names in adjacency-matrix order.
#' @param top_n Number of top nodes selected from each changepoint; the plot
#'   shows their union.
#' @param summary Use posterior median or mean of sample-wise row L2 norms.
#' @param mar Four plot margins; increase the second for long node labels.
#' @returns Invisibly, a list with l2_rank, active_nodes and L2_mat. Draws on
#'   the current base-graphics device; does not return a ggplot.
#' @details Multi-panel base-graphics layout is used and then reset to one
#'   panel. Use a separate device rather than nesting in another layout.
#' @export
plot_node_involvement_L2 = function(heat_out,
                                     node_labels,
                                     top_n = 8,
                                     summary = c("median", "mean"),
                                     mar = c(5, 10, 2, 1)) {
  summary = match.arg(summary)
  K = length(heat_out$B_diff)
  # Normalize node labels, and allow configurable label margins.
  node_names = unname(as.character(node_labels))
  .bcmlr_number(top_n, "top_n", lower = 1, integer = TRUE)
  if (!is.numeric(mar) || is.complex(mar) || length(mar) != 4L || any(!is.finite(mar)) ||
      any(mar < 0)) stop("mar must contain four nonnegative margins.")
  
  # Sanity checks
  if (K < 1) stop("heat_out$B_diff is empty.")
  if (length(node_names) != heat_out$B_diff[[1]]$n) {
    stop("node_labels must have length equal to the number of nodes.")
  }
  # top_n was validated above; clamp before integer conversion.
  top_n = as.integer(min(top_n, length(node_names)))
  
  # Posterior L2 norms for each changepoint
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
  top_nodes = lapply(l2_ls, function(x) head(x$node, top_n))
  active_nodes = sort(unique(unlist(top_nodes)))
  
  # L2 matrix for selected nodes across all changepoints
  L2_mat = sapply(seq_len(K), function(k) {
    x = l2_ls[[k]]
    x$l2_norm[match(active_nodes, x$node)]
  })
  # Keep a matrix when only one node or one changepoint is displayed.
  L2_mat = matrix(L2_mat, nrow = length(active_nodes), ncol = K)
  
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
    mar = mar,
    family = "serif",
    font.lab = 1,
    font.axis = 1,
    font.main = 2,
    cex.lab = 1.5,
    cex.main = 1.5
  )
  
  graphics::image(
    # Explicit cell boundaries support one-node/one-changepoint images.
    x = seq(0.5, ncol(L2_mat) + 0.5, by = 1),
    y = seq(0.5, nrow(L2_mat) + 0.5, by = 1),
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
  
  # Shared color bar
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
