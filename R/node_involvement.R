# Convert labels such as "B2 - B1" to plotmath strings; internal only.
as_Bdiff_plotmath_label <- function(x) {
  hit <- grepl("^\\s*B[0-9]+\\s*-\\s*B[0-9]+\\s*$", x)
  x[hit] <- sub("^\\s*B([0-9]+)\\s*-\\s*B([0-9]+)\\s*$",
                "B[\\1] - B[\\2]", x[hit])
  x
}

#' Rank nodes by posterior row-wise L2 norms
#'
#' Compute each node's row norm within each posterior draw of a coefficient
#' difference, summarize those norms across draws, and rank nodes in descending
#' order. This is not the norm of an entry-wise posterior summary matrix.
#' @param B_diff One element of `heat_out$B_diff`, as returned by
#'   [compute_B_diffs_from_out()] or [plot_B_diff_heatmaps()]. It contains
#'   `draws_upper` (draws by unordered dyads), `n` (node count), `upper_idx`
#'   (the logical strict upper-triangle mask), and `label`. The node count and
#'   mask can be inferred when omitted. A numeric vector represents one draw.
#' @param node_labels Optional node names in adjacency-matrix order. NULL uses
#'   node indices as labels.
#' @param summary Posterior median (default) or mean of each node's draw-wise
#'   row norms.
#' @param output Return a `"data.frame"` or a `"matrix"`. A matrix is character
#'   because it also contains node names and difference labels.
#' @param label_format Use plotmath strings (for example, `"B[2] - B[1]"`) or
#'   plain difference labels (for example, `"B2 - B1"`).
#' @returns A data frame, or character matrix, with columns `diff`, `rank`,
#'   `node`, `node_name`, and `l2_norm`. Ties are ordered by original node index.
#' @details Only off-diagonal dyads contribute to the norms. Each unordered
#'   dyad contributes once to each of its two endpoint nodes.
#' @examples
#' d <- list(draws_upper = rbind(c(3, 4, 0), c(3, 4, 0)), n = 3,
#'           label = "B2 - B1")
#' rank_nodes_by_Bdiff_L2(d, c("A", "B", "C"))
#' @seealso [plot_node_involvement_L2()], [compute_B_diffs_from_out()]
#' @export
rank_nodes_by_Bdiff_L2 <- function(B_diff, node_labels = NULL,
                                  summary = c("median", "mean"),
                                  output = c("data.frame", "matrix"),
                                  label_format = c("plotmath", "plain")) {
  summary <- match.arg(summary)
  output <- match.arg(output)
  label_format <- match.arg(label_format)
  if (!is.list(B_diff) || is.null(B_diff$draws_upper)) {
    stop("B_diff must contain draws_upper; use heat_out$B_diff[[k]].",
         call. = FALSE)
  }
  draws_upper <- B_diff$draws_upper
  if (is.null(dim(draws_upper)) && is.numeric(draws_upper)) {
    draws_upper <- matrix(draws_upper, nrow = 1L)
  }
  draws_upper <- .bcmlr_matrix(draws_upper, "B_diff$draws_upper")
  n_edge <- ncol(draws_upper)
  n <- B_diff$n
  if (is.null(n)) n <- (1 + sqrt(1 + 8 * n_edge)) / 2
  .bcmlr_number(n, "B_diff$n", lower = 2, integer = TRUE)
  if (n_edge != n * (n - 1) / 2) {
    stop("draws_upper must have n * (n - 1) / 2 columns.", call. = FALSE)
  }
  n <- as.integer(n)
  expected_idx <- upper.tri(matrix(FALSE, n, n), diag = FALSE)
  upper_idx <- B_diff$upper_idx
  if (is.null(upper_idx)) upper_idx <- expected_idx
  if (!is.matrix(upper_idx) || !is.logical(upper_idx) ||
      !identical(dim(upper_idx), c(n, n)) || anyNA(upper_idx) ||
      !all(upper_idx == expected_idx)) {
    stop("upper_idx must be the logical strict upper-triangle mask.",
         call. = FALSE)
  }
  if (is.null(node_labels)) node_labels <- as.character(seq_len(n))
  node_labels <- unname(as.character(node_labels))
  .bcmlr_labels(node_labels, n, "node_labels")
  diff_label <- if (is.null(B_diff$label)) "Bdiff" else B_diff$label
  .bcmlr_labels(diff_label, 1L, "B_diff$label")
  if (label_format == "plotmath") diff_label <- as_Bdiff_plotmath_label(diff_label)

  # Preserve the paper's order of operations: norms first, posterior summary second.
  upper_ij <- which(upper_idx, arr.ind = TRUE)
  E <- matrix(0, nrow = n_edge, ncol = n)
  E[cbind(seq_len(n_edge), upper_ij[, 1])] <- 1
  E[cbind(seq_len(n_edge), upper_ij[, 2])] <- 1
  L <- sqrt((draws_upper^2) %*% E)
  if (any(!is.finite(L))) stop("Non-finite row norms; check coefficient magnitudes.")
  l2 <- if (summary == "median") apply(L, 2, stats::median) else colMeans(L)
  ans <- data.frame(diff = diff_label, node = seq_len(n),
                    node_name = node_labels, l2_norm = as.numeric(l2),
                    row.names = NULL)
  ans <- ans[order(-ans$l2_norm, ans$node), , drop = FALSE]
  ans$rank <- seq_len(n)
  ans <- ans[, c("diff", "rank", "node", "node_name", "l2_norm")]
  rownames(ans) <- NULL
  if (output == "matrix") as.matrix(ans) else ans
}

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
#'   The displayed nodes are the union of the top `top_n` nodes at each
#'   changepoint. Their scores use norms computed within draws, followed by
#'   the requested posterior summary. `top_n` is capped at the node count.
#' @examples
#' d <- list(draws_upper = rbind(c(3, 4, 0), c(3, 4, 0)), n = 3,
#'           label = "B2 - B1")
#' plot_node_involvement_L2(list(B_diff = list(d)), c("A", "B", "C"))
#' @seealso [rank_nodes_by_Bdiff_L2()], [plot_B_diff_heatmaps()]
#' @export
plot_node_involvement_L2 = function(heat_out,
                                     node_labels,
                                     top_n = 8,
                                     summary = c("median", "mean"),
                                     mar = c(5, 10, 2, 1)) {
  summary = match.arg(summary)
  if (!is.list(heat_out) || !is.list(heat_out$B_diff) ||
      !length(heat_out$B_diff)) {
    stop("heat_out must contain a nonempty B_diff list.", call. = FALSE)
  }
  K = length(heat_out$B_diff)
  # Normalize node labels, and allow configurable label margins.
  node_names = unname(as.character(node_labels))
  .bcmlr_number(top_n, "top_n", lower = 1, integer = TRUE)
  if (!is.numeric(mar) || is.complex(mar) || length(mar) != 4L || any(!is.finite(mar)) ||
      any(mar < 0)) stop("mar must contain four nonnegative margins.")
  
  # Sanity checks
  if (K < 1) stop("heat_out$B_diff is empty.")
  if (!length(node_names)) stop("node_labels must not be empty.", call. = FALSE)
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
