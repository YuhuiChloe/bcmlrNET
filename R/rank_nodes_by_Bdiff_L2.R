as_Bdiff_plotmath_label <- function(x) {
  x <- as.character(x)

  hit <- grepl("^\\s*B[0-9]+\\s*-\\s*B[0-9]+\\s*$", x)

  x[hit] <- sub(
    "^\\s*B([0-9]+)\\s*-\\s*B([0-9]+)\\s*$",
    "B[\\1] - B[\\2]",
    x[hit]
  )

  x
}

#' Rank nodes by posterior row-wise L2 norms
#'
#' @param B_diff One element of `heat_out$B_diff`, containing `draws_upper`
#'   (draws by upper-triangle dyads), `n`, `upper_idx`, and an optional `label`.
#'   The node count and mask are inferred when omitted, as in the paper script.
#' @param node_labels Node names in adjacency-matrix order; NULL uses node indices.
#' @param summary Summarize draw-wise row norms by their posterior median or mean.
#' @param output Return a data frame or a matrix. The matrix is character because
#'   it includes text labels.
#' @param label_format Use plotmath strings or plain difference labels.
#' @returns Nodes in decreasing order of their summarized row norms, with columns
#'   `diff`, `rank`, `node`, `node_name`, and `l2_norm`.
#' @details Each row norm is computed within a posterior draw, before taking its
#'   median or mean across draws. Only off-diagonal dyads contribute.
#' @examples
#' d <- list(draws_upper = rbind(c(3, 4, 0), c(-3, -4, 0)), n = 3,
#'           label = "B2 - B1")
#' rank_nodes_by_Bdiff_L2(d, c("A", "B", "C"))
#' @export
rank_nodes_by_Bdiff_L2 <- function(B_diff, node_labels = NULL,
                                   summary = c("median", "mean"),
                                   output = c("data.frame", "matrix"),
                                   label_format = c("plotmath", "plain")) {
  label_format = match.arg(label_format)
  summary = match.arg(summary)
  output = match.arg(output)
  #--------------#
  # Checks ####
  #--------------#
  if (is.null(B_diff$draws_upper)) stop("B_diff$draws_upper is missing. Use heat_out$B_diff[[k]], not heat_out$B_diff_sum[[k]].")
  draws_upper = B_diff$draws_upper
  if (is.null(dim(draws_upper))) draws_upper = matrix(draws_upper, nrow = 1)
  draws_upper = as.matrix(draws_upper)
  storage.mode(draws_upper) = "double"
  n_edge = ncol(draws_upper)
  if (is.null(B_diff$n)) {
    n = (1 + sqrt(1 + 8 * n_edge)) / 2
    if (n != as.integer(n)) stop("Cannot infer n from ncol(B_diff$draws_upper).")
    n = as.integer(n)
  } else {
    n = as.integer(B_diff$n)
  }
  if (n < 2L) stop("At least two nodes are needed.")
  if (is.null(B_diff$upper_idx)) {
    upper_idx = upper.tri(matrix(FALSE, n, n), diag = FALSE)
  } else {
    upper_idx = B_diff$upper_idx
  }
  if (!all(dim(upper_idx) == c(n, n))) stop("B_diff$upper_idx must be n x n.")
  if (sum(upper_idx) != n_edge) stop("ncol(B_diff$draws_upper) must equal sum(B_diff$upper_idx).")
  if (is.null(node_labels)) node_labels = as.character(seq_len(n))
  if (length(node_labels) != n) stop("node_labels must have length n.")
  node_labels = as.character(node_labels)
  diff_label = if (is.null(B_diff$label)) "Bdiff" else B_diff$label
  if (label_format == "plotmath") {
    diff_label = as_Bdiff_plotmath_label(diff_label)
  }
  #-------------------------------#
  # L2 norm per posterior draw ####
  #-------------------------------#
  upper_ij = which(upper_idx, arr.ind = TRUE)
  E = matrix(0, nrow = n_edge, ncol = n)
  E[cbind(seq_len(n_edge), upper_ij[, 1])] = 1
  E[cbind(seq_len(n_edge), upper_ij[, 2])] = 1
  L = sqrt((draws_upper^2) %*% E)
  l2 = if (summary == "median") apply(L, 2, stats::median) else colMeans(L)
  #------------#
  # Output ####
  #------------#
  out = data.frame(diff = diff_label, node = seq_len(n), node_name = node_labels, l2_norm = as.numeric(l2), row.names = NULL)
  out = out[order(-out$l2_norm), ]
  out$rank = seq_len(nrow(out))
  out = out[, c("diff", "rank", "node", "node_name", "l2_norm")]
  rownames(out) = NULL
  return(if (output == "matrix") as.matrix(out) else out)
}
