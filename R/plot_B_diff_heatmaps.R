# ============================================================
# Colors and color-scale helpers
# ============================================================

B_heatmap_cols = function(n_col = 301) {
  grDevices::colorRampPalette(
    c("#053061", "#2166AC", "#F7F7F7", "#B2182B", "#67001F")
  )(n_col)
}

get_Bdiff_zlim = function(summary_list, q = 0.99, nonzero_only = TRUE) {
  vals = unlist(
    lapply(
      summary_list,
      function(obj) obj$point[upper.tri(obj$point, diag = FALSE)]
    ),
    use.names = FALSE
  )

  vals = abs(vals[is.finite(vals)])
  if (nonzero_only) vals = vals[vals > 0]
  if (length(vals) == 0L) return(1)

  zlim = as.numeric(
    stats::quantile(vals, probs = q, names = FALSE, type = 8)
  )

  if (!is.finite(zlim) || zlim <= 0) {
    zlim = max(vals, na.rm = TRUE)
  }

  if (!is.finite(zlim) || zlim <= 0) {
    zlim = 1
  }

  zlim
}


# ============================================================
# Color bars
# ============================================================

plot_B_colorbar_bottom = function(zlim, cols, n_ticks = 5) {
  if (!is.finite(zlim) || zlim <= 0) {
    stop("zlim must be positive.")
  }

  graphics::par(
    mar = c(2.6, 0.3, 2.2, 3.2),
    pty = "m",
    xaxs = "i",
    yaxs = "i",
    family = "Times",
    font.axis = 1,
    font.lab = 1,
    cex.axis = 1.10
  )

  graphics::plot.new()
  graphics::plot.window(
    xlim = c(-zlim, zlim),
    ylim = c(0, 1)
  )

  brks = seq(
    -zlim,
    zlim,
    length.out = length(cols) + 1L
  )

  for (ii in seq_along(cols)) {
    graphics::rect(
      brks[ii],
      0.35,
      brks[ii + 1L],
      0.75,
      col = cols[ii],
      border = NA
    )
  }

  ticks = pretty(
    c(-zlim, zlim),
    n = n_ticks
  )

  ticks = ticks[
    ticks >= -zlim &
      ticks <= zlim
  ]

  graphics::axis(
    1,
    at = ticks,
    labels = formatC(
      ticks,
      digits = 2,
      format = "fg"
    ),
    cex.axis = 1.10
  )

  graphics::box()
}


plot_B_colorbar = function(zlim, cols, n_ticks = 7) {
  if (!is.finite(zlim) || zlim <= 0) {
    stop("zlim must be positive.")
  }

  graphics::par(
    mar = c(2.6, 0.4, 2.2, 4.5),
    pty = "m",
    xaxs = "i",
    yaxs = "i",
    family = "Times",
    font.axis = 1,
    font.lab = 1,
    cex.axis = 1.15
  )

  graphics::plot.new()
  graphics::plot.window(
    xlim = c(0, 1),
    ylim = c(-zlim, zlim)
  )

  brks = seq(
    -zlim,
    zlim,
    length.out = length(cols) + 1L
  )

  for (ii in seq_along(cols)) {
    graphics::rect(
      0.20,
      brks[ii],
      0.55,
      brks[ii + 1L],
      col = cols[ii],
      border = NA
    )
  }

  ticks = pretty(
    c(-zlim, zlim),
    n = n_ticks
  )

  ticks = ticks[
    ticks >= -zlim &
      ticks <= zlim
  ]

  graphics::axis(
    4,
    at = ticks,
    labels = formatC(
      ticks,
      digits = 2,
      format = "fg"
    ),
    las = 1,
    cex.axis = 1.15
  )

  graphics::box()
}


# ============================================================
# Construct B matrices from posterior draws
# ============================================================

B_upper_from_U_lambda = function(U, lambda, upper_idx) {
  active = which(lambda != 0)

  if (length(active) == 0L) {
    return(numeric(sum(upper_idx)))
  }

  U_active = U[, active, drop = FALSE]
  lambda_active = lambda[active]

  B = tcrossprod(
    sweep(
      U_active,
      2,
      lambda_active,
      `*`
    ),
    U_active
  )

  B[upper_idx]
}


B_upper_draws_from_out = function(
    out,
    k,
    U_name = "U",
    Lambda_name = "Lambda") {

  Uarr = out[[U_name]]
  Larr = out[[Lambda_name]]

  if (is.null(Larr) && Lambda_name == "Lambda") {
    Larr = out[["lambda"]]
  }

  if (is.null(Uarr)) {
    stop("out$U is missing.")
  }

  if (is.null(Larr)) {
    stop("out$Lambda, or out$lambda, is missing.")
  }

  dU = dim(Uarr)
  dL = dim(Larr)

  if (length(dU) != 4L) {
    stop("out$U must have dimension num_save x K x n x D.")
  }

  if (length(dL) != 3L) {
    stop("out$Lambda must have dimension num_save x K x D.")
  }

  num_save = dU[1]
  K = dU[2]
  n = dU[3]
  D = dU[4]

  if (!all(dL == c(num_save, K, D))) {
    stop("Dimensions of out$U and out$Lambda are incompatible.")
  }

  if (k < 1 || k > K) {
    stop("k must be between 1 and K.")
  }

  upper_idx = upper.tri(
    matrix(FALSE, n, n),
    diag = FALSE
  )

  draws_upper = matrix(
    NA_real_,
    nrow = num_save,
    ncol = sum(upper_idx)
  )

  for (s in seq_len(num_save)) {
    U_sk = matrix(
      Uarr[s, k, , ],
      nrow = n,
      ncol = D
    )

    lambda_sk = as.numeric(
      Larr[s, k, ]
    )

    draws_upper[s, ] = B_upper_from_U_lambda(
      U_sk,
      lambda_sk,
      upper_idx
    )
  }

  list(
    draws_upper = draws_upper,
    upper_idx = upper_idx,
    n = n,
    k = k,
    label = paste0("B", k)
  )
}


# ============================================================
# Posterior summaries
# ============================================================

summarize_upper_draws = function(
    obj,
    point_summary = c("median", "mean"),
    ci_probs = c(0.025, 0.975)) {

  point_summary = match.arg(point_summary)

  draws_upper = obj$draws_upper
  upper_idx = obj$upper_idx
  n = obj$n

  make_mat = function(v) {
    M = matrix(0, n, n)
    M[upper_idx] = v
    M = M + t(M)
    diag(M) = 0
    M
  }

  if (point_summary == "median") {
    if (requireNamespace("matrixStats", quietly = TRUE)) {
      point_vec = matrixStats::colMedians(
        draws_upper
      )
    } else {
      point_vec = apply(
        draws_upper,
        2,
        stats::median
      )
    }
  } else {
    point_vec = colMeans(
      draws_upper
    )
  }

  point_mat = make_mat(
    point_vec
  )

  if (is.null(ci_probs)) {
    return(
      list(
        point = point_mat,
        mean = point_mat,
        lower = NULL,
        upper = NULL,
        sig = NULL,
        point_summary = point_summary,
        label = obj$label,
        n = n
      )
    )
  }

  if (requireNamespace("matrixStats", quietly = TRUE)) {
    qmat = matrixStats::colQuantiles(
      draws_upper,
      probs = ci_probs,
      drop = FALSE
    )
  } else {
    qmat = t(
      apply(
        draws_upper,
        2,
        stats::quantile,
        probs = ci_probs,
        names = FALSE
      )
    )
  }

  lower_vec = qmat[, 1]
  upper_vec = qmat[, 2]

  sig_vec = (
    lower_vec > 0
  ) | (
    upper_vec < 0
  )

  sig_mat = matrix(
    FALSE,
    n,
    n
  )

  sig_mat[upper_idx] = sig_vec
  sig_mat = sig_mat | t(sig_mat)
  diag(sig_mat) = FALSE

  list(
    point = point_mat,
    mean = point_mat,
    lower = make_mat(lower_vec),
    upper = make_mat(upper_vec),
    sig = sig_mat,
    point_summary = point_summary,
    ci_probs = ci_probs,
    label = obj$label,
    n = n
  )
}


# ============================================================
# Plot one heatmap
# ============================================================

plot_B_heatmap = function(
    summary_obj,
    main = NULL,
    zlim = NULL,
    cols = B_heatmap_cols(),
    clip_to_zlim = FALSE,
    show_sig = TRUE,
    moved_nodes = NULL,
    grid_every = NULL,
    axis_n_breaks = 10,
    show_y_axis = TRUE) {

  M = summary_obj$point
  sig = summary_obj$sig
  n = nrow(M)

  if (is.null(main)) {
    main = paste0(
      "Posterior ",
      summary_obj$point_summary,
      " of ",
      summary_obj$label
    )
  }

  diag(M) = 0

  if (is.null(zlim)) {
    vals = M[
      upper.tri(
        M,
        diag = FALSE
      )
    ]

    vals = vals[
      is.finite(vals)
    ]

    zlim = max(
      abs(vals),
      na.rm = TRUE
    )

    if (!is.finite(zlim) || zlim <= 0) {
      zlim = 1
    }
  }

  if (length(zlim) == 2L) {
    zlim = max(
      abs(zlim)
    )
  }

  if (!is.finite(zlim) || zlim <= 0) {
    zlim = 1
  }

  if (clip_to_zlim) {
    M = pmax(
      pmin(
        M,
        zlim
      ),
      -zlim
    )
  }

  graphics::image(
    x = seq_len(n),
    y = seq_len(n),
    z = t(
      M[
        n:1,
        ,
        drop = FALSE
      ]
    ),
    zlim = c(-zlim, zlim),
    col = cols,
    xlab = "",
    ylab = "",
    main = main,
    axes = FALSE,
    asp = 1,
    useRaster = TRUE
  )

  br = unique(
    round(
      seq(
        1,
        n,
        length.out = axis_n_breaks
      )
    )
  )

  graphics::axis(
    1,
    at = br,
    labels = br
  )

  if (show_y_axis) {
    graphics::axis(
      2,
      at = br,
      labels = n - br + 1,
      las = 2,
      line = 0
    )
  }

  graphics::box()

  if (!is.null(grid_every)) {
    gr = seq(
      0.5,
      n + 0.5,
      by = grid_every
    )

    graphics::abline(
      v = gr,
      h = gr,
      col = "grey90",
      lwd = 0.3
    )
  }

  if (!is.null(moved_nodes)) {
    moved_nodes = sort(
      unique(
        as.integer(
          moved_nodes
        )
      )
    )

    moved_nodes = moved_nodes[
      moved_nodes >= 1 &
        moved_nodes <= n
    ]

    graphics::abline(
      v = moved_nodes,
      lty = 2,
      lwd = 1
    )

    graphics::abline(
      h = n - moved_nodes + 1,
      lty = 2,
      lwd = 1
    )
  }

  if (
    show_sig &&
    !is.null(sig) &&
    any(sig, na.rm = TRUE)
  ) {
    idx = which(
      sig,
      arr.ind = TRUE
    )

    graphics::points(
      idx[, 2],
      n - idx[, 1] + 1,
      pch = 15,
      cex = 0.15
    )
  }

  invisible(zlim)
}


# ============================================================
# Compute adjacent B differences
# ============================================================

#' Compute adjacent dyad-coefficient differences from posterior draws
#'
#' @param out Fitted output containing `U` (draws by segments by nodes by latent
#'   dimensions) and `Lambda` (draws by segments by latent dimensions).
#'   The legacy field name `lambda` is also accepted.
#' @param point_summary Posterior median or mean for the displayed summaries.
#' @param ci_probs Two posterior quantile probabilities; NULL omits intervals.
#' @returns A list with `B_upper`, `B_diff`, and `B_diff_sum`. Draw objects contain
#'   `draws_upper`, `upper_idx`, and `n`. Summary objects contain `point`, `mean`,
#'   `lower`, `upper`, and `sig`, plus labels and metadata.
#' @details Reconstructs segment coefficients from `U` and `Lambda` and subtracts
#'   adjacent segments within each draw. Diagonals in summary matrices are zero.
#'   The original `mean` field is retained as an alias of `point`: it contains
#'   medians when `point_summary = "median"`. Use `point` for the selected summary.
#'   All segment and difference draws are retained, which can require substantial
#'   memory for a large network.
#' @examples
#' fit <- list(U = array(0, c(3, 2, 3, 1)), Lambda = array(1, c(3, 2, 1)))
#' for (s in 1:3) fit$U[s, 2, , 1] <- c(1, 2, 3) * s
#' heat_out <- compute_B_diffs_from_out(fit)
#' heat_out$B_diff_sum[[1]]$point
#' @export
compute_B_diffs_from_out = function(
    out,
    point_summary = c("median", "mean"),
    ci_probs = c(0.025, 0.975)) {

  point_summary = match.arg(
    point_summary
  )

  if (is.null(out$U)) {
    stop("out$U is missing.")
  }

  K = dim(out$U)[2]

  if (is.null(K) || K < 2L) {
    stop("Need at least two classes in out$U.")
  }

  B_upper = lapply(
    seq_len(K),
    function(kk) {
      B_upper_draws_from_out(
        out,
        k = kk
      )
    }
  )

  names(B_upper) = paste0(
    "B",
    seq_len(K)
  )

  B_diff = vector(
    "list",
    K - 1L
  )

  for (kk in seq_len(K - 1L)) {
    B_diff[[kk]] = B_upper[[kk + 1L]]

    B_diff[[kk]]$draws_upper =
      B_upper[[kk + 1L]]$draws_upper -
      B_upper[[kk]]$draws_upper

    B_diff[[kk]]$label = paste0(
      "B",
      kk + 1L,
      " - B",
      kk
    )
  }

  names(B_diff) = paste0(
    "B",
    2:K,
    "_minus_B",
    seq_len(K - 1L)
  )

  B_diff_sum = lapply(
    B_diff,
    summarize_upper_draws,
    point_summary = point_summary,
    ci_probs = ci_probs
  )

  names(B_diff_sum) = names(
    B_diff
  )

  list(
    B_upper = B_upper,
    B_diff = B_diff,
    B_diff_sum = B_diff_sum
  )
}


# ============================================================
# Plot B-difference heatmaps
# ============================================================

#' Plot adjacent dyad-coefficient difference heatmaps
#'
#' @inheritParams compute_B_diffs_from_out
#' @param layout Optional c(rows, columns) for the heatmaps.
#' @param output_B_diff Return the computed differences and summaries?
#' @param plot_heatmaps Draw the heatmaps? FALSE computes results without plotting.
#' @param which_heatmaps Optional indices of adjacent differences to display.
#' @param zlim Optional color-limit magnitudes: one value, one per displayed
#'   heatmap, or one per difference. Use a common magnitude with a shared colorbar.
#' @param zlim_q Quantile used to choose an automatic common color limit.
#' @param common_zlim Use a common scale? FALSE retains the paper-script default.
#' @param cols Vector of colors; defaults to the original diverging palette.
#' @param add_colorbar Include a shared or per-panel colorbar?
#' @param show_sig Mark dyads whose posterior interval excludes zero?
#' @param moved_nodes Optional node indices to mark, or one vector per difference.
#' @param grid_every Optional grid-line spacing.
#' @param axis_n_breaks Number of requested node-index axis breaks.
#' @param mar Four base-graphics margins for each heatmap.
#' @returns Invisibly, the output of [compute_B_diffs_from_out()] with plotting
#'   settings when plotting is performed; NULL when `output_B_diff = FALSE`.
#' @details `which_heatmaps` affects the display only; all differences remain in
#'   the returned object. Original fonts, colors, scales, and clipping behavior
#'   are retained. Use `common_zlim = TRUE, zlim_q = 1` to cover the full range
#'   of displayed point summaries. Values beyond narrower limits can be blank,
#'   as in the supplied paper script. This function uses base graphics.
#' @examples
#' fit <- list(U = array(0, c(3, 2, 3, 1)), Lambda = array(1, c(3, 2, 1)))
#' for (s in 1:3) fit$U[s, 2, , 1] <- c(1, 2, 3) * s
#' heat_out <- plot_B_diff_heatmaps(fit, common_zlim = TRUE, zlim_q = 1,
#'                                 show_sig = FALSE)
#' @export
plot_B_diff_heatmaps = function(
    out,
    layout = NULL,
    output_B_diff = TRUE,
    plot_heatmaps = TRUE,
    which_heatmaps = NULL,
    point_summary = c("median", "mean"),
    ci_probs = c(0.025, 0.975),
    zlim = NULL,
    zlim_q = 0.99,
    common_zlim = FALSE,
    cols = B_heatmap_cols(),
    add_colorbar = TRUE,
    show_sig = TRUE,
    moved_nodes = NULL,
    grid_every = NULL,
    axis_n_breaks = 10,
    mar = c(2.2, 1.5, 2.2, 0.05)) {

  point_summary = match.arg(
    point_summary
  )

  res = compute_B_diffs_from_out(
    out,
    point_summary = point_summary,
    ci_probs = ci_probs
  )

  if (
    !is.logical(plot_heatmaps) ||
    length(plot_heatmaps) != 1L
  ) {
    stop(
      "plot_heatmaps must be TRUE or FALSE."
    )
  }

  if (!plot_heatmaps) {
    if (output_B_diff) {
      return(
        invisible(res)
      )
    }

    return(
      invisible(NULL)
    )
  }

  all_idx = seq_along(
    res$B_diff_sum
  )

  if (is.null(which_heatmaps)) {
    plot_idx = all_idx
  } else {
    plot_idx = sort(
      unique(
        as.integer(
          which_heatmaps
        )
      )
    )

    if (
      length(plot_idx) == 0L ||
      any(!is.finite(plot_idx)) ||
      any(plot_idx < 1L) ||
      any(plot_idx > length(all_idx))
    ) {
      stop(
        "which_heatmaps must contain valid B-difference indices."
      )
    }
  }

  B_diff_sum_plot =
    res$B_diff_sum[
      plot_idx
    ]

  m = length(
    B_diff_sum_plot
  )

  if (is.null(layout)) {
    ncol_plot = ceiling(
      sqrt(m)
    )

    nrow_plot = ceiling(
      m / ncol_plot
    )

    layout = c(
      nrow_plot,
      ncol_plot
    )
  }

  if (length(layout) != 2L) {
    stop(
      "layout must be c(nrow, ncol)."
    )
  }

  if (prod(layout) < m) {
    stop(
      "layout does not have enough panels."
    )
  }

  if (!is.null(zlim)) {

    if (length(zlim) == 1L) {

      zlim_each = rep(
        abs(zlim),
        m
      )

    } else if (length(zlim) == m) {

      zlim_each = abs(
        zlim
      )

    } else if (
      length(zlim) ==
      length(res$B_diff_sum)
    ) {

      zlim_each = abs(
        zlim[
          plot_idx
        ]
      )

    } else {

      stop(
        paste(
          "zlim must be NULL, a scalar,",
          "one value per plotted heatmap,",
          "or one value per B difference."
        )
      )
    }

  } else if (common_zlim) {

    zlim_common = get_Bdiff_zlim(
      B_diff_sum_plot,
      q = zlim_q,
      nonzero_only = TRUE
    )

    zlim_each = rep(
      zlim_common,
      m
    )

  } else {

    zlim_each = vapply(
      B_diff_sum_plot,
      function(obj) {
        z = max(
          abs(obj$point),
          na.rm = TRUE
        )

        if (!is.finite(z)) {
          NA_real_
        } else {
          z
        }
      },
      numeric(1)
    )

    zlim_each[
      zlim_each == 0
    ] = 0.01
  }

  zlim_each[
    !is.finite(zlim_each) |
      zlim_each <= 0
  ] = 1

  old_par = graphics::par(
    no.readonly = TRUE
  )
  graphics::par(oma = c(0, 0.7, 0, 0))

  on.exit(
    {
      graphics::layout(1)
      graphics::par(old_par)
    },
    add = TRUE
  )

  nrow_plot = layout[1]
  ncol_plot = layout[2]
  n_panel = prod(layout)

  if (
    add_colorbar &&
    common_zlim
  ) {

    panel_ids = matrix(
      seq_len(n_panel),
      nrow = nrow_plot,
      byrow = TRUE
    )

    layout_mat = cbind(
      panel_ids,
      rep(
        n_panel + 1L,
        nrow_plot
      )
    )

    graphics::layout(
      layout_mat,
      widths = c(
        rep(
          1,
          ncol_plot
        ),
        0.28
      ),
      heights = rep(
        1,
        nrow_plot
      )
    )

  } else if (
    add_colorbar &&
    !common_zlim
  ) {

    layout_mat = matrix(
      NA_integer_,
      nrow = 2L * nrow_plot,
      ncol = ncol_plot
    )

    id = 1L

    for (rr in seq_len(nrow_plot)) {
      for (cc in seq_len(ncol_plot)) {

        layout_mat[
          2L * rr - 1L,
          cc
        ] = id

        layout_mat[
          2L * rr,
          cc
        ] = id + 1L

        id = id + 2L
      }
    }

    graphics::layout(
      layout_mat,
      widths = rep(
        1,
        ncol_plot
      ),
      heights = rep(
        c(1, 0.16),
        nrow_plot
      )
    )

  } else {

    graphics::par(
      mfrow = layout
    )
  }


  # ==========================================================
  # Plot heatmaps
  #
  # IMPORTANT:
  # The exact same margin is used for every heatmap.
  # This keeps all heatmaps exactly the same size.
  #
  # Only the leftmost heatmap in each row displays y labels.
  # ==========================================================

  for (jj in seq_len(m)) {

    kk = plot_idx[jj]

    moved_kk =
      if (is.list(moved_nodes)) {
        moved_nodes[[kk]]
      } else {
        moved_nodes
      }

    show_y =
      ((jj - 1L) %% ncol_plot == 0L)

    graphics::par(
      mar = mar,
      pty = "s",
      xaxs = "i",
      yaxs = "i",
      mgp = c(1.5, 0.35, 0),
      tcl = -0.1,
      family = "Times",
      font.axis = 1,
      font.lab = 1,
      font.main = 2,
      cex.axis = 1.15,
      cex.lab = 1.20,
      cex.main = 1.5
    )

    plot_B_heatmap(
      summary_obj =
        res$B_diff_sum[[kk]],

      main = bquote(
        "Posterior" ~
          .(point_summary) ~
          "of" ~
          B[.(kk + 1)] -
          B[.(kk)]
      ),

      zlim = zlim_each[jj],
      cols = cols,
      clip_to_zlim = FALSE,
      show_sig = show_sig,
      moved_nodes = moved_kk,
      grid_every = grid_every,
      axis_n_breaks = axis_n_breaks,
      show_y_axis = show_y
    )

    if (
      add_colorbar &&
      !common_zlim
    ) {

      plot_B_colorbar_bottom(
        zlim_each[jj],
        cols
      )
    }
  }


  # ==========================================================
  # Blank unused panels
  # ==========================================================

  if (n_panel > m) {

    for (
      ii in seq_len(
        n_panel - m
      )
    ) {

      graphics::plot.new()

      if (
        add_colorbar &&
        !common_zlim
      ) {
        graphics::plot.new()
      }
    }
  }


  # ==========================================================
  # Shared colorbar
  # ==========================================================

  if (
    add_colorbar &&
    common_zlim
  ) {

    plot_B_colorbar(
      zlim_each[1],
      cols
    )
  }


  # ==========================================================
  # Output
  # ==========================================================

  names(zlim_each) =
    names(
      res$B_diff_sum
    )[plot_idx]

  res$zlim_each =
    zlim_each

  res$which_heatmaps =
    plot_idx

  res$common_zlim =
    common_zlim

  res$zlim_q =
    zlim_q

  res$cols =
    cols

  res$layout =
    layout

  if (output_B_diff) {
    return(
      invisible(res)
    )
  }

  invisible(NULL)
}


# ============================================================
# Enron example
# ============================================================

# heat_out = plot_B_diff_heatmaps(
#   out = bcmlrNet_Enron,
#   which_heatmaps = 1:3,
#   layout = c(1, 3),
#   point_summary = "median",
#   zlim_q = 1,
#   common_zlim = TRUE,
#   add_colorbar = TRUE,
#   output_B_diff = TRUE
# )
