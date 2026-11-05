# POLECAT workflow from Chapter 4. Run after installing bcmlrNET.
# These are analysis settings, not a short software test.
library(bcmlrNET)
data("POLECAT", package = "bcmlrNET")
A_list <- POLECAT$A_list
node_labels <- POLECAT$node_labels
time_labels <- POLECAT$time_labels
date_list <- POLECAT$date_list
S <- compute_netstats(A_list, network_stats = c("edges", "triangles", "isolates"))

# Select the number of changepoints and refit if it is positive.
fit <- bcmlrNET_model_select(
  A = A_list, S = S, max_num_cp = 10, min_size = 10,
  eta_dir = 1.01, b_tau = 0.5,
  num_iter = 10000, num_warmup = 5000, seed = 2026
)
fit$num_cp
fit$probs_num_cp
print(plot_histogram(fit$num_cp_dist, max_num = 10))

if (fit$num_cp > 0) {
  print(time_labels[fit$Kappa_mode])
  print(plot_changepoints(fit$Kappa, S = S, dates = date_list))
  print(plot_beta_diff_CIs(
    fit$Beta, mfrow = c(fit$num_cp, 1),
    coef_labels = c("Edges", "Triangles", "Isolates")
  ))
  heat_out <- plot_B_diff_heatmaps(
    fit, point_summary = "median", common_zlim = TRUE,
    zlim_q = 1, show_sig = FALSE
  )
  plot_node_involvement_L2(heat_out, node_labels)
}

# Optional alternative: fit a prespecified number of changepoints.
# fixed_fit <- bcmlrNET(
#   A = A_list, S = S, num_CP = 3, min_size = 10,
#   num_iter = 10000, num_warmup = 5000, seed = 2026
# )
# fixed_fit$Kappa_mode
