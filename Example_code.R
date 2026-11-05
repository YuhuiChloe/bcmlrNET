# Real-data examples: POLECAT and Enron
# Run after installing bcmlrNET, including both packaged datasets.
# Run either dataset section independently after loading the package below.
# Each example uses 10000 iterations with 5000 warmup iterations. These are
# analysis settings, not a short software test; assess convergence before
# interpreting posterior summaries. The two fits use separate object names.
# B-difference heatmaps retain posterior dyad draws and can require substantial
# memory, especially for the 197-node POLECAT series.

library(bcmlrNET)

# ============================================================
# 1. POLECAT: 184 weekly networks on 197 countries
# ============================================================

data("POLECAT", package = "bcmlrNET")
A_POLECAT <- POLECAT$A_list
node_labels_POLECAT <- POLECAT$node_labels
time_labels_POLECAT <- POLECAT$time_labels
date_list_POLECAT <- POLECAT$date_list

length(A_POLECAT)           # 184
dim(A_POLECAT[[1]])         # 197 197
length(node_labels_POLECAT) # 197

S_POLECAT <- compute_netstats(
  A_POLECAT,
  network_stats = c("edges", "triangles", "isolates")
)
dim(S_POLECAT)              # 184 3

# Select the number of changepoints and refit if it is positive.
# These preserve the settings in the original POLECAT example.
fit_POLECAT <- bcmlrNET_model_select(
  A = A_POLECAT, S = S_POLECAT,
  max_num_cp = 10, min_size = 10,
  eta_dir = 1.01, b_tau = 0.5,
  num_iter = 10000, num_warmup = 5000,
  seed = 2026, print_progress = TRUE
)
print(fit_POLECAT$num_cp)
print(fit_POLECAT$probs_num_cp)
print(fit_POLECAT$auc_diagnostics)

# Posterior distribution of the number of supported changepoints.
print(plot_histogram(fit_POLECAT$num_cp_dist, max_num = 10))

# Location and coefficient plots require a positive selected number of CPs.
if (fit_POLECAT$num_cp > 0) {
  print(time_labels_POLECAT[fit_POLECAT$Kappa_mode])

  # Network statistics with changepoint locations and credible intervals.
  # plot_changepoints() calls plot_netstats_cp_cis(); list() means that no
  # competing-method estimates are supplied.
  print(plot_changepoints(
    fit_POLECAT$Kappa,
    S = S_POLECAT, dates = date_list_POLECAT,
    kappa_mode = fit_POLECAT$Kappa_mode,
    other_cps = list(),
    stat_names = c("Edges", "Triangles", "Isolates"),
    title = "POLECAT"
  ))

  # Credible intervals for changes in network-statistic coefficients.
  print(plot_beta_diff_CIs(
    fit_POLECAT$Beta, mfrow = c(fit_POLECAT$num_cp, 1),
    coef_labels = c("Edges", "Triangles", "Isolates")
  ))

  # Dyad-level differences and node involvement.
  heat_POLECAT <- plot_B_diff_heatmaps(
    fit_POLECAT, point_summary = "median", common_zlim = TRUE,
    zlim_q = 1, show_sig = FALSE
  )
  plot_node_involvement_L2(heat_POLECAT, node_labels_POLECAT)
} else {
  message("POLECAT: zero changepoints selected; location and difference plots are skipped.")
}

# Optional alternative: fit a prespecified number of changepoints.
# Uncomment to request three changepoints instead of selecting the number.
# fixed_fit_POLECAT <- bcmlrNET(
#   A = A_POLECAT, S = S_POLECAT, num_CP = 3, min_size = 10,
#   prior_for_U = "Spike_Slab",
#   num_iter = 10000, num_warmup = 5000, seed = 2026
# )
# time_labels_POLECAT[fixed_fit_POLECAT$Kappa_mode]

# ============================================================
# 2. Enron: 100 weekly networks on 89 employees
# ============================================================

data("Enron", package = "bcmlrNET")
A_Enron <- Enron$A_list
node_labels_Enron <- Enron$node_labels
time_labels_Enron <- Enron$time_labels
date_list_Enron <- Enron$date_list

length(A_Enron)           # 100
dim(A_Enron[[1]])         # 89 89
length(node_labels_Enron) # 89
time_labels_Enron[c(1, 100)]
# "2000-06-05/2000-06-11" "2002-04-29/2002-05-05"
# Labels include both the Monday start and Sunday end; date_list_Enron
# contains the Monday dates used on the plots.

S_Enron <- compute_netstats(
  A_Enron,
  network_stats = c("edges", "triangles", "isolates")
)
dim(S_Enron)              # 100 3

# Select the number of changepoints using the supplied Enron analysis settings.
# With min_size = 6 and thinning = 5, some segments may have too few evaluation
# observations for AUC intervals. Inspect the printed auc_diagnostics below.
fit_Enron <- bcmlrNET_model_select(
  A = A_Enron, S = S_Enron,
  max_num_cp = 8, min_size = 6,
  latent_dim = 10, eta_dir = 1.01, b_tau = 0.5,
  thinning = 5, prior_kappa = "default",
  alpha_f = 0.05, threshold = 0.5, sd_beta = 3,
  num_iter = 10000, num_warmup = 5000,
  prior_for_U = "Spike_Slab", lambda_init = "equal",
  sparse_networks = FALSE, center_std_A = TRUE,
  seed = 2026, print_progress = TRUE, print_outputs = FALSE
)
print(fit_Enron$num_cp)
print(fit_Enron$probs_num_cp)
print(fit_Enron$auc_diagnostics)

# Posterior distribution of the number of supported changepoints.
print(plot_histogram(fit_Enron$num_cp_dist, max_num = 8))

if (fit_Enron$num_cp > 0) {
  print(time_labels_Enron[fit_Enron$Kappa_mode])

  # Network statistics and changepoint intervals, without competing methods.
  print(plot_changepoints(
    fit_Enron$Kappa,
    S = S_Enron, dates = date_list_Enron,
    kappa_mode = fit_Enron$Kappa_mode,
    other_cps = list(),
    stat_names = c("Edges", "Triangles", "Isolates"),
    title = "Enron"
  ))

  # Credible intervals for changes in network-statistic coefficients.
  print(plot_beta_diff_CIs(
    fit_Enron$Beta, mfrow = c(fit_Enron$num_cp, 1),
    coef_labels = c("Edges", "Triangles", "Isolates")
  ))

  # Dyad-level differences and node involvement.
  heat_Enron <- plot_B_diff_heatmaps(
    fit_Enron, point_summary = "median", common_zlim = TRUE,
    zlim_q = 1, show_sig = FALSE
  )
  plot_node_involvement_L2(heat_Enron, node_labels_Enron)
} else {
  message("Enron: zero changepoints selected; location and difference plots are skipped.")
}

# Optional alternative: fit a prespecified number of changepoints.
# For illustration, uncomment to request four changepoints.
# fixed_fit_Enron <- bcmlrNET(
#   A = A_Enron, S = S_Enron, num_CP = 4, min_size = 6,
#   prior_for_U = "Spike_Slab",
#   num_iter = 10000, num_warmup = 5000, seed = 2026
# )
# time_labels_Enron[fixed_fit_Enron$Kappa_mode]

# Dataset documentation:
# help("POLECAT", package = "bcmlrNET")
# help("Enron", package = "bcmlrNET")
