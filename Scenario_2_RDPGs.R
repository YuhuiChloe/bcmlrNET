# ============================================================
# Simulation Scenario 2:
# Multi-D0 random dot product graph with piecewise-constant
# latent positions
# ============================================================

sim_scenario2_rdpg <- function(
    n = 50,
    Tt = 300,
    D0 = 4,
    cp = c(75, 225),
    alpha_x = rep(1, D0),
    delta_x = 0.12,
    seed = NULL
) {
  # --------------#
  # Sanity checks #
  # --------------#
  if (!is.null(seed)) set.seed(seed)
  if (length(cp) != 2) stop("cp must be a length-2 vector.")
  if (D0 < 2) stop("D0 must be at least 2.")
  if (length(alpha_x) != D0) stop("alpha_x must have length D0.")
  if (any(alpha_x <= 0)) stop("All entries of alpha_x must be positive.")
  if (delta_x <= 0 || delta_x >= 1) stop("delta_x must satisfy 0 < delta_x < 1.")
  if (cp[1] < 1 || cp[1] >= cp[2] || cp[2] >= Tt) {
    stop("Require 1 <= cp[1] < cp[2] < Tt.")
  }
  # ------------------------------------------------#
  # Simple Dirichlet sampler using Gamma draws
  # ------------------------------------------------#
  rdirichlet_simple = function(n_draws, alpha) {
    out = matrix(0, nrow = n_draws, ncol = length(alpha))
    for (r in 1:n_draws) {
      g = rgamma(length(alpha), shape = alpha, rate = 1)
      out[r, ] = g / sum(g)
    }
    return(out)
  }
  # ------------------------------------------------#
  # Baseline latent positions x_i^*
  # ------------------------------------------------#
  X_base = rdirichlet_simple(n, alpha_x)
  # Standard basis vectors e1, e2 in R^D0
  e1 = c(1, rep(0, D0 - 1))
  e2 = c(0, 1, rep(0, D0 - 2))
  # ------------------------------------------------#
  # Segment-specific latent positions
  # ------------------------------------------------#
  X_seg1 = X_base
  
  X_seg2 = X_base
  for (i in 1:5) {
    X_seg2[i, ] = (1 - delta_x) * X_base[i, ] + delta_x * e1
  }
  for (i in (n - 4):n) {
    X_seg2[i, ] = (1 - delta_x) * X_base[i, ] + delta_x * e2
  }
  
  X_seg3 = X_base
  # ------------------------------------------------#
  # Helper to get X_t by segment
  # ------------------------------------------------#
  get_X_t = function(t) {
    if (t <= cp[1]) {
      X_seg1
    } else if (t <= cp[2]) {
      X_seg2
    } else {
      X_seg3
    }
  }
  # ------------------------------------------------#
  # Simulate adjacency matrices:
  # A_{ij,t} ~ Bernoulli(x_{i,t}^T x_{j,t})
  # ------------------------------------------------#
  A = vector("list", Tt)
  for (t in 1:Tt) {
    X_t = get_X_t(t)
    A_t = matrix(0, n, n)
    
    for (i in 1:(n - 1)) {
      for (j in (i + 1):n) {
        prob_ij = sum(X_t[i, ] * X_t[j, ])
        prob_ij = min(max(prob_ij, 0), 1)
        A_t[i, j] = rbinom(1, size = 1, prob = prob_ij)
        A_t[j, i] = A_t[i, j]
      }
    }
    diag(A_t) = 0
    A[[t]] = A_t
  }
  return(list(
    A = A,
    cp = cp,
    D0 = D0,
    delta_x = delta_x,
    alpha_x = alpha_x,
    basis = list(e1 = e1, e2 = e2),
    X = list(
      baseline = X_base,
      segment1 = X_seg1,
      segment2 = X_seg2,
      segment3 = X_seg3
    )
  ))
}
