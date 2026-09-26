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
    alpha_x = rep(1, D0 + 1),
    delta_x = 0.1,
    num_shift = 5,
    seed = NULL
) {
  # --------------#
  # Sanity checks #
  # --------------#
  if (!is.null(seed)) set.seed(seed)
  if (length(cp) != 2) stop("cp must be a length-2 vector.")
  if (D0 < 2) stop("D0 must be at least 2.")
  if (length(alpha_x) != D0 + 1) {
    stop("alpha_x must have length D0 + 1 as the D0th coordinate will be discarded for the retained latent vectors remain in the sub-simplex.")
  }
  if (any(alpha_x <= 0)) stop("All entries of alpha_x must be positive.")

  if (num_shift > n / 2) {
    stop("The number of nodes with shifts in latent dimensions cannot be larger than half of the total number of nodes.")
  }
  if (num_shift < 0) stop("num_shift must be nonnegative.")
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

  # Helper for checking that retained latent positions
  # lie in the D0-dimensional sub-simplex:
  # x_id >= 0 and sum_d x_id <= 1.
  check_subsimplex = function(X, name, tol = 1e-12) {
    if (any(!is.finite(X))) stop(paste0(name, " contains non-finite entries."))
    if (any(X < -tol)) stop(paste0(name, " contains negative entries."))
    if (any(rowSums(X) > 1 + tol)) {
      stop(paste0(name, " has row sums greater than 1; latent dot products may be invalid."))
    }
  }

  # ------------------------------------------------#
  # Baseline latent positions x_i^*
  # ------------------------------------------------#
  # Draw a full (D0 + 1)-dimensional Dirichlet vector.
  X_full_base = rdirichlet_simple(n, alpha_x)

  # Keep only the first D0 coordinates. Do NOT renormalize.
  # The discarded last coordinate is slack mass.
  X_base = X_full_base[, 1:D0, drop = FALSE]
  X_slack = X_full_base[, D0 + 1]
  check_subsimplex(X_base, "X_base")

  # Standard basis vectors e1, e2 in R^D0
  e1 = c(1, rep(0, D0 - 1))
  e2 = c(0, 1, rep(0, D0 - 2))

  # ------------------------------------------------#
  # Segment-specific latent positions
  # ------------------------------------------------#
  X_seg1 = X_base

  X_seg2 = X_base
  if (num_shift > 0) {
    for (i in seq_len(num_shift)) {
      X_seg2[i, ] = (1 - delta_x) * X_base[i, ] + delta_x * e1
    }
    for (i in (n - num_shift + 1):n) {
      X_seg2[i, ] = (1 - delta_x) * X_base[i, ] + delta_x * e2
    }
  }
  check_subsimplex(X_seg2, "X_seg2")

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

        # >>>  This should already be in [0, 1] because X_t
        # >>>  is in the sub-simplex. Stop if there is a real
        # >>>  violation; only clamp tiny numerical error.
        if (!is.finite(prob_ij) || prob_ij < -1e-12 || prob_ij > 1 + 1e-12) {
          stop("Invalid Bernoulli probability from x_i^T x_j. Check alpha_x, D0, and transformations of X.")
        }
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
    # Return the full Dirichlet draws and discarded slack mass
    # so the construction can be inspected.
    X_full = X_full_base,
    X_slack = X_slack,
    basis = list(e1 = e1, e2 = e2),
    X = list(
      baseline = X_base,
      segment1 = X_seg1,
      segment2 = X_seg2,
      segment3 = X_seg3
    )
  ))
}
