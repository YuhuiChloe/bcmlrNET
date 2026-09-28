# ============================================================#
# Simulation Scenario 4:
# Mixed membership stochastic block model with fixed
# mixed memberships, changing community connection probabilities,
# and dyad-level temporal dependence
# ============================================================#

sim_scenario4_mmsbm <- function(n = 50, 
                                Tt = 300, 
                                C = 4, 
                                cp = c(75, 225), 
                                rho = 0.1, 
                                p_range = c(0.4, 0.5), 
                                q_range = c(0.2, 0.3), 
                                group_sizes = c(10, 10, 15, 15),
                                alpha_main = 8,
                                alpha_other = 1,
                                delta_within = 0.2,
                                delta_between = 0.2,
                                seed = NULL) {
  # --------------#
  # Sanity checks #
  # --------------#
  if (!is.null(seed)) set.seed(seed)
  if (length(cp) != 2) stop("cp must be a length-2 vector.")
  if (C < 2) stop("C must be at least 2.")
  if (length(group_sizes) != C) stop("group_sizes must have length C.")
  if (any(group_sizes <= 0) || any(group_sizes != floor(group_sizes))) stop("group_sizes must contain positive integers.")
  if (sum(group_sizes) != n) stop("The entries of group_sizes must sum to n.")
  if (alpha_main <= 0 || alpha_other <= 0) stop("alpha_main and alpha_other must be positive.")
  if (delta_within < 0 || delta_between < 0) stop("delta_within and delta_between must be nonnegative.")
  if (rho < 0 || rho >= 1) stop("rho must satisfy 0 <= rho < 1.")
  if (cp[1] < 1 || cp[1] >= cp[2] || cp[2] >= Tt) {
    stop("Require 1 <= cp[1] < cp[2] < Tt.")
  }
  if (any(p_range < 0) || any(p_range > 1) || any(q_range < 0) || any(q_range > 1)) {
    stop("p_range and q_range must be contained in [0, 1].")
  }
  if (p_range[1] > p_range[2] || q_range[1] > q_range[2]) {
    stop("The lower endpoint of each probability range must be no larger than the upper endpoint.")
  }
  if (p_range[1] - delta_within < 0) stop("p_range[1] - delta_within must be nonnegative.")
  if (q_range[2] + delta_between > 1) stop("q_range[2] + delta_between must be no larger than 1.")
  # ---------------------------------------------------------#
  # Draw p and q once, then change B across segments
  # ---------------------------------------------------------#
  p_in = runif(1, min = p_range[1], max = p_range[2]) # within-community connection probability
  q_out = runif(1, min = q_range[1], max = q_range[2]) # between-community connection probability
  
  # Helper to construct a C by C community connection probability matrix
  make_B = function(within_prob, between_prob) {
    B = matrix(between_prob, nrow = C, ncol = C)
    diag(B) = within_prob
    return(B)
  }
  
  B1 = make_B(within_prob = p_in,
              between_prob = q_out)
  B2 = make_B(within_prob = p_in - delta_within,
              between_prob = q_out)
  B3 = make_B(within_prob = p_in - delta_within,
              between_prob = q_out + delta_between)
  
  if (any(B1 < 0) || any(B1 > 1) ||
      any(B2 < 0) || any(B2 > 1) ||
      any(B3 < 0) || any(B3 > 1)) {
    stop("All entries of B1, B2, and B3 must be contained in [0, 1].")
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
  # --------------------------------------------------#
  # Helper to keep each mixed membership vector on the simplex
  # --------------------------------------------------#
  normalize_rows = function(Z) {
    Z = Z / rowSums(Z)
    return(Z)
  }
  # -----------------------------------------#
  # Build one fixed mixed membership matrix Z
  # -----------------------------------------#
  # nodes in the same group mostly belong to the same community, but memberships are not absolute
  z_primary = rep(1:C, times = group_sizes)
  Z = matrix(0, nrow = n, ncol = C)
  for (i in 1:n) {
    alpha_i = rep(alpha_other, C)
    # if node i is assigned to primary community c, its alpha_i has alpha_main in entry c
    alpha_i[z_primary[i]] = alpha_main
    Z[i, ] = rdirichlet_simple(n_draws = 1, alpha = alpha_i)
  }
  Z = normalize_rows(Z)
  
  Z1 = Z
  Z2 = Z # community memberships remain fixed after the first changepoint
  Z3 = Z # community memberships remain fixed after the second changepoint
  # ------------------------------------------#
  # Segment-specific probability matrices P^(s)
  # ------------------------------------------#
  P1 = Z %*% B1 %*% t(Z)
  P2 = Z %*% B2 %*% t(Z)
  P3 = Z %*% B3 %*% t(Z)
  
  diag(P1) = 0
  diag(P2) = 0
  diag(P3) = 0
  # ---------------------------------------------------------#
  # Helper to pick P_t by segment
  # ---------------------------------------------------------#
  # P1, P2, P3 are globally defined
  get_P_t = function(t) {
    if (t <= cp[1]) {
      return (P1)
    } else if (t <= cp[2]) {
      return (P2)
    } else {
      return (P3)
    }
  }
  # ---------------------------------------------------------#
  # Simulate undirected adjacency matrices with dyad-level
  # Markov dependence:
  # A_{ij,1} ~ Bernoulli(p_ij^(1))
  # A_{ij,t+1} | A_{ij,t}
  # ---------------------------------------------------------#
  A = vector("list", Tt)
  # t = 1, initialize A_{ij,1}
  A1 = matrix(0, n, n)
  P_init = get_P_t(1)
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      A1[i, j] = rbinom(1, size = 1, prob = P_init[i, j]) # a Bernoulli draw 
      A1[j, i] = A1[i, j] # copy to the corresponding off-diagonal entry 
    }
  }
  diag(A1) = 0 # set all diagonal entries to zero
  A[[1]] = A1
  
  # t = 2, ..., Tt
  for (t in 2:Tt) {
    A_prev = A[[t - 1]]
    A_curr = matrix(0, n, n)
    P_curr = get_P_t(t)
    
    for (i in 1:(n - 1)) {
      for (j in (i + 1):n) {
        if (A_prev[i, j] == 1) {
          prob_ij = rho * (1 - P_curr[i, j]) + P_curr[i, j]
        } else {
          prob_ij = (1 - rho) * P_curr[i, j]
        }
        A_curr[i, j] = rbinom(1, size = 1, prob = prob_ij)
        A_curr[j, i] = A_curr[i, j] 
      }
    }
    diag(A_curr) = 0
    A[[t]] = A_curr
  }
  
  return(list(
    A = A, 
    cp = cp,
    C = C,
    group_sizes = group_sizes,
    B = list(B1 = B1, B2 = B2, B3 = B3),
    p_in = p_in,
    q_out = q_out,
    alpha = list(alpha_main = alpha_main, alpha_other = alpha_other),
    delta = list(delta_within = delta_within, delta_between = delta_between),
    z_primary = z_primary,
    Z = list(Z1 = Z1, Z2 = Z2, Z3 = Z3),
    P = list(P1 = P1, P2 = P2, P3 = P3)
  ))
}
