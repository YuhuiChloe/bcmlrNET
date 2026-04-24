# ============================================================#
# Simulation Scenario 4:
# Mixed membership stochastic block model with piecewise-constant
# mixed memberships and dyad-level temporal dependence
# ============================================================#
sim_scenario4_mmsbm <- function(n = 50, 
                                Tt = 300, 
                                C = 2, 
                                cp = c(75, 225), 
                                rho = 0.5, 
                                p_range = c(0.5, 1.0), 
                                q_range = c(0.0, 0.4), 
                                num_swap = 5,
                                alpha_main = 8,
                                alpha_other = 1,
                                delta_z = 0.85,
                                seed = NULL) {
  # --------------#
  # Sanity checks #
  # --------------#
  if (!is.null(seed)) set.seed(seed)
  if (length(cp) != 2) stop("cp must be a length-2 vector.")
  if (C != 2) stop("This implementation assumes C = 2.")
  if (n %% 2 != 0) stop("n must be even for the default baseline split.")
  if (num_swap > n/2) stop("The number of nodes swapping cannot be larger than half of the total number of nodes.")
  if (alpha_main <= 0 || alpha_other <= 0) stop("alpha_main and alpha_other must be positive.")
  if (delta_z <= 0 || delta_z >= 1) stop("delta_z must satisfy 0 < delta_z < 1.")
  if (rho < 0 || rho >= 1) stop("rho must satisfy 0 <= rho < 1.")
  if (cp[1] < 1 || cp[1] >= cp[2] || cp[2] >= Tt) {
    stop("Require 1 <= cp[1] < cp[2] < Tt.")
  }
  # ------------------------------#
  # Draw B once and keep it fixed #
  # ------------------------------#
  p_in = runif(1, min = p_range[1], max = p_range[2]) # within-community connection probability
  q_out = runif(1, min = q_range[1], max = q_range[2]) # between-community connection probability
  B = matrix(c(
    p_in, q_out,
    q_out, p_in
  ), nrow = 2, byrow = TRUE)
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
  # Build segment-specific mixed membership matrices
  # Z^(1), Z^(2), Z^(3)
  # -----------------------------------------#
  # the first half nodes belong to community 1, the last half nodes belong to community 2
  Z1 = matrix(0, nrow = n, ncol = C)
  for (i in 1:n) {
    alpha_i = rep(alpha_other, C)
    if (i <= n / 2) {
      # if node i is among the first half of all the nodes, it's alpha_i = (alpha_main, alpha_other)
      alpha_i[1] = alpha_main
    } else {
      # if node i is among the first half of all the nodes, it's alpha_i = (alpha_other, alpha_main)
      alpha_i[2] = alpha_main
    }
    Z1[i, ] = rdirichlet_simple(n_draws = 1, alpha = alpha_i)
  }
  
  # Standard C-dimensional basis vectors e1, e2
  e1 = c(1, rep(0, C - 1)) 
  e2 = c(0, 1, rep(0, C - 2))
  
  Z2 = Z1
  for (i in 1:num_swap) {
    # the first num_swap nodes move toward community 2
    Z2[i, ] = (1 - delta_z) * Z1[i, ] + delta_z * e2 
  }
  for (i in (n - num_swap + 1):n) {
    # the last num_swap nodes move toward community 1
    Z2[i, ] = (1 - delta_z) * Z1[i, ] + delta_z * e1 
  }
  Z2 = normalize_rows(Z2)
  
  Z3 = Z1 # all moving nodes move back to where they belong at the beginning
  
  # ------------------------------------------#
  # Segment-specific probability matrices P^(s)
  # ------------------------------------------#
  P1 = Z1 %*% B %*% t(Z1)
  P2 = Z2 %*% B %*% t(Z2)
  P3 = Z3 %*% B %*% t(Z3)
  
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
  #
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
    B = B,
    alpha = list(alpha_main = alpha_main, alpha_other = alpha_other),
    delta_z = delta_z,
    basis = list(e1 = e1, e2 = e2),
    z = list(z1 = Z1, z2 = Z2, z3 = Z3),
    Z = list(Z1 = Z1, Z2 = Z2, Z3 = Z3),
    P = list(P1 = P1, P2 = P2, P3 = P3)
  ))
}
