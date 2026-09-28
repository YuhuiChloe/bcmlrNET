# ============================================================#
# Simulation Scenario 1:
# Stochastic block model with piecewise-constant memberships
# and dyad-level temporal dependence
# ============================================================#

sim_scenario1_sbm <- function(n = 50, 
                              Tt = 300, 
                              C = 2, 
                              cp = c(75, 225), 
                              rho = 0.1, 
                              p_range = c(0.5, 1.0), 
                              q_range = c(0.0, 0.4), 
                              num_swap = 5,
                              seed = NULL) {
  # --------------#
  # Sanity checks #
  # --------------#
  if (!is.null(seed)) set.seed(seed)
  if (length(cp) != 2) stop("cp must be a length-2 vector.")
  if (C != 2) stop("This implementation assumes C = 2.")
  if (n %% 2 != 0) stop("n must be even for the default baseline split.")
  if (num_swap >n/2) stop("The number of nodes swapping cannot be larger than half of the total number of nodes.")
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
  # -----------------------------------------#
  # Build segment-specific membership vectors
  # Z^(1), Z^(2), Z^(3)
  # -----------------------------------------#
  # the first half nodes are in community 1, the last half nodes are in community 2
  z1 = c(rep(1, n / 2), rep(2, n / 2)) 
  
  z2 = z1
  z2[1:num_swap] = 2 # the first num_swap nodes move to community 2
  z2[(n-num_swap+1):n] = 1 # the last num_swap nodes move to community 1

  z3 = z1 # all moving nodes move back to where they belong at the beginning
  
  # construct membership matrices
  membership_to_Z = function(z, C) {
    Z = matrix(0, nrow = length(z), ncol = C)
    for (i in 1:length(z)) {
      Z[i, z[i]] = 1
    }
    return(Z)
  }
  Z1 = membership_to_Z(z1, C)
  Z2 = membership_to_Z(z2, C)
  Z3 = membership_to_Z(z3, C)
  
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
    swapped_nodes = c(1:num_swap,(n-num_swap+1):n),
    B = B,
    z = list(z1 = z1, z2 = z2, z3 = z3),
    Z = list(Z1 = Z1, Z2 = Z2, Z3 = Z3),
    P = list(P1 = P1, P2 = P2, P3 = P3)
  ))
}
