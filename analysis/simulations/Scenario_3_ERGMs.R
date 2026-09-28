# ============================================================#
# Simulation Scenario 3:
# Independent ERGMs with piecewise-constant parameters
# ============================================================#

sim_scenario3_ergm_independent <- function(
    n = 50,
    Tt = 300,
    cp = c(75, 225),
    network_stats = c("edges", "gwesp(0.5, fixed = TRUE)"),
    theta = NULL,
    node_attr = NULL,
    seed = NULL,
    max_tries = 1,
    MCMC.burnin = 5000,
    MCMC.interval = 1000
) {
  # --------------#
  # Sanity checks #
  # --------------#
  if (!is.null(seed)) set.seed(seed)
  
  if (!exists("compute_netstats", mode = "function", inherits = TRUE)) {
    stop("compute_netstats must be defined outside sim_scenario3_ergm_independent().")
  }
  
  if (length(cp) != 2) stop("cp must be a length-2 vector.")
  if (cp[1] < 1 || cp[1] >= cp[2] || cp[2] >= Tt) {
    stop("Require 1 <= cp[1] < cp[2] < Tt.")
  }
  if (n < 2) stop("n must be at least 2.")
  if (max_tries < 1) stop("max_tries must be at least 1.")
  if (!is.null(node_attr) && length(node_attr) != n) {
    stop("node_attr must be NULL or a vector of length n.")
  }
  
  required_pkgs = c("network", "ergm")
  missing_pkgs = required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
  if (length(missing_pkgs) > 0) {
    stop(paste0(
      "The following packages must be installed before running independent ERGMs: ",
      paste(missing_pkgs, collapse = ", "), "."
    ))
  }
  
  # ergm uses "triangle", not "triangles"
  ergm_stats = network_stats
  ergm_stats[ergm_stats == "triangles"] = "triangle"
  
  num_seg = length(cp) + 1
  num_stats = length(ergm_stats)
  
  # ------------------------------------------------------------#
  # Default segment-specific ERGM coefficients.
  # Rows correspond to network_stats.
  # Columns correspond to the three segments.
  #
  # Default is for:
  #   edges + gwesp(0.5, fixed = TRUE)
  # ------------------------------------------------------------#
  if (is.null(theta)) {
    theta = matrix(c(
      -2.10, -1.80, -2.10,   # edges
      0.20,  0.12,  0.20    # gwesp
    ), nrow = 2, byrow = TRUE)
  }
  
  if (!is.matrix(theta)) {
    theta = as.matrix(theta)
  }
  
  if (!all(dim(theta) == c(num_stats, num_seg))) {
    stop("theta must have length(network_stats) rows and length(cp) + 1 columns.")
  }
  
  rownames(theta) = ergm_stats
  colnames(theta) = paste0("segment", 1:num_seg)
  
  # ------------------------------------------------#
  # Helper: convert a network object to a symmetric
  # zero-diagonal adjacency matrix.
  # ------------------------------------------------#
  network_to_adjacency = function(net_obj) {
    if (inherits(net_obj, "network.list")) {
      net_obj = net_obj[[1]]
    }
    if (is.list(net_obj) && !inherits(net_obj, "network")) {
      net_obj = net_obj[[1]]
    }

    A_tmp = as.matrix(net_obj, matrix.type = "adjacency")
    A_tmp = matrix(as.numeric(A_tmp), nrow = nrow(A_tmp), ncol = ncol(A_tmp))
    A_tmp = 1 * ((A_tmp + t(A_tmp)) > 0)
    diag(A_tmp) = 0

    return(A_tmp)
  }
  
  # -----------------------------------------#
  # Segment starts and ends for observations
  # -----------------------------------------#
  starts = c(1, cp + 1)
  ends = c(cp, Tt)
  
  # -----------------------------------------#
  # Empty basis network for ERGM simulation
  # -----------------------------------------#
  base_net = network::network(
    matrix(0, nrow = n, ncol = n),
    directed = FALSE,
    loops = FALSE,
    matrix.type = "adjacency",
    ignore.eval = TRUE
  )
  
  if (!is.null(node_attr)) {
    network::set.vertex.attribute(base_net, "node_attr", node_attr)
  }
  
  ergm_formula = stats::as.formula(
    paste("base_net ~", paste(ergm_stats, collapse = " + "))
  )
  
  # -----------------------------------------#
  # Simulate independent ERGMs
  # -----------------------------------------#
  A = vector("list", Tt)
  segment_id = integer(Tt)
  
  for (seg in 1:num_seg) {
    
    for (t in starts[seg]:ends[seg]) {
      sim_net = NULL
      last_error = NULL
      
      for (attempt in 1:max_tries) {
        sim_net = tryCatch(
          stats::simulate(
            ergm_formula,
            coef = as.numeric(theta[, seg]),
            nsim = 1,
            output = "network",
            control = ergm::control.simulate.formula(
              MCMC.burnin = MCMC.burnin,
              MCMC.interval = MCMC.interval
            )
          ),
          error = function(e) {
            last_error <<- e
            return(NULL)
          }
        )
        
        if (!is.null(sim_net)) break
      }
      
      if (is.null(sim_net)) {
        stop(paste0(
          "Independent ERGM simulation failed at time ", t,
          " in segment ", seg,
          " after ", max_tries,
          " attempts. Last error: ",
          conditionMessage(last_error)
        ))
      }
      
      A[[t]] = network_to_adjacency(sim_net)
      segment_id[t] = seg
    }
  }
  
  # ------------------------------------------------------------#
  # Use your externally defined compute_netstats(A_list).
  # This returns only edges and triangles.
  # ------------------------------------------------------------#
  S = compute_netstats(A)
  
  return(list(
    A = A,
    S = S,
    cp = cp,
    segment_id = segment_id,
    network_stats = network_stats,
    ergm_stats = ergm_stats,
    theta = theta,
    node_attr = node_attr
  ))
}
