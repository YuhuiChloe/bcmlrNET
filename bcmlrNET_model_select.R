bcmlrNET_model_select <- function(A, S, latent_dim = 10, eta_dir = 1.5, b_tau = 1, thinning = 5, prior_kappa = "default", alpha_f = 0.1, 
                                 threshold = 0.5, min_size = 30, max_num_cp = 5, sd_beta = 3,
                                 num_iter = 5000, num_warmup = 2500, print_progress = FALSE, print_outputs = FALSE,
                                 prior_for_U = c("Spike_Slab", "Gaussian"), lambda_init = c("equal", "prior", "zero"),
                                 sparse_networks = FALSE, center_std_A = TRUE, seed = NULL){
  prior_for_U = match.arg(prior_for_U)
  lambda_init = match.arg(lambda_init) 
  if (is.null(A) && is.null(S)){
    stop("Please enter a list of adjacency matrices or/and a matrix of network statistics.")
  }
  Tt_A = if (!is.null(A)) length(A) else -Inf
  Tt_S = if (!is.null(S)) dim(S)[1] else -Inf
  if (max_num_cp > max(Tt_A/min_size -1, Tt_S/min_size -1)){
    stop("The value of max_num_cp exceeded the maximum possible number of changepoints under the specified min_size.")
  }
  # Select the number of changepoints
  old_out = bcmlrNET(A, S, num_CP = max_num_cp, latent_dim = latent_dim, eta_dir = eta_dir, b_tau = b_tau, prior_kappa = prior_kappa, 
                    alpha_f = alpha_f, min_size = min_size, sd_beta = sd_beta, 
                    thinning = thinning, num_iter = num_iter, num_warmup = num_warmup,
                    print_outputs = FALSE, print_progress = print_progress, model_selection = TRUE,
                    prior_for_U = prior_for_U, sparse_networks = sparse_networks, center_std_A = center_std_A,
                    seed = seed) 
  boolean = (old_out$AUC <= threshold) # TRUE if less than threshold, FALSE otherwise
  boolean[is.na(boolean)] = TRUE # Treat non-informative/failed AUC as no rejection.
  rej = rowSums(!as.matrix(boolean)) # Count number of FALSE in each row
  tb = table(factor(rej, levels = seq_len(max_num_cp+1)-1))/length(rej)  # The levels start from 0 to max_num_cp
  num_cp = as.integer(names(tb)[which.max(tb)]) # the posterior mode is taken as the estimated num of CP
  probs = as.numeric(tb) # the posterior probabilities of number of CPs from 0 to max_num_cp
  names(probs) = names(tb) 
  print(num_cp)
  if (num_cp == 0){
    print("There is no change.")
    out = list()
    out$num_cp = num_cp # estimate of the number of CPs 
    out$max_cps = old_out$Kappa_mode # the max number of posterior modes from the "old output" before refitting on all the data
    out$num_cp_dist = rej # posterior distribution of the number of CPs
    out$probs_num_cp = probs 
  }else{
    # Refit with all the data: 
    out = bcmlrNET(A, S, num_CP = num_cp, latent_dim = latent_dim, eta_dir = eta_dir, b_tau = b_tau, prior_kappa = prior_kappa,
                   alpha_f = alpha_f,  min_size = min_size, sd_beta = sd_beta, 
                   thinning = 1, num_iter = num_iter, num_warmup = num_warmup,
                   print_outputs = FALSE, print_progress = print_progress, model_selection = FALSE,
                   prior_for_U = prior_for_U, sparse_networks = sparse_networks, center_std_A = center_std_A,
                   seed = seed) 
    out$num_cp = num_cp # estimate of the number of CPs 
    out$max_cps = old_out$Kappa_mode # the max number of posterior modes from the "old output" before refitting on all the data
    out$num_cp_dist = rej # posterior distribution of the number of CPs
    out$probs_num_cp = probs 
  }
  if (print_outputs){
    cat(str(out))
  }
  return(out)
}
