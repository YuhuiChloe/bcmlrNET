bcmlrNET_model_select <- function(A, S, latent_dim = 5, init = "even", thinning = 5, prior_kappa = "default", alpha_f = 0.1, 
                                  threshold = 0.5, min_size = 30, max_num_cp = 5, 
                                  num_iter = 5000, num_warmup = 2500, print_progress = FALSE, print_outputs = FALSE){
  if (!is.matrix(S)){
    stop("Please convert your summary statistics to a matrix. If you have factor variables, apply a dummification or one-hot encoding.")
  }
  if (!is.array(A)){
    stop("Please input a 3D array of adjacency matrices. The first dimension is the number of observation.")
  }
  if (max_num_cp > (dim(A)[1]/min_size -1)){
    stop("The value of max_num_cp exceeded the maximum possible number of changepoints under the specified min_size.")
  }
  # Select the number of changepoints
  old_out = bcmlrNET(A, S, num_CP = max_num_cp, latent_dim = latent_dim, init = init, prior_kappa = prior_kappa, 
                     alpha_f = alpha_f, thinning = thinning, min_size = min_size, num_iter = num_iter, num_warmup = num_warmup,
                  print_outputs = FALSE, print_progress = print_progress, model_selection = TRUE)
  boolean = (old_out$AUC < threshold) # TRUE if less than threshold, FALSE otherwise
  rej = rowSums(!as.matrix(boolean)) # Count number of FALSE in each row
  tb = table(factor(rej, levels = seq_len(max_num_cp+1)-1))/length(rej)  # The levels start from 0 to max_num_cp
  num_cp = as.integer(names(tb)[which.max(tb)]) # the posterior mode is taken as the estimated num of CP
  probs = as.numeric(tb) # the posterior probabilities of number of CPs from 0 to max_num_cp
  names(probs) = names(tb) 
  if (num_cp == 0){
    print("There is no change.")
    out = list()
    out$num_cp = num_cp # estimate of the number of CPs 
    out$max_cps = old_out$Kappa_mode # the max number of posterior modes from the "old output" before refitting on all the data
    out$num_cp_dist = rej # posterior distribution of the number of CPs
    out$probs_num_cp = probs 
  }else{
    # Refit with all the data: 
    out = bcmlrNET(A, S, init = init, num_CP = num_cp, latent_dim = latent_dim, prior_kappa = prior_kappa,
                   thinning = 1, min_size = min_size, num_iter = num_iter, num_warmup = num_warmup,
                print_outputs = FALSE, print_progress = print_progress, model_selection = FALSE)
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
