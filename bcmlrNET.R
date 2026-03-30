bcmlrNET <- function(A, S, num_CP, latent_dim = 5, init = "even", prior_kappa = "default", alpha_f = 0.1, 
                     min_size = 30, sd_beta = 3, thinning = 1, 
                     num_iter = 5000, num_warmup = 2500, print_outputs = FALSE, 
                     print_progress = FALSE, model_selection = FALSE){
  # S is a matrix with Tt rows and p columns of summary statistics.
  # A is an input series with Tt matrices in dimensions n \times n. 
  n = dim(A)[2] # number of nodes
  p = dim(S)[2] # number of summary statistics
  Tt_all = dim(A)[1]  # total number of observations
  D = latent_dim # dimensions of latent variables per node 
  #===============#
  # Sanity checks ####
  #===============#
  if (dim(A)[1] != dim(S)[1]){
    stop("The series of adjacency matrices (A) and summary statistics (S) must have the same length. ")
  }
  if (dim(A)[2]!=dim(A)[3]){
    stop("Adjacency matrices must be square.")
  }
  if (latent_dim > n){
    stop("The number of latent variables (D) cannot exceed the number of nodes (n).")
  }
  #=================#
  # Data processing ####
  #=================#
  # A function to center & standardize a matrix x
  center_standardize <- function(x, m, sdv){
    for (i in 1:ncol(x)){
      x[,i] = x[,i] - m[i]
      x[,i] = x[,i]/sdv[i]
    }
    return(x)
  }
  
  # Use all observations (thinning =1) when num_CP is known
  if (thinning ==1){ 
    Tt = Tt_all
    all_idx = seq_len(Tt) # all indices from the input data
    train_idx = all_idx
    holdout_idx = all_idx
    
    # Standardize and center the summary statistcis S
    if (dim(S)[2] == 1){
      m = mean(S, na.rm = TRUE)
    }else{
      m = colMeans(S, na.rm = TRUE) # column means in the training set
    }
    sdv = apply(S, 2, FUN = sd)
    sdv[sdv == 0] = 1 # Avoid divide-by-zero in standardization: sdv = 0 could occur if the input data has a column of the same constant. 
    S = center_standardize(S, m, sdv)
    S_auc = S
    S_all = S
    
    A_all = A
    A_auc = A
  }
  # Set thinning >1 when selecting the number of changepoints (model_selection = TRUE)
  else if (thinning > 1){ 
    Tt = dim(S)[1]
    all_idx = seq_len(Tt) # all indices from the input data
    holdout_idx = all_idx[all_idx%%thinning == 0] # the indices of the held-out samples
    train_idx = all_idx[-holdout_idx] # the indices for fitting the model (the trainng set indices)
    
    S_all = S
    S_auc = S_all[holdout_idx, , drop = FALSE] # rows of S selected to calculate AUC
    S = S_all[-holdout_idx, , drop = FALSE] # rows of S selected to fit the model  
    
    A_all = A
    A_auc = A_all[holdout_idx, , , drop = FALSE]
    A = A_all[-holdout_idx, , ,drop = FALSE]
    
    # Standardize and center the summary statistcis S
    if (ncol(S) == 1){
      m = mean(S, na.rm = TRUE)
    }else{
      m = colMeans(S, na.rm = TRUE) # column means in the training set
    }
    sdv = apply(S, 2, FUN = sd, na.rm = TRUE) # standard deviation of the training set
    sdv[sdv == 0] = 1 # Avoid divide-by-zero in standardization: sdv = 0 could occur if the input data has a column of the same constant. 
    S = center_standardize(S, m, sdv) # center & standardize the subseries used to fit the model with its column means
    S_auc = center_standardize(S_auc, m, sdv) # center & standardize the held-out samples with the same column means
  }
  
  Tt_auc = dim(S_auc)[1]
  Tt = dim(S)[1] # number of observations
  K = num_CP + 1 # number of classes
  L = num_CP # number of changepoints
  #===================#
  # Utility functions ####
  #===================#
  # https://arxiv.org/pdf/2205.04997
  # https://gregorygundersen.com/blog/2020/02/09/log-sum-exp/
  # A function that mitigates computational issues in massive division or multiplication of extreme values.
  logsumexp <- function(x){
    c <- max(x)
    return(c + log(sum(exp(x - c))))
  }
  # A function that converts a change point vector to a response matrix with the change point structure.
  kappa_to_Y <- function(kappa, N, L){
    class <- rep(1, kappa[1])
    for(j in 2:L){
      class <- c(class, rep(j, kappa[j] - kappa[j-1]))
    }
    class <- c(class, rep(L+1, N - kappa[L]))
    class <- as.factor(class)
    Y <- model.matrix(~ -1 + class, data=class)
    return(Y)
  }
  #=======================#
  # Create an output list ####
  #=======================#
  out = list()
  out$tau_sq = rep(0, num_iter-num_warmup)
  #out$Q = array(0, dim = c(num_iter-num_warmup, Tt, K))
  out$Lambda = array(0, dim = c(num_iter-num_warmup, K, D))
  out$U = array(0, dim = c(num_iter-num_warmup, K, n, D))
  out$Z = array(0, dim = c(num_iter-num_warmup, n, K))
  out$Beta = array(0, dim = c(num_iter-num_warmup, p, K))
  out$W = array(0, dim = c(num_iter-num_warmup, n, K))
  out$pi = matrix(0, nrow = num_iter-num_warmup, ncol = 3)
  #out$omega = array(0, dim = c(num_iter-num_warmup, Tt, K))
  #out$phi = array(0, dim = c(num_iter-num_warmup, Tt, K))
  if (num_CP==1){
    out$Kappa = rep(0,L)
  }else if (num_CP>1){
    out$Kappa = matrix(data = 0, nrow = num_iter-num_warmup, ncol = L)
  }
  if (model_selection){
    if (num_CP>1){
      out$AUC = matrix(0, nrow = num_iter-num_warmup, ncol = L)
    }
  }else{
    out$AUC = NULL
  }
  #============================#
  # hyperparameters for priors ####
  #============================#
  # Prior quantities for beta (coefficients of summary statistics):
  m0 <- array(data=0, dim=c(p,K-1)) # m0[,k] is the prior mean of the kth coef vec
  V0 <- array(data=0, dim=c(p,p,K-1)) # V0[,,k] is the prior cov. of the kth coef vec
  V0inv <- array(data=0, dim=c(p,p,K-1)) # V0inv[,,k] is the inverse of V0[,,j]
  sd0 <- sd_beta # prior variance for beta (coefficients of summary statistics)
  for(k in 1:(K-1)){
    V0[,,k] <- sd0*diag(p)
    V0inv[,,k] <- 1/sd0*diag(p)
  }
  # priors for tau_sq
  a_tau = 1
  b_tau = 1
  # priors for w_ik ~ Beta(a_w, b_w)
  a_w = 1
  b_w = 1
  # prior for pi (Dirichlet)
  eta_dir = 2 # must be >1
  #===============================#
  # Initialize working quantities ####
  #===============================#
  # Q = matrix(0, nrow = Tt, ncol = K) # Frobenius inner products
  tau_sq = rigamma(n = 1, a = a_tau, b = b_tau) # inverse gamma variable
  Lambda = matrix(0, nrow = K, ncol = D) # Impose sparsity on the latent variables
  Z = matrix(rbinom(n = n*K, size = 1, prob = 0.5), nrow = n, ncol = K) # Bernoulli variables
  U = array(0, dim = c(K, n, D)) # latent variables
  for (i in 1:n){
    for (k in 1:K){
      if (Z[i,k] ==1){
        U[k,i,] = mvrnorm(n =1, mu = rep(0,D), Sigma = tau_sq*diag(D))
      }else if(Z[i,k] ==0){
        U[k,i,] = rep(0,D)
      }
    }
  }
  Beta = matrix(0, nrow = p, ncol = K) # coefficients for summary statistics
  W = matrix(rbeta(n*K, shape1 = a_w, shape2 = b_w), nrow = n, ncol = K) # Beta variables
  pi = rep(1/3,3)
  omega = matrix(0, nrow = Tt, ncol = K)
  if (init == "even"){
    if(num_CP==1){
      kappa = round(Tt/2)
    }else if (num_CP>1){
      kappa = rep(0, L)
      for (l in 1:L){
        kappa[l] = round(Tt/K) *l
      }
    }
  }
  original_kappa = train_idx[kappa]
  #=================================#
  # Begin Gibbs sampling (iter)  ####
  #=================================#
  if (print_progress){
    pb = txtProgressBar(min = 0, max = num_iter, style = 3)
  }
  start = Sys.time()
  for (iter in 1: num_iter){
    if (num_CP == 1){
      Y = c(rep(0, times = kappa), rep(1, times = Tt - kappa))
      delta = matrix(Y-1/2, nrow = Tt, ncol = 1)
    }else if (num_CP > 1){
      Y = kappa_to_Y(kappa, Tt, L)
      delta = Y-1/2
    }
    # Update Q_tk per iteration
    Q_iter = matrix(0, nrow = Tt, ncol = K)
    # update quantities in class k
    # Begin for-loop of k = 1, ..., (K-1) ####
    for (k in 1:(K-1)){
      #============================#
      # After each parameter in class k is updated (in each iteration of k), 
      # get the latest phi_{tk}, c_{tk} and eta_{tk} based on all classes.
      #============================#
      #=============#
      # Compute phi ####
      #=============#
      phi = matrix(0, nrow = Tt, ncol = K)
      for (kk in 1:(K-1)){
        # We must use all observed matrices A_t, 
        # just like we use all vectors x_i in a non-network data set. 
        for (t in 1:Tt){
          A_t = A[t,,]
          # the i_th row of R_tmp is r_{ti,kk}^{\top} (1*D)
          R_t_kk = A_t %*% U[kk,,] %*% diag(Lambda[kk,]) # n*n n*D D*D -> n*D
          # Q_{tk} = 1/2 * tr(U_k^{\top} A_t U_k Lambda_k) 
          # * represents element-wise multiplication
          # trace = sum(diag()), but here we have a diagonal matrix: sum will do. 
          Q_t_kk  = 0.5*sum(U[kk,,] * R_t_kk)
          # inner product of summary statistics and coefficients
          StBkk = as.numeric(crossprod(S[t,], Beta[,kk]))
          # phi is T * K with the Kth column = 0.  
          phi[t,kk] = Q_t_kk + StBkk
        }
      }# End for-loop of kk = 1, ..., (K-1) 
      # Now, we have phi. 
      #===========================#
      # Compute c_{k} and eta_{k} ####
      #===========================#
      if(num_CP ==1){
        c_k = rep(0, Tt)
      }else if (num_CP>1){
        c_k = log(rowSums(exp(phi[,-k]))) # Tt * 1
      }
      eta_k = phi[,k] - c_k
      #============================#
      # Update PG variables, omega ####
      #============================#
      for (t in 1:Tt){
        omega[t,k] <- rpg(num=1, h = 1, z = eta_k[t]) # update omega
      }
      # (The current loop is for class k.)
      #===============================================#
      # Update z_{ik}, u_{ik}, w_{ik} for each node i ####
      #===============================================#
      # Begin for-loop of i = 1, ...n ####
      for (i in 1:n){
        #====================# 
        # Compute H_ik, g_ik ####
        #====================#
        # Not as interesting as the fitted values phi_{tk}. 
        # So, we don't need to track their MCMC samples. 
        H_ik = matrix(0, nrow = D, ncol = D)
        g_ik = rep(0, D)
        for (t in 1:Tt){
          A_t = A[t,,] # n*n
          R_tk = A_t %*% U[k,,] %*% diag(Lambda[k,]) #n*D
          r_tik = R_tk[i,] # D*1
          Q_tk = 0.5*sum(U[k,,]*R_tk) *1*1
          StBk = as.numeric(crossprod(S[t,], Beta[,k])) 
          b_tik = Q_tk - t(U[k,i,])%*%r_tik + StBk - c_k[t]
          # Accumulate H_ik and g_ik over t = 1,...,T
          H_ik = H_ik + omega[t,k]*r_tik%*%t(r_tik)
          g_ik = g_ik + r_tik*as.numeric(delta[t,k] - b_tik*omega[t,k])
        } # End for-loop of t = 1, ..., Tt
        # Two working quantities, H_ik, g_ik, are updated. 
        #============================================# 
        # Compute M_ik and related working quantities####
        #============================================#
        M_ik = H_ik + diag(D)*(1/tau_sq)
        # https://math.stackexchange.com/questions/3158303/using-cholesky-decomposition-to-compute-covariance-matrix-determinant
        # If A is a symmetric p.d. matrix, then 
        # det(A) = det(LL^{\top}) = (det(L))^2 = (prod of diagnoal entries of L)^2
        # and log(det(A)) = 2 sum(log(diag(L)))
        # L is the lower triangular matrix of the Cholesky decomposition of A. 
        # The same equation also holds for the upper triangular matrix
        Uchol_M = chol(M_ik) # Upper triangular matrix of the Choleskey decomposition of M_ik
        logdet_M = 2 * sum(log(diag(Uchol_M))) 
        Minv_g = backsolve(Uchol_M, forwardsolve(t(Uchol_M), g_ik)) # (U^{\top}U)^{-1} g
        g_Minv_g = t(g_ik) %*% Minv_g # g^{\top} (U^{\top}U)^{-1} g 
        #=========================# 
        # Update Z_ik (Bernoulli) ####
        #=========================#
        # Now we can compute logW_0 and logW_1 
        # (where we omitted a term that will be cancelled out in the posterior probabilities of z_ik)
        logW0 = log(W[i,k])
        logW1 = log(1-W[i,k]) - (D / 2) * log(tau_sq) - 0.5 * logdet_M + 0.5 * g_Minv_g 
        # Probability of Z_ik = 1
        prob_z1 = 1/(1+ exp(logW0-logW1))
        # use the max trick for stability?:  
        # maxlog = max(logW0, logW1)
        # W0_shift = exp(logW0 - maxlog)
        # W1_shift = exp(logW1 - maxlog)
        # prob_z1 = W1_shift / (W0_shift + W1_shift)
        Z[i,k] = rbinom(1, size = 1, prob = prob_z1)
        #===================# 
        # Update u_ik (D*1) ####
        #===================#
        if (Z[i,k] == 0){
          U[k,i,] = rep(0,D)
        }else{
          # Reusse Uchol_M = chol(M_ik), where M_ik = V_ik^{-1} is the precision matrix
          m_ik = backsolve(Uchol_M, forwardsolve(t(Uchol_M), g_ik))
          stdGaussians = matrix(rnorm(D), D, 1) 
          U[k,i,] = m_ik + backsolve(Uchol_M, stdGaussians) 
        }
        #======================# 
        # Update w_ik (scalar) ####
        #======================#
        W[i,k] = rbeta(1, a_w + 1 - Z[i, k], b_w + Z[i, k])
      }# End for-loop of i = 1, ...n ####
      #=======================# 
      # Update Lambda_k (D*D) ####
      #=======================#
      # Need to recompute Q_tk because, in the current iteration, 
      # it appeared in the local for-loop for computing H_ik and g_ik 
      # before u_ik was updated 
      Q_k = rep(0, Tt) # a vector of Q_1k, ..., Q_Tk
      q_k = matrix(0, nrow = Tt, ncol = D) # the (t,d)th entry is q_tdk
      for (t in 1:Tt){
        A_t = A[t,,]
        for (d in 1:D) {
          
          # q_tdk = sum_{i<j} a_ijt u_ikd u_jkd = 0.5 * u_kd^T A_t u_kd, where u_kd is in n*1
          q_k[t, d] = 0.5 * as.numeric(crossprod(U[k,,d], A_t %*% U[k,,d]))
        }
        # Q_k = \sum_{d=1}^D \Lambda_k q_{kt}
        Q_k[t] = sum( diag(Lambda[k,]) %*% q_k[t,] ) # 1*1, where q_k[t,] is D*1
      } # Done updating Q_k and q_k
      # Now, we can update Lambda_dk for each d 
      Lambda_k = Lambda[k,]
      for (d in 1:D){
        lambda_dk_old = Lambda_k[d]
        alpha_k_minus_d = Q_k - lambda_dk_old*q_k[d] + as.numeric(S %*% Beta[,k]) - c_k   # Tt*1
        cand_vals = c(-1, 0, 1) # support of lambda_dk, (3 candidate values) 
        logweights = rep(0, length(cand_vals))
        for (r in 1:length(cand_vals)){
          v = cand_vals[r]
          eta_cand_k = alpha_k_minus_d + v*q_k[,d]   # Tt*1
          # element-wise multiplication between Tt*1 vectors in the sum():
          logweights[r] = log(pi[r]) + sum(delta[,k]*eta_cand_k - 0.5*omega[,k]*eta_cand_k^2) 
        } 
        # Use the maximum shift trick 
        max_logweights = max(logweights)
        shifted_weights = exp(logweights - max_logweights)
        # Get posterior probabilities for lambda's
        probs_lambda = shifted_weights/(sum(shifted_weights)) # 3*1 vector that should sum up to 1
        # Update lambda_dk
        lambda_dk_new = sample(cand_vals, size = 1, prob = probs_lambda)
        Lambda_k[d] = lambda_dk_new
        # Update Q_k as lambda_dk is updated
        Q_k = Q_k + (lambda_dk_new - lambda_dk_old)*q_k[,d] # Q_t is Tt*1 and q_t[,d] is also Tt*1
      }# End for-loop of d
      Lambda[k,] = Lambda_k
      # Done updating lambda_dk for each d (k is still fixed in the current iteration)
      #=====================# 
      # Update beta_k (p*1) ####
      #=====================#
      normals = matrix(rnorm(p), p, 1) # A p-vector of standard normals
      precision_beta = as.numeric(crossprod(x = S, y = diag(omega[,k]) %*% S)) + V0inv[,,k]
      Uchol_beta = chol(precision_beta) # Upper Cholesky factor
      # crossprod(S, delta[, k] + omega[,k] * (c_k - Q_k)) + V0inv[,,k] %*% m0[, k]
      Beta[,k] = chol2inv(Uchol_beta) %*% 
        (t(S) %*% (delta[,k] + diag(omega[,k]) %*% (c_k - Q_k)) + V0inv[,,k] %*% m0[, k]) 
      + backsolve(Uchol_beta, normals)
      # Save FINAL Q_tk for this class after U and Lambda updates
      Q_iter[, k] = Q_k
    } # End of for loop for k = 1, ..., (K-1) ####
      #=======================================# 
      # Update pi (shared across all classes) ####
      #=======================================#
      # accumulate the numbers of lambda's = -1, 0, or 1 over all classes
      n_neg1 = 0
      n_0 = 0 
      n_pos1 = 0
      # recall: we set Lambda is in dimensions K*d*d (D*D diagonal matrix per class)
      for (k in 1:(K-1)){
        n_neg1 = n_neg1 + sum( diag(Lambda[k,]) == -1 )
        n_0 = n_0 + sum( diag(Lambda[k,]) == 0 )
        n_pos1 = n_pos1 + sum( diag(Lambda[k,]) == 1 )
      }
      # eta_dir is the hyperparameter in Dirichlet(D^eta_dir, 1, 1)
      alpha_pi = c(D^eta_dir + n_neg1, 1 + n_0, 1 + n_pos1) 
      # use the DirichletReg package to draw Dirichlet samples. Another package is MCMCpack. 
      # pi = DirichletReg::rdirichlet(n =1, alpha = alpha_pi) 
      # Dirichlet draw via Gammas
      gamma_draw = rgamma(3, shape = alpha_pi, rate = 1)
      pi = gamma_draw / sum(gamma_draw)
      #===========================================# 
      # Update tau_sq (shared across all classes) ####
      #===========================================#
      sum_m_k = 0
      sum_ss_k = 0
      for (k in 1:(K-1)){
        for (i in 1:n){
          # Z: n*K
          if (Z[i,k] == 1){ 
            sum_m_k = sum_m_k + 1
            sum_ss_k = sum_ss_k + sum(U[k,i,]^2) # U: K*n*D, sum(U[k,i,]^2) = (u_ik)^T u_ik
          }
        }
      }
      #tau_sq = 1/rgamma(n=1, shape = a_tau + (D/2) *sum_m_k, rate =b_tau + 0.5*sum_ss_k)
      tau_sq = rigamma(n = 1, a = a_tau + (D/2) *sum_m_k, b = b_tau + 0.5*sum_ss_k) 
    #} 
    # Done updating all other parameters except the changepoints
    
    #=============================#
    # Update kappa (changepoints) ####
    #=============================#
    if (num_CP == 1){
      Phi_exp <- exp(phi[,1]) 
      P = Phi_exp/(1+Phi_exp)
      # range <- 1:(N-1)
      range <- (1+min_size):(Tt-1-min_size)
      lr <- length(range)
      log_prob_vec_prop <- rep(0, lr) 
      if (is.character(prior_kappa) && prior_kappa == "default"){
        for(l in 1:lr){
          log_prob_vec_prop[l] = sum(log(1-P[1:range[l]])) + sum(log(P[(range[l]+1):Tt]))
          - length(1:range[l])*log(length(1:range[l]))
          - length((range[l]+1):Tt)*log(length((range[l]+1):Tt))
        }
      }else if (is.function(prior_kappa)){
        jt_pmf_kappa = prior_kappa(kappa)
        for(l in 1:lr){
          log_prob_vec_prop[l] = sum(log(1-P[1:range[l]])) + sum(log(P[(range[l]+1):Tt])) + log(jt_pmf_kappa)
        }
      }
      prob_vec <- exp(log_prob_vec_prop - logsumexp(log_prob_vec_prop)) 
      if(lr > 1){
        kappa = sample(x=range, size=1, prob=prob_vec)
      }
      #=============================================#
      # Select num of CPs: Is there a changepoint ? ####
      #=============================================#
      original_kappa = train_idx[kappa]
      # 4. Compute AUC
      if (model_selection){
        Y_all = c(rep(0, times = original_kappa), rep(1, times = Tt_all - original_kappa))
        Y_auc = Y_all[holdout_idx]
        #=============#
        # Compute phi ####
        #=============#
        phi_all = matrix(0, nrow = Tt_all, ncol = K)
        for (kk in 1:(K-1)){
          # We must use all observed matrices A_t, 
          # just like we use all vectors x_i in a non-network data set. 
          for (t in 1:Tt){
            # A_t_auc = A_all[t,,]
            # the i_th row of R_tmp is r_{ti,kk}^{\top} (1*D)
            R_t_kk_all = A_all[t,,] %*% U[kk,,] %*% diag(Lambda[kk,]) # n*n n*D D*D -> n*D
            # Q_{tk} = 1/2 * tr(U_k^{\top} A_t U_k Lambda_k) 
            # * represents element-wise multiplication
            # trace = sum(diag()), but here we have a diagonal matrix: sum will do. 
            Q_t_kk_all  = 0.5*sum(U[kk,,] * R_t_kk_all)
            # inner product of summary statistics and coefficients
            StBkk_all = as.numeric(crossprod(S_all[t,], Beta[,kk]))
            # phi is T * K with the Kth column = 0.  
            phi_all[t,kk] = Q_t_kk_all + StBkk_all
          }
        }
        # Update the matrix of probabilities of success based on the fitted coefficients
        Phi_exp_all <- exp(phi_all[,1]) 
        P_all = Phi_exp_all/(1+Phi_exp_all)
        P_auc = P_all[holdout_idx]
        ci = as.numeric(suppressWarnings(ci.auc(as.factor(Y_auc), as.numeric(P_auc), conf.level= 1-alpha_f, quiet = TRUE)))
      }
    }else if(num_CP > 1){
      # Create the P matrix
      Phi_exp <- exp(phi) # Tt * K
      P <- t(apply(Phi_exp, 1, function(x) x/sum(x))) # Normalize the rows of Phi_exp to get P
      if (is.character(prior_kappa) && prior_kappa == "default"){
        for(j in 1:(K-1)){
          # Now calculate the range of values that kappa[j] can take under its full
          # conditional distribution and the vector prob_vec, which is defined
          # such that prob_vec[k] is the log of the conditional probability that
          # kappa[j] = range[k] given all the other data and parameters.
          if(j==1){
            #range <- 1:(kappa[j+1]-1)
            range <- (1 + min_size):(kappa[j+1] - 1 - min_size)
            lr <- length(range)
            log_prob_vec_prop <- rep(0, lr)
            for(l in 1:lr){
              log_prob_vec_prop[l] <- sum(log(P[1:range[l],j])) + sum(log(P[(range[l]+1):kappa[j+1],j+1]))
              - length(1:range[l])*log(length(1:range[l]))
              - length((range[l]+1):kappa[j+1])*log(length((range[l]+1):kappa[j+1]))
            }
            #prob_vec <- exp(log_prob_vec_prop)/sum(exp(log_prob_vec_prop))
            prob_vec <- exp(log_prob_vec_prop - logsumexp(log_prob_vec_prop))
          } else if(j==L){
            #range <- (kappa[j-1]+1):(Tt-1)
            range <- (kappa[j-1]+1 + min_size):(Tt - 1 - min_size)
            lr <- length(range)
            log_prob_vec_prop <- rep(0, lr)
            for(l in 1:lr){
              log_prob_vec_prop[l] <- sum(log(P[(kappa[j-1]+1):range[l],j])) + sum(log(P[(range[l]+1):Tt,j+1]))
              - length((kappa[j-1]+1):range[l])*log(length((kappa[j-1]+1):range[l]))
              - length((range[l]+1):Tt)*log(length((range[l]+1):Tt))
            }
            #prob_vec <- exp(log_prob_vec_prop)/sum(exp(log_prob_vec_prop))
            prob_vec <- exp(log_prob_vec_prop - logsumexp(log_prob_vec_prop))
          } else{
            #range <- (kappa[j-1]+1):(kappa[j+1]-1)
            range <- (kappa[j-1]+1 + min_size):(kappa[j+1]-1 - min_size)
            lr <- length(range)
            log_prob_vec_prop <- rep(0, lr)
            for(l in 1:lr){
              log_prob_vec_prop[l] <- sum(log(P[(kappa[j-1]+1):range[l],j])) + sum(log(P[(range[l]+1):kappa[j+1],j+1]))
              - length((kappa[j-1]+1):range[l])*log(length((kappa[j-1]+1):range[l]))
              - length((range[l]+1):kappa[j+1])*log(length((range[l]+1):kappa[j+1]))
            }
            #prob_vec <- exp(log_prob_vec_prop)/sum(exp(log_prob_vec_prop))
            prob_vec <- exp(log_prob_vec_prop - logsumexp(log_prob_vec_prop))
          }
          if(lr > 1){
            kappa[j] <- sample(x=range, size=1, prob=prob_vec) # update kappa
          }
        }
      }# End of the case when we use default prior for kappa (based on segment lengths) as in the paper.
      else if (is.function(prior_kappa)){ # else if the user have a joint PMF for the kappas 
        jt_pmf_kappa = prior_kappa(kappa) # Get joint PMF of the prior distribution of kappa
        for(j in 1:(K-1)){
          if(j==1){
            range <- (1 + min_size):(kappa[j+1] - 1 - min_size)
            lr <- length(range)
            log_prob_vec_prop <- rep(0, lr)
            for(l in 1:lr){
              log_prob_vec_prop[l] <- sum(log(P[1:range[l],j])) + sum(log(P[(range[l]+1):kappa[j+1],j+1])) + log(jt_pmf_kappa)
            }
            prob_vec <- exp(log_prob_vec_prop - logsumexp(log_prob_vec_prop))
          } else if(j==L){
            range <- (kappa[j-1]+1 + min_size):(Tt - 1 - min_size)
            lr <- length(range)
            log_prob_vec_prop <- rep(0, lr)
            for(l in 1:lr){
              log_prob_vec_prop[l] <- sum(log(P[(kappa[j-1]+1):range[l],j])) + sum(log(P[(range[l]+1):Tt,j+1])) + log(jt_pmf_kappa)
            }
            prob_vec <- exp(log_prob_vec_prop - logsumexp(log_prob_vec_prop))
          } else{
            range <- (kappa[j-1]+1 + min_size):(kappa[j+1]-1 - min_size)
            lr <- length(range)
            log_prob_vec_prop <- rep(0, lr)
            for(l in 1:lr){
              log_prob_vec_prop[l] <- sum(log(P[(kappa[j-1]+1):range[l],j])) + sum(log(P[(range[l]+1):kappa[j+1],j+1])) + log(jt_pmf_kappa)
            }
            prob_vec <- exp(log_prob_vec_prop - logsumexp(log_prob_vec_prop))
          }
          if(lr > 1){
            kappa[j] <- sample(x=range, size=1, prob=prob_vec) # update kappa
          }
        }
      } # End of the case if the user has a joint PMF of kappa
      # Done updating the changepoints (kappas)
      #============================================#
      # Select the number of multiple changepoints: ####
      # Which of the fitted changepoints correspond to true changeppints?  #
      #============================================#
      # Update kappa in the scale of the entire series
      original_kappa = train_idx[kappa]
      #=================================#
      # Compute AUC for model selection ####
      #=================================#
      if (model_selection){
        #=============#
        # Compute phi ####
        #=============#
        phi_all = matrix(0, nrow = Tt_all, ncol = K) # Fitted values
        for (kk in 1:(K-1)){
          # We must use all observed matrices A_t,
          # just like we use all vectors x_i in a non-network data set.
          for (t in 1:Tt_all){
            # A_t = A_all[t,,]
            # the i_th row of R_tmp is r_{ti,kk}^{\top} (1*D)
            R_t_kk_all = A_all[t,,] %*% U[kk,,] %*% diag(Lambda[kk,]) # n*n n*D D*D -> n*D
            # Q_{tk} = 1/2 * tr(U_k^{\top} A_t U_k Lambda_k)
            # * represents element-wise multiplication
            # trace = sum(diag()), but here we have a diagonal matrix: sum will do.
            Q_t_kk_all  = 0.5*sum(U[kk,,] * R_t_kk_all)
            # inner product of summary statistics and coefficients
            StBkk_all = as.numeric(crossprod(S_all[t,], Beta[,kk]))
            # phi is T * K with the Kth column = 0.
            phi_all[t,kk] = Q_t_kk_all + StBkk_all
          }
        }
        # "ratios" are the probabilities of a subject belonging to
        # class j+1 other than class j
        Phi_exp_all = exp(phi_all) # Element-wise exponent of Phi matrix from supplement
        P_all = t(apply(Phi_exp_all, 1, function(x) x/sum(x))) # Normalize the rows of Phi_exp to get P
        Y_all = kappa_to_Y(original_kappa, Tt_all, L) # Y in Tt_all *K
        ci = matrix(0, nrow = K-1, ncol = 3) # each row has a LB, point estimate, and UB
        for (j in 1:(K-1)){
          if(j == 1){
            rg = 1:original_kappa[j+1]
            lg = length(rg)
            ratios = P_all[rg,j+1]/(P_all[rg,j] + P_all[rg,j+1])
            Y_auc = Y_all[rg, j + 1]
            idx_test = which(rg %in% holdout_idx)
            P_auc = ratios[idx_test]
            Y_auc = Y_auc[idx_test]
            ci[j,] = as.numeric(suppressWarnings(ci.auc(as.factor(Y_auc), as.numeric(P_auc), conf.level= 1-alpha_f, quiet = TRUE)))
          }else if (j == K-1){
            rg = (original_kappa[j-1]+1):Tt_all
            lg = length(rg)
            ratios = P_all[rg,j+1]/(P_all[rg,j] + P_all[rg, j+1])
            Y_auc = Y_all[rg, j + 1]
            idx_test = which(rg %in% holdout_idx)
            P_auc = ratios[idx_test]
            Y_auc = Y_auc[idx_test]
            ci[j,] = as.numeric(suppressWarnings(ci.auc(as.factor(Y_auc), as.numeric(P_auc), conf.level= 1-alpha_f, quiet = TRUE)))
          }else{
            rg = (original_kappa[j-1]+1):original_kappa[j+1]
            lg = length(rg)
            ratios = P_all[rg,j+1]/(P_all[rg,j] + P_all[rg, j+1])
            Y_auc = Y_all[rg, j + 1]
            idx_test = which(rg %in% holdout_idx)
            P_auc = ratios[idx_test]
            Y_auc = Y_auc[idx_test]
            ci[j,] = as.numeric(suppressWarnings(ci.auc(as.factor(Y_auc), as.numeric(P_auc), conf.level= 1-alpha_f, quiet = TRUE)))
          }
        }
      } # End of selecting the number of changepoints given that num_CP >1
    } # END of the case "num_CP > 1"
    #=================================#
    # Store updates to the output list ####
    #=================================#
    if ( iter > num_warmup){
      #out$Q[iter-num_warmup,,] = Q_iter
      out$tau_sq[iter-num_warmup] = tau_sq
      out$Lambda[iter-num_warmup,,] = Lambda
      out$U[iter-num_warmup,,,] = U
      out$Z[iter-num_warmup,,] = Z
      out$Beta[iter-num_warmup,,] = Beta
      out$W[iter-num_warmup,,] = W
      out$pi[iter-num_warmup, ] = pi
      #out$omega[iter-num_warmup,,] = omega
      #out$phi[iter-num_warmup,,] = phi
      if (num_CP==1){
        out$Kappa[iter-num_warmup] = original_kappa
        # Use posterior modes as changepoint estimates:
        cp = table(out$Kappa)
        out$Kappa_mode = as.integer(names(cp)[which.max(cp)]) # posterior mode 
        if (model_selection){
          out$AUC[iter-num_warmup] = ci[1] # Take the CI lower bound
        }
      }else if (num_CP >1){
        out$Kappa[iter-num_warmup,] = original_kappa
        # Use posterior modes as changepoint estimates:
        cps = c()
        for (l in 1:L){
          cp = table(out$Kappa[,l])
          cp_mode = as.integer(names(cp)[which.max(cp)])
          cps = c(cps, cp_mode)
        }
        out$Kappa_mode = cps # posterior modes
        if (model_selection){
          out$AUC[iter-num_warmup,] = ci[,1] # Take the CI lower bound
        }
      }
    }
    if (print_progress){ setTxtProgressBar(pb, iter)}
  } # End of the Gibbs sampler (iter = 1, ..., num_iter)
  end = Sys.time()
  if (print_progress){ close(pb) }
  #================================#
  # END Gibbs sampling (iter)   ####
  #================================#
  # out$runTime = end-start
  return(out)
} 
