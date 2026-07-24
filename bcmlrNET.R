library(Rcpp)
library(RcppArmadillo)
library(BayesLogit) # for rpg()
library(LearnBayes) # for computing inverse gamma: rigamma(), alternative to 1/rgamma()
library(MASS) # for mvrnorm()
library(pROC) # for ROC curves and AUC values
library(DirichletReg)

options(pROCProgress = list(name = "none"))# Disable pROC progress bar globally

Rcpp::sourceCpp(code = '
#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

inline double plogis_stable_cpp(const double x) {
  if (x >= 0.0) {
    double e = std::exp(-x);
    return 1.0 / (1.0 + e);
  } else {
    double e = std::exp(x);
    return e / (1.0 + e);
  }
}

// [[Rcpp::export]]
arma::vec compute_Q_vec_cube_cpp(const arma::cube& A,
                                 const arma::mat& U_k,
                                 const arma::vec& lambda_k) {
  const int TT = A.n_slices;
  arma::vec Q(TT, arma::fill::zeros);
  arma::mat U_lambda = U_k.each_row() % lambda_k.t();

  for (int t = 0; t < TT; ++t) {
    arma::mat R_t = A.slice(t) * U_lambda;
    Q[t] = 0.5 * arma::accu(U_k % R_t);
  }
  return Q;
}

// [[Rcpp::export]]
Rcpp::List update_class_cpp(const arma::cube& A,
                            arma::mat U_k,
                            arma::ivec Z_k,
                            arma::vec W_k,
                            const arma::vec& lambda_k,
                            const arma::vec& StBk,
                            const arma::vec& c_k,
                            const arma::vec& delta_k,
                            const arma::vec& omega_k,
                            const double tau_sq,
                            const double a_w,
                            const double b_w,
                            const bool use_spike_slab) {
  const int n  = U_k.n_rows;
  const int D  = U_k.n_cols;
  const int TT = A.n_slices;

  std::vector<arma::mat> AU_list(TT);
  arma::mat q_k(TT, D, arma::fill::zeros);
  arma::vec Q_cache(TT, arma::fill::zeros);

  const arma::mat I_D = arma::eye(D, D);
  const double prior_prec = 1.0 / tau_sq;
  const double log_tau_sq = std::log(tau_sq);

  for (int t = 0; t < TT; ++t) {
    AU_list[t] = A.slice(t) * U_k;
    q_k.row(t) = 0.5 * arma::sum(U_k % AU_list[t], 0);
  }
  Q_cache = q_k * lambda_k;

  for (int i = 0; i < n; ++i) {
    arma::rowvec old_u = U_k.row(i);
    arma::mat H_ik(D, D, arma::fill::zeros);
    arma::vec g_ik(D, arma::fill::zeros);

    for (int t = 0; t < TT; ++t) {
      arma::rowvec r_tik = AU_list[t].row(i) % lambda_k.t();
      double b_tik = Q_cache[t] - arma::dot(old_u, r_tik) + StBk[t] - c_k[t];
      H_ik += omega_k[t] * (r_tik.t() * r_tik);
      g_ik += r_tik.t() * (delta_k[t] - b_tik * omega_k[t]);
    }

    arma::mat M_ik = H_ik + prior_prec * I_D; 
    arma::mat Uchol_M;
    bool chol_ok = arma::chol(Uchol_M, M_ik, "upper"); // explicit Cholesky status check.
    if (!chol_ok) {
      Rcpp::stop("Cholesky failed in update_class_cpp for the A-block U update."); 
    }
    arma::vec Minv_g =
      arma::solve(arma::trimatu(Uchol_M),
                  arma::solve(arma::trimatl(Uchol_M.t()), g_ik));

    int z_new = 1;
    if (use_spike_slab) {
      double logdet_M = 2.0 * arma::sum(arma::log(Uchol_M.diag()));
      double w_i = std::min(std::max(W_k[i], 1.0e-12), 1.0 - 1.0e-12); // avoid log(0) in spike/slab odds.

      double logW0 = std::log(w_i); 
      double logW1 = std::log(1.0 - w_i) - 
        0.5 * static_cast<double>(D) * log_tau_sq -
        0.5 * logdet_M +
        0.5 * arma::dot(g_ik, Minv_g);

      z_new = static_cast<int>(R::rbinom(1.0, plogis_stable_cpp(logW1 - logW0)));
    }
    Z_k[i] = z_new;

    arma::rowvec new_u(D, arma::fill::zeros);
    if (!use_spike_slab || z_new == 1) {
      arma::vec z_std(D);
      for (int d = 0; d < D; ++d) z_std[d] = R::rnorm(0.0, 1.0);
      new_u = (Minv_g + arma::solve(arma::trimatu(Uchol_M), z_std)).t();
    }

    arma::rowvec delta_u = new_u - old_u;
    U_k.row(i) = new_u;

    if (arma::accu(delta_u % delta_u) > 0.0) {
      for (int t = 0; t < TT; ++t) {
        arma::rowvec row_i_AU = AU_list[t].row(i);
        arma::rowvec dq = delta_u % row_i_AU;
        q_k.row(t) += dq;
        Q_cache[t] += arma::dot(lambda_k.t(), dq);
        AU_list[t] += A.slice(t).col(i) * delta_u;
      }
    }

    if (use_spike_slab) {
      W_k[i] = R::rbeta(a_w + 1.0 - static_cast<double>(Z_k[i]),
                        b_w + static_cast<double>(Z_k[i]));
    } else {
      W_k[i] = 0.0;
    }
  }

  return Rcpp::List::create(
    Rcpp::Named("U") = U_k,
    Rcpp::Named("Z") = Z_k,
    Rcpp::Named("W") = W_k,
    Rcpp::Named("q") = q_k,
    Rcpp::Named("Q") = Q_cache
  );
}

// ============================================================
// Cached class update for the hybrid R/C version.
// This keeps AU = A_t * U_k, q_k, and Q_k inside C++ across
// calls using an external pointer.  The goal is to avoid
// rebuilding the large AU cache for every Spike_Slab class update.
// ============================================================

struct BcmlrClassCache {
  arma::cube AU;  // n x D x TT, AU.slice(t) = A.slice(t) * U_k
  arma::mat q;   // TT x D, q[t,d] = 0.5 * sum_i U[i,d] * AU[i,d,t]
  arma::vec Q;   // TT, Q = q * lambda_k
};

inline void rebuild_class_cache_inplace(BcmlrClassCache* cache,
                                        const arma::cube& A,
                                        const arma::mat& U_k,
                                        const arma::vec& lambda_k) {
  const int n  = U_k.n_rows;
  const int D  = U_k.n_cols;
  const int TT = A.n_slices;

  cache->AU.set_size(n, D, TT);
  cache->q.set_size(TT, D);
  cache->q.zeros();

  for (int t = 0; t < TT; ++t) {
    cache->AU.slice(t) = A.slice(t) * U_k;
    cache->q.row(t) = 0.5 * arma::sum(U_k % cache->AU.slice(t), 0);
  }

  cache->Q = cache->q * lambda_k;
}

// [[Rcpp::export]]
SEXP make_class_cache_xptr_cpp(const arma::cube& A,
                               const arma::mat& U_k,
                               const arma::vec& lambda_k) {
  BcmlrClassCache* cache = new BcmlrClassCache();
  rebuild_class_cache_inplace(cache, A, U_k, lambda_k);
  Rcpp::XPtr<BcmlrClassCache> ptr(cache, true);
  return ptr;
}

// [[Rcpp::export]]
arma::vec get_class_cache_Q_cpp(SEXP cache_sexp) {
  Rcpp::XPtr<BcmlrClassCache> cache(cache_sexp);
  return cache->Q;
}

// [[Rcpp::export]]
arma::mat get_class_cache_q_cpp(SEXP cache_sexp) {
  Rcpp::XPtr<BcmlrClassCache> cache(cache_sexp);
  return cache->q;
}

// [[Rcpp::export]]
void set_class_cache_Q_cpp(SEXP cache_sexp,
                           const arma::vec& Q_new) {
  Rcpp::XPtr<BcmlrClassCache> cache(cache_sexp);
  if (cache->Q.n_elem != Q_new.n_elem) {
    Rcpp::stop("set_class_cache_Q_cpp: Q_new has incompatible length.");
  }
  cache->Q = Q_new;
}

// Optional exact rebuild. Use only for debugging or occasional drift checks.
// Do not call every iteration if speed is the goal.
// [[Rcpp::export]]
void rebuild_class_cache_xptr_cpp(SEXP cache_sexp,
                                  const arma::cube& A,
                                  const arma::mat& U_k,
                                  const arma::vec& lambda_k) {
  Rcpp::XPtr<BcmlrClassCache> cache(cache_sexp);
  rebuild_class_cache_inplace(cache.get(), A, U_k, lambda_k);
}

// [[Rcpp::export]]
Rcpp::List update_class_spikeslab_cached_xptr_cpp(const arma::cube& A,
                                                  arma::mat U_k,
                                                  arma::ivec Z_k,
                                                  arma::vec W_k,
                                                  const arma::vec& lambda_k,
                                                  const arma::vec& StBk,
                                                  const arma::vec& c_k,
                                                  const arma::vec& delta_k,
                                                  const arma::vec& omega_k,
                                                  const double tau_sq,
                                                  const double a_w,
                                                  const double b_w,
                                                  const bool use_spike_slab,
                                                  SEXP cache_sexp) {
  const int n  = U_k.n_rows;
  const int D  = U_k.n_cols;
  const int TT = A.n_slices;

  Rcpp::XPtr<BcmlrClassCache> cache(cache_sexp);

  // Rebuild only if dimensions are inconsistent. This should not happen in a
  // normal run, but it protects the external pointer path from stale cache shape.
  if (cache->AU.n_rows != static_cast<arma::uword>(n) ||
      cache->AU.n_cols != static_cast<arma::uword>(D) ||
      cache->AU.n_slices != static_cast<arma::uword>(TT) ||
      cache->q.n_rows != static_cast<arma::uword>(TT) ||
      cache->q.n_cols != static_cast<arma::uword>(D) ||
      cache->Q.n_elem != static_cast<arma::uword>(TT)) {
    rebuild_class_cache_inplace(cache.get(), A, U_k, lambda_k);
  }

  arma::cube& AU_cube = cache->AU;
  arma::mat& q_k = cache->q;
  arma::vec& Q_cache = cache->Q;

  const arma::mat I_D = arma::eye(D, D);
  const double prior_prec = 1.0 / tau_sq;
  const double log_tau_sq = std::log(tau_sq);
  const arma::rowvec lambda_row = lambda_k.t();

  for (int i = 0; i < n; ++i) {
    arma::rowvec old_u = U_k.row(i);
    arma::mat H_ik(D, D, arma::fill::zeros);
    arma::vec g_ik(D, arma::fill::zeros);

    for (int t = 0; t < TT; ++t) {
      arma::rowvec r_tik = AU_cube.slice(t).row(i) % lambda_row;
      double b_tik = Q_cache[t] - arma::dot(old_u, r_tik) + StBk[t] - c_k[t];

      H_ik += omega_k[t] * (r_tik.t() * r_tik);
      g_ik += r_tik.t() * (delta_k[t] - b_tik * omega_k[t]);
    }

    if (!H_ik.is_finite() || !g_ik.is_finite()) {
      Rcpp::stop("Non-finite H_ik or g_ik in update_class_spikeslab_cached_xptr_cpp.");
    }

    arma::mat M_ik = H_ik + prior_prec * I_D;
    arma::mat Uchol_M;

    bool chol_ok = arma::chol(Uchol_M, M_ik, "upper");
    if (!chol_ok) {
      Rcpp::stop("Cholesky failed in update_class_spikeslab_cached_xptr_cpp.");
    }

    arma::vec Minv_g =
      arma::solve(arma::trimatu(Uchol_M),
                  arma::solve(arma::trimatl(Uchol_M.t()), g_ik));

    int z_new = 1;
    if (use_spike_slab) {
      double logdet_M = 2.0 * arma::sum(arma::log(Uchol_M.diag()));
      double w_i = std::min(std::max(W_k[i], 1.0e-12), 1.0 - 1.0e-12);

      double logW0 = std::log(w_i);
      double logW1 = std::log(1.0 - w_i) -
        0.5 * static_cast<double>(D) * log_tau_sq -
        0.5 * logdet_M +
        0.5 * arma::dot(g_ik, Minv_g);

      z_new = static_cast<int>(R::rbinom(1.0, plogis_stable_cpp(logW1 - logW0)));
    }
    Z_k[i] = z_new;

    arma::rowvec new_u(D, arma::fill::zeros);
    if (!use_spike_slab || z_new == 1) {
      arma::vec z_std(D);
      for (int d = 0; d < D; ++d) z_std[d] = R::rnorm(0.0, 1.0);
      new_u = (Minv_g + arma::solve(arma::trimatu(Uchol_M), z_std)).t();
    }

    arma::rowvec delta_u = new_u - old_u;
    U_k.row(i) = new_u;

    if (arma::accu(delta_u % delta_u) > 0.0) {
      for (int t = 0; t < TT; ++t) {
        arma::rowvec row_i_AU = AU_cube.slice(t).row(i);
        arma::rowvec dq = delta_u % row_i_AU;

        q_k.row(t) += dq;
        Q_cache[t] += arma::dot(lambda_row, dq);

        // Rank-one cache update:
        // A_t * U_new = A_t * U_old + A_t[,i] * delta_u.
        AU_cube.slice(t) += A.slice(t).col(i) * delta_u;
      }
    }

    if (use_spike_slab) {
      W_k[i] = R::rbeta(a_w + 1.0 - static_cast<double>(Z_k[i]),
                        b_w + static_cast<double>(Z_k[i]));
    } else {
      W_k[i] = 0.0;
    }
  }

  return Rcpp::List::create(
    Rcpp::Named("U") = U_k,
    Rcpp::Named("Z") = Z_k,
    Rcpp::Named("W") = W_k,
    Rcpp::Named("q") = q_k,
    Rcpp::Named("Q") = Q_cache
  );
}

')

bcmlrNET = function(A, S, num_CP, latent_dim = 5, eta_dir = 1.5, b_tau = 1, 
                    prior_kappa = "default", alpha_f = 0.1, 
                    min_size = 30, sd_beta = 3, thinning = 1, 
                    num_iter = 5000, num_warmup = 2500, print_outputs = FALSE, 
                    print_progress = FALSE, model_selection = FALSE,
                    prior_for_U = c("Gaussian", "Spike_Slab"), sparse_networks = FALSE, 
                    center_std_A = TRUE, seed  = NULL,
                    lambda_init = c("equal", "prior", "zero")){

  if (!is.null(seed)) set.seed(seed)
  # S is a matrix with Tt rows and p columns of summary statistics.
  # A is an input series with Tt matrices in dimensions n \times n.
  #===============#
  # Sanity checks ####
  #===============#
  # A function that returns nothing and does sanity checks 
  # on the input list of adjacency matrices
  sanitycheck_on_A = function(A){
    # sanity checks on A: 
    if (!is.list(A)){ 
      stop("Please input a list of adjacency matrices.") 
    }else if (is.list(A)){
      if (length(A) == 1){
        stop("Please include more than one adjacency matrices in your list.") 
      }else if (length(A) > 1){
        # check if all adjacency matrices have the same dimensions (i.e., same number of nodes)
        first_dims = dim(A[[1]])
        all_dims_match = all(sapply(A, function (m) identical(dim(m), first_dims)))
        if (!all_dims_match){
          stop("All adjacency matrices in your list must have the same dimensions.")
        }
        all_finite = all(sapply(A, function (m) all(is.finite(as.matrix(m)))))
        if (!all_finite){
          stop("All adjacency matrices in your list must contain only finite values.") # prevent NA/Inf from entering the A-block updates.
        }
        all_zero_diag = all(sapply(A, function (m) all(abs(diag(as.matrix(m))) < 1e-8)))
        if (!all_zero_diag){
          stop("All adjacency matrices in your list must have zero diagonal entries.") #  the 1/2 <A, B> calculation assumes no self-loops.
        }
        all_symmetric = all(sapply(A, function (m) isTRUE(all.equal(as.matrix(m), t(as.matrix(m)), tolerance = 1e-8, check.attributes = FALSE))))
        if (!all_symmetric){
          stop("All adjacency matrices in your list must be symmetric.") # the A-block update assumes undirected networks.
        }
      } # end of the case if A_list contains more than one matrices
    } # end of the case if A is a list
  }
  # Track whether A or S was provided by the user before placeholder inputs are created.
  has_A = !is.null(A)
  has_S = !is.null(S)
  if (is.null(A) && is.null(S)){
    stop("Please enter a list of adjacency matrices or/and a matrix of network statistics.")
  }else if (is.null(A) && !is.null(S)){
    Tt_all = dim(S)[1]  # total number of observations
    A = list() # set A to zero
    for (t in 1:Tt_all){
      A[[t]] = matrix(0,1,1)
    }
  }else if (!is.null(A) && is.null(S)){
    Tt_all = length(A)  # total number of observations
    S = matrix(0,nrow = Tt_all, ncol = 1)  # set S to zero
    sanitycheck_on_A(A)
  }# end of the case if A is not NULL but S is NULL
  else if (!is.null(A) && !is.null(S)){
    Tt_all = length(A)  # total number of observations
    sanitycheck_on_A(A) 
    if (length(A) != dim(S)[1]){
      stop("The series of adjacency matrices (A) and summary statistics (S) must have the same length. ")
    }
  }
  n = dim(A[[1]])[2] # number of nodes
  p = dim(S)[2] # number of summary statistics
  # Tt_all = length(A)  # total number of observations
  D = latent_dim # dimensions of latent variables per node 
  prior_for_U = match.arg(prior_for_U)
  lambda_init = match.arg(lambda_init) 
  if (!is.logical(sparse_networks) || length(sparse_networks) != 1){
    stop("sparse_networks must be TRUE or FALSE.")
  }
  if (!is.logical(center_std_A) || length(center_std_A) != 1){
    stop("center_std_A must be TRUE or FALSE.")
  }
  # if (latent_dim > n){
  #   stop("The number of latent variables (D) cannot exceed the number of nodes (n).")
  # }
  #=================#
  # Data processing ####
  #=================#
  # A function to center & standardize a matrix x
  center_standardize = function(x, m, sdv){
    for (i in 1:ncol(x)){
      x[,i] = x[,i] - m[i]
      x[,i] = x[,i]/sdv[i]
    }
    return(x)
  }
  
  # A function that converts adjacency matrices to sparse matrices.
  to_sparse_list = function(A_ls){
    out_ls = vector(mode = "list", length = length(A_ls))
    for (tt in 1:length(A_ls)){
      if (methods::is(A_ls[[tt]], "Matrix")){
        out_ls[[tt]] = A_ls[[tt]]
      }else{
        out_ls[[tt]] = Matrix::Matrix(A_ls[[tt]], sparse = TRUE)
      }
    }
    return(out_ls)
  }
  
  # A function that converts a list of adjacency matrices to a matrix of upper-triangular entries.
  A_list_to_upper = function(A_ls, upper_idx){
    TT = length(A_ls)
    num_edge = sum(upper_idx)
    A_upper = matrix(0, nrow = TT, ncol = num_edge)
    if (num_edge == 0){
      return(A_upper)
    }
    for (tt in 1:TT){
      A_upper[tt,] = as.numeric(A_ls[[tt]][upper_idx])
    }
    return(A_upper)
  }
  
  # A function that converts a matrix of upper-triangular entries back to a list of symmetric matrices.
  upper_to_A_list = function(A_upper, upper_idx, n_node){
    TT = nrow(A_upper)
    A_ls = vector(mode = "list", length = TT)
    for (tt in 1:TT){
      A_t = matrix(0, nrow = n_node, ncol = n_node)
      if (ncol(A_upper) > 0){
        A_t[upper_idx] = A_upper[tt,]
        A_t = A_t + t(A_t)
      }
      diag(A_t) = 0
      A_ls[[tt]] = A_t
    }
    return(A_ls)
  }
  
  # A function that centers & standardizes adjacency matrices using training moments of the upper triangle.
  center_standardize_A = function(A_ls, m_A, sd_A, upper_idx, n_node){
    A_upper = A_list_to_upper(A_ls, upper_idx)
    if (ncol(A_upper) > 0){
      for (ee in 1:ncol(A_upper)){
        A_upper[,ee] = A_upper[,ee] - m_A[ee]
        A_upper[,ee] = A_upper[,ee] / sd_A[ee]
      }
    }
    return(upper_to_A_list(A_upper, upper_idx, n_node))
  }
  
  pack_A_cube = function(A_ls){
    n_local = nrow(A_ls[[1]])
    TT_local = length(A_ls)
    out = array(0, dim = c(n_local, n_local, TT_local))
    for (tt in 1:TT_local){
      out[,,tt] = as.matrix(A_ls[[tt]])
    }
    storage.mode(out) = "double"
    out
  }
  
  upper_idx_A = upper.tri(matrix(0, nrow = n, ncol = n), diag = FALSE)
  
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
    
    if (sparse_networks){
      A_all_raw = to_sparse_list(A)
    }else{
      A_all_raw = A
    }
    if (center_std_A){
      A_upper_train = A_list_to_upper(A_all_raw, upper_idx_A)
      if (ncol(A_upper_train) == 0){
        m_A = numeric(0)
        sd_A = numeric(0)
      }else{
        m_A = colMeans(A_upper_train, na.rm = TRUE)
        sd_A = apply(A_upper_train, 2, FUN = sd, na.rm = TRUE)
        sd_A[sd_A == 0] = 1
      }
      A_all = center_standardize_A(A_all_raw, m_A, sd_A, upper_idx_A, n)
      A_auc = A_all
      A = A_all
    }else{
      A_all = A_all_raw
      A_auc = A_all
      A = A_all
    }
  }
  # Set thinning >1 when selecting the number of changepoints (model_selection = TRUE)
  else if (thinning > 1){ 
    min_size = round(min_size*(1 - 1/thinning))
    Tt = dim(S)[1]
    all_idx = seq_len(Tt) # all indices from the input data
    holdout_idx = all_idx[all_idx%%thinning == 0] # the indices of the held-out samples
    train_idx = all_idx[-holdout_idx] # the indices for fitting the model (the trainng set indices)
    
    S_all = S
    S_auc = S_all[holdout_idx, , drop = FALSE] # rows of S selected to calculate AUC
    S = S_all[-holdout_idx, , drop = FALSE] # rows of S selected to fit the model  
    
    if (sparse_networks){
      A_all_raw = to_sparse_list(A)
    }else{
      A_all_raw = A
    }
    A_auc_raw = A_all_raw[holdout_idx]
    A_train_raw = A_all_raw[-holdout_idx]
    
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
    S_all = center_standardize(S_all, m, sdv)
    
    if (center_std_A){
      A_upper_train = A_list_to_upper(A_train_raw, upper_idx_A)
      if (ncol(A_upper_train) == 0){
        m_A = numeric(0)
        sd_A = numeric(0)
      }else{
        m_A = colMeans(A_upper_train, na.rm = TRUE)
        sd_A = apply(A_upper_train, 2, FUN = sd, na.rm = TRUE)
        sd_A[sd_A == 0] = 1
      }
      A = center_standardize_A(A_train_raw, m_A, sd_A, upper_idx_A, n)
      A_auc = center_standardize_A(A_auc_raw, m_A, sd_A, upper_idx_A, n)
      A_all = center_standardize_A(A_all_raw, m_A, sd_A, upper_idx_A, n)
    }else{
      A = A_train_raw
      A_auc = A_auc_raw
      A_all = A_all_raw
    }
  }
  
  Tt_auc = dim(S_auc)[1]
  Tt = dim(S)[1] # number of observations
  K = num_CP + 1 # number of classes
  L = num_CP # number of changepoints
  if (has_A){
    A_cpp = pack_A_cube(A)
    if (model_selection){
      A_auc_cpp = pack_A_cube(A_auc)
    }
  }
  
  #===================#
  # Utility functions ####
  #===================#
  # https://arxiv.org/pdf/2205.04997
  # https://gregorygundersen.com/blog/2020/02/09/log-sum-exp/
  # A function that mitigates computational issues in massive division or multiplication of extreme values.
  logsumexp = function(x){
    # Return -Inf when all entries are non-finite to avoid NaN in probability normalization.
    x_finite = x[is.finite(x)]
    if (length(x_finite) == 0){
      return(-Inf)
    }
    c = max(x_finite)
    return(c + log(sum(exp(x_finite - c))))
  }
  # Normalize log-weights safely for sample(); 
  # fall back to a uniform distribution if all weights are invalid.
  normalize_logweights = function(logweights){
    if (length(logweights) == 1){
      return(1)
    }
    if (all(!is.finite(logweights))){
      return(rep(1/length(logweights), length(logweights)))
    }
    max_logweights = max(logweights[is.finite(logweights)])
    shifted_weights = exp(logweights - max_logweights)
    shifted_weights[!is.finite(shifted_weights)] = 0
    sum_shifted = sum(shifted_weights)
    if (!is.finite(sum_shifted) || sum_shifted <= 0){
      return(rep(1/length(logweights), length(logweights)))
    }
    return(shifted_weights/sum_shifted)
  }
  # Robust AUC CI helper for model selection.
  # Undefined AUCs are stored as NA so final AUC summaries are not inflated.
  # Constant finite scores with both classes are recorded as non-informative AUC = 0.5.
  # Fix ROC direction because larger scores mean class 1; this permits AUC < 0.5.
  safe_ci_auc = function(Y_auc, scores, alpha_f){
    keep = is.finite(scores) & !is.na(Y_auc)
    Y_auc = as.numeric(Y_auc[keep])
    scores = as.numeric(scores[keep])
    if (length(Y_auc) < 2 || length(unique(Y_auc)) < 2){
      return(c(NA_real_, NA_real_, NA_real_)) 
    }
    if (length(unique(scores)) < 2){
      return(c(0.5, 0.5, 0.5)) 
    }
    ci_tmp = tryCatch(
      as.numeric(suppressWarnings(ci.auc(Y_auc, scores, levels = c(0, 1), direction = "<", conf.level = 1-alpha_f, quiet = TRUE))), 
      error = function(e) c(NA_real_, NA_real_, NA_real_) 
    )
    if (length(ci_tmp) != 3 || any(!is.finite(ci_tmp))){
      ci_tmp = c(NA_real_, NA_real_, NA_real_) 
    }
    return(ci_tmp)
  }
  # A function that converts a change point vector to a response matrix with the change point structure.
  kappa_to_Y = function(kappa, N, L){
    class = rep(1, kappa[1])
    for(j in 2:L){
      class = c(class, rep(j, kappa[j] - kappa[j-1]))
    }
    class = c(class, rep(L+1, N - kappa[L]))
    class = as.factor(class)
    Y = model.matrix(~ -1 + class, data=class)
    return(Y)
  }
  
  # A function that caches A_t %*% U_k, q_tdk, Q_tk, and S_t^T beta_k for a fixed class.
  build_class_cache = function(A_ls, S_mat, U_k, Lambda_k, Beta_k){
    TT = length(A_ls)
    AU_list = vector(mode = "list", length = TT)
    q_k = matrix(0, nrow = TT, ncol = D)
    for (t in 1:TT){
      AU_t = as.matrix(A_ls[[t]] %*% U_k)
      AU_list[[t]] = AU_t
      q_k[t,] = 0.5 * colSums(U_k * AU_t)
    }
    Q_k = as.numeric(q_k %*% Lambda_k)
    SB_k = as.numeric(S_mat %*% Beta_k)
    return(list(AU = AU_list, q = q_k, Q = Q_k, SB = SB_k))
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
  out$pi = array(0, dim = c(num_iter-num_warmup, D, 3)) # store dimension-specific pi_d = (Pr(0), Pr(1), Pr(-1)) for each latent dimension d.
  dimnames(out$pi) = list(NULL, paste0("dim_", 1:D), c("zero", "pos1", "neg1")) # make the lambda/probability ordering explicit.
  #out$omega = array(0, dim = c(num_iter-num_warmup, Tt, K))
  #out$phi = array(0, dim = c(num_iter-num_warmup, Tt, K))
  if (num_CP==1){
    out$Kappa = rep(0, num_iter-num_warmup)
  }else if (num_CP>1){
    out$Kappa = matrix(data = 0, nrow = num_iter-num_warmup, ncol = L)
  }
  if (model_selection){
    if (num_CP==1){
      out$AUC = rep(0, num_iter-num_warmup)
    }else if (num_CP>1){
      out$AUC = matrix(0, nrow = num_iter-num_warmup, ncol = L)
    }
  }else{
    out$AUC = NULL
  }
  #============================#
  # hyperparameters for priors ####
  #============================#
  # Prior quantities for beta (coefficients of summary statistics):
  m0 = array(data=0, dim=c(p,K-1)) # m0[,k] is the prior mean of the kth coef vec
  V0 = array(data=0, dim=c(p,p,K-1)) # V0[,,k] is the prior cov. of the kth coef vec
  V0inv = array(data=0, dim=c(p,p,K-1)) # V0inv[,,k] is the inverse of V0[,,j]
  sd0 = sd_beta # prior variance for beta (coefficients of summary statistics)
  for(k in 1:(K-1)){
    V0[,,k] = sd0*diag(p)
    V0inv[,,k] = 1/sd0*diag(p)
  }
  # priors for tau_sq
  a_tau = 4
  b_tau = b_tau
  # priors for w_ik ~ Beta(a_w, b_w)
  a_w = 1
  b_w = 1
  # prior for pi (Dirichlet)
  eta_dir = eta_dir # must be >1
  lambda_support = c(0, 1, -1) # : order pi as (Pr(lambda = 0), Pr(lambda = 1), Pr(lambda = -1)); zero is the favored value.
  #===============================#
  # Initialize working quantities ####
  #===============================#
  # Q = matrix(0, nrow = Tt, ncol = K) # Frobenius inner products
  tau_sq = rigamma(n = 1, a = a_tau, b = b_tau) # inverse gamma variable
  Z = matrix(0, nrow = n, ncol = K) # Bernoulli variables
  Z[,K] = rep(0, n)
  U = array(0, dim = c(K, n, D)) # latent variables
  if (prior_for_U == "Gaussian"){
    for (i in 1:n){
      for (k in 1:(K-1)){
        Z[i,k] = 1
        U[k,i,] = mvrnorm(n =1, mu = rep(0,D), Sigma = tau_sq*diag(D))
      }
    }
    W = matrix(0, nrow = n, ncol = K) # Beta variables are unused under the Gaussian prior for U
  }else if (prior_for_U == "Spike_Slab"){
    Z = matrix(rbinom(n = n*K, size = 1, prob = 0.5), nrow = n, ncol = K) # Bernoulli variables
    Z[,K] = rep(0, n)
    for (i in 1:n){
      for (k in 1:(K-1)){
        if (Z[i,k] ==1){
          U[k,i,] = mvrnorm(n =1, mu = rep(0,D), Sigma = tau_sq*diag(D))
        }else if(Z[i,k] ==0){
          U[k,i,] = rep(0,D)
        }
      }
    }
    W = matrix(rbeta(n*K, shape1 = a_w, shape2 = b_w), nrow = n, ncol = K) # Beta variables
  }
  Beta = matrix(0, nrow = p, ncol = K) # coefficients for summary statistics
  pi = matrix(0, nrow = D, ncol = 3) # pi[d,] is dimension-specific, not shared across dimensions.
  colnames(pi) = c("zero", "pos1", "neg1") # explicit ordering matches lambda_support = c(0, 1, -1).
  for (d in 1:D){
    pi[d,] = as.numeric(DirichletReg::rdirichlet(n = 1, alpha = c(d^eta_dir, 1, 1))) # Dirichlet(d^eta, 1, 1) increasingly favors lambda = 0 as d increases.
  }
  # Lambda imposes sparsity on the latent variables. 
  Lambda = matrix(0, nrow = K, ncol = D) # initialize dimension-by-dimension so lambda_dk can use pi[d,] when requested.
  for (k in 1:(K-1)){
    for (d in 1:D){
      # We use equal-probability initialization by default, but make alternatives available for sensitivity checks.
      if (lambda_init == "equal"){
        # Initialize with an equal probability for each value
        Lambda[k,d] = sample(x = lambda_support, size = 1, prob = rep(1/3, 3))
      }else if (lambda_init == "prior"){
        # Initialize with its prior: sample from {0, 1, -1}, with dimension-specific pi_d.
        Lambda[k,d] = sample(x = lambda_support, size = 1, prob = pi[d,])
      }else if (lambda_init == "zero"){
        # Initialize from the fully sparse state.
        Lambda[k,d] = 0
      }
    }
  }
  Lambda[K,] = rep(0, D) # make the last reference class be zero. 
  omega = matrix(0, nrow = Tt, ncol = K)
  #if (init == "even"){
  # Initialize CPs at time points that equally split the series 
  if(num_CP==1){
    kappa = round(Tt/2)
  }else if (num_CP>1){
    kappa = rep(0, L)
    for (l in 1:L){
      kappa[l] = round(Tt/K) *l
    }
  }
  #}
  original_kappa = train_idx[kappa]
  #=============================#
  # Initialize class caches ####
  #=============================#
  Q_current = matrix(0, nrow = Tt, ncol = K)
  SB_current = matrix(0, nrow = Tt, ncol = K)
  # C++ external-pointer caches for the network block.
  # These are used only for the Spike_Slab path so AU = A_t %*% U_k
  # does not need to be rebuilt from scratch at every class update.
  class_cache = vector(mode = "list", length = K - 1)
  for (k in 1:(K-1)){
    if (has_A){
      U_k_init = matrix(U[k,,], nrow = n, ncol = D)
      if (prior_for_U == "Spike_Slab"){
        class_cache[[k]] = make_class_cache_xptr_cpp(A_cpp, U_k_init, Lambda[k,])
        Q_current[,k] = get_class_cache_Q_cpp(class_cache[[k]])
      }else{
        Q_current[,k] = compute_Q_vec_cube_cpp(A_cpp, U_k_init, Lambda[k,])
      }
    }
    if (has_S){
      SB_current[,k] = as.numeric(S %*% Beta[,k])
    }
  }
  if (model_selection){
    Q_auc_current = matrix(0, nrow = Tt_auc, ncol = K)
    if (has_A){
      for (k in 1:(K-1)){
        Q_auc_current[,k] = compute_Q_vec_cube_cpp(A_auc_cpp, matrix(U[k,,], nrow = n, ncol = D), Lambda[k,])
      }
    }
  }
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
    #===========================#
    # Update Q_tk per iteration #
    #===========================#
    Q_iter = Q_current
    StB = SB_current
    phi = Q_iter + StB
    # update quantities in class k
    # Begin for-loop of k = 1, ..., (K-1)####
    for (k in 1:(K-1)){
      # After each parameter in class k is updated (in each iteration of k), 
      # get the latest phi_{tk}, c_{tk} and eta_{tk} based on all classes.
      #===========================#
      # Compute c_{k} and eta_{k} ####
      #===========================#
      if(num_CP ==1){
        c_k = rep(0, Tt)
      }else if (num_CP>1){
        c_k = apply(phi[, -k, drop = FALSE], 1, logsumexp) # Tt * 1 , less stable: log(rowSums(exp(phi[,-k]))) 
      }
      eta_k = phi[,k] - c_k
      #============================#
      # Update PG variables, omega ####
      #============================#
      omega[,k] = rpg(num = Tt, h = 1, z = eta_k) # update omega using BayesLogit::rpg
      # (The current loop is for class k.)
      U_k = matrix(U[k,,], nrow = n, ncol = D)
      Lambda_k = Lambda[k,]
      StBk_vec = StB[,k]
      #===============================================#
      # Update z_{ik}, u_{ik}, w_{ik} for each node i ####
      #===============================================#
      if (has_A){
        if (prior_for_U == "Spike_Slab"){
          # Cached Spike_Slab class update.
          # The C++ external pointer stores AU, q, and Q for this class.
          # This avoids rebuilding A_t %*% U_k from scratch at every class update.
          class_upd = update_class_spikeslab_cached_xptr_cpp(
            A = A_cpp,
            U_k = U_k,
            Z_k = as.integer(Z[,k]),
            W_k = W[,k],
            lambda_k = Lambda_k,
            StBk = StBk_vec,
            c_k = c_k,
            delta_k = delta[,k],
            omega_k = omega[,k],
            tau_sq = tau_sq,
            a_w = a_w,
            b_w = b_w,
            use_spike_slab = TRUE,
            cache_sexp = class_cache[[k]]
          )
        }else{
          # Keep the Gaussian path exactly as before.
          class_upd = update_class_cpp(
            A = A_cpp,
            U_k = U_k,
            Z_k = as.integer(Z[,k]),
            W_k = W[,k],
            lambda_k = Lambda_k,
            StBk = StBk_vec,
            c_k = c_k,
            delta_k = delta[,k],
            omega_k = omega[,k],
            tau_sq = tau_sq,
            a_w = a_w,
            b_w = b_w,
            use_spike_slab = FALSE
          )
        }
        U_k = class_upd$U
        U[k,,] = U_k
        Z[,k] = class_upd$Z
        W[,k] = class_upd$W
        q_k = class_upd$q
        Q_k = class_upd$Q
        #=======================# 
        # Update Lambda_k (D*D) ####
        #=======================#
        Lambda_k = Lambda[k,]
        for (d in 1:D){
          lambda_dk_old = Lambda_k[d]
          alpha_k_minus_d = Q_k - lambda_dk_old*q_k[,d] + StBk_vec - c_k   # Tt*1 
          cand_vals = lambda_support # support is ordered as c(0, 1, -1), matching pi[d,] = (Pr(0), Pr(1), Pr(-1)).
          logweights = rep(0, length(cand_vals))
          for (r in 1:length(cand_vals)){
            v = cand_vals[r] # evaluate each candidate value
            eta_cand_k = alpha_k_minus_d + v*q_k[,d]   # Tt*1  # use candidate v instead of lambda in the math expression
            # element-wise multiplication between Tt*1 vectors in the sum():
            logweights[r] = log(pi[d,r]) + sum(delta[,k]*eta_cand_k - 0.5*omega[,k]*eta_cand_k^2) # use dimension-specific prior mass for lambda_dk.
          } 
          # Normalize log-weights safely so sample() never receives NA probabilities.
          probs_lambda = normalize_logweights(logweights) # 3*1 vector that should sum up to 1
          #probs_lambda = exp(logweights)/(sum(exp(logweights)))
          # Update lambda_dk
          lambda_dk_new = sample(cand_vals, size = 1, prob = probs_lambda)
          Lambda_k[d] = lambda_dk_new
          # Update Q_k as lambda_dk is updated
          Q_k = Q_k + (lambda_dk_new - lambda_dk_old)*q_k[,d] # Q_t is Tt*1 and q_t[,d] is also Tt*1
        }# End for-loop of d
        # Lambda updates change Q_k but not AU or q_k.
        # Therefore, only Q_k must be synchronized back into the C++ cache.
        if (prior_for_U == "Spike_Slab"){
          set_class_cache_Q_cpp(class_cache[[k]], Q_k)
        }
      }else{
        # If A is NULL, keep the network block fixed at zero.
        Z[,k] = rep(0, n)
        U_k = matrix(0, nrow = n, ncol = D)
        U[k,,] = U_k
        W[,k] = rep(0, n)
        q_k = matrix(0, nrow = Tt, ncol = D)
        Q_k = rep(0, Tt)
        Lambda_k = rep(0, D)
      }
      Lambda[k,] = Lambda_k
      # Done updating lambda_dk for each d (k is still fixed in the current iteration)
      #=====================# 
      # Update beta_k (p*1) ####
      #=====================#
      if (has_S){
        normals = matrix(rnorm(p), p, 1) # A p-vector of standard normals
        precision_beta = crossprod(x = S, y = S * omega[,k]) + V0inv[,,k] # no need for as.numeric(), keep matrix precision, no vector recycling
        Uchol_beta = chol(precision_beta) # Upper Cholesky factor
        # crossprod(S, delta[, k] + omega[,k] * (c_k - Q_k)) + V0inv[,,k] %*% m0[, k]
        # Beta[,k] = chol2inv(Uchol_beta) %*% 
        #   (t(S) %*% (delta[,k] + omega[,k] * (c_k - Q_k)) + V0inv[,,k] %*% m0[, k]) + # used * directly instead of using diag(omega[,k]) %*% vector
        #   backsolve(Uchol_beta, normals)
        beta_mean = backsolve(Uchol_beta, forwardsolve(t(Uchol_beta), t(S) %*% (delta[,k] + omega[,k] * (c_k - Q_k)) + V0inv[,,k] %*% m0[, k]))
        Beta[,k] = as.numeric(beta_mean + backsolve(Uchol_beta, normals))
        StBk_vec = as.numeric(S %*% Beta[,k])
      }else{
        # If S is NULL, keep Beta fixed at zero.
        Beta[,k] = rep(0, p)
        StBk_vec = rep(0, Tt)
      }
      if (model_selection && (iter > num_warmup) ){
        if (has_A){
          Q_auc_current[,k] = compute_Q_vec_cube_cpp(A_auc_cpp, U_k, Lambda_k)
        }
      }
      # Save FINAL Q_tk for the current class k after U and Lambda updates
      Q_iter[, k] = Q_k
      StB[, k] = StBk_vec
      phi[, k] = Q_k + StBk_vec
    } # End of for loop for k = 1, ..., (K-1) ####
    #===========================================# 
    # Update pi_d for each latent dimension d ####
    #===========================================#
    if (has_A){
      # pi is dimension-specific.  For each d, pi[d,] = (Pr(lambda_dk = 0), Pr(lambda_dk = 1), Pr(lambda_dk = -1)).
      # the reference class K has Lambda[K,] fixed at zero for identifiability, so it is not counted as a random draw from pi_d.
      for (d in 1:D){
        n_0 = sum(Lambda[1:(K-1), d] == 0)
        n_pos1 = sum(Lambda[1:(K-1), d] == 1)
        n_neg1 = sum(Lambda[1:(K-1), d] == -1)
        alpha_pi = c(d^eta_dir + n_0, 1 + n_pos1, 1 + n_neg1) # Dirichlet(d^eta, 1, 1) favors zero more strongly as d increases.
        pi[d,] = as.numeric(DirichletReg::rdirichlet(n = 1, alpha = alpha_pi)) # update pi_d separately for each latent dimension.
      }
      # Dirichlet draw via Gammas for one dimension d, if needed:
      # gamma_draw = rgamma(3, shape = alpha_pi, rate = 1)
      # pi[d,] = gamma_draw / sum(gamma_draw)
      #===========================================# 
      # Update tau_sq (shared across all classes) ####
      #===========================================#
      if (prior_for_U == "Spike_Slab"){
        sum_m_k = sum(Z[,1:(K-1), drop = FALSE])
      }else if (prior_for_U == "Gaussian"){
        sum_m_k = n * (K-1)
      }
      sum_ss_k = sum(U[1:(K-1), , , drop = FALSE]^2)
      #tau_sq = 1/rgamma(n=1, shape = a_tau + (D/2) *sum_m_k, rate =b_tau + 0.5*sum_ss_k)
      tau_sq = rigamma(n = 1, a = a_tau + (D/2) *sum_m_k, b = b_tau + 0.5*sum_ss_k)
    }else{
      # If A is NULL, keep pi and tau_sq fixed at zero.
      pi = matrix(0, nrow = D, ncol = 3) # keep the same dimension-specific pi shape when A is absent.
      colnames(pi) = c("zero", "pos1", "neg1") # explicit ordering is retained for downstream output.
      tau_sq = 0
    }
    #} # The updates of tau_sq and pi should not be included in the class update for-loop 
    # Done updating all other parameters except the changepoints
    
    #=============================#
    # Update kappa (changepoints) ####
    #=============================#
    # recompute phi using current U, Lambda, Beta before updating kappa
    phi = matrix(0, nrow = Tt, ncol = K) 
    for (kk in 1:(K-1)){ 
      phi[,kk] = Q_iter[,kk] + as.numeric(S %*% Beta[,kk]) 
    } 
    if (num_CP == 1){
      # Phi_exp = exp(phi[,1]) 
      # P = Phi_exp/(1+Phi_exp)
      P = plogis(phi[,1]) # numerically stable logistic probabilities
      # range = 1:(N-1)
      range = (1 + min_size):(Tt - 1 - min_size) 
      lr = length(range)
      log_prob_vec_prop = rep(0, lr) 
      if (is.character(prior_kappa) && prior_kappa == "default"){
        for(l in 1:lr){
          log_prob_vec_prop[l] = sum(log(1-P[1:range[l]])) + sum(log(P[(range[l]+1):Tt])) - 
            length(1:range[l])*log(length(1:range[l])) - 
            length((range[l]+1):Tt)*log(length((range[l]+1):Tt)) 
        }
      }else if (is.function(prior_kappa)){
        for(l in 1:lr){
          jt_pmf_kappa = prior_kappa(range[l]) # evaluate custom prior at candidate kappa
          log_prob_vec_prop[l] = sum(log(1-P[1:range[l]])) + sum(log(P[(range[l]+1):Tt])) + log(jt_pmf_kappa)
        }
      }
      prob_vec = normalize_logweights(log_prob_vec_prop) 
      if(lr > 1){
        kappa = sample(x=range, size=1, prob=prob_vec)
      }
      #=============================================#
      # Select num of CPs: Is there a changepoint ? ####
      #=============================================#
      original_kappa = train_idx[kappa]
      # 4. Compute AUC
      if (model_selection && (iter > num_warmup) ){
        Y_all = c(rep(0, times = original_kappa), rep(1, times = Tt_all - original_kappa))
        Y_auc = Y_all[holdout_idx]
        phi_auc = matrix(0, nrow = Tt_auc, ncol = K)
        for (kk in 1:(K-1)){
          phi_auc[,kk] = Q_auc_current[,kk] + as.numeric(S_auc %*% Beta[,kk])
        }
        # Update the matrix of probabilities of success based on the fitted coefficients
        P_auc = plogis(phi_auc[,1])
        ci = safe_ci_auc(Y_auc, P_auc, alpha_f) 
      }
    }else if(num_CP > 1){
      # Create the P matrix
      # Phi_exp = exp(phi) # Tt * K
      # P = t(apply(Phi_exp, 1, function(x) x/sum(x))) # Normalize the rows of Phi_exp to get P
      P = t(apply(phi, 1, function(x) normalize_logweights(x))) # stable softmax with a safe fallback when a row is non-finite
      if (is.character(prior_kappa) && prior_kappa == "default"){
        for(j in 1:(K-1)){
          # Now calculate the range of values that kappa[j] can take under its full
          # conditional distribution and the vector prob_vec, which is defined
          # such that prob_vec[k] is the log of the conditional probability that
          # kappa[j] = range[k] given all the other data and parameters.
          if(j==1){
            #range = 1:(kappa[j+1]-1)
            range = (1 + min_size):(kappa[j+1] - 1 - min_size) 
            lr = length(range)
            log_prob_vec_prop = rep(0, lr)
            for(l in 1:lr){
              log_prob_vec_prop[l] = sum(log(P[1:range[l],j])) + sum(log(P[(range[l]+1):kappa[j+1],j+1])) - 
                length(1:range[l])*log(length(1:range[l])) - 
                length((range[l]+1):kappa[j+1])*log(length((range[l]+1):kappa[j+1])) 
            }
            #prob_vec = exp(log_prob_vec_prop)/sum(exp(log_prob_vec_prop))
            prob_vec = normalize_logweights(log_prob_vec_prop)
          } else if(j==L){
            #range = (kappa[j-1]+1):(Tt-1)
            range = (kappa[j-1] + 1 + min_size):(Tt - 1 - min_size) 
            lr = length(range)
            log_prob_vec_prop = rep(0, lr)
            for(l in 1:lr){
              log_prob_vec_prop[l] = sum(log(P[(kappa[j-1]+1):range[l],j])) + sum(log(P[(range[l]+1):Tt,j+1])) - 
                length((kappa[j-1]+1):range[l])*log(length((kappa[j-1]+1):range[l])) - 
                length((range[l]+1):Tt)*log(length((range[l]+1):Tt)) 
            }
            #prob_vec = exp(log_prob_vec_prop)/sum(exp(log_prob_vec_prop))
            prob_vec = normalize_logweights(log_prob_vec_prop)
          } else{
            #range = (kappa[j-1]+1):(kappa[j+1]-1)
            range = (kappa[j-1] + 1 + min_size):(kappa[j+1] - 1 - min_size) 
            lr = length(range)
            log_prob_vec_prop = rep(0, lr)
            for(l in 1:lr){
              log_prob_vec_prop[l] = sum(log(P[(kappa[j-1]+1):range[l],j])) + sum(log(P[(range[l]+1):kappa[j+1],j+1])) - 
                length((kappa[j-1]+1):range[l])*log(length((kappa[j-1]+1):range[l])) - 
                length((range[l]+1):kappa[j+1])*log(length((range[l]+1):kappa[j+1])) 
            }
            #prob_vec = exp(log_prob_vec_prop)/sum(exp(log_prob_vec_prop))
            prob_vec = normalize_logweights(log_prob_vec_prop)
          }
          if(lr > 1){
            kappa[j] = sample(x=range, size=1, prob=prob_vec) # update kappa
          }
        }
      }# End of the case when we use default prior for kappa (based on segment lengths) as in the paper.
      else if (is.function(prior_kappa)){ # else if the user have a joint PMF for the kappas 
        for(j in 1:(K-1)){
          if(j==1){
            range = (1 + min_size):(kappa[j+1] - 1 - min_size) 
            lr = length(range)
            log_prob_vec_prop = rep(0, lr)
            for(l in 1:lr){
              kappa_prop = kappa # evaluate custom prior at candidate kappa vector
              kappa_prop[j] = range[l] 
              jt_pmf_kappa = prior_kappa(kappa_prop) 
              log_prob_vec_prop[l] = sum(log(P[1:range[l],j])) + sum(log(P[(range[l]+1):kappa[j+1],j+1])) + log(jt_pmf_kappa)
            }
            prob_vec = normalize_logweights(log_prob_vec_prop)
          } else if(j==L){
            range = (kappa[j-1] + 1 + min_size):(Tt - 1 - min_size) 
            lr = length(range)
            log_prob_vec_prop = rep(0, lr)
            for(l in 1:lr){
              kappa_prop = kappa # evaluate custom prior at candidate kappa vector
              kappa_prop[j] = range[l] 
              jt_pmf_kappa = prior_kappa(kappa_prop) 
              log_prob_vec_prop[l] = sum(log(P[(kappa[j-1]+1):range[l],j])) + sum(log(P[(range[l]+1):Tt,j+1])) + log(jt_pmf_kappa)
            }
            prob_vec = normalize_logweights(log_prob_vec_prop)
          } else{
            range = (kappa[j-1] + 1 + min_size):(kappa[j+1] - 1 - min_size) 
            lr = length(range)
            log_prob_vec_prop = rep(0, lr)
            for(l in 1:lr){
              kappa_prop = kappa 
              kappa_prop[j] = range[l] 
              jt_pmf_kappa = prior_kappa(kappa_prop) 
              log_prob_vec_prop[l] = sum(log(P[(kappa[j-1]+1):range[l],j])) + sum(log(P[(range[l]+1):kappa[j+1],j+1])) + log(jt_pmf_kappa)
            }
            prob_vec = normalize_logweights(log_prob_vec_prop)
          }
          if(lr > 1){
            kappa[j] = sample(x=range, size=1, prob=prob_vec) # update kappa
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
      if (model_selection && (iter > num_warmup) ){
        #=============#
        # Compute phi ####
        #=============#
        phi_auc = matrix(0, nrow = Tt_auc, ncol = K) # Fitted values
        for (kk in 1:(K-1)){
          phi_auc[,kk] = Q_auc_current[,kk] + as.numeric(S_auc %*% Beta[,kk])
        }
        # "ratios" are the probabilities of a subject belonging to
        # class j+1 other than class j
        # Pairwise AUC scores below are computed directly from phi_auc with plogis().
        # This avoids possible 0/0 ratios after full K-class softmax underflow.
        Y_auc_all = kappa_to_Y(original_kappa, Tt_all, L)[holdout_idx, , drop = FALSE] # Y in Tt_auc *K
        ci = matrix(0, nrow = K-1, ncol = 3) # each row has a LB, point estimate, and UB
        for (j in 1:(K-1)){
          if(j == 1){
            rg = 1:original_kappa[j+1]
            idx_test = which(holdout_idx %in% rg)
            ratios = plogis(phi_auc[idx_test, j+1] - phi_auc[idx_test, j]) 
            Y_auc = Y_auc_all[idx_test, j + 1]
            ci[j,] = safe_ci_auc(Y_auc, ratios, alpha_f) 
          }else if (j == K-1){
            rg = (original_kappa[j-1]+1):Tt_all
            idx_test = which(holdout_idx %in% rg)
            ratios = plogis(phi_auc[idx_test, j+1] - phi_auc[idx_test, j]) 
            Y_auc = Y_auc_all[idx_test, j + 1]
            ci[j,] = safe_ci_auc(Y_auc, ratios, alpha_f) 
          }else{
            rg = (original_kappa[j-1]+1):original_kappa[j+1]
            idx_test = which(holdout_idx %in% rg)
            ratios = plogis(phi_auc[idx_test, j+1] - phi_auc[idx_test, j]) 
            Y_auc = Y_auc_all[idx_test, j + 1]
            ci[j,] = safe_ci_auc(Y_auc, ratios, alpha_f) 
          }
        }
      } # End of selecting the number of changepoints given that num_CP >1
    } # END of the case "num_CP > 1"
    Q_current = Q_iter
    SB_current = StB
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
      out$pi[iter-num_warmup, , ] = pi # store the D x 3 matrix of dimension-specific pi_d values.
      #out$omega[iter-num_warmup,,] = omega
      #out$phi[iter-num_warmup,,] = phi
      
      # Store changepoints and AUCs: 
      if (num_CP==1){
        out$Kappa[iter-num_warmup] = original_kappa
        if (model_selection){
          out$AUC[iter-num_warmup] = ci[1] # Take the CI lower bound
        }
      }else if (num_CP >1){
        out$Kappa[iter-num_warmup,] = original_kappa
        if (model_selection){
          out$AUC[iter-num_warmup,] = ci[,1] # Take the CI lower bound
        }
      }
    } # Done storing samples in the current iteration if it's after burn-in
    
    if (print_progress){ setTxtProgressBar(pb, iter)}
  } # End of the Gibbs sampler (iter = 1, ..., num_iter)
  end = Sys.time()
  if (print_progress){ close(pb) }
  #================================#
  # END Gibbs sampling (iter)   ####
  #================================#
  
  # Store posterior mode changepoints: 
  if (num_CP == 1){
    # Use posterior modes as changepoint estimates:
    cp = table(out$Kappa)
    out$Kappa_mode = as.integer(names(cp)[which.max(cp)]) # posterior mode
  }else if (num_CP > 1){
    # Use posterior modes as changepoint estimates:
    cps = c()
    for (l in 1:L){
      cp = table(out$Kappa[,l])
      cp_mode = as.integer(names(cp)[which.max(cp)])
      cps = c(cps, cp_mode)
    }
    out$Kappa_mode = cps # posterior modes
  }
  # out$runTime = end-start
  return(out)
}
