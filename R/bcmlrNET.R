#' Fit changepoints in a dynamic network sequence
#'
#' Fit the generalized Bayesian model conditional on a positive number of
#' changepoints. C++ routines accelerate the network-parameter updates.
#'
#' @param A List of equally sized, finite, symmetric adjacency matrices with
#'   zero diagonals, in chronological order. Supply A, S, or both.
#' @param S Optional numeric matrix (or all-numeric data frame) with one row
#'   per observation and one column per network statistic.
#' @param num_CP Positive integer number of fitted changepoints. Zero is not
#'   a supported direct fit; use [bcmlrNET_model_select()] for selection.
#' @param latent_dim Number of latent dimensions per node.
#' @param eta_dir Finite exponent in the dimension-specific Dirichlet prior.
#' @param b_tau Positive inverse-gamma scale for the shared latent variance.
#' @param prior_kappa Either `"default"` (the segment-length prior) or a function
#'   returning a finite nonnegative prior mass for a changepoint vector.
#' @param alpha_f Significance level for AUC confidence intervals during selection.
#' @param min_size Nonnegative segment-spacing parameter. The existing candidate
#'   bounds require more than `min_size` observations per training segment.
#'   When `thinning > 1`, the working value is
#'   `round(min_size * (1 - 1/thinning))`. Changepoints index the last observation
#'   in the preceding segment. This is not an exact original-time minimum after
#'   holding observations out.
#' @param sd_beta Positive prior variance (not standard deviation) for each
#'   network-statistic coefficient. The historical argument name is retained.
#' @param thinning Positive integer holdout spacing, not posterior-sample
#'   thinning. Values above one hold out every `thinning`-th observation.
#'   With one, all observations are used for fitting and evaluation.
#' @param num_iter Total Gibbs iterations, including warmup.
#' @param num_warmup Number of initial iterations discarded.
#' @param print_outputs Print the structure of the returned list?
#' @param print_progress Display a progress bar?
#' @param model_selection Compute AUC confidence-interval lower bounds for
#'   neighboring-segment comparisons? Usually set by the selection wrapper.
#' @param prior_for_U Prior for latent node parameters: `"Gaussian"` or
#'   `"Spike_Slab"`. This direct-fitting function defaults to `"Gaussian"`.
#' @param sparse_networks Use Matrix objects during preprocessing? The current
#'   C++ implementation still receives a dense adjacency cube.
#' @param center_std_A Center and standardize dyads using fitting observations?
#' @param seed Optional integer seed passed to `set.seed()`.
#' @param lambda_init Latent-weight initialization: `"equal"`, `"prior"`, or `"zero"`.
#' @returns A list containing posterior samples `Kappa`, `Beta`, `U`, `Lambda`,
#'   `Z`, `W`, `pi`, `tau_sq`, and point estimates `Kappa_mode` (marginal
#'   posterior modes). `AUC`, when requested, contains confidence-interval
#'   lower bounds, not AUC point estimates. The number of retained samples is
#'   `num_iter - num_warmup`. `Kappa` is a vector for one changepoint and a
#'   samples-by-changepoints matrix otherwise; indices use the original series.
#'   `Beta` is samples-by-statistics-by-segments; `U` is
#'   samples-by-segments-by-nodes-by-latent-dimensions; `Lambda` is
#'   samples-by-segments-by-latent-dimensions. `Z` and `W` are
#'   samples-by-nodes-by-segments; `pi` is samples-by-latent-dimensions-by-three.
#' @details Segment-indexed samples are returned in chronological order. For
#'   one changepoint the first segment (`Y = 0`) is the reference; for multiple
#'   changepoints the last segment is the reference. Reference-slot `W` values
#'   are unused zero placeholders, not estimated probabilities. Posterior
#'   samples from absent network or statistic blocks are placeholders too.
#'   Separate marginal changepoint modes need not form a valid joint partition.
#'   The original rounded initializer and segment bounds are retained; an
#'   inadmissible initialization stops rather than silently changing settings.
#' @seealso [bcmlrNET_model_select()], [plot_changepoints()]
#' @export
# C++ code lives in src/bcmlrNET.cpp; sourceCpp/library calls removed.
# BayesLogit, LearnBayes and MASS calls use explicit namespaces.
bcmlrNET = function(A = NULL, S = NULL, num_CP, latent_dim = 10, eta_dir = 1.01, b_tau = 0.5, 
                    prior_kappa = "default", alpha_f = 0.05, 
                    min_size = 30, sd_beta = 3, thinning = 1, 
                    num_iter = 5000, num_warmup = 2500, print_outputs = FALSE, 
                    print_progress = FALSE, model_selection = FALSE,
                    prior_for_U = c("Spike_Slab", "Gaussian"), sparse_networks = FALSE, 
                    center_std_A = TRUE, seed  = NULL,
                    lambda_init = c("equal", "prior", "zero")){
  
  # Validate scalar arguments before using them or setting the seed.
  check_number <- function(x, name) {
    if (!is.numeric(x) || is.complex(x) || length(x) != 1L ||
        !is.finite(x)) {
      stop(name, " must be one finite real number.", call. = FALSE)
    }
  }
  
  check_integer <- function(x, name, lower = 1) {
    check_number(x, name)
    if (x != floor(x) || x < lower || x > .Machine$integer.max) {
      stop(
        name, " must be an integer between ", lower, " and ",
        .Machine$integer.max, ".", call. = FALSE
      )
    }
  }
  
  check_flag <- function(x, name) {
    if (!is.logical(x) || length(x) != 1L || is.na(x)) {
      stop(name, " must be TRUE or FALSE.", call. = FALSE)
    }
  }
  
  check_integer(num_CP, "num_CP")
  check_integer(latent_dim, "latent_dim")
  check_integer(num_iter, "num_iter")
  check_integer(num_warmup, "num_warmup", lower = 0)
  check_integer(thinning, "thinning")
  check_integer(min_size, "min_size", lower = 0)
  
  if (num_warmup >= num_iter) {
    stop("num_warmup must be smaller than num_iter.", call. = FALSE)
  }
  
  check_number(b_tau, "b_tau")
  check_number(sd_beta, "sd_beta")
  check_number(alpha_f, "alpha_f")
  check_number(eta_dir, "eta_dir")
  
  if (b_tau <= 0 || sd_beta <= 0) {
    stop("b_tau and sd_beta must be positive.", call. = FALSE)
  }
  if (alpha_f <= 0 || alpha_f >= 1) {
    stop("alpha_f must be strictly between 0 and 1.", call. = FALSE)
  }
  
  check_flag(print_outputs, "print_outputs")
  check_flag(print_progress, "print_progress")
  check_flag(model_selection, "model_selection")
  check_flag(sparse_networks, "sparse_networks")
  check_flag(center_std_A, "center_std_A")
  
  if (!is.null(seed)) {
    check_integer(seed, "seed", lower = -.Machine$integer.max)
  }
  
  valid_default_prior <- is.character(prior_kappa) &&
    length(prior_kappa) == 1L && !is.na(prior_kappa) &&
    prior_kappa == "default"
  if (!valid_default_prior && !is.function(prior_kappa)) {
    stop("prior_kappa must be 'default' or a function.", call. = FALSE)
  }
  
  prior_for_U <- match.arg(prior_for_U)
  lambda_init <- match.arg(lambda_init)
  
  # Record which inputs were supplied before creating placeholders.
  has_A <- !is.null(A)
  has_S <- !is.null(S)
  if (!has_A && !has_S) {
    stop("Supply A, S, or both.", call. = FALSE)
  }
  
  # S: one row per observation and one column per network statistic.
  if (has_S) {
    if (is.data.frame(S)) {
      numeric_columns <- vapply(
        S, function(x) is.numeric(x) && !is.complex(x), logical(1)
      )
      if (!all(numeric_columns)) {
        stop("Every column of S must be numeric.", call. = FALSE)
      }
      S <- as.matrix(S)
    }
    if (!is.matrix(S) || !is.numeric(S) || is.complex(S) ||
        nrow(S) < 2L || ncol(S) < 1L || any(!is.finite(S))) {
      stop(
        "S must be a finite numeric matrix with at least two rows and one column.",
        call. = FALSE
      )
    }
  }
  
  # A: equally sized, symmetric adjacency matrices with zero diagonals.
  if (has_A) {
    if (!is.list(A) || is.data.frame(A) || length(A) < 2L) {
      stop("A must be a list of at least two adjacency matrices.", call. = FALSE)
    }
    first_dims <- dim(A[[1L]])
    for (tt in seq_along(A)) {
      A_t <- A[[tt]]
      if (!(is.matrix(A_t) || methods::is(A_t, "Matrix"))) {
        stop("A[[", tt, "]] must be a matrix.", call. = FALSE)
      }
      dims <- dim(A_t)
      if (length(dims) != 2L || any(dims < 1L) || dims[1L] != dims[2L]) {
        stop("A[[", tt, "]] must be a nonempty square matrix.", call. = FALSE)
      }
      if (!identical(dims, first_dims)) {
        stop("All adjacency matrices must have the same dimensions.", call. = FALSE)
      }
      A_check <- as.matrix(A_t)
      if (!(is.numeric(A_check) || is.logical(A_check)) ||
          is.complex(A_check) || any(!is.finite(A_check))) {
        stop("A[[", tt, "]] must contain finite real values.", call. = FALSE)
      }
      if (any(abs(diag(A_check)) >= 1e-8)) {
        stop("A[[", tt, "]] must have a zero diagonal.", call. = FALSE)
      }
      if (!isTRUE(all.equal(
        A_check, t(A_check), tolerance = 1e-8, check.attributes = FALSE
      ))) {
        stop("A[[", tt, "]] must be symmetric.", call. = FALSE)
      }
    }
  }
  
  if (has_A && has_S && length(A) != nrow(S)) {
    stop("length(A) must equal nrow(S).", call. = FALSE)
  }
  Tt_all <- if (has_A) length(A) else nrow(S)
  
  # Preserve the original placeholder convention for an absent input.
  if (!has_A) {
    A <- lapply(seq_len(Tt_all), function(tt) matrix(0, 1L, 1L))
  }
  if (!has_S) {
    S <- matrix(0, nrow = Tt_all, ncol = 1L)
  }
  n <- nrow(A[[1L]])
  p <- ncol(S)
  D <- latent_dim
  
  # Each changepoint is the last observation in the preceding segment.
  # Internal kappa indexes the training series; original_kappa indexes
  # the original input series.
  if (!is.null(seed)) set.seed(seed)
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
  
  # Fit on all observations when thinning = 1; no data are held out.
  if (thinning ==1){ 
    Tt = Tt_all
    all_idx = seq_len(Tt) # all indices from the input data
    train_idx = all_idx
    holdout_idx = all_idx
    
    # Standardize and center the summary statistics S
    if (dim(S)[2] == 1){
      m = mean(S, na.rm = TRUE)
    }else{
      m = colMeans(S, na.rm = TRUE) # column means in the training set
    }
    sdv = apply(S, 2, FUN = sd)
    sdv[sdv == 0] = 1 # Avoid divide-by-zero in standardization: sdv = 0 could occur if the input data has a column of the same constant. 
    S = center_standardize(S, m, sdv)
    S_auc = S
    
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
  # With thinning > 1, hold out every thinning-th observation from fitting.
  else if (thinning > 1){ 
    # Preserve the existing conversion of min_size to the training scale.
    min_size <- round(min_size * (1 - 1 / thinning))
    Tt <- nrow(S)
    all_idx <- seq_len(Tt)
    holdout_idx <- all_idx[all_idx %% thinning == 0L]
    train_idx <- all_idx[all_idx %% thinning != 0L]
    
    if (length(holdout_idx) == 0L) {
      stop("thinning leaves no held-out observations.", call. = FALSE)
    }
    if (length(train_idx) < 2L) {
      # 2 is the minimum number of training observations to calculate sample standard deviations
      stop("At least two training observations are required.", call. = FALSE)
    }
    
    S_all <- S
    S_auc <- S_all[holdout_idx, , drop = FALSE]
    S <- S_all[train_idx, , drop = FALSE]
    
    A_all_raw <- if (sparse_networks) to_sparse_list(A) else A
    A_auc_raw <- A_all_raw[holdout_idx]
    A_train_raw <- A_all_raw[train_idx]
    
    # Standardize and center the summary statistics S
    if (ncol(S) == 1){
      m = mean(S, na.rm = TRUE)
    }else{
      m = colMeans(S, na.rm = TRUE) # column means in the training set
    }
    sdv = apply(S, 2, FUN = sd, na.rm = TRUE) # standard deviation of the training set
    sdv[sdv == 0] = 1 # Avoid divide-by-zero in standardization: sdv = 0 could occur if the input data has a column of the same constant. 
    S = center_standardize(S, m, sdv) # center & standardize the subseries used to fit the model with its column means
    S_auc = center_standardize(S_auc, m, sdv) # center & standardize the held-out samples with the same column means
    
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
    }else{
      A = A_train_raw
      A_auc = A_auc_raw
    }
  }
  
  Tt_auc = dim(S_auc)[1]
  Tt = dim(S)[1] # number of observations
  K = num_CP + 1 # number of classes
  L = num_CP # number of changepoints
  
  # Each segment must contain at least min_size + 1 training observations.
  if (Tt < K * (min_size + 1)) {
    stop(
      "Too few training observations for the requested number and lengths of segments.",
      call. = FALSE
    )
  }
  if (has_A){
    A_cpp = pack_A_cube(A)
    if (model_selection){
      A_auc_cpp = pack_A_cube(A_auc)
    }
  }
  
  #===================#
  # Utility functions ####
  #===================#
  # Internal linear predictors must be finite.
  logsumexp <- function(x) {
    if (any(!is.finite(x))) {
      stop("Non-finite linear predictors.", call. = FALSE)
    }
    x_max <- max(x)
    x_max + log(sum(exp(x - x_max)))
  }

  # -Inf denotes zero weight; invalid conditionals must not become uniform.
  normalize_logweights <- function(logweights) {
    if (anyNA(logweights) || any(logweights == Inf) ||
        !any(is.finite(logweights))) {
      stop("Invalid conditional log-weights; cannot sample.", call. = FALSE)
    }
    weights <- exp(logweights - max(logweights))
    weights / sum(weights)
  }
  # Robust AUC CI helper for model selection.
  # Unavailable intervals are stored as NA and must be reported by the wrapper.
  # Constant finite scores with both classes are recorded as non-informative AUC = 0.5.
  # Larger scores favor Y_auc = 1 (the later segment); allow AUC < 0.5.
  safe_ci_auc = function(Y_auc, scores, alpha_f){
    if (any(!is.finite(scores))) {
      stop("Non-finite AUC scores.", call. = FALSE)
    }
    Y_auc = as.numeric(Y_auc)
    scores = as.numeric(scores)
    if (length(Y_auc) < 2 || length(unique(Y_auc)) < 2){
      return(c(NA_real_, NA_real_, NA_real_)) 
    }
    if (length(unique(scores)) < 2){
      return(c(0.5, 0.5, 0.5)) 
    }
    ci_tmp = tryCatch(
      as.numeric(suppressWarnings(pROC::ci.auc(Y_auc, scores, levels = c(0, 1), direction = "<", conf.level = 1-alpha_f, quiet = TRUE))),
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
  # Evaluate one custom-prior value; zero excludes a candidate.
  checked_prior_value <- function(kappa_candidate) {
    value <- prior_kappa(kappa_candidate)
    if (!is.numeric(value) || is.complex(value) || length(value) != 1L ||
        !is.finite(value) || value < 0) {
      stop(
        "prior_kappa must return one finite, nonnegative real number.",
        call. = FALSE
      )
    }
    return(value)
  }
  
  
  #=======================#
  # Create an output list ####
  #=======================#
  out = list()
  out$tau_sq = rep(0, num_iter-num_warmup)
  out$Lambda = array(0, dim = c(num_iter-num_warmup, K, D))
  out$U = array(0, dim = c(num_iter-num_warmup, K, n, D))
  out$Z = array(0, dim = c(num_iter-num_warmup, n, K))
  out$Beta = array(0, dim = c(num_iter-num_warmup, p, K))
  out$W = array(0, dim = c(num_iter-num_warmup, n, K))
  out$pi = array(0, dim = c(num_iter-num_warmup, D, 3)) # store dimension-specific pi_d = (Pr(0), Pr(1), Pr(-1)) for each latent dimension d.
  dimnames(out$pi) = list(NULL, paste0("dim_", 1:D), c("zero", "pos1", "neg1")) # make the lambda/probability ordering explicit.
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
  V0inv = array(data=0, dim=c(p,p,K-1)) # Prior precision matrices
  sd0 = sd_beta # prior variance for beta (coefficients of summary statistics)
  for(k in 1:(K-1)){
    V0inv[,,k] = 1/sd0*diag(p)
  }
  # priors for tau_sq
  a_tau = 4
  # priors for w_ik ~ Beta(a_w, b_w)
  a_w = 1
  b_w = 1
  # prior for pi (Dirichlet)
  lambda_support = c(0, 1, -1) # Matches the columns of pi: zero, pos1, neg1.
  #===============================#
  # Initialize working quantities ####
  #===============================#
  # Explicit dependency namespace for installed-package use.
  tau_sq = LearnBayes::rigamma(n = 1, a = a_tau, b = b_tau) # inverse gamma variable
  Z = matrix(0, nrow = n, ncol = K) # Bernoulli variables
  Z[,K] = rep(0, n)
  U = array(0, dim = c(K, n, D)) # latent variables
  if (prior_for_U == "Gaussian"){
    for (i in 1:n){
      for (k in 1:(K-1)){
        Z[i,k] = 1
        # Explicit dependency namespace for installed-package use.
        U[k,i,] = MASS::mvrnorm(n =1, mu = rep(0,D), Sigma = tau_sq*diag(D))
      }
    }
    W = matrix(0, nrow = n, ncol = K) # Beta variables are unused under the Gaussian prior for U
  }else if (prior_for_U == "Spike_Slab"){
    Z = matrix(rbinom(n = n*K, size = 1, prob = 0.5), nrow = n, ncol = K) # Bernoulli variables
    Z[,K] = rep(0, n)
    for (i in 1:n){
      for (k in 1:(K-1)){
        if (Z[i,k] ==1){
          # Explicit dependency namespace for installed-package use.
          U[k,i,] = MASS::mvrnorm(n =1, mu = rep(0,D), Sigma = tau_sq*diag(D))
        }else if(Z[i,k] ==0){
          U[k,i,] = rep(0,D)
        }
      }
    }
    W = matrix(rbeta(n*K, shape1 = a_w, shape2 = b_w), nrow = n, ncol = K) # Beta variables
  }
  # W is unused in the reference slot; zero is a placeholder, not an estimate.
  W[,K] = 0
  Beta = matrix(0, nrow = p, ncol = K) # coefficients for summary statistics
  pi = matrix(0, nrow = D, ncol = 3) # pi[d,] is dimension-specific, not shared across dimensions.
  colnames(pi) = c("zero", "pos1", "neg1") # explicit ordering matches lambda_support = c(0, 1, -1).
  for (d in 1:D){
    pi[d,] = as.numeric(DirichletReg::rdirichlet(n = 1, alpha = c(d^eta_dir, 1, 1))) # Dirichlet concentrations for lambda = 0, 1, and -1.
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
  Lambda[K,] = rep(0, D) # Internal reference slot K is fixed at zero.
  omega = matrix(0, nrow = Tt, ncol = K)
  # Initialize CPs at time points that equally split the series 
  if(num_CP==1){
    kappa = round(Tt/2)
  }else if (num_CP>1){
    kappa = rep(0, L)
    for (l in 1:L){
      kappa[l] = round(Tt/K) *l
    }
  }
  # Keep the initializer unchanged, but reject an inadmissible state.
  segment_lengths <- diff(c(0, kappa, Tt))
  if (length(kappa) != L || any(!is.finite(kappa)) ||
      any(kappa != floor(kappa)) || any(segment_lengths <= min_size)) {
    stop(
      paste0(
        "An initialized segment has too few training observations. ",
        "Review the initialization or fitting settings before sampling."
      ),
      call. = FALSE
    )
  }
  if (is.function(prior_kappa) && checked_prior_value(kappa) == 0) {
    stop("The initialized changepoints have zero mass under prior_kappa.",
         call. = FALSE)
  }
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
    on.exit(close(pb), add = TRUE)
  }
  for (iter in 1: num_iter){
    if (num_CP == 1){
      # The first segment (Y = 0) is the binary reference class.
      # Internally, slot 1 holds fitted Y = 1 coefficients; slot 2 is zero.
      # Segment-indexed samples are reordered chronologically before return.
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
      # Explicit dependency namespace for installed-package use.
      omega[,k] = BayesLogit::rpg(num = Tt, h = 1, z = eta_k) # update omega using BayesLogit::rpg
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
        # Update the D diagonal weights for class k ####
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
      # Internal reference slot K is fixed at zero and excluded from pi updates.
      for (d in 1:D){
        n_0 = sum(Lambda[1:(K-1), d] == 0)
        n_pos1 = sum(Lambda[1:(K-1), d] == 1)
        n_neg1 = sum(Lambda[1:(K-1), d] == -1)
        alpha_pi = c(d^eta_dir + n_0, 1 + n_pos1, 1 + n_neg1) # Prior concentrations plus frequencies among nonreference slots.
        pi[d,] = as.numeric(DirichletReg::rdirichlet(n = 1, alpha = alpha_pi)) # update pi_d separately for each latent dimension.
      }
      #===========================================# 
      # Update tau_sq (shared across all classes) ####
      #===========================================#
      if (prior_for_U == "Spike_Slab"){
        sum_m_k = sum(Z[,1:(K-1), drop = FALSE])
      }else if (prior_for_U == "Gaussian"){
        sum_m_k = n * (K-1)
      }
      sum_ss_k = sum(U[1:(K-1), , , drop = FALSE]^2)
      # Explicit dependency namespace for installed-package use.
      tau_sq = LearnBayes::rigamma(n = 1, a = a_tau + (D/2) *sum_m_k, b = b_tau + 0.5*sum_ss_k)
    }else{
      # If A is NULL, keep pi and tau_sq fixed at zero.
      pi = matrix(0, nrow = D, ncol = 3) # keep the same dimension-specific pi shape when A is absent.
      colnames(pi) = c("zero", "pos1", "neg1") # explicit ordering is retained for downstream output.
      tau_sq = 0
    }
    # Done updating all other parameters except the changepoints
    
    #=============================#
    # Update kappa (changepoints) ####
    #=============================#
    # Use the updated network and summary-statistic contributions.
    phi = Q_iter + StB
    if (any(!is.finite(phi))) {
      stop("Non-finite linear predictors in the changepoint update.",
           call. = FALSE)
    } 
    if (num_CP == 1){
      log_p0 = stats::plogis(phi[,1], lower.tail = FALSE, log.p = TRUE)
      log_p1 = stats::plogis(phi[,1], log.p = TRUE)
      range = (1 + min_size):(Tt - 1 - min_size) 
      lr = length(range)
      log_prob_vec_prop = rep(0, lr) 
      if (is.character(prior_kappa) && prior_kappa == "default"){
        for(l in 1:lr){
          log_prob_vec_prop[l] = sum(log_p0[1:range[l]]) + sum(log_p1[(range[l]+1):Tt]) -
            length(1:range[l])*log(length(1:range[l])) - 
            length((range[l]+1):Tt)*log(length((range[l]+1):Tt)) 
        }
      }else if (is.function(prior_kappa)){
        for(l in 1:lr){
          jt_pmf_kappa <- checked_prior_value(range[l])# evaluate custom prior at candidate kappa
          log_prob_vec_prop[l] = sum(log_p0[1:range[l]]) + sum(log_p1[(range[l]+1):Tt]) + log(jt_pmf_kappa)
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
        # Use log-odds as AUC scores to avoid probability-scale ties.
        ci = safe_ci_auc(Y_auc, phi_auc[,1], alpha_f)
      }
    }else if(num_CP > 1){
      log_P = t(apply(phi, 1, function(x) {
        shifted = x - max(x)
        shifted - log(sum(exp(shifted)))
      }))
      if (is.character(prior_kappa) && prior_kappa == "default"){
        for(j in 1:(K-1)){
          # Evaluate and normalize the conditional log-weights for kappa[j].
          if(j==1){
            range = (1 + min_size):(kappa[j+1] - 1 - min_size) 
            lr = length(range)
            log_prob_vec_prop = rep(0, lr)
            for(l in 1:lr){
              log_prob_vec_prop[l] = sum(log_P[1:range[l],j]) + sum(log_P[(range[l]+1):kappa[j+1],j+1]) -
                length(1:range[l])*log(length(1:range[l])) - 
                length((range[l]+1):kappa[j+1])*log(length((range[l]+1):kappa[j+1])) 
            }
            prob_vec = normalize_logweights(log_prob_vec_prop)
          } else if(j==L){
            range = (kappa[j-1] + 1 + min_size):(Tt - 1 - min_size) 
            lr = length(range)
            log_prob_vec_prop = rep(0, lr)
            for(l in 1:lr){
              log_prob_vec_prop[l] = sum(log_P[(kappa[j-1]+1):range[l],j]) + sum(log_P[(range[l]+1):Tt,j+1]) -
                length((kappa[j-1]+1):range[l])*log(length((kappa[j-1]+1):range[l])) - 
                length((range[l]+1):Tt)*log(length((range[l]+1):Tt)) 
            }
            prob_vec = normalize_logweights(log_prob_vec_prop)
          } else{
            range = (kappa[j-1] + 1 + min_size):(kappa[j+1] - 1 - min_size) 
            lr = length(range)
            log_prob_vec_prop = rep(0, lr)
            for(l in 1:lr){
              log_prob_vec_prop[l] = sum(log_P[(kappa[j-1]+1):range[l],j]) + sum(log_P[(range[l]+1):kappa[j+1],j+1]) -
                length((kappa[j-1]+1):range[l])*log(length((kappa[j-1]+1):range[l])) - 
                length((range[l]+1):kappa[j+1])*log(length((range[l]+1):kappa[j+1])) 
            }
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
              jt_pmf_kappa = checked_prior_value(kappa_prop) 
              log_prob_vec_prop[l] = sum(log_P[1:range[l],j]) + sum(log_P[(range[l]+1):kappa[j+1],j+1]) + log(jt_pmf_kappa)
            }
            prob_vec = normalize_logweights(log_prob_vec_prop)
          } else if(j==L){
            range = (kappa[j-1] + 1 + min_size):(Tt - 1 - min_size) 
            lr = length(range)
            log_prob_vec_prop = rep(0, lr)
            for(l in 1:lr){
              kappa_prop = kappa # evaluate custom prior at candidate kappa vector
              kappa_prop[j] = range[l] 
              jt_pmf_kappa = checked_prior_value(kappa_prop) 
              log_prob_vec_prop[l] = sum(log_P[(kappa[j-1]+1):range[l],j]) + sum(log_P[(range[l]+1):Tt,j+1]) + log(jt_pmf_kappa)
            }
            prob_vec = normalize_logweights(log_prob_vec_prop)
          } else{
            range = (kappa[j-1] + 1 + min_size):(kappa[j+1] - 1 - min_size) 
            lr = length(range)
            log_prob_vec_prop = rep(0, lr)
            for(l in 1:lr){
              kappa_prop = kappa 
              kappa_prop[j] = range[l] 
              jt_pmf_kappa = checked_prior_value(kappa_prop) 
              log_prob_vec_prop[l] = sum(log_P[(kappa[j-1]+1):range[l],j]) + sum(log_P[(range[l]+1):kappa[j+1],j+1]) + log(jt_pmf_kappa)
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
      # Which fitted changepoints appear to correspond to true changepoints?  #
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
        # Pairwise log-odds have the same ordering as conditional
        # probabilities, without saturation at zero or one.
        Y_auc_all = kappa_to_Y(original_kappa, Tt_all, L)[holdout_idx, , drop = FALSE] # Y in Tt_auc *K
        ci = matrix(0, nrow = K-1, ncol = 3) # each row has a LB, point estimate, and UB
        for (j in 1:(K-1)){
          if(j == 1){
            rg = 1:original_kappa[j+1]
            idx_test = which(holdout_idx %in% rg)
            scores = phi_auc[idx_test, j+1] - phi_auc[idx_test, j]
            Y_auc = Y_auc_all[idx_test, j + 1]
            ci[j,] = safe_ci_auc(Y_auc, scores, alpha_f)
          }else if (j == K-1){
            rg = (original_kappa[j-1]+1):Tt_all
            idx_test = which(holdout_idx %in% rg)
            scores = phi_auc[idx_test, j+1] - phi_auc[idx_test, j]
            Y_auc = Y_auc_all[idx_test, j + 1]
            ci[j,] = safe_ci_auc(Y_auc, scores, alpha_f)
          }else{
            rg = (original_kappa[j-1]+1):original_kappa[j+1]
            idx_test = which(holdout_idx %in% rg)
            scores = phi_auc[idx_test, j+1] - phi_auc[idx_test, j]
            Y_auc = Y_auc_all[idx_test, j + 1]
            ci[j,] = safe_ci_auc(Y_auc, scores, alpha_f)
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
      out$tau_sq[iter-num_warmup] = tau_sq
      out$Lambda[iter-num_warmup,,] = Lambda
      out$U[iter-num_warmup,,,] = U
      out$Z[iter-num_warmup,,] = Z
      out$Beta[iter-num_warmup,,] = Beta
      out$W[iter-num_warmup,,] = W
      out$pi[iter-num_warmup, , ] = pi # store the D x 3 matrix of dimension-specific pi_d values.
      
      # Store changepoints and AUC confidence-interval lower bounds:
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
  #================================#
  # END Gibbs sampling (iter)   ####
  #================================#
  
  # Point estimates (marginal posterior modes) of changepoints:
  if (num_CP == 1){
    cp = table(out$Kappa)
    out$Kappa_mode = as.integer(names(cp)[which.max(cp)]) # posterior mode
  }else if (num_CP > 1){
    cps = c()
    for (l in 1:L){
      cp = table(out$Kappa[,l])
      cp_mode = as.integer(names(cp)[which.max(cp)])
      cps = c(cps, cp_mode)
    }
    out$Kappa_mode = cps # posterior modes
  }
  # Return segment-indexed samples in chronological order.
  # Reference: first segment (Y = 0) for one changepoint;
  # last chronological segment for multiple changepoints.
  if (num_CP == 1) {
    segment_order <- c(2L, 1L)
    out$Beta <- out$Beta[, , segment_order, drop = FALSE]
    out$U <- out$U[, segment_order, , , drop = FALSE]
    out$Lambda <- out$Lambda[, segment_order, , drop = FALSE]
    out$Z <- out$Z[, , segment_order, drop = FALSE]
    out$W <- out$W[, , segment_order, drop = FALSE]
  }
  
  if (print_outputs) {
    utils::str(out)
  }
  return(out)
}
