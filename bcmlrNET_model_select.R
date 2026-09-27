# Select the number of changepoints using AUC confidence-interval lower bounds.
# Refit on the full series when the selected number is positive.
#
# Selection outputs:
#   num_cp       : selected number (posterior mode) of changepoints.
#   max_cps      : the number of changepoint point estimates from the initial fit.
#   num_cp_dist  : posterior samples used to select the number of changepoints.
#   probs_num_cp : probabilities for 0, ..., max_num_cp changepoints.
#   auc_diagnostics : first-stage AUC availability summaries (see below).
# A positive num_cp also returns the full-series fit from bcmlrNET().
# With num_cp = 0, only the selection outputs above are returned.
#
# warn_small_auc controls the small-group warning only. Diagnostics are always
# returned; other interval failures still warn, and all-unavailable AUCs stop.
bcmlrNET_model_select <- function(
    A = NULL, S = NULL, latent_dim = 10, eta_dir = 1.5, b_tau = 1,
    thinning = 5, prior_kappa = "default", alpha_f = 0.1,
    threshold = 0.5, min_size = 30, max_num_cp = 5, sd_beta = 3,
    num_iter = 5000, num_warmup = 2500,
    print_progress = FALSE, print_outputs = FALSE,
    prior_for_U = c("Spike_Slab", "Gaussian"),
    lambda_init = c("equal", "prior", "zero"),
    sparse_networks = FALSE, center_std_A = TRUE,
    seed = NULL, warn_small_auc = TRUE) {

  #-----------------#
  # Validate inputs #
  #-----------------#
  check_integer <- function(x, name, lower = 1) {
    if (!is.numeric(x) || is.complex(x) || length(x) != 1L ||
        !is.finite(x) || x != floor(x) || x < lower ||
        x > .Machine$integer.max) {
      stop(name, " must be a valid integer of at least ", lower, ".",
           call. = FALSE)
    }
  }
  check_flag <- function(x, name) {
    if (!is.logical(x) || length(x) != 1L || is.na(x)) {
      stop(name, " must be TRUE or FALSE.", call. = FALSE)
    }
  }

  check_integer(max_num_cp, "max_num_cp")
  check_integer(min_size, "min_size", lower = 0)
  check_integer(num_iter, "num_iter")
  check_integer(num_warmup, "num_warmup", lower = 0)
  if (num_warmup >= num_iter) {
    stop("num_warmup must be smaller than num_iter.", call. = FALSE)
  }
  if (!is.numeric(threshold) || is.complex(threshold) ||
      length(threshold) != 1L || !is.finite(threshold) ||
      threshold < 0 || threshold > 1) {
    stop("threshold must be one finite number in [0, 1].", call. = FALSE)
  }
  check_flag(print_progress, "print_progress")
  check_flag(print_outputs, "print_outputs")
  check_flag(warn_small_auc, "warn_small_auc")
  prior_for_U <- match.arg(prior_for_U)
  lambda_init <- match.arg(lambda_init)

  if (is.null(A) && is.null(S)) {
    stop("Supply A, S, or both.", call. = FALSE)
  }
  if (!is.null(A) &&
      (!is.list(A) || is.data.frame(A) || length(A) < 2L)) {
    stop("A must be a list of at least two adjacency matrices.", call. = FALSE)
  }
  if (!is.null(S) &&
      (!(is.matrix(S) || is.data.frame(S)) ||
       nrow(S) < 2L || ncol(S) < 1L)) {
    stop("S must have at least two rows and one column.", call. = FALSE)
  }
  if (!is.null(A) && !is.null(S) && length(A) != nrow(S)) {
    stop("length(A) must equal nrow(S).", call. = FALSE)
  }

  # Preliminary feasibility check; bcmlrNET() also checks the training series.
  Tt <- if (!is.null(A)) length(A) else nrow(S)
  if (min_size > 0 && max_num_cp > Tt / min_size - 1) {
    stop("max_num_cp exceeds the preliminary series-length limit.",
         call. = FALSE)
  }

  check_integer(thinning, "thinning")

  #--------------------------------#
  # Stage 1: Fit the initial model #
  #--------------------------------#
  if (print_progress || print_outputs) {
    message(
      "Stage 1: Initial fit (number of changepoints = ", max_num_cp, ").\n",
      if (thinning > 1L) {
        paste0(
          "Holding out one observation every ", thinning,
          " observations for AUC evaluation; fitting to the remaining observations."
        )
      } else {
        "Using the entire series (thinning = 1); no observations are held out."
      }
    )
  }

  old_out <- bcmlrNET(
    A = A, S = S, num_CP = max_num_cp,
    latent_dim = latent_dim, eta_dir = eta_dir, b_tau = b_tau,
    prior_kappa = prior_kappa, alpha_f = alpha_f,
    min_size = min_size, sd_beta = sd_beta, lambda_init = lambda_init,
    thinning = thinning, num_iter = num_iter, num_warmup = num_warmup,
    print_outputs = FALSE, print_progress = print_progress,
    model_selection = TRUE, prior_for_U = prior_for_U,
    sparse_networks = sparse_networks, center_std_A = center_std_A,
    seed = seed
  )

  #---------------------------------#
  # Validate AUC interval summaries #
  #---------------------------------#
  if (!is.list(old_out) || is.null(old_out$AUC)) {
    stop("The initial fit did not return AUC confidence-interval lower bounds.",
         call. = FALSE)
  }
  auc_lb <- as.matrix(old_out$AUC)
  if (!is.numeric(auc_lb) || is.complex(auc_lb) ||
      nrow(auc_lb) != num_iter - num_warmup ||
      ncol(auc_lb) != max_num_cp) {
    stop("Unexpected AUC output from the initial bcmlrNET fit.", call. = FALSE)
  }

  valid_auc <- is.finite(auc_lb)
  if (!any(valid_auc)) {
    stop("All AUC intervals are unavailable; selection cannot proceed.",
         call. = FALSE)
  }
  if (any(auc_lb[valid_auc] < 0 | auc_lb[valid_auc] > 1)) {
    stop("AUC confidence-interval lower bounds must lie in [0, 1].",
         call. = FALSE)
  }

  #-------------------------#
  # Compute AUC diagnostics #
  #-------------------------#
  # Reconstruct evaluation-set sizes from the first-stage changepoint samples.
  # Tt is the original series length; kappa samples use original indices.
  cuts <- cbind(0, as.matrix(old_out$Kappa), Tt)
  cumulative_eval <- floor(cuts / thinning)

  n_eval <- cumulative_eval[, -1, drop = FALSE] -
    cumulative_eval[, -ncol(cumulative_eval), drop = FALSE]

  # Each comparison uses two neighboring segments. Limited samples means
  # fewer than two evaluation observations in at least one segment.
  # With thinning = 1, all observations are used for evaluation.
  small_group <- n_eval[, -ncol(n_eval), drop = FALSE] < 2 |
    n_eval[, -1, drop = FALSE] < 2

  # Summarize availability once after Stage 1; do not recompute AUC intervals.
  # The two unavailable categories partition n_unavailable.
  # One row per first-stage changepoint, not per final selected changepoint.
  # All n_* columns count interval evaluations, not observations or segments.
  auc_diagnostics <- data.frame(
    fitted_changepoint = seq_len(max_num_cp),
    n_evaluations = nrow(auc_lb),
    n_unavailable = colSums(!valid_auc),
    pct_unavailable = 100 * colMeans(!valid_auc),
    n_unavailable_limited_samples = colSums(!valid_auc & small_group),
    n_unavailable_other = colSums(!valid_auc & !small_group)
  )

  if (sum(auc_diagnostics$n_unavailable_other) > 0L) {
    warning(
      "Some AUC intervals are unavailable despite at least two evaluation ",
      "observations in each neighboring segment. See auc_diagnostics.",
      call. = FALSE
    )
  } else if (any(!valid_auc) && warn_small_auc) {
    warning(
      "Some AUC intervals are unavailable because a neighboring segment has ",
      "fewer than two evaluation observations; consider increasing min_size. ",
      "See auc_diagnostics.",
      call. = FALSE
    )
  }

  #-----------------------------------#
  # Select the number of changepoints #
  #-----------------------------------#
  # A fitted changepoint is supported when its AUC lower bound exceeds threshold.
  # Unavailable intervals contribute no support; equality does not exceed the threshold.
  boolean <- auc_lb <= threshold
  boolean[!valid_auc] <- TRUE
  rej <- rowSums(!boolean)
  tb <- table(factor(rej, levels = seq_len(max_num_cp + 1) - 1)) / length(rej)
  num_cp <- as.integer(names(tb)[which.max(tb)])
  probs <- as.numeric(tb)
  names(probs) <- names(tb)

  #-----------------------------------#
  # Stage 2: Refit on the full series #
  #-----------------------------------#
  if (print_progress || print_outputs) {
    message(
      "Selected number of changepoints: ", num_cp, ".\n",
      if (num_cp > 0L) {
        "Stage 2: Refitting with this number on the entire series."
      } else {
        "No refitting is performed because zero changepoints were selected."
      }
    )
  }

  if (num_cp == 0L) {
    out <- list()
  } else {
    out <- bcmlrNET(
      A = A, S = S, num_CP = num_cp,
      latent_dim = latent_dim, eta_dir = eta_dir, b_tau = b_tau,
      prior_kappa = prior_kappa, alpha_f = alpha_f,
      min_size = min_size, sd_beta = sd_beta, lambda_init = lambda_init,
      thinning = 1, num_iter = num_iter, num_warmup = num_warmup,
      print_outputs = FALSE, print_progress = print_progress,
      model_selection = FALSE, prior_for_U = prior_for_U,
      sparse_networks = sparse_networks, center_std_A = center_std_A,
      seed = seed
    )
  }

  #---------------------#
  # Assemble the output #
  #---------------------#
  out$num_cp <- num_cp
  # Point estimates from the initial fit with max_num_cp changepoints.
  out$max_cps <- old_out$Kappa_mode
  # Posterior samples of the number of fitted changepoints that appear
  # to correspond to true changepoints under the selection procedure.
  out$num_cp_dist <- rej
  out$probs_num_cp <- probs
  out$auc_diagnostics <- auc_diagnostics
  if (print_outputs) utils::str(out)
  return(out)
}
