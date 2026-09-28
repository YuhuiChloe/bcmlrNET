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
// Cached class update for the R/C++ implementation.
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
