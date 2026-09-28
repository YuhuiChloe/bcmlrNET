# Package registration and namespace declarations.
#' Generalized Bayesian Changepoint Analysis for Dynamic Networks
#'
#' Fixed-dimensional Gibbs fitting, selection of the number of changepoints,
#' posterior coefficient summaries, and network visualizations.
#' @useDynLib bcmlrNET, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom rlang .data
#' @importFrom stats model.matrix rbeta rbinom rnorm sd
#' @importFrom utils head tail setTxtProgressBar txtProgressBar
#' @keywords internal
"_PACKAGE"
