#' Compute simple statistics for undirected binary networks
#' @param A_list Nonempty list of binary, symmetric, zero-diagonal adjacency
#'   matrices on the same set of at least two nodes.
#' @param network_stats Statistics to return, in order: `edges`, `density`,
#'   `triangles`, `isolates`. The alias `triangle` is accepted.
#' @returns A numeric matrix with one row per network and named columns.
#' @examples
#' A <- matrix(c(0, 1, 1, 1, 0, 1, 1, 1, 0), 3, 3)
#' compute_netstats(list(A), c("edges", "triangles", "isolates"))
#' @export
compute_netstats <- function(A_list, network_stats = c("density", "triangles")) {
  if (!is.list(A_list) || is.data.frame(A_list) || !length(A_list)) {
    stop("A_list must be a nonempty list of adjacency matrices.", call. = FALSE)
  }
  if (!is.character(network_stats) || !length(network_stats) || anyNA(network_stats)) {
    stop("network_stats must be a nonempty character vector.", call. = FALSE)
  }
  network_stats[network_stats == "triangle"] <- "triangles"
  allowed <- c("edges", "density", "triangles", "isolates")
  if (any(!network_stats %in% allowed) || anyDuplicated(network_stats)) {
    stop("Choose distinct statistics from edges, density, triangles, isolates.",
         call. = FALSE)
  }
  first_dims <- dim(A_list[[1L]])
  S <- matrix(0, nrow = length(A_list), ncol = length(network_stats),
              dimnames = list(NULL, network_stats))
  for (tt in seq_along(A_list)) {
    A <- A_list[[tt]]
    if (!(is.matrix(A) || methods::is(A, "Matrix")) ||
        length(dim(A)) != 2L || nrow(A) < 2L || nrow(A) != ncol(A) ||
        !identical(dim(A), first_dims)) {
      stop("All networks must be equally sized square matrices with at least two nodes.",
           call. = FALSE)
    }
    A <- as.matrix(A)
    if (!(is.numeric(A) || is.logical(A)) || is.complex(A) ||
        anyNA(A) || any(!A %in% c(0, 1)) || any(diag(A) != 0) ||
        !isTRUE(all.equal(A, t(A), check.attributes = FALSE))) {
      stop("Network statistics require binary, symmetric matrices with zero diagonals.",
           call. = FALSE)
    }
    n <- nrow(A)
    # The statistical definitions are unchanged from the research helper.
    if ("edges" %in% network_stats) S[tt, "edges"] <- sum(A) / 2
    if ("density" %in% network_stats) S[tt, "density"] <- sum(A) / (n * (n - 1))
    if ("triangles" %in% network_stats) {
      S[tt, "triangles"] <- sum(diag(A %*% A %*% A)) / 6
    }
    if ("isolates" %in% network_stats) S[tt, "isolates"] <- sum(rowSums(A) == 0)
  }
  return(S)
}
