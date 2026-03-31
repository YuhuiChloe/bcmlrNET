library(ergm)
library(network)

simulate_sbm_series_array = function(T_len, n, cp, z, B_list) {
  K = length(cp) + 1
  if (length(z) != n) stop("length(z) must equal n.")
  if (length(B_list) != length(cp) + 1) {
    stop("Need length(cp) + 1 matrices in B_list.")
  }
  if (any(sort(cp) != cp)) stop("cp must be sorted in increasing order.")
  if (any(cp < 1 | cp >= T_len)) stop("cp values must be between 1 and T_len - 1.")
  
  for (k in seq_along(B_list)) {
    if (!all(dim(B_list[[k]]) == c(K, K))) {
      stop(paste("B_list[[", k, "]] must be a ", K, " x ", K, " matrix.", sep = ""))
    }
  }
  simulate_sbm = function(n, z, B) {
    A = matrix(0, n, n)
    for (i in 1:(n - 1)) {
      for (j in (i + 1):n) {
        p = B[z[i], z[j]]
        A[i, j] = rbinom(1, 1, p)
        A[j, i] = A[i, j]
      }
    }
    diag(A) = 0
    return(A)
  }
  compute_ergm_stats = function(A) {
    g = network(A, directed = FALSE, matrix.type = "adjacency")
    s = summary(g ~ edges + triangle + isolates)
    c(
      edges = unname(s["edges"]),
      triangles = unname(s["triangle"]),
      isolates = unname(s["isolates"])
    )
  }
  A_array = array(0, dim = c(T_len, n, n))
  stats_mat = matrix(NA, nrow = T_len, ncol = 3)
  colnames(stats_mat) = c("edges", "triangles", "isolates")
  for (t in 1:T_len) {
    seg = findInterval(t, cp) + 1
    B = B_list[[seg]]
    
    A_t = simulate_sbm(n, z, B)
    A_array[t, , ] = A_t
    stats_mat[t, ] = compute_ergm_stats(A_t)
  }
  return(list(
    A_array = A_array,
    stats_mat = stats_mat
  ))
}

# =============== #
# Example code ####
# =============== #
T_len = 300
n = 15
K = 3
cp = c(50, 250)
z = c(rep(1, 5), rep(2, 5), rep(3, 5))
B1 = matrix(c(
  0.50, 0.05, 0.10,
  0.05, 0.50, 0.04,
  0.10, 0.04, 0.50
), nrow = K, byrow = TRUE)

B2 = matrix(c(
  0.40, 0.07, 0.03,
  0.07, 0.40, 0.00,
  0.03, 0.00, 0.40
), nrow = K, byrow = TRUE)

B3 = matrix(c(
  0.50, 0.05, 0.10,
  0.05, 0.50, 0.04,
  0.10, 0.04, 0.50
), nrow = K, byrow = TRUE)

out = simulate_sbm_series_array(T_len = T_len,n = n,cp = cp,z = z, B_list = list(B1, B2, B3))
A_array = out$A_array
stats_mat = out$stats_mat
dim(A_array)
dim(stats_mat)

#==========#
# Plots ####
#==========#
par(mfrow = c(3,1))
for (d in 1:3){
  plot(stats_mat[,d])
}

plot_adj = function(A, z) {
  n = nrow(A)
  ord = order(z)
  A_ord = A[ord, ord]
  image(
    1:n, 1:n,
    t(A_ord[n:1, ]),
    col = c("white", "black"),
    xlab = "Node",
    ylab = "Node",
    main = "Adjacency matrix"
  )
  block_sizes = table(z[ord])
  cuts = cumsum(block_sizes)
  for (b in cuts[-length(cuts)]) {
    abline(v = b + 0.5, col = "red", lwd = 2)
    abline(h = n - b + 0.5, col = "red", lwd = 2)
  }
}
plot_adj(as.matrix(A_array[1, ,]), z)
