# --------------------------------------------------
# Simulate a series of adjacency matrices from a
# piecewise SBM with changepoints at 100 and 500,
# and compute edges / triangles / isolates using ergm
# --------------------------------------------------
set.seed(123)
library(ergm)
library(network)

T_len <- 600          # number of time points
n <- 27            # number of nodes
K <- 3                # number of blocks
cp <- c(100, 500)     # changepoints

# Fixed block memberships
# Example: 15 nodes in block 1, 15 in block 2, 10 in block 3
z <- c(rep(1, 9), rep(2, 9), rep(3, 9))

B1 <- matrix(c(
  0.25, 0.05, 0.03,
  0.05, 0.20, 0.04,
  0.03, 0.04, 0.15
), nrow = K, byrow = TRUE)

B2 <- matrix(c(
  0.25, 0.05, 0.03,
  0.10, 0.40, 0.08,
  0.03, 0.04, 0.15
), nrow = K, byrow = TRUE)

B3 <- matrix(c(
  0.30, 0.02, 0.01,
  0.02, 0.12, 0.06,
  0.01, 0.06, 0.25
), nrow = K, byrow = TRUE)

# ---------------------------
# 2. Function to simulate one SBM adjacency matrix
# ---------------------------
simulate_sbm <- function(n, z, B) {
  A <- matrix(0, n, n)
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      p <- B[z[i], z[j]]
      A[i, j] <- rbinom(1, 1, p)
      A[j, i] <- A[i, j]
    }
  }
  diag(A) <- 0
  return(A)
}
# ---------------------------
# 3. Function to compute summary statistics using ergm
# ---------------------------
compute_ergm_stats <- function(A) {
  # Convert adjacency matrix to network object
  g <- network(A, directed = FALSE, matrix.type = "adjacency")
  # Compute ERGM statistics
  s <- summary(g ~ edges + triangle + isolates)
  # Return as named numeric vector
  c(
    edges = unname(s["edges"]),
    triangles = unname(s["triangle"]),
    isolates = unname(s["isolates"])
  )
}
# ---------------------------
# 4. Simulate the time series
# ---------------------------
A_list <- vector("list", T_len)
stats_mat <- matrix(NA, nrow = T_len, ncol = 3)
colnames(stats_mat) <- c("edges", "triangles", "isolates")

for (t in 1:T_len) {
  # Choose the SBM matrix according to segment
  if (1<= t && t<= cp[1]) {
    B <- B1
  } else if (cp[1] < t  && t<= cp[2]) {
    B <- B2
  } else {
    B <- B3
  }
  # Simulate adjacency matrix
  A_t <- simulate_sbm(n, z, B)
  A_list[[t]] <- A_t
  # Compute ergm-based summaries
  stats_mat[t, ] <- compute_ergm_stats(A_t)
}

head(stats_mat)
dim(stats_mat)
length(A_list)

# Make the list of adjacency matrices as a 3D array
A_array <- array(unlist(A_list), dim = c(length(A_list), n, n))
is.array(A_array)

#================#
# Plot the data  #
#================#
par(mfrow = c(3,1))
for (d in 1:3){
  plot(stats_mat[,d])
}

# A is an n x n adjacency matrix
# z is the true block membership vector
plot_adj <- function(A, z) {
  n <- nrow(A)
  ord <- order(z[1:n])
  A_ord <- A[ord, ord]
  
  image(
    1:n, 1:n,
    t(A_ord[n:1, ]),
    col = c("white", "black"),
    xlab = "Node",
    ylab = "Node",
    main = "Adjacency matrix"
  )
  
  block_sizes <- table(z[1:n][ord])
  cuts <- cumsum(block_sizes)
  
  for (b in cuts[-length(cuts)]) {
    abline(v = b + 0.5, col = "red", lwd = 2)
    abline(h = n - b + 0.5, col = "red", lwd = 2)
  }
}

plot_adj(A_list[[1]], z)
