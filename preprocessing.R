find_outliers_sd <- function(X, threshold = 3, na.rm = TRUE) {
  if (!is.matrix(X) && !is.data.frame(X)) {
    stop("X must be a matrix or data.frame")
  }

  X <- as.matrix(X)

  # compute z-scores columnwise
  Z <- scale(X, center = TRUE, scale = TRUE)

  # identify outliers
  outlier_matrix <- abs(Z) > threshold

  # indices of outlying entries
  outlier_indices <- which(outlier_matrix, arr.ind = TRUE)

  # corresponding values
  outlier_values <- X[outlier_matrix]

  list(
    outlier_matrix = outlier_matrix,
    outlier_indices = outlier_indices,
    outlier_values = outlier_values,
    z_scores = Z
  )
}

select_proteins <- function(cor_mat, auc_values, r2_threshold = 0.64) {

  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("Package 'igraph' is required.")
  }

  if (!is.matrix(cor_mat) || nrow(cor_mat) != ncol(cor_mat)) {
    stop("cor_mat must be a square correlation matrix.")
  }

  p <- ncol(cor_mat)

  if (length(auc_values) != p) {
    stop("Length of auc_values must match cor_mat dimensions.")
  }

  if (r2_threshold < 0 || r2_threshold > 1) {
    stop("r2_threshold must be between 0 and 1.")
  }

  # --- Construct adjacency matrix ---
  r2_mat <- cor_mat^2

  adj_mat <- r2_mat >= r2_threshold
  diag(adj_mat) <- FALSE

  # --- Build graph ---
  g <- igraph::graph_from_adjacency_matrix(
    adj_mat,
    mode = "undirected",
    diag = FALSE
  )

  comps <- igraph::components(g)
  membership <- comps$membership

  # --- Select highest AUC per component ---
  keep <- integer(0)
  representatives <- list()

  for (cid in unique(membership)) {
    members <- which(membership == cid)

    if (length(members) == 1) {
      best <- members
    } else {
      best <- members[which.max(auc_values[members])]
    }

    keep <- c(keep, best)
    representatives[[as.character(cid)]] <- best
  }

  list(
    keep_indices = sort(keep),
    components = membership,
    representatives = representatives,
    graph = g
  )
}
