#' Perform cluster analysis with the K-means clustering algorithm.
#'
#' @param  object  A `Cellustering` instance.
#' @param  k  The number of clusters.
#' @param  dimensions  The number of principal components used for clustering.
#'
#' @return  The `Cellustering` instance with clustering results and
#' visualization added to the `clustering` slot.
#' @export
#'
#' @importFrom stats kmeans
#'
#' @examples
#' # Run the necessary to catch up the progress
#' data("pbmc_small")
#' pbmc_small_9 <- Cellustering(pbmc_small@data)
#' pbmc_small_9 <- qc_plot(pbmc_small_9)
#' pbmc_small_9 <- qc_filter(pbmc_small_9)
#' pbmc_small_9 <- normalize(pbmc_small_9, scale_factor = 1e6)
#' pbmc_small_9 <- find_HVG(pbmc_small_9, n_feature = 20)
#' pbmc_small_9 <- scale_data(pbmc_small_9)
#' pbmc_small_9 <- principal_component_analysis(pbmc_small_9)
#' pbmc_small_9 <- select_proper_dimension(pbmc_small_9,
#'   lower_bound = 1, upper_bound = 5
#' )
#'
#' pbmc_small_9 <- kmeans(pbmc_small_9, 2)
#'
kmeans <- function(object,
                   k,
                   dimensions = 10) {
  # Check if the object belongs to "Cellustering class"
  umap <- FALSE
  if (!is(object, "Cellustering")) {
    stop("Please input a Cellustering object.")
  }

  # Check if the user have run `principal_component_analysis`` first
  if (!object@progress["PCA"]) {
    stop("Please run principal_component_analysis() first.")
  }
  features <- object@reduced_dimension$scaled_data
  PC1 <- object@reduced_dimension$PC1
  PC2 <- object@reduced_dimension$PC2
  N <- nrow(features)

  # Check if other parameters are in correct data types
  proper_type <- all(
    is.numeric(dimensions), is.logical(umap)
  )
  proper_length <- all(
    length(dimensions) == 1, length(umap) == 1
  )
  proper_range <- all(
    dimensions > 0, identical(dimensions, floor(dimensions)),
    dimensions <= N
  )

  if (!(proper_type && proper_length && proper_range)) {
    stop("Invalid input. Please ensure your input is meaningful.")
  }

  # Project the data to selected dimensions.
  eigenvectors <- object@reduced_dimension$eigenvectors
  feature_matrix <- as.matrix(features)
  filtered_eigenvectors <- eigenvectors[, 1:dimensions]
  projected_data <- data.frame(t(feature_matrix) %*% filtered_eigenvectors)
  object@clustering$data <- projected_data
  fast_kmeans <- stats::kmeans(projected_data, k)
  cluster_assignment <- fast_kmeans$cluster
  centroids <- fast_kmeans$centers
  cluster_result <- list(clusters = cluster_assignment, centroids = centroids)

  #  Randomly choose k data points as initial centroids.
  # set.seed(0)
  # initial_index <- sample(1:N, k)
  # distance <- function(x, y) {
  #   sqrt(sum((x - y)^2))
  # }
  # initial_centroids <- projected_data[initial_index, ]
  # centroids <- initial_centroids
  # cluster_assignment <- rep(0, N)
  # max_iterations <- 500
  # distance_matrix <- matrix(0, nrow = N, ncol = k)
  # for (time in 1:max_iterations) {
  #   for (i in 1:N) {
  #     for (j in 1:k) {
  #       distance_matrix[i, j] <- distance(projected_data[i, ], centroids[j, ])
  #     }
  #   }
  #   # Assign points to closest centroid
  #   new_cluster_assignment <- apply(distance_matrix, 1, which.min)
  #   # If no change in assignments, break the loop
  #   if (all(cluster_assignment == new_cluster_assignment)) {
  #     break
  #   }

  #   # Update cluster assignments
  #   cluster_assignment <- new_cluster_assignment

  #   # Update centroids
  #   for (i in 1:k) {
  #     centroids[i, ] <- colMeans(projected_data[cluster_assignment == i, ])
  #   }
  # }

  # Visualize
  visualization_base <- object@reduced_dimension$PCA_plot
  if (umap) {
    cluster_plot <- visualization_base
  } else {
    k_means_plot <- visualization_base +
      # geom_point(
      #   mapping = aes(color = cluster_assignment),
      #   shape = 16, size = 1
      # ) +
      geom_point(
        color = cluster_assignment,
        shape = 16, size = 1
      ) +
      labs(
        title = "K Means Plot",
        x = paste("PC1 ( =", PC1, ")"),
        y = paste("PC2 ( =", PC2, ")")
      )
    cluster_plot <- k_means_plot
  }
  show(cluster_plot)
  object@clustering$clustering_result <- cluster_result
  object@clustering$clustering_plot <- cluster_plot
  object@progress["clustering"] <- TRUE
  invisible(object)
}

#' Calculate the average silhouette coefficient for the clustering results.
#'
#' @param object  A `Cellustering` instance.
#'
#' @return The `Cellustering` instance with the average silhouette coefficient
#' added to the `clustering` slot.
#' @export
#'
#' @importFrom stats dist
#'
#' @examples
#' # Run the necessary to catch up the progress
#' data("pbmc_small")
#' pbmc_small_0 <- Cellustering(pbmc_small@data)
#' pbmc_small_0 <- qc_plot(pbmc_small_0)
#' pbmc_small_0 <- qc_filter(pbmc_small_0)
#' pbmc_small_0 <- normalize(pbmc_small_0, scale_factor = 1e6)
#' pbmc_small_0 <- find_HVG(pbmc_small_0, n_feature = 20)
#' pbmc_small_0 <- scale_data(pbmc_small_0)
#' pbmc_small_0 <- principal_component_analysis(pbmc_small_0)
#' pbmc_small_0 <- select_proper_dimension(pbmc_small_0,
#'   lower_bound = 1, upper_bound = 5
#' )
#' pbmc_small_0 <- kmeans(pbmc_small_0, 2)
#'
#' pbmc_small_0 <- compute_silhouette(pbmc_small_0)
compute_silhouette <- function(object) {
  # Check if the object belongs to "Cellustering class"
  if (!is(object, "Cellustering")) {
    stop("Please input a Cellustering object.")
  }
  # Compute pairwise distances
  data <- object@clustering$data
  distances <- as.matrix(dist(data))
  clusters <- object@clustering$clustering_result$clusters

  # Number of data points
  N <- nrow(data)

  # Compute silhouette values for each data point
  silhouette_values <- vector(mode = "double", length = N)
  for (i in 1:N) {
    # Indices of data points in the same cluster
    same_cluster <- which(clusters == clusters[i])

    # Average distance to data points in the same cluster
    a_i <- mean(distances[i, same_cluster])

    # Indices of data points in other clusters
    other_clusters <- unique(clusters[-same_cluster])

    b_i <- min(sapply(other_clusters, function(j) {
      mean(distances[i, clusters == j])
    }))

    # Silhouette value
    silhouette_values[i] <- (b_i - a_i) / max(a_i, b_i)
  }
  # Average silhouette value
  mean_silhouette <- mean(silhouette_values)
  print(paste("The mean value of silhouette coefficients:", mean_silhouette))
  object@progress["evaluate"] <- TRUE
  object@clustering$silhouette_coefficient <- mean_silhouette
  invisible(object)
}
