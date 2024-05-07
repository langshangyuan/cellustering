#' Function Description
#'
#' @param  parameter_1  Parameter one Description.
#' @param  parameter_2  Parameter two Description.
#'
#' @return  Parameter two Description.
#' @export
#'
#' @importFrom stats kmeans
#'
#' @examples
# library(umap)
kmeans <- function(object,
                   k,
                   dimensions = 10,
                   umap = FALSE) {
  # Check if the object belongs to "Cellustering class"
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
  # print(cluster_assignment)
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
  #   # new_cluster_assignment <- apply(distance_matrix, 2, which.min)
  #   print(time)
  #   new_cluster_assignment <- apply(distance_matrix, 1, which.min)
  #   print(new_cluster_assignment)
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
  # while (TRUE) {
  #   for (i in 1:N) {
  #     for (j in 1:k) {
  #       distance_matrix[i, j] <- distance(projected_data[i, ], centroids[j, ])
  #     }
  #   }
  #   # Assign points to closest centroid
  #   # new_cluster_assignment <- apply(distance_matrix, 2, which.min)
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
  # print(cluster_assignment)
  # cluster_result <- list(clusters = cluster_assignment, centroids = centroids)

  # Visualize
  visualization_base <- object@reduced_dimension$PCA_plot
  ## (如果 umap == FALSE，用 pca visualize)
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

#' Function Description
#'
#' @param object  Parameter one Description.
#'
#' @return  Parameter two Description.
#' @export
#'
#' @importFrom stats dist
#'
#' @examples
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
  mean(silhouette_values)
}
