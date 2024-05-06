#' Function Description
#' @param  object  Parameter one Description.
#' @param  object  Parameter one Description.
#'
#' @return  Parameter two Description.
#' @export
#'
#' @examples
#' # principal_component_analysis(123) # report error
#' # principal_component_analysis(pbmc, PC1 = 0.1, PC2 = 3000)  # report error
#' pbmc <- principal_component_analysis(pbmc, PC1 = 3, PC2 = 4)
#' @importFrom ggplot2 theme element_text
principal_component_analysis <- function(object,
                                         PC1 = 1,
                                         PC2 = 2) {
  # Check if the object belongs to "Cellustering class"
  if (!is(object, "Cellustering")) {
    stop("Please input a Cellustering object.")
  }

  # Check if the user have run `scale_data` first
  if (!object@progress["scale_data"]) {
    stop("Please run scale_data() first.")
  }
  features <- object@reduced_dimension$scaled_data

  # Check if other parameters are in correct data types
  proper_type <- all(
    is.numeric(PC1), is.numeric(PC2)
  )
  proper_length <- all(
    length(PC1) == 1, length(PC2) == 1
  )
  proper_range <- all(
    PC1 <= nrow(features), PC2 <= nrow(features)
  )


  if (!(proper_type && proper_length && proper_range)) {
    stop("Invalid input. Please ensure your input is meaningful.")
  }

  feature_matrix <- as.matrix(features)
  covariance_matrix <- 1 / ncol(features) * feature_matrix %*% t(feature_matrix)
  eigen_decomposition <- eigen(covariance_matrix, symmetric = TRUE)
  eigenvalues <- eigen_decomposition$values
  eigenvectors <- eigen_decomposition$vectors
  object@reduced_dimension$eigenvalues <- eigenvalues
  object@reduced_dimension$eigenvectors <- eigenvectors

  # Plot the two-dimensional reconstruction
  filtered_eigenvectors <- eigenvectors[, c(PC1, PC2)]
  projected_data <- data.frame(t(feature_matrix) %*% filtered_eigenvectors)
  colnames(projected_data) <- c("PC1", "PC2")
  PCA_plot <- ggplot(projected_data, aes(x = PC1, y = PC2)) +
    geom_point(shape = 16, size = 1, color = "#1f78b4") +
    labs(
      title = "PCA Plot",
      x = paste("PC1 ( =", PC1, ")"),
      y = paste("PC2 ( =", PC2, ")")
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 16, face = "bold"),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12),
      legend.position = "none"
    )
  object@reduced_dimension$PC1 <- PC1
  object@reduced_dimension$PC2 <- PC2
  object@reduced_dimension$PCA_plot <- PCA_plot
  print(PCA_plot)

  # Return the object
  object@progress["PCA"] <- TRUE
  invisible(object)
}

# Examples

#' Function Description
#'
#' @param  parameter_1  Parameter one Description.
#' @param  parameter_2  Parameter two Description.
#'
#' @return  Parameter two Description.
#' @export
#' @importFrom utils head
#' @importFrom stats dnorm sd
#'
#' @examples
select_proper_dimension <- function(object,
                                    lower_bound = 1,
                                    upper_bound = 30) {
  # Check if the object belongs to "Cellustering class"
  if (!is(object, "Cellustering")) {
    stop("Please input a Cellustering object.")
  }

  # Check if the user have run principal_component_analysis first
  if (!object@progress["PCA"]) {
    stop("Please run principal_component_analysis() first.")
  }

  features <- object@reduced_dimension$scaled_data
  N <- nrow(features)

  # Check if other parameters are in correct data types
  proper_type <- all(
    is.numeric(lower_bound), is.numeric(upper_bound)
  )
  proper_length <- all(
    length(lower_bound) == 1, length(upper_bound) == 1
  )
  proper_range <- all(
    lower_bound > 0, identical(lower_bound, floor(lower_bound)),
    upper_bound > 0, identical(upper_bound, floor(upper_bound)),
    lower_bound <= N, upper_bound <= N, lower_bound <= upper_bound
  )


  if (!(proper_type && proper_length && proper_range)) {
    stop("Invalid input. Please ensure your input is meaningful.")
  }
  # Do the Scree plot
  eigenvalues <- object@reduced_dimension$eigenvalues
  if (upper_bound < N) {
    tail_eigenvalue_sum <- c(sum(eigenvalues[(upper_bound + 1):N]))
  } else {
    tail_eigenvalue_sum <- c(0)
  }
  if (upper_bound > lower_bound) {
    for (i in upper_bound:(lower_bound + 1)) {
      tail_eigenvalue_sum <- c(
        tail_eigenvalue_sum, tail(tail_eigenvalue_sum, 1) + eigenvalues[i]
      )
    }
  }
  tail_eigenvalue_sum <- rev(tail_eigenvalue_sum)
  scree_plot_data <- data.frame(
    dim = lower_bound:upper_bound,
    eigenvalue_sum = tail_eigenvalue_sum
  )
  scree_plot <- ggplot(scree_plot_data, aes(x = dim, y = eigenvalue_sum)) +
    geom_point(size = 1, color = "#1f78b4") +
    geom_line(color = "#1f78b4") +
    labs(
      title = "Scree Plot",
      x = paste("Number of PC"),
      y = paste("Eigenvalue Sums")
    ) +
    theme_minimal()
  print(scree_plot)

  # Do the profile-likelihood plot
  eigenvalue_sum <- sum(eigenvalues)
  squared_eigenvalue_sum <- sum(eigenvalues^2)
  ## derive two means
  head_eigenvalue_sum <- eigenvalue_sum - tail_eigenvalue_sum
  tail_eigenvalue_mean <- tail_eigenvalue_sum / (N - (lower_bound:upper_bound))
  head_eigenvalue_mean <- head_eigenvalue_sum / (lower_bound:upper_bound)
  ## derive variance
  if (upper_bound < N) {
    squared_tail_eigenvalue_sum <- c(sum(eigenvalues[(upper_bound + 1):N]^2))
  } else {
    squared_tail_eigenvalue_sum <- c(0)
  }
  if (upper_bound > lower_bound) {
    for (i in upper_bound:(lower_bound + 1)) {
      squared_tail_eigenvalue_sum <- c(
        squared_tail_eigenvalue_sum,
        tail(squared_tail_eigenvalue_sum, 1) + eigenvalues[i]^2
      )
    }
  }
  squared_tail_eigenvalue_sum <- rev(squared_tail_eigenvalue_sum)
  squared_head_eigenvalue_sum <- squared_eigenvalue_sum -
    squared_tail_eigenvalue_sum

  vars <- (squared_head_eigenvalue_sum -
    head_eigenvalue_mean^2 * (lower_bound:upper_bound) +
    squared_tail_eigenvalue_sum - tail_eigenvalue_mean^2 *
      (N - (lower_bound:upper_bound))
  ) / N
  if (upper_bound == N) {
    vars[length(vars)] <- squared_head_eigenvalue_sum[length(vars)] / N -
      head_eigenvalue_mean[length(vars)]^2
  }
  ## derive profile log likelihood
  log.likelihoods <- rep(0, length(lower_bound:upper_bound))
  for (i in 1:length(lower_bound:upper_bound)) {
    L <- (lower_bound:upper_bound)[i]
    if (L < N) {
      log.likelihood <-
        sum(dnorm(eigenvalues[1:L],
          mean = head_eigenvalue_mean[i], sd = sqrt(vars[i]),
          log = TRUE
        )) +
        sum(dnorm(eigenvalues[(L + 1):N],
          mean = tail_eigenvalue_mean[i], sd = sqrt(vars[i]),
          log = TRUE
        ))
    } else {
      log.likelihood <-
        sum(dnorm(eigenvalues[1:N],
          mean = head_eigenvalue_mean[i], sd = sqrt(vars[i]),
          log = TRUE
        ))
    }
    log.likelihoods[i] <- log.likelihood
  }
  ## do plot
  likelihood_plot_data <- data.frame(
    dim = lower_bound:upper_bound,
    l = log.likelihoods
  )
  likelihood_plot <- ggplot(likelihood_plot_data, aes(x = dim, y = l)) +
    geom_point(size = 1, color = "#33a02c") +
    labs(
      title = "Profile Log Likelihood Plot",
      x = paste("Number of PC"),
      y = paste("Profile Log Likelihood")
    ) +
    theme_minimal()
  if (lower_bound < upper_bound) {
    likelihood_plot <- likelihood_plot +
      geom_line(
        mapping = aes(x = dim, y = l),
        data = likelihood_plot_data,
        color = "#33a02c"
      )
  }
  ## suggest dimension
  suggested_dimension <- (lower_bound:upper_bound)[which.max(log.likelihoods)]
  print(paste("The suggested dimension to reduce to:", suggested_dimension))
  object@reduced_dimension$suggested_dimension <- suggested_dimension

  # Combine, save, and show plot
  dimension_selection_plot <- cowplot::plot_grid(
    scree_plot, likelihood_plot,
    ncol = 2
  )
  print(dimension_selection_plot)
  object@reduced_dimension$dimension_selection_plot <- dimension_selection_plot

  # Return the object
  object@progress["dimension_selection"] <- TRUE
  invisible(object)
}
