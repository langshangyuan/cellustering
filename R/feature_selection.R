#' Identify, visualize, and keep features called highly variable genes (HVGs).
#'
#' @param  object  A `Cellustering` instance.
#' @param  n_feature  The number of features to select.
#' @param  loess_span  The span parameter for the local regression model,
#' controlling smoothness of the fitted curve.
#'
#' @return  The `Cellusetering` instance with HVG names and plot added to the
#' `HVG` slot.
#' @export
#'
#' @importFrom ggplot2 ggplot geom_point aes geom_point labs scale_x_continuous
#' scale_y_continuous
#' @importFrom stats loess predict var
#'
#' @examples
#' # find_HVG(123) # report error
#' # scale_data(pbmc) # report error
#'
#' # Run the necessary to catch up the progress
#' data("pbmc_small")
#' pbmc_small_5 <- Cellustering(pbmc_small@data)
#' pbmc_small_5 <- qc_plot(pbmc_small_5)
#' pbmc_small_5 <- qc_filter(pbmc_small_5)
#' pbmc_small_5 <- normalize(pbmc_small_5, scale_factor = 1e6)
#'
#' pbmc_small_5 <- find_HVG(pbmc_small_5, n_feature = 20)
#' pbmc_small_5@HVG$HVG_plot
find_HVG <- function(object,
                     n_feature = 2000,
                     loess_span = 0.5) {
  # Check if the object belongs to "Cellustering class"
  if (!is(object, "Cellustering")) {
    stop("Please input a Cellustering object.")
  }

  # Check if other parameters are correctly input
  proper_type <- all(
    is.numeric(n_feature), is.numeric(loess_span)
  )
  proper_length <- all(
    length(n_feature) == 1, length(loess_span) == 1
  )
  proper_range <- all(
    n_feature > 0, identical(n_feature, floor(n_feature)),
    loess_span > 0, loess_span < 1
  )

  if (!(proper_type && proper_length && proper_range)) {
    stop("Invalid input. Please ensure your input is meaningful.")
  }

  # Check if the user have run normalize first
  if (!object@progress["normalize"]) {
    stop("Please run normalize() first.")
  }

  # Compute the mean and variance of each gene's expression (each row)
  gene_means <- apply(object@data, 1, mean)
  gene_vars <- apply(object@data, 1, var)

  # Fit a local regression model using `loess()` function loess.span is a
  # parameter for this function
  fit <- loess(gene_vars ~ gene_means, span = loess_span)

  # Predict each variance using the model and each mean value
  predicted_vars <- predict(fit, newdata = data.frame(gene_means = gene_means))

  # Compute clip.max, it's the square root of the number of cells
  clip_max <- sqrt(ncol(object@data))

  # Scale each count with the row mean and the predicted variance
  scaled_counts <- (object@data - gene_means) / sqrt(predicted_vars)

  # Clip the scaled counts if it's larger than clip.max
  scaled_counts[scaled_counts > clip_max] <- clip_max

  # Recompute each gene's variance
  recomputed_vars <- apply(scaled_counts, 1, var)

  # Select top-n_feature variances as HVG
  HVG_index <- order(recomputed_vars, decreasing = TRUE)[1:n_feature]

  # Add the names of HVGs into the "HVG" slot
  object@HVG$HVG_names <- rownames(object@data)[HVG_index]

  # Prepare data for plotting
  plot_data <- data.frame(
    Mean = gene_means,
    Variance = recomputed_vars,
    HVG = rep(FALSE, length(gene_means))
  )
  # Identify highly variable genes
  plot_data$HVG[HVG_index] <- TRUE

  # Plot without labels
  HVG_plot <- ggplot(
    plot_data,
    aes(x = .data$Mean, y = .data$Variance, color = .data$HVG)
  ) +
    geom_point(size = 0.5) +
    ggplot2::scale_color_manual(
      values = c("#1f78b4", "#33a02c"),
      name = "Gene Counts",
      labels = c(
        paste("Low variance:", sum(!plot_data$HVG)),
        paste("High variance:", sum(plot_data$HVG))
      )
    ) +
    labs(
      title = "Feature Selection Plot",
      x = "Average Expression",
      y = "Standardized Variance"
    ) +
    theme_minimal() +
    scale_x_continuous(
      trans = "sqrt", breaks = c(0.1, 1, 10), limits = c(0, 10)
    ) +
    scale_y_continuous(
      trans = "log", breaks = c(0.1, 1, 10), limits = c(0.1, 10)
    )

  # Save the two plots into the object's HVG slot
  object@HVG$HVG_plot <- HVG_plot

  # Return the object
  object@progress["find_HVG"] <- TRUE
  invisible(object)
}
