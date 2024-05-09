# globalVariables(
#   names = c("count_depth"),
#   package = "Cellustering",
#   add = TRUE
# )

#' Function to plot distribution of total count and feature count or violin plot
#'
#' @param object A `Cellustering` instance.
#'
#' @return  A `Cellusetering` instance that ready to generate plot and for later
#' process.
#' @export
#'
#' @examples
#' pbmc <- read_10x("hg19")
#' # qc_plot(123) # report Error
#' # qc_filter(pbmc) # report Error
#' pbmc_small <- qc_plot(pbmc_small)
#' head(pbmc@quality$cell_data)
#' head(pbmc@quality$feature_data)
qc_plot <- function(object) {
  # Check if the object belongs to "Cellustering class"
  if (!is(object, "Cellustering")) {
    stop("Please input a Cellustering object.")
  }

  # Calculate quantities related to QC
  object <- object |>
    qc_summarize() |>
    qc_plot_cell() |>
    qc_plot_gene() |>
    qc_plot_violin()
  object@progress["qc_plot"] <- TRUE
  # show(object@progress)
  show(object@quality$cell_plot)
  invisible(object)
}

qc_summarize <- function(object) {
  # Calculate quantities related to QC
  original_data <- object@data
  cell_total_counts <- colSums(original_data)
  cell_unique_counts <- colSums(original_data != 0)

  feature_total_counts <- rowSums(original_data)
  feature_unique_counts <- rowSums(original_data != 0)

  mito_genes <- grepl("^(MT|mt)", rownames(original_data))
  mito_counts <- colSums(original_data[mito_genes, , drop = FALSE])
  mito_percent <- mito_counts / cell_total_counts * 100

  # Save the derived quantities into the "QC" slot
  cell_quality <- data.frame(
    count_depth = cell_total_counts,
    feature_count = cell_unique_counts,
    mito_percent = mito_percent
  )

  feature_quality <- data.frame(
    total_expression = feature_total_counts,
    cell_count = feature_unique_counts
  )

  object@quality$cell <- cell_quality
  object@quality$feature <- feature_quality
  object
}

#' @importFrom ggplot2 ggplot aes geom_histogram geom_line geom_point labs
#' theme_minimal
#' @importFrom cowplot plot_grid
qc_plot_cell <- function(object) {
  cell_data <- object@quality$cell
  cell_plot_total <- ggplot(cell_data, aes(x = count_depth)) +
    geom_histogram(
      binwidth = 200, fill = "#1f78b4", color = "white", lwd = 0.2
    ) +
    labs(
      title = "Distribution of Count Depth",
      x = "Count Depth", y = "Frequency"
    ) +
    theme_minimal()

  cell_plot_unique <- ggplot(cell_data, aes(x = feature_count)) +
    geom_histogram(
      binwidth = 50, fill = "#33a02c", color = "white", lwd = 0.2
    ) +
    labs(
      title = "Distribution of Feature Counts",
      x = "Feature Counts", y = "Frequency"
    ) +
    theme_minimal()

  barcode_rank <- rank(-cell_data$count_depth)
  cell_plot_rank <- ggplot(
    cell_data,
    aes(x = barcode_rank, y = count_depth)
  ) +
    geom_line(color = "#6a3d9a") +
    labs(
      title = "Count Depth vs. Barcode Rank",
      x = "Barcode Rank", y = "Count Depth"
    ) +
    theme_minimal()

  cell_plot_joint <- ggplot(cell_data, aes(
    x = count_depth, y = feature_count,
    color = mito_percent
  )) +
    geom_point(position = "jitter", size = 0.001) +
    ggplot2::scale_color_gradientn(
      colors = c(ggplot2::alpha("lightgrey", 0.3), "#6a3d9a")
    ) +
    labs(
      title = "Joint Scatter Plot",
      x = "Count Depth", y = "Feature Counts"
    ) +
    theme_minimal()

  cell_plot <- plot_grid(
    cell_plot_total, cell_plot_unique, cell_plot_rank, cell_plot_joint,
    ncol = 2
  )
  object@quality$cell_plot <- cell_plot
  object
}

#' @importFrom ggplot2 ggplot aes geom_histogram labs theme_minimal
qc_plot_gene <- function(object) {
  # Create QC plots for genes
  feature_data <- object@quality$feature

  trunc_condition <- feature_data$total_expression <= 25
  trunc_feature_total_count <- data.frame(
    total_count = feature_data$total_expression[trunc_condition]
  )

  plot_gene_total <- ggplot2::ggplot(
    trunc_feature_total_count,
    aes(x = total_count)
  ) +
    ggplot2::geom_histogram(
      binwidth = 1, fill = "#1f78b4", color = "white", lwd = 0.2
    ) +
    ggplot2::labs(
      title = "Truncated Distribution of Total Expression",
      x = "Total expression", y = "Frequency"
    ) +
    ggplot2::coord_cartesian(xlim = c(0, 25)) +
    ggplot2::theme_minimal()


  plot_gene_unique <- ggplot2::ggplot(feature_data, aes(x = cell_count)) +
    ggplot2::geom_histogram(
      binwidth = 1, fill = "#33a02c", color = "white", lwd = 0.2
    ) +
    ggplot2::labs(
      title = "Truncated Distribution of Cell Counts",
      x = "Number of Cells", y = "Frequency"
    ) +
    ggplot2::coord_cartesian(xlim = c(0, 25)) +
    ggplot2::theme_minimal()

  feature_plot <- plot_grid(
    plot_gene_total, plot_gene_unique,
    ncol = 2
  )
  object@quality$feature_plot <- feature_plot
  object
}

#' @importFrom ggplot2 ggplot aes geom_violin labs theme_minimal
qc_plot_violin <- function(object) {
  cell_data <- object@quality$cell
  violin_cell_total <- ggplot(cell_data, aes(x = 1, y = count_depth)) +
    geom_violin(fill = "#1f78b4", color = "black", lwd = 0.2) +
    labs(title = "Count Depth", x = "", y = "") +
    theme_minimal()

  violin_cell_unique <- ggplot(
    cell_data,
    aes(x = 1, y = feature_count)
  ) +
    geom_violin(fill = "#33a02c", color = "black", lwd = 0.2) +
    labs(title = "Feature Counts", x = "", y = "") +
    theme_minimal()

  violin_cell_mito <- ggplot(cell_data, aes(x = 1, y = mito_percent)) +
    geom_violin(fill = "#6a3d9a", color = "black", lwd = 0.2) +
    labs(title = "Mitochondrial Percent", x = "", y = "") +
    theme_minimal()
  violin_plot <- plot_grid(
    violin_cell_total, violin_cell_unique, violin_cell_mito,
    ncol = 3
  )

  object@quality$violin_plot <- violin_plot
  object
}

#' Filter out low-quality cells and genes based on user-defined conditions.
#'
#' @param  object  A `Cellustering` instance.
#' @param  min_count_depth  The minimum count depth allowed for a cell.
#' @param  max_count_depth  The maximum count depth allowed for a cell.
#' @param  min_feature  The minimum number of features allowed for a cell.
#' @param  max_feature  The maximum number of features allowed for a cell.
#' @param  max_mito_percent  The maximum mitochondrial gene percentage allowed
#' for a cell.
#' @param  min_total_expression  The minimum total expression allowed for a
#' gene.
#' @param  min_cell  The minimum number of cells in which a gene must be
#'                   expressed.

#' @return  A `Cellustering` instance with cells and genes satisfying the
#' filtering conditions.
#' @export
#'
#' @examples
#' # qc_filter(123) # report Error
#' # qc_filter(pbmc, min.count.depth = -1) # report Error
#' # qc_filter(pbmc, min.cells = "abc") # report Error
#' # Run the necessary to catch up the progress
#' pbmc_small <- qc_plot(pbmc_small)
#'
#' pbmc_small <- qc_filter(pbmc_small)
#' dim(pbmc_small@data)
#'
#' pbmc_small <- qc_filter(pbmc_small,
#'   min_feature = 200, max_feature = 2500,
#'   max_mito_percent = 5, min_cell = 3
#' )
#' dim(pbmc_small@data)
#' # Then we check the plots again
#' pbmc_small <- qc_plot(pbmc_small)
qc_filter <- function(object,
                      min_count_depth = 0,
                      max_count_depth = Inf,
                      min_feature = 0,
                      max_feature = Inf,
                      max_mito_percent = 100,
                      min_total_expression = 0,
                      min_cell = 0) {
  # Check if the object belongs to "Cellustering class"
  if (!is(object, "Cellustering")) {
    stop("Please input a Cellustering object.")
  }

  object <- object |>
    qc_summarize()

  # Check if other parameters are in correct data types
  proper_type <- all(
    is.numeric(min_count_depth), is.numeric(max_count_depth),
    is.numeric(min_feature), is.numeric(max_feature),
    is.numeric(max_mito_percent), is.numeric(min_total_expression),
    is.numeric(min_cell)
  )
  proper_length <- all(
    length(min_count_depth) == 1, length(max_count_depth) == 1,
    length(min_feature) == 1, length(max_feature) == 1,
    length(max_mito_percent), length(min_total_expression) == 1,
    length(min_cell) == 1
  )
  proper_range <- all(
    min_count_depth >= 0, max_count_depth >= 0,
    min_feature >= 0, max_feature >= 0,
    max_mito_percent >= 0, max_mito_percent <= 100,
    min_total_expression >= 0, min_cell >= 0
  )
  if (!(proper_type && proper_length && proper_range)) {
    stop("Invalid input. Please ensure your input is meaningful.")
  }
  if (!object@progress["qc_plot"]) {
    stop("Please run qc_plot() first.")
  }
  # Filter cells and genes satisfying given constraints
  gene_index <- which(
    object@quality$feature$total_expression >= min_total_expression &
      object@quality$feature$cell_count >= min_cell
  )
  cell_index <- which(
    object@quality$cell$count_depth >= min_count_depth &
      object@quality$cell$count_depth <= max_count_depth &
      object@quality$cell$feature_count >= min_feature &
      object@quality$cell$feature_count <= max_feature &
      object@quality$cell$mito_percent <= max_mito_percent
  )
  object@data <- object@data[gene_index, cell_index]
  # Return the object
  object@progress["qc_filter"] <- TRUE
  invisible(object)
}
