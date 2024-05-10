# Define the "Cellustering" class
#' @importFrom methods is new show
setClass("Cellustering",
  slots = c(
    data = "data.frame",
    quality = "list",
    HVG = "list",
    reduced_dimension = "list",
    clustering = "list",
    progress = "logical"
  ),
  prototype = list(
    data = data.frame(),
    quality = list(),
    HVG = list(),
    reduced_dimension = list(),
    clustering = list(),
    progress = vector(mode = "logical", length = 9)
  )
)

#' Initialize a Cellustering instance
#'
#' @param data The data storing 10X single-cell RNA data
#'
#' @return A fresh Cellustering instance with data describing the cell names
#' and their features.
#' @export
Cellustering <- function(data) {
  progress <- c(
    qc_plot = FALSE,
    qc_filter = FALSE,
    normalize = FALSE,
    find_HVG = FALSE,
    scale_data = FALSE,
    PCA = FALSE,
    dimension_selection = FALSE,
    kmeans = FALSE,
    evaluate = FALSE
  )

  new("Cellustering",
    data = data,
    quality = list(),
    HVG = list(),
    reduced_dimension = list(),
    clustering = list(),
    progress = progress
  )
}
