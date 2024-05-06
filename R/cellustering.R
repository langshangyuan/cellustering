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
    progress = vector(mode = "logical", length = 8)
  )
)

# Initialize a Cellustering instance
Cellustering <- function(data,
                         quality = list(),
                         HVG = list(),
                         reduced_dimension = list(),
                         clustering = list()) {
  progress <- c(
    qc_plot = FALSE,
    qc_filter = FALSE,
    normalize = FALSE,
    find_HVG = FALSE,
    scale_data = FALSE,
    PCA = FALSE,
    dimension_selection = FALSE,
    kmeans = FALSE,
    run.evaluate = FALSE
  )

  new("Cellustering",
    data = data,
    quality = quality,
    HVG = HVG,
    reduced_dimension = reduced_dimension,
    clustering = clustering,
    progress = progress
  )
}
