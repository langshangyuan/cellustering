# Define the "Cellustering" class
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

# setValidity("Cellustering", function(object) {
#   if (length(object@data) == 0) {

#   }
# })

# Initialize a Cellustering instance
Cellustering <- function(data,
                         quality = list(),
                         HVG = list(),
                         reduced_dimension = list(),
                         clustering = list()) {
  progress <- c(
    qc.plot = FALSE,
    qc.filter = FALSE,
    run.normalize = FALSE,
    find.HVG = FALSE,
    scale.data = FALSE,
    run.PCA = FALSE,
    dim.select = FALSE,
    run.kmeans = FALSE,
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
