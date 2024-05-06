#' Function Description
#'
#' @param  object  Parameter one Description.
#'
#' @return  Parameter two Description.
#' @export
#'
#' @examples
#' x <- test_value
#' function(test_value) {
#'   # Examples
#'   pbmc <- scale_data(pbmc)
#' }
#' pbmc@reduced_dimension$scaled_data[1:5, 1:10]
#'
#' # Check if the row sums equal to 0
#' row.sums <- rowSums(pbmc@reduced_dimension$scaled_data)
#' all.equal(row.sums, rep(0, 2000), check.attributes = FALSE)
#'
#' # Check if the row sd equal to 1
#' sds <- apply(pbmc@reduced_dimension$scaled_data, 1, sd)
#' all.equal(sds, rep(1, 2000), check.attributes = FALSE)
scale_data <- function(object) {
  # Check if the object belongs to "Cellustering class"
  if (!is(object, "Cellustering")) {
    stop("Please input a Cellustering object.")
  }

  # Check if the user have run `find_HVG`` first
  if (!object@progress["find_HVG"]) {
    stop("Please run find_HVG() first.")
  }

  # Create a new data.frame, containing only the HVG data
  HVG_data <- object@data[object@HVG$HVG_names, , drop = FALSE]

  # Do scaling on this data.frame: minus row mean, then devided by the row sd
  HVG_means <- rowMeans(HVG_data)
  HVG_sds <- apply(HVG_data, 1, sd)
  scaled_data <- (HVG_data - HVG_means) / HVG_sds

  # Save this dataframe into the "reduced_dimension" slot, with name
  # "scaled_data"
  object@reduced_dimension$scaled_data <- scaled_data

  # Return the object
  object@progress["scale_data"] <- TRUE
  return(object)
}