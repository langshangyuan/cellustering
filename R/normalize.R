#' Normalize the data using the Counts per Million (CPM) normalization
#' protocol. The purpose of normalization is to remove the variations in count
#' depths among different cells, and to place the gene expression levels of
#' different cells on a comparable scale.
#'
#' @param object A `Cellustering` instance.
#' @param scale_factor The level to which each cell's total count depth is
#' normalized.
#' @param log_transformation Determine whether to perform log-transformation
#' after the CPM normalization.
#'
#' @return  The `Cellustering` instance with normalized data stored in the
#' `data` slot.
#' @export
#'
#' @examples
#' # normalize(123) # report error
#' # normalize(pbmc, scale.factor = -100) # report error
#' # normalize(pbmc, log.transformation = "abc") # report error
#' # find_HVG(pbmc) # report error
#'
#' # Run the necessary to catch up the progress
#' data("pbmc_small")
#' pbmc_small_4 <- Cellustering(pbmc_small@data)
#' pbmc_small_4 <- qc_plot(pbmc_small_4)
#' pbmc_small_4 <- qc_filter(pbmc_small_4)
#'
#' col_sum_before <- sum(pbmc_small_4@data[, 1])
#' element_before <- pbmc_small_4@data[30, 1]
#' pbmc_small_4 <- normalize(pbmc_small_4, scale_factor = 1e6)
#' col_sum_after <- sum(pbmc_small_4@data[, 1])
#' element_after <- pbmc_small_4@data[30, 1]
#' all.equal(element_after, log(element_before / col_sum_before * 1e6 + 1))
#' all.equal(sum(exp(pbmc_small_4@data[, 1]) - 1), 1e6)
normalize <- function(object,
                      scale_factor = 1e6,
                      log_transformation = TRUE) {
  # Check if the object belongs to "Cellustering class"
  if (!is(object, "Cellustering")) {
    stop("Please input a Cellustering object.")
  }

  # Check if other parameters are in correct data types
  proper_type <- all(
    is.numeric(scale_factor), is.logical(log_transformation)
  )
  proper_length <- all(
    length(scale_factor) == 1, length(log_transformation) == 1
  )
  proper_range <- scale_factor > 0

  if (!(proper_type && proper_length && proper_range)) {
    stop("Invalid input. Please ensure your input is meaningful.")
  }

  # Perform CPM transformation
  object@data <- object@data |>
    dplyr::mutate_at(
      names(object@data),
      function(col) col / sum(col) * scale_factor
    )

  # Perform log.transformation
  if (log_transformation) {
    object@data <- log(object@data + 1)
  }

  # Return the object
  object@progress["normalize"] <- TRUE
  invisible(object)
}
