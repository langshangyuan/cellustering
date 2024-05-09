#' Get path to Cellustering example
#'
#' cellustering comes bundled with a number of sample files in its
#' `inst/extdata` directory. This function make them easy to access
#'
#' @param file Name of file. If `NULL`, the example files will be listed.
#' @export
#' @examples
#' cellustering_example()
#' cellustering_example("hg19")
cellustering_example <- function(file = NULL) {
  if (is.null(file)) {
    dir(system.file("extdata", package = "cellustering"))
  } else {
    system.file("extdata", file, package = "cellustering", mustWork = TRUE)
  }
}

#' PBMC Small
#'
#' The `pbmc_small` is a tiny dataset tailored from the `pbmc` dataset for
#' illustration the entire workflow of the `cellustering` package.
#'
#' @format ` pbmc_small`
#' A `Cellustering` instance whose data is a dataframe with 230 rows and 80
#' columns:
#' @source <https://github.com/satijalab/seurat>
"pbmc_small"
