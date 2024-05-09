#' Initialize a `Cellustering` instance from a directory containing files
#' storing 10X single-cell RNA data in Market Exchange (MEX) file format.
#'
#' @param directory A directory containing single cell RNA data.
#'
#' @return  A `Cellustering` instance.
#' @export
#'
#' @importFrom utils read.table
#'
#' @examples
#' # pbmc <- read_10x(cellustering_example("hg19"))
#' # pbmc@data[c("CD3D", "TCL1A", "MS4A1"), 1:30]
read_10x <- function(directory) {
  # Check if the directory exists
  if (!(is.character(directory) && dir.exists(directory))) {
    stop(
      "Please input an valid directory",
      call. = FALSE
    )
  }

  # Check if the directory contains required files
  valid_matrix <- check_file(directory, "matrix.mtx")
  valid_barcodes <- check_file(directory, "barcodes.tsv")
  valid_genes <- check_file(directory, "genes.tsv")
  qualified <- all(valid_matrix, valid_barcodes, valid_genes)
  if (!qualified) {
    stop("Please input a directory path containing required files.")
  }

  # Read in the data
  matrix <- as.matrix(
    Matrix::readMM(paste0(directory, "/matrix.mtx"))
  )
  barcodes <- readLines(paste0(directory, "/barcodes.tsv"))
  genes <- read.table(
    paste0(directory, "/genes.tsv"),
    header = FALSE, stringsAsFactors = FALSE
  )[, 2]

  # Create a Cellustering instance
  colnames(matrix) <- barcodes
  rownames(matrix) <- genes
  data <- data.frame(matrix)
  invisible(Cellustering(data = data))
}

check_file <- function(path, pattern) {
  list.files(path = path, pattern = pattern, full.names = TRUE)
  length(list.files) != 0
}
