test_that(
  "A single RNA dataset stored in the format of 10x can be successfully read",
  {
    pbmc <- read_10x(cellustering_example("hg19"))
    expect_equal(dim(pbmc@data), c(32738, 2700))
  }
)

test_that("Try to read a non-existence directory or non-invalid input", {
  expect_error(read_10x("pbmc"))
  expect_error(read_10x(pbmc))
})
