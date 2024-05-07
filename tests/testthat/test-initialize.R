test_that(
  "A single RNA dataset stored in the format of 10x can be successfully read",
  {
    pbmc <- read_10x("hg19")
    expect_equal(dim(pbmc@data), c(32738, 2700))
  }
)

test_that("Try to read a non-existence directory", {
  expect_error(read_10x("pbmc"))
})
