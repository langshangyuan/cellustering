test_that(
  "Check the Cellustering will report incorrect operation order",
  {
    data("pbmc_small")
    pbmc_small_raw <- Cellustering(pbmc_small@data)
    pbmc_small_plot <- qc_plot(pbmc_small_raw)
    pbmc_small_filter <- qc_filter(pbmc_small_plot)
    pbmc_small_normalized <- normalize(pbmc_small_filter)
    expect_error(scale_data(pbmc_small_normalized))
    pbmc_small_HVG <- find_HVG(pbmc_small_normalized, n_feature = 20)
    pbmc_small_scaled <- scale_data(pbmc_small_HVG)
    expect_equal(
      pbmc_small_scaled@progress["qc_plot"],
      c(qc_plot = TRUE)
    )
    expect_equal(
      pbmc_small_scaled@progress["qc_filter"],
      c(qc_filter = TRUE)
    )
    expect_equal(
      pbmc_small_scaled@progress["normalize"],
      c(normalize = TRUE)
    )
    expect_equal(
      pbmc_small_scaled@progress["find_HVG"],
      c(find_HVG = TRUE)
    )
    expect_equal(
      pbmc_small_scaled@progress["scale_data"],
      c(scale_data = TRUE)
    )
  }
)

test_that("The scale_data ensures the the mean of HVG is 0 and variance is 1", {
  data("pbmc_small")
  pbmc_small_raw <- Cellustering(pbmc_small@data)
  pbmc_small_plot <- qc_plot(pbmc_small_raw)
  pbmc_small_filter <- qc_filter(pbmc_small_plot)
  pbmc_small_normalized <- normalize(pbmc_small_filter)
  pbmc_small_HVG <- find_HVG(pbmc_small_normalized, n_feature = 20)
  pbmc_small_scaled <- scale_data(pbmc_small_HVG)
  row_sums <- rowSums(pbmc_small_scaled@reduced_dimension$scaled_data)
  expect_equal(row_sums, rep(0, 20), ignore_attr = TRUE)
  # Check if the row sd equal to 1
  sds <- apply(pbmc_small_scaled@reduced_dimension$scaled_data, 1, sd)
  expect_equal(sds, rep(1, 20), ignore_attr = TRUE)
})
