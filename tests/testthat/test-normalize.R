test_that(
  "Check the Cellustering will report incorrect operation order",
  {
    data("pbmc_small")
    pbmc_small_raw <- Cellustering(pbmc_small@data)
    pbmc_small_plot <- qc_plot(pbmc_small_raw)
    pbmc_small_filter <- qc_filter(pbmc_small_plot)
    pbmc_small_normalized <- normalize(pbmc_small_filter)
    expect_equal(
      pbmc_small_normalized@progress["qc_plot"],
      c(qc_plot = TRUE)
    )
    expect_equal(
      pbmc_small_normalized@progress["qc_filter"],
      c(qc_filter = TRUE)
    )
    expect_equal(
      pbmc_small_normalized@progress["normalize"],
      c(normalize = TRUE)
    )
  }
)

test_that("The normalize only remove the variations in count depths", {
  data("pbmc_small")
  pbmc_small_raw <- Cellustering(pbmc_small@data)
  pbmc_small_plot <- qc_plot(pbmc_small_raw)
  pbmc_small_filter <- qc_filter(pbmc_small_plot)
  expect_error(normalize(123))
  expect_error(normalize(pbmc_small_filter, scale_factor = -100))
  col_sum_before <- sum(pbmc_small_filter@data[, 1])
  element_before <- pbmc_small_filter@data[30, 1]
  pbmc_small_normalized <- normalize(pbmc_small_filter)
  col_sum_after <- sum(pbmc_small_normalized@data[, 1])
  element_after <- pbmc_small_normalized@data[30, 1]
  expect_equal(element_after, log(element_before / col_sum_before * 1e6 + 1))
  expect_equal(sum(exp(pbmc_small_normalized@data[, 1]) - 1), 1e6)
})
