test_that(
  "Check the Cellustering will report incorrect operation order",
  {
    data("pbmc_small")
    pbmc_small_raw <- Cellustering(pbmc_small@data)
    pbmc_small_plot <- qc_plot(pbmc_small_raw)
    pbmc_small_filter <- qc_filter(pbmc_small_plot)
    expect_error(find_HVG(pbmc_small_filter))
    pbmc_small_normalized <- normalize(pbmc_small_filter)
    pbmc_small_HVG <- find_HVG(pbmc_small_normalized, n_feature = 20)
    expect_equal(
      pbmc_small_HVG@progress["qc_plot"],
      c(qc_plot = TRUE)
    )
    expect_equal(
      pbmc_small_HVG@progress["qc_filter"],
      c(qc_filter = TRUE)
    )
    expect_equal(
      pbmc_small_HVG@progress["normalize"],
      c(normalize = TRUE)
    )
    expect_equal(
      pbmc_small_HVG@progress["find_HVG"],
      c(find_HVG = TRUE)
    )
  }
)

test_that(
  "The find_HVG only saves the highly variable genes",
  {
    data("pbmc_small")
    pbmc_small_raw <- Cellustering(pbmc_small@data)
    pbmc_small_plot <- qc_plot(pbmc_small_raw)
    pbmc_small_filter <- qc_filter(pbmc_small_plot)
    pbmc_small_normalized <- normalize(pbmc_small_filter)
    pbmc_small_HVG <- find_HVG(pbmc_small_normalized, n_feature = 20)
    expect_equal(
      length(pbmc_small_HVG@HVG$HVG_name),
      20
    )
    expect_equal(
      class(pbmc_small_HVG@HVG$HVG_plot),
      c("gg", "ggplot")
    )
  }
)
