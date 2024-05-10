test_that(
  "The qc_plot function will summarize the data, generate nice visualizations",
  {
    data("pbmc_small")
    expect_error(qc_plot(pbmc_small@data))
    pbmc_small_plot <- Cellustering(pbmc_small@data)
    pbmc_small_plot <- qc_plot(pbmc_small_plot)
    expect_equal(
      dim(pbmc_small_plot@quality$cell),
      c(dim(pbmc_small@data)[2], 3)
    )
    expect_equal(
      dim(pbmc_small_plot@quality$feature),
      c(dim(pbmc_small@data)[1], 2)
    )
    expect_equal(
      class(pbmc_small_plot@quality$cell_plot),
      c("gg", "ggplot")
    )
    expect_equal(
      class(pbmc_small_plot@quality$feature_plot),
      c("gg", "ggplot")
    )
    expect_equal(
      class(pbmc_small_plot@quality$violin_plot),
      c("gg", "ggplot")
    )
    expect_equal(
      pbmc_small_plot@progress["qc_plot"],
      c(qc_plot = TRUE)
    )
  }
)

test_that(
  "Check the Cellustering will report incorrect operation order",
  {
    data("pbmc_small")
    pbmc_small_raw <- Cellustering(pbmc_small@data)
    expect_error(qc_filter(pbmc_small_raw))
    pbmc_small_plot <- qc_plot(pbmc_small_raw)
    pbmc_small_filter <- qc_filter(pbmc_small_plot)
    expect_equal(
      pbmc_small_filter@progress["qc_plot"],
      c(qc_plot = TRUE)
    )
    expect_equal(
      pbmc_small_filter@progress["qc_filter"],
      c(qc_filter = TRUE)
    )
  }
)

test_that(
  "The qc_filter function should have the idempotent property",
  {
    data("pbmc_small")
    pbmc_small_raw <- Cellustering(pbmc_small@data)
    pbmc_small_plot <- qc_plot(pbmc_small_raw)
    pbmc_small_filter <- qc_filter(pbmc_small_plot)
    expect_equal(
      dim(pbmc_small_plot@data),
      dim(pbmc_small_filter@data)
    )
    pbmc_small_filter_2 <- qc_filter(pbmc_small_filter)
    expect_equal(
      dim(pbmc_small_filter@data),
      dim(pbmc_small_filter_2@data)
    )
    expect_error(
      qc_filter(
        pbmc_small_plot,
        min_count_depth = -1,
        max_count_depth = -2, min_feature = -1,
        max_mito_percent = 101
      )
    )
    pbmc_small_filter_3 <- qc_filter(pbmc_small_filter_2,
      min_feature = 200, max_feature = 2500,
      max_mito_percent = 5
    )
    pbmc_small_filter_4 <- qc_filter(pbmc_small_filter_3,
      min_feature = 200, max_feature = 2500,
      max_mito_percent = 5
    )
    expect_equal(
      dim(pbmc_small_filter_3@data),
      dim(pbmc_small_filter_4@data)
    )
  }
)
