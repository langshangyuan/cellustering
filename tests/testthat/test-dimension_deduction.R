test_that(
  "Check the Cellustering will report incorrect operation order",
  {
    data("pbmc_small")
    pbmc_small_raw <- Cellustering(pbmc_small@data)
    pbmc_small_plot <- qc_plot(pbmc_small_raw)
    pbmc_small_filter <- qc_filter(pbmc_small_plot)
    pbmc_small_normalized <- normalize(pbmc_small_filter)
    pbmc_small_HVG <- find_HVG(pbmc_small_normalized, n_feature = 20)
    expect_error(principal_component_analysis(pbmc_small_HVG))
    pbmc_small_scaled <- scale_data(pbmc_small_HVG)
    pbmc_small_PCA <- principal_component_analysis(pbmc_small_scaled)
    expect_equal(
      pbmc_small_PCA@progress["qc_plot"],
      c(qc_plot = TRUE)
    )
    expect_equal(
      pbmc_small_PCA@progress["qc_filter"],
      c(qc_filter = TRUE)
    )
    expect_equal(
      pbmc_small_PCA@progress["normalize"],
      c(normalize = TRUE)
    )
    expect_equal(
      pbmc_small_PCA@progress["find_HVG"],
      c(find_HVG = TRUE)
    )
    expect_equal(
      pbmc_small_PCA@progress["scale_data"],
      c(scale_data = TRUE)
    )
    expect_equal(
      pbmc_small_PCA@progress["PCA"],
      c(PCA = TRUE)
    )
  }
)

test_that(
  "Check the PCA will get the desired form",
  {
    data("pbmc_small")
    pbmc_small_raw <- Cellustering(pbmc_small@data)
    pbmc_small_plot <- qc_plot(pbmc_small_raw)
    pbmc_small_filter <- qc_filter(pbmc_small_plot)
    pbmc_small_normalized <- normalize(pbmc_small_filter)
    pbmc_small_HVG <- find_HVG(pbmc_small_normalized, n_feature = 20)
    pbmc_small_scaled <- scale_data(pbmc_small_HVG)
    pbmc_small_PCA <- principal_component_analysis(pbmc_small_scaled)
    expect_equal(
      length(pbmc_small_PCA@reduced_dimension$eigenvalues),
      20
    )
    expect_equal(
      dim(pbmc_small_PCA@reduced_dimension$eigenvectors),
      c(20, 20)
    )
    expect_equal(
      class(pbmc_small_PCA@reduced_dimension$PCA_plot),
      c("gg", "ggplot")
    )
  }
)

test_that(
  "Check the select_dimension will produce the desired plot",
  {
    data("pbmc_small")
    pbmc_small_raw <- Cellustering(pbmc_small@data)
    pbmc_small_plot <- qc_plot(pbmc_small_raw)
    pbmc_small_filter <- qc_filter(pbmc_small_plot)
    pbmc_small_normalized <- normalize(pbmc_small_filter)
    pbmc_small_HVG <- find_HVG(pbmc_small_normalized, n_feature = 20)
    pbmc_small_scaled <- scale_data(pbmc_small_HVG)
    pbmc_small_PCA <- principal_component_analysis(pbmc_small_scaled)
    pbmc_small_select <- select_proper_dimension(pbmc_small_PCA,
      lower_bound = 1, upper_bound = 5
    )
    expect_equal(
      pbmc_small_select@reduced_dimension$suggested_dimension,
      2
    )
    expect_equal(
      class(pbmc_small_select@reduced_dimension$dimension_selection_plot),
      c("gg", "ggplot")
    )
  }
)
