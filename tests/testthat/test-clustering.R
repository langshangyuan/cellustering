test_that(
  "Check the Cellustering will report incorrect operation order",
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
    expect_error(kmeans(pbmc_small_scaled))
    pbmc_small_kmeans <- kmeans(pbmc_small_select, 2)
    expect_equal(
      pbmc_small_kmeans@progress["qc_plot"],
      c(qc_plot = TRUE)
    )
    expect_equal(
      pbmc_small_kmeans@progress["qc_filter"],
      c(qc_filter = TRUE)
    )
    expect_equal(
      pbmc_small_kmeans@progress["normalize"],
      c(normalize = TRUE)
    )
    expect_equal(
      pbmc_small_kmeans@progress["find_HVG"],
      c(find_HVG = TRUE)
    )
    expect_equal(
      pbmc_small_kmeans@progress["scale_data"],
      c(scale_data = TRUE)
    )
    expect_equal(
      pbmc_small_kmeans@progress["PCA"],
      c(PCA = TRUE)
    )
    expect_equal(
      pbmc_small_kmeans@progress["clustering"],
      c(clustering = TRUE)
    )
  }
)

test_that(
  "Check the k-means gets reasonable results",
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
    pbmc_small_kmeans <- kmeans(pbmc_small_select, 2)
    expect_equal(
      sort(unique(pbmc_small_kmeans@clustering$clustering_result$clusters)),
      c(1, 2)
    )
    expect_equal(
      class(pbmc_small_kmeans@clustering$clustering_plot),
      c("gg", "ggplot")
    )
  }
)
