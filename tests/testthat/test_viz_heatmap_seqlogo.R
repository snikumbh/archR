context("Test heatmap and seqlogo")
library(ggplot2)
library(ggseqlogo)
library(gridExtra)

test_that("Given object is matrix", {
  testFeaturesMat <- rnorm(10000) # err
  testPositionLabels <- seq(25)
  expect_error(
    viz_all_factors_in_combined_heatmaps_seqlogos(testFeaturesMat, position_labels = testPositionLabels),
    "not of type matrix"
  )
})

test_that("Handling empty matrix", {
  testFeaturesMat <- matrix()
  testPositionLabels <- seq(25)
  expect_error(
    viz_all_factors_in_combined_heatmaps_seqlogos(testFeaturesMat,
      position_labels = testPositionLabels
    ),
    "Empty"
  )
})

test_that("Position labels inadequate", {
  testFeaturesMat <- matrix(rnorm(10000), nrow = 200)
  testPositionLabels <- seq(20)
  expect_error(
    viz_all_factors_in_combined_heatmaps_seqlogos(testFeaturesMat, position_labels = testPositionLabels),
    "Inadequate"
  )
})

test_that("Position labels over-abundant", {
  testFeaturesMat <- matrix(rnorm(10000), nrow = 200)
  testPositionLabels <- seq(60)
  expect_error(
    viz_all_factors_in_combined_heatmaps_seqlogos(testFeaturesMat, position_labels = testPositionLabels),
    "Overabundant"
  )
})

test_that("Combined heatmap and seqlogo works", {
  # setting seed enables proper comparison between ggplot objects since we use
  # rnorm
  set.seed(11223344)
  # test variables
  testPositionLabels <- seq(25)
  testFeaturesMat <- matrix(rnorm(100), ncol = 1)
  # test plot, include the function directly inside because it returns nothing.
  vdiffr::expect_doppelganger(
    "heatmap seqlogo plot ex",
    viz_all_factors_in_combined_heatmaps_seqlogos(featuresMatrix = testFeaturesMat, position_labels = testPositionLabels)
  )
})
