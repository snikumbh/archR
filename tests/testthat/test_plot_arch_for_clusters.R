context("Test plotting architectures")

test_that("Given samples matrix object is matrix", {

  testSamplesMat <- rnorm(10000) # err
  testFeaturesMat <- matrix(rnorm(10000), nrow = 200)
  testPositionLabels <- seq(50)
  nClusterTest <- 5
  clusterSolTest <- list(clust_sol = c(), reordering_idx = c())
  seqsTest <- NULL
  expect_error(plot_arch_for_clusters(testSamplesMat, testFeaturesMat,
                                      nCluster = nClusterTest,
                                      clustering_sol = clusterSolTest,
                                      seqs = seqsTest,
                                      position_labels = testPositionLabels)
               , "not of type matrix")

})

test_that("Given features matrix object is matrix", {

  testSamplesMat <- matrix(rnorm(10000), nrow = 200)
  testFeaturesMat <- rnorm(10000) # err
  testPositionLabels <- seq(50)
  nClusterTest <- 5
  clusterSolTest <- list(clust_sol = c(), reordering_idx = c())
  seqsTest <- NULL
  expect_error(plot_arch_for_clusters(testSamplesMat, testFeaturesMat,
                                      nCluster = nClusterTest,
                                      clustering_sol = clusterSolTest,
                                      seqs = seqsTest,
                                      position_labels = testPositionLabels)
               , "not of type matrix")

})

test_that("Handling empty samples matrix", {

  testSamplesMat <- matrix()
  testFeaturesMat <- matrix(rnorm(10000), nrow = 200)
  testPositionLabels <- seq(50)
  nClusterTest <- 5
  clusterSolTest <- list(clust_sol = c(), reordering_idx = c())
    #get_clusters(t(samplesMatrix), clustMethod = "kmeans", final_nClust)
  seqsTest <- NULL
  expect_error(plot_arch_for_clusters(testSamplesMat, testFeaturesMat,
                                      nCluster = nClusterTest,
                                      clustering_sol = clusterSolTest,
                                      seqs = seqsTest,
                                      position_labels = testPositionLabels)
               , "Empty")

})

test_that("Handling empty features matrix", {

  testSamplesMat <- matrix(rnorm(10000), nrow = 200)
  testFeaturesMat <- matrix()
  testPositionLabels <- seq(50)
  nClusterTest <- 5
  clusterSolTest <- list(clust_sol = c(), reordering_idx = c())
  seqsTest <- NULL
  expect_error(plot_arch_for_clusters(testSamplesMat, testFeaturesMat,
                                      nCluster = nClusterTest,
                                      clustering_sol = clusterSolTest,
                                      seqs = seqsTest,
                                      position_labels = testPositionLabels)
               , "Empty")

})


test_that("Handling negative nCluster value", {

  testSamplesMat <- matrix(rnorm(10000), nrow = 200)
  testFeaturesMat <- matrix()
  testPositionLabels <- seq(50)
  nClusterTest <- 5
  clusterSolTest <- list(clust_sol = c(), reordering_idx = c())
  seqsTest <- NULL
  expect_error(plot_arch_for_clusters(testSamplesMat, testFeaturesMat,
                                      nCluster = nClusterTest,
                                      clustering_sol = clusterSolTest,
                                      seqs = seqsTest,
                                      position_labels = testPositionLabels)
               , "Empty")

})

test_that("Handling list as nCluster value", {

  testSamplesMat <- matrix(rnorm(10000), nrow = 200)
  testFeaturesMat <- matrix(rnorm(10000), nrow = 200)
  testPositionLabels <- seq(50)
  nClusterTest <- seq(5)
  clusterSolTest <- list(clust_sol = c(), reordering_idx = c())
  seqsTest <- NULL
  expect_error(plot_arch_for_clusters(testSamplesMat, testFeaturesMat,
                                      nCluster = nClusterTest,
                                      clustering_sol = clusterSolTest,
                                      seqs = seqsTest,
                                      position_labels = testPositionLabels)
               , "Expecting only one value")

})

test_that("Handling negative nCluster value", {

  testSamplesMat <- matrix(rnorm(10000), nrow = 200)
  testFeaturesMat <- matrix(rnorm(10000), nrow = 200)
  testPositionLabels <- seq(50)
  nClusterTest <- -5
  clusterSolTest <- list(clust_sol = c(), reordering_idx = c())
  seqsTest <- NULL
  expect_error(plot_arch_for_clusters(testSamplesMat, testFeaturesMat,
                                      nCluster = nClusterTest,
                                      clustering_sol = clusterSolTest,
                                      seqs = seqsTest,
                                      position_labels = testPositionLabels)
               , "non-negative")

})


test_that("Handling different nCluster value than in clustering_sol", {

  testSamplesMat <- matrix(rnorm(10000), nrow = 50)
  testFeaturesMat <- matrix(rnorm(10000), nrow = 200, ncol = 50)
  testPositionLabels <- seq(50)
  nClusterTest <- 5
  clusterSolTest <- get_clusters(testSamplesMat, clustMethod = "kmeans",
                                 nCluster = nClusterTest+2)
  seqsTest <- NULL
  expect_error(plot_arch_for_clusters(testSamplesMat, testFeaturesMat,
                                      nCluster = nClusterTest,
                                      clustering_sol = clusterSolTest,
                                      seqs = seqsTest,
                                      position_labels = testPositionLabels)
               , "mismatch")

})


