context("Clustering functions - Get nClusters")

test_that("Handles empty matrix", {
  testMat <- matrix() # err
  clustMethodTest <- "kmeans"
  # distMethodTest <- "euclidean"
  nClustersTest <- 5
  expect_error(
    get_clusters(
      givenMat = testMat, clustMethod = clustMethodTest,
      nCluster = nClustersTest
    ),
    "Empty"
  )
})

test_that("Handles non-positive #clusters", {
  testMat <- matrix(rnorm(10000), nrow = 200)
  clustMethodTest <- "kmeans"
  # distMethodTest <- "euclidean"
  nClustersTest <- -5
  expect_error(
    get_clusters(
      givenMat = testMat, clustMethod = clustMethodTest,
      nCluster = nClustersTest
    ),
    "Ask for at least 1 clusters"
  )
})

test_that("Handles #clusters > #sequences", {
  testMat <- matrix(rnorm(10000), nrow = 200)
  clustMethodTest <- "kmeans"
  # distMethodTest <- "euclidean"
  nClustersTest <- ncol(testMat) + 5
  expect_error(
    get_clusters(
      givenMat = testMat, clustMethod = clustMethodTest,
      nCluster = nClustersTest
    ),
    "nClusters more than #sequences"
  )
})

test_that("Handles a list as nCluster value", {
  testMat <- matrix(rnorm(10000), nrow = 200)
  clustMethodTest <- "kmeans"
  # distMethodTest <- "euclidean"
  nClustersTest <- seq(3, 5)
  expect_error(
    get_clusters(
      givenMat = testMat, clustMethod = clustMethodTest,
      nCluster = nClustersTest
    ),
    "Expecting only one value"
  )
})

test_that("Handles wrongly given clustMethod", {
  testMat <- matrix(rnorm(10000), nrow = 200)
  clustMethodTest <- "kmns" # err
  # distMethodTest <- "euclidean"
  nClustersTest <- 5
  expect_error(
    get_clusters(
      givenMat = testMat, clustMethod = clustMethodTest,
      nCluster = nClustersTest
    ),
    "Wrong clustMethod"
  )
})

# test_that("Handles wrongly given distMethod with hclust", {
#   testMat <- matrix(rnorm(10000), nrow = 200)
#   clustMethodTest <- "hclust"
#   distMethodTest <- "" # err, with hclust
#   nClustersTest <- 5
#   expect_error(
#     get_clusters(
#       givenMat = testMat, clustMethod = clustMethodTest,
#       nCluster = nClustersTest,
#       distMethod.hclust = distMethodTest
#     ),
#     "Wrong distMethod.hclust"
#   )
# })


test_that("Returned List element names are right.", {
  #
  # with kmeans
  testMat <- matrix(rnorm(10000), nrow = 200)
  clustMethodTest <- "kmeans"
  nClustersTest <- 5
  givenAns <- get_clusters(
    givenMat = testMat, clustMethod = clustMethodTest,
    nCluster = nClustersTest
  )
  testAnsNames <- c("clustType", "clust_sol", "reordering_idx")
  expect_equal(names(givenAns), testAnsNames)
  #
  # # with hclust
  # testMat <- matrix(rnorm(10000), nrow = 200)
  # clustMethodTest <- "hclust"
  # distMethodTest <- "euclidean"
  # nClustersTest <- 5
  # givenAns <- get_clusters(
  #   givenMat = testMat, clustMethod = clustMethodTest,
  #   nCluster = nClustersTest,
  #   distMethod.hclust = distMethodTest
  # )
  # testAnsNames <- c("clust_sol", "reordering_idx")
  # expect_equal(names(givenAns), testAnsNames)
})
