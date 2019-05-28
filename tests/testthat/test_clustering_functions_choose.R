context("Clustering functions - Choosing nClusters")

test_that("Choosing clusters handles wrong clustMethod", {

  testMat <- matrix(rnorm(10000), nrow = 200)
  clustMethodTest <- "kmeas" # err
  distMethodTest <- "euclidean"
  expect_error(choose_clusters(testMat, distMethod = distMethodTest,
                               clustMethod = clustMethodTest,
                               nCluster_vals_test = seq(3,5))
               , "Wrong clustMethod"
               )

})
#
test_that("Choosing clusters handles wrong distMethod", {
  testMat <- matrix(rnorm(10000), nrow = 200)
  clustMethodTest <- "kmeans"
  distMethodTest <- "eulcidean" # err
  expect_error(choose_clusters(testMat, distMethod = distMethodTest,
                               clustMethod = clustMethodTest,
                               nCluster_vals_test = seq(3,5))
               , "Wrong distMethod"
              )

})
#
test_that("Choosing clusters handles empty matrix", {
  testMat <- matrix() # err
  clustMethodTest <- "kmeans"
  distMethodTest <- "euclidean"
  expect_error(choose_clusters(testMat, distMethod = distMethodTest,
                               clustMethod = clustMethodTest,
                               nCluster_vals_test = seq(3,5))
               , "Empty"
  )
})

test_that("Choosing clusters handles non-positive #clusters", {
  testMat <- matrix(rnorm(10000), nrow = 200)
  clustMethodTest <- "kmeans"
  distMethodTest <- "euclidean"
  nClustersTest <- -5
  expect_error(choose_clusters(testMat, distMethod = distMethodTest,
                               clustMethod = clustMethodTest,
                               nCluster_vals_test = nClustersTest)
                , "at least 2 clusters"
  )
})

test_that("Choosing clusters handles #clusters > #sequences", {
  testMat <- matrix(rnorm(10000), nrow = 200)
  clustMethodTest <- "kmeans"
  distMethodTest <- "euclidean"
  nClustersTest <- ncol(testMat) + 5
  expect_error(choose_clusters(testMat, distMethod = distMethodTest,
                               clustMethod = clustMethodTest,
                               nCluster_vals_test = seq(nClustersTest-10,
                                                        nClustersTest))
               , "nClusters more than #sequences"
  )
})




# test_that("Choosing clusters handles kmeans well", {
#   testMat <- matrix(rnorm(10000), nrow = 200)
#   clustMethodTest <- "kmeans"
#   distMethodTest <- "euclidean"
#   nClustersTest <- seq(3,5)
#   #
#   #
#   expect_error(choose_clusters(testMat, distMethod = distMethodTest,
#                                clustMethod = clustMethodTest,
#                                nCluster_vals_test = nClustersTest)
#                , "nClusters more than #sequences"
#   )
# })
#
# test_that("Choosing clusters handles hclust well", {
#   testMat <- matrix(rnorm(10000), nrow = 200)
#   clustMethodTest <- "hclust"
#   distMethodTest <- "euclidean"
#   nClustersTest <- seq(3,5)
#   #
#   #
#   expect_error(choose_clusters(testMat, distMethod = distMethodTest,
#                                clustMethod = clustMethodTest,
#                                nCluster_vals_test = nClustersTest)
#                , "nClusters more than #sequences"
#   )
# })
